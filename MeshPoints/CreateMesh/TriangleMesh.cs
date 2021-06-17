using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using MeshPoints.Classes;
using Rhino.Geometry.Collections;

namespace MeshPoints.CreateMesh
{


    public class TriangleMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateTriangleMesh class.
        /// </summary>
        public TriangleMesh()
          : base("Triangle Mesh", "TriMesh",
              "Creates a triangle mesh on a planar brep. Constrained Laplacian Smooth are performed to improve mesh quality.",
              "SmartMesh", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Brep", "brep", "Planar brep.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Edge Node Count", "count", "Wanted amount of edge nodes. Note that this is only a target value.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Inner Nodes", "points", "A list with points inside the brep.", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "mesh", "Mesh (triangle-elements).", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Inputs
            Brep meshSurface = new Brep();
            double totalEdgeNodeCount = 0;
            List<Point3d> innerNodeList = new List<Point3d>();
            DA.GetData(0, ref meshSurface);
            DA.GetData(1, ref totalEdgeNodeCount);
            DA.GetDataList(2, innerNodeList);

            // 1. Create nodes along the edges of the surface and flatten.
            List<List<Point3d>> edgeNodesSurface = CreateEdgePointsByCount(meshSurface, totalEdgeNodeCount);

            // 2. Flatten the list here, as the ListList-structure of CreateEdgePointsByCount is useful elsewhere.
            List<Point3d> flattenedEdgeNodes = new List<Point3d>();
            foreach(List<Point3d> edge in edgeNodesSurface)
            {
                // Do not add endnode of edge as this is duplicate of startnode of next edge
                for ( int i = 0; i<edge.Count()-1; i++)
                {
                    flattenedEdgeNodes.Add(edge[i]);
                }
            }

            // 3. Create points inside the surface by creating a bounding box, populating it with points, and culling all points not inside the surface.
            Brep boundingBoxSurface = CreateBoundingBoxFromBrep(meshSurface);
            if (boundingBoxSurface == null) { return; }
            List<Point3d> nodeGridBoundingBox = innerNodeList;
            List<Point3d> nodesInsideSurface = CullPointsOutsideSurface(nodeGridBoundingBox, flattenedEdgeNodes, meshSurface);


            // 4. Collect flat list of all points to use for triangle meshing and cast to compatible data structure (Node2List) for Delaunay method.
            List<Point3d> nodeCollection = new List<Point3d>();
            nodeCollection.AddRange(flattenedEdgeNodes);
            nodeCollection.AddRange(nodesInsideSurface);
            Point3d[] culledNodes = Point3d.CullDuplicates(nodeCollection, 0.001);
            var meshNodes = new Grasshopper.Kernel.Geometry.Node2List(culledNodes);

            // 5. Throw all our points into the Delaunay mesher. Adjust jitter_amount as needed.
            var meshFaces = new List<Grasshopper.Kernel.Geometry.Delaunay.Face>();
            var triangleMesh = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Mesh(meshNodes, 0.01, ref meshFaces); // todo: what is "double jitter_amount"?
            
            // 6. Sometimes the mesh acts up; in these cases it is necessary to cull mesh faces that are outside the surface.
            Mesh culledTriangleMesh = CullMeshFacesOutsideSurface(triangleMesh, meshSurface);
            
            if (!culledTriangleMesh.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The mesh is invalid. Check for duplicate points.");
            }

            // 7. Get initial edges and elements of mesh using mesh topology properties
            var initialEdgeAndElementList = GetInitialEdgesAndElements(culledTriangleMesh);
            List<qEdge> globalEdgeList = initialEdgeAndElementList.Item1;
            List<qElement> globalElementList = initialEdgeAndElementList.Item2;
            foreach (qElement element in globalElementList)
            {
                element.FixEdgeOrder();
            }

            // 8. Do global smoothing
            DoGlobalSmoothing(meshSurface, globalEdgeList, globalElementList);

            // 9. Convert mesh to SmartMesh Class and create final mesh
            var meshProperties = ConvertToMainMeshClasses(globalElementList);
            Mesh newMesh = new Mesh();
            foreach (Element e in meshProperties.Item2)
            {
                Mesh mesh = new Mesh();
                mesh.Vertices.Add(e.Nodes[0].Coordinate); 
                mesh.Vertices.Add(e.Nodes[1].Coordinate); 
                mesh.Vertices.Add(e.Nodes[2].Coordinate); 
                mesh.Faces.AddFace(0, 1, 2);
                mesh.Normals.ComputeNormals();  //Control if needed
                mesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                mesh.Compact(); //to ensure that it calculate

                //create global mesh
                newMesh.Append(mesh);
            }
            newMesh.Weld(0.1);

            // Outputs
            DA.SetData(0, newMesh);
        }
        #region Methods
        /// <summary>
        /// Cull unwanted mesh faces by checking if their center points are outside the actual surface of the mesh.
        /// </summary>
        /// <returns>A <see cref="Mesh"/> with (hopefully) no outside mesh faces.</returns>
        private Mesh CullMeshFacesOutsideSurface(Mesh meshSurface, Brep brep)
        {
            Mesh insideFaces = meshSurface.DuplicateMesh();
            for (int i = meshSurface.Faces.Count-1; i>0; i--) // reverse iteration to maintain indices
            {
                if (!IsPointOnBrepSurface(meshSurface.Faces.GetFaceCenter(i), brep))
                {
                    insideFaces.Faces.RemoveAt(i);
                }
            }
            return insideFaces;
        }

        /// <summary>
        /// Takes an input point and a Brep surface. If the distance between input point 
        /// and the closest point on the Brep ~ 0, the point is deemed on the surface.
        /// </summary>
        /// <returns>True if point is on Brep.</returns>
        private bool IsPointOnBrepSurface(Point3d point, Brep brep)
        {
            var testPointSurfaceDistance = point.DistanceTo(brep.ClosestPoint(point));

            if (testPointSurfaceDistance < RhinoMath.SqrtEpsilon) 
            { 
                return true; 
            }
            else 
            { 
                return false; 
            }
        }

        /// <summary>
        /// Takes a list of <see cref="Point3d"/> and checks if points are inside an input <see cref="Brep"/> surface.
        /// </summary>
        /// <returns>A list of <see cref="Point3d"/> containing points inside the input <see cref="Brep"/> surface.</returns>
        private List<Point3d> CullPointsOutsideSurface(List<Point3d> pointGrid, List<Point3d> edgeNodes, Brep meshSurface)
        {
            var insidePoints = new List<Point3d>();
            HashSet<Point3d> edgeNodesSet = new HashSet<Point3d>(edgeNodes);
            foreach (Point3d point in pointGrid)
            {
                if (IsPointOnBrepSurface(point, meshSurface) && !edgeNodesSet.Contains(point))
                {

                    insidePoints.Add(point);
                }
            }
            return insidePoints;
        }

        /// <summary>
        /// Creates a <see cref="Brep"/> rectangle around an arbitrary plane <see cref="Brep"/> geometry and fills it with points.
        /// </summary>
        private Brep CreateBoundingBoxFromBrep(Brep meshSurface)
        {
            double zAxis = 0.0; // double representing the points' placement on the z-axis 
            BoundingBox boundingBox = meshSurface.GetBoundingBox(false);
            Brep boundingBoxBrep = Brep.CreateFromCornerPoints(
                new Point3d(boundingBox.Min.X, boundingBox.Min.Y, zAxis),
                new Point3d(boundingBox.Max.X, boundingBox.Min.Y, zAxis),
                new Point3d(boundingBox.Max.X, boundingBox.Max.Y, zAxis),
                new Point3d(boundingBox.Min.X, boundingBox.Max.Y, zAxis),
                RhinoMath.ZeroTolerance
                );

            return boundingBoxBrep;
        }

        /// <summary>
        /// Creates points on the <see cref="BrepEdge"/>s of a <see cref="Brep"/> based on a given total edge node count.
        /// </summary>
        /// <returns>Returnds a list of <see cref="Point3d"/> along the edge of a <see cref="Brep"/>.</returns>
        private List<List<Point3d>> CreateEdgePointsByCount(Brep meshSurface, double totalEdgeNodeCount)
        {
            var edgePoints = new List<List<Point3d>>();

            double totalEdgeLength = 0;
            foreach (Curve edge in meshSurface.Edges)
            {
                totalEdgeLength += edge.GetLength();
            }
            foreach (Curve edge in meshSurface.Edges)
            {
                double[] tValues;
                var innerEdgePoints = new List<Point3d>();
                double edgeLength = edge.GetLength();
                var edgeNodeCount = Convert.ToInt32(totalEdgeNodeCount * (edgeLength / totalEdgeLength) + 1);
                
                tValues = edge.DivideByCount(edgeNodeCount, true);
                foreach ( var t in tValues)
                {
                    innerEdgePoints.Add(edge.PointAt(t));
                }
                edgePoints.Add(innerEdgePoints);
            }
            return edgePoints;
        }

        /// <summary>
        /// Creates a grid of points in a bounding box <see cref="Brep"/> based on a given number of wanted internal nodes in a surface.
        /// </summary>
        /// <returns>A List of <see cref="Point3d"/> describing a grid of points in a rectangle.</returns>
        private Tuple<List<Node>, List<Element>> ConvertToMainMeshClasses(List<qElement> globalqElementList)
        {
            // Create global qNode list
            List<qNode> globalqNodeList = new List<qNode>();
            foreach (qElement qrElement in globalqElementList)
            {
                List<qNode> qrElementNodes = qrElement.GetNodesOfElement();
                foreach (qNode qrNode in qrElementNodes)
                {
                    if (!globalqNodeList.Contains(qrNode))
                    {
                        globalqNodeList.Add(qrNode);
                    }
                }
            }

            // Get node list
            List<Node> nodes = new List<Node>();
            Node newNode = new Node();
            for (int i = 0; i < globalqNodeList.Count; i++)
            {
                qNode qrNode = globalqNodeList[i];

                // Assign properties
                if (qrNode.BoundaryNode)
                {
                    newNode = new Node(i, qrNode.Coordinate, true, true);
                }
                else
                {
                    newNode = new Node(i, qrNode.Coordinate, false, false);
                }
                nodes.Add(newNode);
            }

            // Get element list
            List<Element> elements = new List<Element>();
            for (int i = 0; i < globalqElementList.Count; i++)
            {

                List<qNode> qrElementNodes = globalqElementList[i].GetNodesOfElement();
                List<int> connectivity = new List<int>();

                // Get element connectivity:
                foreach (qNode qrNode in qrElementNodes)
                {
                    int id = globalqNodeList.IndexOf(qrNode);
                    connectivity.Add(id);
                }

                // Create elements
                List<Node> elementNodes = new List<Node>();
                Mesh elementMesh = new Mesh();

                foreach (int index in connectivity)
                {
                    elementNodes.Add(nodes[index]);
                    elementMesh.Vertices.Add(nodes[index].Coordinate);
                }
                Element element = new Element(i, elementNodes, connectivity);

                // Assign mesh
                if (element.Type == "Triangle")
                {
                    elementMesh.Faces.AddFace(0, 1, 2);
                }
                else
                {
                    elementMesh.Faces.AddFace(0, 1, 2, 3);
                }
                elementMesh.Faces.AddFace(0, 1, 2);
                element.Mesh = elementMesh;

                // Add element to element list
                elements.Add(element);
            }
            return Tuple.Create(nodes, elements);
        }
        private Tuple<List<qEdge>, List<qElement>> GetInitialEdgesAndElements(Mesh mesh)
        {
            // summary: create edges and elements from the initial mesh topology
            List<qEdge> edgeList = CreateInitialEdges(mesh);
            List<qElement> elementList = CreateInitialElements(mesh, edgeList);
            SetNeighborElements(mesh, elementList, edgeList);
            /* to do: fixslett
            foreach (qElement element in elementList)
            {
                element.FixElementEdgeAndAngle();
            }*/
            return Tuple.Create(edgeList, elementList);
        }
        private List<qNode> CreateInitalNodes(Mesh mesh)
        {
            // summary: create nodes from initial mesh
            List<qNode> nodeList = new List<qNode>();
            MeshTopologyVertexList topologyVertexList = mesh.TopologyVertices;
            for (int i = 0; i < topologyVertexList.Count; i++)
            {
                qNode node = new qNode();
                int meshVertexIndex = topologyVertexList.MeshVertexIndices(i)[0];
                node.Coordinate = mesh.Vertices.Point3dAt(meshVertexIndex);
                nodeList.Add(node);
            }

            List<Polyline> nakedEdges = mesh.GetNakedEdges().ToList();
            for (int i = 0; i < nodeList.Count; i++)
            {
                Point3d point = nodeList[i].Coordinate;
                bool isBoundaryNode = false;
                foreach (Polyline line in nakedEdges)
                {
                    NurbsCurve edge = line.ToNurbsCurve();
                    edge.ClosestPoint(nodeList[i].Coordinate, out double parameter);
                    if (point.DistanceTo(edge.PointAt(parameter)) <= 0.001)
                    {
                        isBoundaryNode = true;
                    }
                    else { nodeList[i].BoundaryNode = false; }
                }
                if (isBoundaryNode) { nodeList[i].BoundaryNode = true; }
                else { nodeList[i].BoundaryNode = false; }
            }
            return nodeList;
        }
        private List<qEdge> CreateInitialEdges(Mesh mesh)
        {
            // summary: create edges from intial mesh
            List<qEdge> edgeList = new List<qEdge>();
            MeshTopologyEdgeList topologyEdgeList = mesh.TopologyEdges;
            List<qNode> nodeList = CreateInitalNodes(mesh);

            for (int i = 0; i < topologyEdgeList.Count; i++)
            {
                Rhino.IndexPair edgeTopoVerticesIndex = topologyEdgeList.GetTopologyVertices(i);
                qNode startNode = nodeList[edgeTopoVerticesIndex[0]];
                qNode endNode = nodeList[edgeTopoVerticesIndex[1]];

                qEdge edge = new qEdge(startNode, endNode);
                edgeList.Add(edge);
            }
            return edgeList;
        }
        private List<qElement> CreateInitialElements(Mesh mesh, List<qEdge> globalEdgeList)
        {
            // summary: create elements from initial mesh
            qElement element = new qElement();
            List<qElement> elementList = new List<qElement>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                List<qEdge> edgeListElement = new List<qEdge>();

                int[] elementEdgesIndex = mesh.TopologyEdges.GetEdgesForFace(i);
                foreach (int n in elementEdgesIndex)
                {
                    edgeListElement.Add(globalEdgeList[n]);
                }
                element = new qElement(edgeListElement);
                elementList.Add(element);
            }
            return elementList;
        }
        private void SetNeighborElements(Mesh mesh, List<qElement> globalElementList, List<qEdge> globalEdgeList)
        {
            // summary: assign neighbor elements to edges. Limited to initial mesh.
            for (int i = 0; i < globalEdgeList.Count; i++)
            {
                int[] connectedElementsIndex = mesh.TopologyEdges.GetConnectedFaces(i);

                // warning message if topology is may contain an error
                if (connectedElementsIndex.Length == 3)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "SetNeighborElements: Mesh topology may be wrong. Check corner edges if connected elements are correct.");
                }

                globalEdgeList[i].Element1 = globalElementList[connectedElementsIndex[0]];

                if (connectedElementsIndex.Length == 2)
                {
                    globalEdgeList[i].Element2 = globalElementList[connectedElementsIndex[1]];
                }
            }
            return;
        }
        private Tuple<List<qEdge>, List<qEdge>> UpdateGlobalEdgeList_NodePosition(qNode smoothNode, Point3d newCoordinate, List<qEdge> globalEdgeList)
        {
            // summary: modify node location of oldNode to smoothNode location and update connected edges
            // silje comment: sjekk naboelementer. Blir deres edgelist oppdatert?
            List<qEdge> oldEdges = new List<qEdge>();
            List<qEdge> newEdges = new List<qEdge>();
            List<qEdge> globalEdgeListCopy = new List<qEdge>(globalEdgeList);
            List<qEdge> connectedEdges = smoothNode.GetConnectedEdges(globalEdgeListCopy);

            //for (int i = 0; i < connectedEdges.Count; i++) // silje comment: fjerne denne loopen?
            //{
            foreach (qEdge edge in connectedEdges)
            {
                int id = globalEdgeListCopy.IndexOf(edge);

                if (globalEdgeListCopy[id].StartNode == smoothNode)
                {
                    globalEdgeList[id].StartNode.Coordinate = newCoordinate;
                }
                else if (globalEdgeListCopy[id].EndNode == smoothNode)
                {
                    globalEdgeList[id].EndNode.Coordinate = newCoordinate;
                }

                globalEdgeList[id].Length = globalEdgeList[id].CalculateLength(globalEdgeList[id].StartNode, globalEdgeList[id].EndNode);
                globalEdgeList[id].EdgeLine = globalEdgeList[id].VisualizeLine(globalEdgeList[id].StartNode, globalEdgeList[id].EndNode);

                int edgeId1 = globalEdgeListCopy[id].Element1.EdgeList.IndexOf(edge);
                globalEdgeList[id].Element1.EdgeList[edgeId1] = globalEdgeList[id];
                globalEdgeList[id].Element1.GetContourOfElement();
                globalEdgeList[id].Element1.CalculateAngles(); // to do:isquad?

                if (!smoothNode.BoundaryNode)
                {
                    int edgeId2 = globalEdgeListCopy[id].Element2.EdgeList.IndexOf(edge);
                    globalEdgeList[id].Element2.EdgeList[edgeId2] = globalEdgeList[id];
                    globalEdgeList[id].Element2.GetContourOfElement();
                    globalEdgeList[id].Element2.CalculateAngles(); // to do: isquad?
                }

                newEdges.Add(globalEdgeList[id]);
                oldEdges.Add(globalEdgeListCopy[id]);
                //todo: else { add runtimemessage }
            }
            //}
            return Tuple.Create(oldEdges, newEdges);
        }
        private void UpdateGlobalElementList_ChangedEdges(List<qEdge> newEdges, List<qEdge> oldEdges, List<qElement> globalElementList)
        {
            // make list with old and new edges at index where edges are changed.
            List<qElement> globalElementListOld = new List<qElement>(globalElementList);

            foreach (qElement oldElement in globalElementListOld)
            {
                for (int i = 0; i < oldEdges.Count; i++)
                {
                    int elementId = globalElementListOld.IndexOf(oldElement);
                    qEdge oldEdge = oldEdges[i];
                    if (oldElement.EdgeList.Contains(oldEdge))
                    {
                        int edgeId = oldElement.EdgeList.IndexOf(oldEdge);

                        globalElementList[elementId].EdgeList[edgeId] = newEdges[i];

                        globalElementList[elementId].GetContourOfElement();
                        globalElementList[elementId].CalculateAngles(); // to do: isquad?

                    }
                }
            }
        }
        private void UpdateSmoothableNodeList(qNode node, List<qNode> smoothableNodes, List<qEdge> globalEdgeList)
        {
            List<qEdge> connectedEdges = node.GetConnectedEdges(globalEdgeList);
            foreach (qEdge edge in connectedEdges)
            {
                if (!smoothableNodes.Contains(edge.GetOppositeNode(node)))
                {
                    smoothableNodes.Add(edge.GetOppositeNode(node));
                    edge.GetOppositeNode(node).OBS = false; //todo: er dette rett??
                }
            }
        }
        private void UpdateDistortionMetric(qNode node, List<qEdge> globalEdgeList)
        {
            List<qElement> connectedElements = node.GetConnectedElements(globalEdgeList);
            foreach (qElement element in connectedElements)
            {
                element.DistortionMetric = CalculateDistortionMetric(element);
            }
        }
        private List<qNode> GetGlobalNodeList(List<qEdge> globalEdgeList)
        {
            List<qNode> globalNodeList = new List<qNode>();
            foreach (qEdge edge in globalEdgeList)
            {
                if (!globalNodeList.Contains(edge.StartNode)) { globalNodeList.Add(edge.StartNode); }
                if (!globalNodeList.Contains(edge.EndNode)) { globalNodeList.Add(edge.EndNode); }
            }
            return globalNodeList;
        }

        //
        // The implementation of Global Smoothing are based on the paper:
        // "An approach to Combined Laplacian and Optimization-Based Smoothing for Triangular, Quadrilateral and
        // Quad-Dominant Meshes" (1998) by Cannan, Tristano, and Staten
        // 
        // In addition, modification and assumtions by the paper by Karl Erik Levik are implementet:
        // "Q-Morph - Implementing a Quadrilateral Meshing Algorithm" (2002) by Levik, Karl Erik
        // 
        private void DoGlobalSmoothing(Brep brep, List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            List<qNode> globalNodeList = GetGlobalNodeList(globalEdgeList);
            List<qNode> smoothableNodes = new List<qNode>(globalNodeList);
            double moveTolerance = 0.01; // constant: fix: todo.
            double maxModelDimension = 0;

            // 0. Make list with brep vertex coordinates
            List<Point3d> vertexPoints = new List<Point3d>();
            foreach (BrepVertex vertex in brep.Vertices)
            {
                vertexPoints.Add(vertex.Location);
            }

            // 1. Calculate initial distortion metrics for all elements
            foreach (qElement element in globalElementList)
            {
                element.DistortionMetric = CalculateDistortionMetric(element);
            }

            // 2. Calculate maximum model dimension (assume qmorph <- todo: skriv bedre)
            foreach (qElement element in globalElementList)
            {
                foreach (qEdge edge in element.EdgeList)
                {
                    double edgeLength = edge.Length;
                    if (edgeLength > maxModelDimension) { maxModelDimension = edgeLength; }
                }
            }

            // 3. Set iteration to 1
            int iterations = 1;

            // 4. Start to smooth nodes
            int n = 0;

            bool continueSmooth = true;
            while (continueSmooth & n < 100)
            {
                bool nodeIsMoved = false;
                double maxDistanceMoved = 0;
                List<qNode> smoothableNodesCopy = new List<qNode>(smoothableNodes);
                foreach (qNode node in smoothableNodesCopy)
                {
                    if (!smoothableNodesCopy.Contains(node) | vertexPoints.Contains(node.Coordinate)) { smoothableNodes.Remove(node); continue; } // todo: forlag: ta BC-edge og flytt vektor.
                    if (iterations == 1) { node.OBS = false; }


                    // 4.1. Perform Constrained Laplacian Smooth
                    if (!node.OBS)
                    {
                        qNode movedNode = new qNode();
                        // 4.1.0. Perform Laplacian Smooth on Boundary Nodes or perform Constrained Laplacian Smooth
                        if (node.BoundaryNode)
                        {
                            List<qEdge> boundaryEdge = new List<qEdge>();
                            Vector3d laplacianVector = LaplacianSmooth(node, globalEdgeList);
                            Point3d newPoint = new Point3d(node.Coordinate.X + laplacianVector.X, node.Coordinate.Y + laplacianVector.Y, node.Coordinate.Z + laplacianVector.Z);
                            BrepEdge edge = PointOnEdge(node, brep);
                            edge.ClosestPoint(newPoint, out double parameter);
                            movedNode = new qNode(edge.PointAt(parameter), node.BoundaryNode);
                        }
                        else
                        {
                            // 4.1.1. Move node
                            movedNode = ConstrainedLaplacianSmooth(node, globalEdgeList, globalElementList);
                        }

                        // 4.1.2. Check if move shall be allowed;
                        double distanceMoved = (movedNode.Coordinate - node.Coordinate).Length;
                        if (distanceMoved < moveTolerance)
                        {
                            movedNode = new qNode(node.Coordinate, node.BoundaryNode); // reset node to old node
                            smoothableNodes.Remove(node);
                        }
                        else
                        {
                            var update = UpdateGlobalEdgeList_NodePosition(node, movedNode.Coordinate, globalEdgeList);
                            UpdateGlobalElementList_ChangedEdges(update.Item1, update.Item2, globalElementList);
                            UpdateSmoothableNodeList(node, smoothableNodes, globalEdgeList);
                            UpdateDistortionMetric(node, globalEdgeList);
                            if (distanceMoved > maxDistanceMoved) { maxDistanceMoved = distanceMoved; }
                            nodeIsMoved = true;
                        }
                    }
                    /*
                    // 4.2 Perform Optimization-based Smoothing
                    if (iterations >= 2)
                    {
                        // 4.2.1 Find minimum distortion metric from adjacent elements to node:
                        List<qElement> connectedElements = node.GetConnectedElements(globalEdgeList);
                        double minDistortion = 100;
                        foreach (qElement element in connectedElements)
                        {
                            if (element.DistortionMetric < minDistortion) { minDistortion = element.DistortionMetric; }
                        }

                        // 4.2.2. Check if OBS shall be performed:
                        if (minDistortion <= 0.1) // constant proposed by paper
                        {
                            qNode movedNode = OptimizationBasedSmoothing(node, maxModelDimension, globalEdgeList);
                            double distanceMoved = (movedNode.Coordinate - node.Coordinate).Length;

                            if (movedNode.OBS)
                            {
                                node.OBS = true;
                                var update = UpdateGlobalEdgeList_NodePosition(node, movedNode.Coordinate, globalEdgeList);
                                UpdateGlobalElementList_ChangedEdges(update.Item1, update.Item2, globalElementList);
                                UpdateSmoothableNodeList(node, smoothableNodes, globalEdgeList);
                                UpdateDistortionMetric(node, globalEdgeList);
                                if (distanceMoved > maxDistanceMoved) { maxDistanceMoved = distanceMoved; }
                                nodeIsMoved = true;
                            }
                            else { node.OBS = false; }
                        }
                    }*/
                }
                iterations++;
                if (!nodeIsMoved | maxDistanceMoved <= 1.75 * moveTolerance) { continueSmooth = false; }
                n++;
            }
        }
        private double CalculateDistortionMetric(qElement element)
        {
            double my = 0; // distortion metric
            if (!element.IsQuad)
            {
                my = DistorionMetricTriangle(0, 1, 2, element);
            }
            else
            {
                List<double> alphas = new List<double>();
                int invertedTriangles = 0;
                for (int i = 0; i < 4; i++)
                {
                    if (i == 0)
                    {
                        double a = DistorionMetricTriangle(0, 1, 3, element);
                        alphas.Add(a);
                        if (a < 0) { invertedTriangles++; }
                    }
                    else if (i == 1)
                    {
                        double a = DistorionMetricTriangle(1, 2, 3, element);
                        alphas.Add(a);
                        if (a < 0) { invertedTriangles++; }
                    }
                    else if (i == 2)
                    {
                        double a = DistorionMetricTriangle(0, 1, 2, element);
                        alphas.Add(a);
                        if (a < 0) { invertedTriangles++; }
                    }
                    else
                    {
                        double a = DistorionMetricTriangle(0, 2, 3, element);
                        alphas.Add(a);
                        if (a < 0) { invertedTriangles++; }
                    }
                }

                // Set negval
                double negval = 0;
                double coincidentTolerance = 0.01; // Constant: no background for selecting this value. todo: fix tolerance
                List<double> angles = element.AngleList;

                if (angles.Min() < (double)6 * Math.PI / (double)180 | CoincidentNodes(element) < coincidentTolerance | invertedTriangles == 2) { negval = 1; }
                else if (invertedTriangles == 3) { negval = 2; }
                else if (invertedTriangles == 4) { negval = 3; }
                else { negval = 0; }

                // Calculat Distortion metric:
                my = alphas.Min() - negval;
            }

            return my;
        }
        private double DistorionMetricTriangle(int node1, int node2, int node3, qElement element)
        {
            List<qNode> elementNodes = element.GetNodesOfElement();
            double alpha = 0;
            Vector3d AB = elementNodes[node2].Coordinate - elementNodes[node1].Coordinate;
            Vector3d CB = elementNodes[node2].Coordinate - elementNodes[node3].Coordinate;
            Vector3d CA = elementNodes[node1].Coordinate - elementNodes[node3].Coordinate;

            if (!element.IsQuad)
            {
                // Calculate alpha
                if (!element.IsInverted())
                {
                    alpha = 2 * Math.Sqrt(3) * Vector3d.CrossProduct(CA, CB).Length / (Math.Pow(CA.Length, 2) + Math.Pow(AB.Length, 2) + Math.Pow(CB.Length, 2));
                }
                else
                {
                    alpha = -2 * Math.Sqrt(3) * Vector3d.CrossProduct(CA, CB).Length / (Math.Pow(CA.Length, 2) + Math.Pow(AB.Length, 2) + Math.Pow(CB.Length, 2));
                }
            }
            else
            {
                // Check if triangle from quads is inverted
                Point3d A = elementNodes[node1].Coordinate;
                Point3d B = elementNodes[node2].Coordinate;
                Point3d C = elementNodes[node3].Coordinate;
                bool isInverted = false;
                double area = 0.5 * (A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y));
                if (area <= 0) { isInverted = true; }

                // Calculate alpha
                if (!isInverted)
                {
                    alpha = 4 * Vector3d.CrossProduct(CA, CB).Length / (Math.Pow(CA.Length, 2) + Math.Pow(AB.Length, 2) + Math.Pow(CB.Length, 2));
                }
                else
                {
                    alpha = -4 * Vector3d.CrossProduct(CA, CB).Length / (Math.Pow(CA.Length, 2) + Math.Pow(AB.Length, 2) + Math.Pow(CB.Length, 2));
                }
            }
            return alpha;
        }
        private double CoincidentNodes(qElement element)
        {
            List<qNode> elementNodes = element.GetNodesOfElement();
            List<double> nodeDistance = new List<double>();
            nodeDistance.Add((elementNodes[1].Coordinate - elementNodes[0].Coordinate).Length);
            nodeDistance.Add((elementNodes[2].Coordinate - elementNodes[0].Coordinate).Length);
            nodeDistance.Add((elementNodes[3].Coordinate - elementNodes[0].Coordinate).Length);
            nodeDistance.Add((elementNodes[2].Coordinate - elementNodes[1].Coordinate).Length);
            nodeDistance.Add((elementNodes[3].Coordinate - elementNodes[1].Coordinate).Length);
            nodeDistance.Add((elementNodes[3].Coordinate - elementNodes[2].Coordinate).Length);
            return nodeDistance.Min();
        }
        private qNode ConstrainedLaplacianSmooth(qNode node, List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            // 1. Move node with Laplacian smooth
            Vector3d laplacianVector = LaplacianSmooth(node, globalEdgeList);
            Point3d newPoint = new Point3d(node.Coordinate.X + laplacianVector.X, node.Coordinate.Y + laplacianVector.Y, node.Coordinate.Z + laplacianVector.Z);
            qNode newNode = new qNode(newPoint, node.BoundaryNode);

            // 2. Loop to find final node position
            for (int i = 0; i < 20; i++) // constant proposed in paper: 20
            {
                // 2.1 Defining variables used to check if move is acceptable:
                int posN = 0;
                int negN = 0;
                int upN = 0;
                int downN = 0;
                int invN = 0;
                double theta = 0;
                double thetaMax = 200 * Math.PI / 180; // constant chosen
                double deltaDistMetric = 0;

                // 2.1. Update edge and element list
                List<qElement> connectedElements = node.GetConnectedElements(globalEdgeList);
                int N = connectedElements.Count;

                // 2.2. Calculate Acceptance Criteria for each element 
                foreach (qElement element in connectedElements)
                {
                    List<qNode> elementNodes = element.GetNodesOfElement();
                    elementNodes[elementNodes.IndexOf(node)] = newNode;
                    qElement newElement = CreateElementFromNodes(elementNodes);

                    double newDisMetric = CalculateDistortionMetric(newElement);
                    double oldDisMetric = element.DistortionMetric;

                    if (newDisMetric > oldDisMetric)
                    {
                        posN++;
                    }
                    else if (newDisMetric < oldDisMetric)
                    {
                        negN++;
                    }

                    if ((oldDisMetric < 0 & newDisMetric >= 0) | (oldDisMetric < 0 & newDisMetric > oldDisMetric) | oldDisMetric < 0.05 & newDisMetric >= 0.05) // constant of 0.05 as proposed in paper 
                    {
                        upN++;
                    }
                    else if ((oldDisMetric >= 0 & newDisMetric < 0) | (oldDisMetric < 0 & newDisMetric < oldDisMetric) | oldDisMetric >= 0.05 & newDisMetric < 0.05) // constant of 0.05 as proposed in paper 
                    {
                        downN++;
                    }
                    if (newElement.IsInverted())
                    {
                        invN++;
                    }
                    double testAngle = newElement.AngleList.Max();
                    if (testAngle > theta) { theta = testAngle; }

                    deltaDistMetric = (newDisMetric - oldDisMetric) + deltaDistMetric;
                }

                // 2.4. Check node location
                deltaDistMetric = deltaDistMetric / N;
                if (negN == N | invN > 0 | downN > upN | deltaDistMetric < -0.05 | theta > thetaMax) // constant as proposed in paper, 
                {
                    laplacianVector = laplacianVector * 0.5;
                    newPoint = new Point3d(node.Coordinate.X + laplacianVector.X, node.Coordinate.Y + laplacianVector.Y, node.Coordinate.Z + laplacianVector.Z);
                    newNode = new qNode(newPoint, newNode.BoundaryNode);
                    if (i == 19) { newNode = node; }
                }
                else if (posN == N | (upN > 0 & downN == 0) | (upN >= downN & deltaDistMetric > -0.05)) { break; }
            }
            return newNode;
        }
        private Vector3d LaplacianSmooth(qNode node, List<qEdge> globalEdgeList)
        {
            Vector3d vectorSum = Vector3d.Zero;
            List<qEdge> connectedEdges = node.GetConnectedEdges(globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                Vector3d vector = Vector3d.Zero;
                vector = edge.GetOppositeNode(node).Coordinate - node.Coordinate;
                vectorSum = vectorSum + vector;
            }
            Vector3d laplacian = vectorSum / (double)connectedEdges.Count;
            return laplacian;
        } 
        private double CalculateGradient(qElement element, qNode node, double delta, string direction)
        {
            if (direction == "x")
            {
                Point3d point = new Point3d(node.Coordinate.X, node.Coordinate.Y, node.Coordinate.Z);
                List<qNode> elementNodes = element.GetNodesOfElement();
                point = new Point3d(point.X + delta, point.Y, point.Z);
                elementNodes[elementNodes.IndexOf(node)] = new qNode(point, node.BoundaryNode);
                qElement newElement = CreateElementFromNodes(elementNodes);
                double myPertubed = CalculateDistortionMetric(newElement); // pertubed distortion metric. 
                double gi = (myPertubed - element.DistortionMetric) / delta;
                return gi;
            }
            else if (direction == "y")
            {
                Point3d point = new Point3d(node.Coordinate.X, node.Coordinate.Y, node.Coordinate.Z);
                List<qNode> elementNodes = element.GetNodesOfElement();
                point = new Point3d(point.X, point.Y + delta, point.Z);
                elementNodes[elementNodes.IndexOf(node)] = new qNode(point, node.BoundaryNode);
                qElement newElement = CreateElementFromNodes(elementNodes);
                double myPertubed = CalculateDistortionMetric(newElement); // pertubed distortion metric. 
                double gi = (myPertubed - element.DistortionMetric) / delta;
                return gi;
            }
            else if (direction == "z")
            {
                Point3d point = new Point3d(node.Coordinate.X, node.Coordinate.Y, node.Coordinate.Z);
                List<qNode> elementNodes = element.GetNodesOfElement();
                point = new Point3d(point.X, point.Y, point.Z + delta);
                elementNodes[elementNodes.IndexOf(node)] = new qNode(point, node.BoundaryNode);
                qElement newElement = CreateElementFromNodes(elementNodes);
                double myPertubed = CalculateDistortionMetric(newElement); // pertubed distortion metric. 
                double gi = (myPertubed - element.DistortionMetric) / delta;
                return gi;
            }
            else { return 0; }
        }
        private qElement CreateElementFromNodes(List<qNode> nodes)
        {
            if (nodes.Count == 4)
            {
                qEdge edge1 = new qEdge(nodes[0], nodes[1]);
                qEdge edge2 = new qEdge(nodes[1], nodes[2]);
                qEdge edge3 = new qEdge(nodes[3], nodes[0]);
                qEdge edge4 = new qEdge(nodes[3], nodes[2]);
                List<qEdge> edgeList = new List<qEdge>() { edge1, edge2, edge3, edge4 };
                qElement element = new qElement(edgeList);
                return element;
            }
            else
            {
                qEdge edge1 = new qEdge(nodes[0], nodes[1]);
                qEdge edge2 = new qEdge(nodes[2], nodes[1]);
                qEdge edge3 = new qEdge(nodes[2], nodes[0]);
                List<qEdge> edgeList = new List<qEdge>() { edge1, edge2, edge3 };
                qElement element = new qElement(edgeList);
                return element;
            }
        }
        private BrepEdge PointOnEdge(qNode node, Brep brep)
        {
            BrepEdge edge = null;
            BrepEdgeList brepEdge = brep.Edges;

            foreach (BrepEdge bEdge in brepEdge) // check if node is on edge
            {
                bool IsOnEdge = node.IsOnEdge(bEdge);
                if (IsOnEdge)
                {
                    edge = bEdge;
                    return edge;
                }
            }
            return edge;
        }

        #endregion

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Icon_Triangle;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("a07a01a6-a771-4d75-9f96-e87ece274885"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;

namespace MeshPoints.Tools
{
    public class MergeSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MergeSmartMesh class.
        /// </summary>
        public MergeSmartMesh()
          : base("Merge SmartMesh", "merge",
              "Merge multiple SmartMesh to one SmartMesh.",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("List of SmartMesh", "smartMeshes", "List of SmartMeshes to merge.", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "smartMesh", "SmartMesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Grid Information", "grids", "Grid information is build as: grid information -> grid groups -> grids -> nodes.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Assumptions: merge only surface edges, edges to merge must share start/edge point, must have equal number of nodes. For solid, sweep the merged SmartMesh.
            // Assumption: U and V direction match

            // Input

            List<SmartMesh> smartMeshList = new List<SmartMesh>(); 
            DA.GetDataList(0,  smartMeshList);

            // Code
            SmartMesh mergedSmartMesh = new SmartMesh();
            List<List<List<Point3d>>> gridInformation = new List<List<List<Point3d>>>();

            SmartMesh mesh1 = smartMeshList[0];
            CreateGridInformation(mesh1, gridInformation, null, null);

            // Loop the meshes to be merged
            for (int i = 1; i < smartMeshList.Count; i++) 
            {
                SmartMesh mesh2 = smartMeshList[i];

                // Clone nodes and elements
                List<Node> newNodes = new List<Node>();
                List<Element> newElements = new List<Element>();
                var mesh1Prop = CloneNodesAndElements(mesh1);
                newNodes.AddRange(mesh1Prop.Item1);
                newElements.AddRange(mesh1Prop.Item2);
                var mesh2Prop = CloneNodesAndElements(mesh2);
                newNodes.AddRange(mesh2Prop.Item1);
                newElements.AddRange(mesh2Prop.Item2);


                // Find nodes to merge
                var edges = FindEdgeToMerge(mesh1, mesh2);
                BrepEdge edge1 = edges.Item1;
                BrepEdge edge2 = edges.Item2;

                var nodesToMerge = FindNodesToMerge(mesh1, mesh2, edge1, edge2);
                List<Node> nodesOnEdge1 = nodesToMerge.Item1;
                List<Node> nodesOnEdge2 = nodesToMerge.Item2;
                if (nodesOnEdge1[0] == null) { return; }

                // Update grid information
                CreateGridInformation(mesh2, gridInformation, nodesOnEdge1, nodesOnEdge2);

                // Merge nodes
                int numNodesToMerge = nodesOnEdge1.Count; 
                int startIdNewNodes = mesh1.Nodes.Count + mesh2.Nodes.Count - numNodesToMerge * 2; 

                for (int j = 0; j < numNodesToMerge; j++)
                {
                    // Find node from mesh2 that is equal to node from mesh1, and delete it. Update connected elements.
                    Node node1 = nodesOnEdge1[j];
                    Node node2 = FindClosestNode(node1, nodesOnEdge2);

                    MergeNodes(node1, node2, startIdNewNodes + j, newNodes, newElements);
                }
                    
                // Update global id of nodes
                List<List<int>> oldToNewId = new List<List<int>>();
                for (int j = 0; j < (newNodes.Count); j++) // loop nodes from mesh2
                {
                    List<int> list = new List<int>() { newNodes[j].GlobalId, j }; // fix denne
                    oldToNewId.Add(list);
                    newNodes[j].GlobalId = j;
                }

                // Update elements and global mesh
                int counter1 = 0;
                int counter2 = 0;

                for (int j = 0; j < newElements.Count; j++) // loop elements from mesh1
                {
                    for (int k = 0; k < 4; k++)
                    {
                        // loop through - old to new id of nodes to get the updated id
                        for (int n = 0; n < mesh1.Nodes.Count + mesh2.Nodes.Count - numNodesToMerge * 2; n++)
                        {
                            bool perform = false;
                            List<int> oldToNew = oldToNewId[n];
                            if (j < mesh1.Elements.Count && n < mesh1.Nodes.Count - numNodesToMerge)
                            {
                                if (oldToNew[0] == mesh1.Elements[counter1].Connectivity[k]) {perform = true;}
                            }
                            else if (j >= mesh1.Elements.Count && n >= mesh1.Nodes.Count - numNodesToMerge)
                            {
                                if (oldToNew[0] == mesh2.Elements[counter2].Connectivity[k]) {perform = true;}
                            }

                            if (perform)
                            {
                                newElements[j].Nodes[k] = newNodes[oldToNew[1]];
                                newElements[j].Connectivity[k] = oldToNew[1];
                                newElements[j].Id = j;
                            }
                        }
                    }

                    if (j < mesh1.Elements.Count) { counter1++; }
                    else { counter2++; }
                }

                // Create merged SmartMesh
                mergedSmartMesh = new SmartMesh(newNodes, newElements, "Surface"); 
                UpdateGeometry(mesh1, mesh2, mergedSmartMesh);
                mergedSmartMesh.GridInformation = gridInformation;

                // Repeat
                mesh1 = mergedSmartMesh;
            }

            // Fix gridInformation
            List<List<List<Node>>> gridInfoNodes = new List<List<List<Node>>>();
            for (int i = 0; i < gridInformation.Count; i++)
            {
                List<List<Node>> gridGroupsNodes = new List<List<Node>>();
                for (int j = 0; j < gridInformation[i].Count; j++)
                {
                    List<Node> gridsNodes = new List<Node>();
                    for (int k = 0; k < gridInformation[i][j].Count; k++)
                    {
                        foreach (Node node in mergedSmartMesh.Nodes)
                        {
                            if ((gridInformation[i][j][k]- node.Coordinate).Length < 0.001)
                            {
                                gridsNodes.Add(node);
                            }
                        }
                    }
                    gridGroupsNodes.Add(gridsNodes);
                }
                gridInfoNodes.Add(gridGroupsNodes);
            }

            DA.SetData(0, mergedSmartMesh);
            DA.SetDataList(1, gridInfoNodes); // to do: slett denne
        }


        // Methods
        private Tuple<List<Node>, List<Element>> CloneNodesAndElements(SmartMesh smartMesh)
        {
            // Create new node list
            List<Node> newNodes = new List<Node>();
            List<Node> oldNodes = new List<Node>();
            oldNodes.AddRange(smartMesh.Nodes);
            foreach (Node node in oldNodes)
            {
                Node newNode = new Node(node.GlobalId, node.Coordinate, node.BC_U, node.BC_V);
                newNode.Type = node.Type;
                newNodes.Add(newNode);
            }

            // Create new element list
            List<Element> newElements = new List<Element>();
            List<Element> oldElements = new List<Element>();
            oldElements.AddRange(smartMesh.Elements);
            foreach (Element element in oldElements)
            {
                List<Node> elementNodes = new List<Node>();
                List<int> elementConnectivity = new List<int>();
                int id = element.Id;
                for (int j = 0; j < 4; j++)
                {
                    Node node = element.Nodes[j];
                    Node newNode = new Node(node.GlobalId, node.Coordinate, node.BC_U, node.BC_V);
                    newNode.Type = node.Type;
                    elementNodes.Add(newNode);
                    elementConnectivity.Add(element.Connectivity[j]);
                }

                Element newElement = new Element(id, elementNodes, elementConnectivity);
                newElement.GetElementMesh();
                newElements.Add(newElement);
            }
            return Tuple.Create(newNodes, newElements);
        }
        private void MergeNodes(Node node1, Node node2, int id, List<Node> newNodes, List<Element> newElements)
        {
            Point3d newNodeCoordinate = (node1.Coordinate + node2.Coordinate) / (double)2;
            Node newNode = new Node(id, newNodeCoordinate);

            if (node1.Type == "Corner") // keep BC if corner node
            {
                newNode.BC_U = true;
                newNode.BC_V = true;
            }
            else
            {
                newNode.Type = "Merged";
            }

            // Find elements to change and update the nodes and connectivity
            foreach (Element element in newElements)
            {
                for (int k = 0; k < 4; k++)
                {
                    if (element.Nodes[k].Coordinate == node1.Coordinate)
                    {
                        element.Nodes[k] = newNode;
                        element.Connectivity[k] = newNode.GlobalId;
                        foreach (Node node in newNodes)
                        {
                            if (node.Coordinate == node1.Coordinate)
                            {
                                newNodes.Remove(node);
                                break;
                            }
                        }
                    }
                    else if (element.Nodes[k].Coordinate == node2.Coordinate)
                    {
                        element.Nodes[k] = newNode;
                        element.Connectivity[k] = newNode.GlobalId;
                        foreach (Node node in newNodes)
                        {
                            if (node.Coordinate == node2.Coordinate)
                            {
                                newNodes.Remove(node);
                                break;
                            }
                        }
                    }
                }
                UpdateElementMesh(element); // to do: check if needed
            }
            newNodes.Add(newNode);

        }
        private Node FindClosestNode(Node nodeToCheck, List<Node> nodes)
        {
            double minDistance = 1000000;
            Node closestNode = null;
            foreach (Node node in nodes)
            {
                double distance = nodeToCheck.Coordinate.DistanceTo(node.Coordinate);
                if (distance < minDistance)
                {
                    minDistance = distance;
                    closestNode = node;
                }
            }
            return closestNode;
        }
        private Tuple<BrepEdge, BrepEdge> FindEdgeToMerge(SmartMesh mesh1, SmartMesh mesh2)
        {
            BrepEdge edge1 = null;
            BrepEdge edge2 = null;

            // Find edges to merge
            for (int j = 0; j < mesh1.Geometry.Edges.Count; j++)
            {
                BrepEdge edgeToTest1 = mesh1.Geometry.Edges[j];
                for (int k = 0; k < mesh2.Geometry.Edges.Count; k++)
                {
                    BrepEdge edgeToTest2 = mesh2.Geometry.Edges[k];
                    List<Point3d> edgePoints1 = new List<Point3d>() { edgeToTest1.PointAtStart, edgeToTest1.PointAtEnd };
                    List<Point3d> edgePoints2 = new List<Point3d>() { edgeToTest2.PointAtStart, edgeToTest2.PointAtEnd };

                    if (edgePoints1.Contains(edgePoints2[0]) & edgePoints1.Contains(edgePoints2[1])) // if edges share start and end points
                    {
                        edge1 = edgeToTest1;
                        edge2 = edgeToTest2;
                        j = mesh1.Geometry.Edges.Count; // break
                        k = mesh2.Geometry.Edges.Count; // break
                    }
                }
            }
            return Tuple.Create(edge1, edge2);
        }
        private Tuple<List<Node>, List<Node>> FindNodesToMerge(SmartMesh mesh1, SmartMesh mesh2, BrepEdge edge1, BrepEdge edge2)
        {
            List<Node> nodesOnEdge1 = new List<Node>();
            List<Node> nodesOnEdge2 = new List<Node>();

            // Find nodes to merge
            foreach (Node node1 in mesh1.Nodes)
            {
                if (node1.IsOnEdge(edge1))
                {
                    nodesOnEdge1.Add(node1);
                }
            }

            foreach (Node node2 in mesh2.Nodes)
            {
                if (node2.IsOnEdge(edge2))
                {
                    nodesOnEdge2.Add(node2);
                }
            }

            if (nodesOnEdge1.Count != nodesOnEdge2.Count)
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Number of nodes on edges to merge do not match."); return Tuple.Create(nodesOnEdge1, nodesOnEdge2); }

            return Tuple.Create(nodesOnEdge1, nodesOnEdge2);
        }
        private void UpdateElementMesh(Element element)
        {
            Mesh mesh = new Mesh();
            foreach (Node node in element.Nodes)
            {
                mesh.Vertices.Add(node.Coordinate);
            };

            mesh.Faces.AddFace(0, 1, 2, 3);
            mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh

            element.Mesh = mesh; // replace old mesh
        }
        private void CreateGridInformation(SmartMesh mesh, List<List<List<Point3d>>> gridInformation, List<Node> nodesOnEdge1, List<Node> nodesOnEdge2)
        {
            // Create grid grops, each node, grid[i], can move relative to neighbor grid nodes, grid[i+-1]. For nodes grid[0] and grid[endIndex] no translation.

            // Assmume all mesh input are structured
            List<List<Point3d>> gridGroupU = new List<List<Point3d>>();
            List<List<Point3d>> gridGroupV = new List<List<Point3d>>();

            // Create grid info for mesh
            if (mesh.GridInformation == null)
            {
                // Create grids
                for (int i = 0; i < mesh.nu; i++) // u grid
                {
                    List<Point3d> gridTemp = new List<Point3d>();
                    for (int j = 0; j < mesh.nv; j++)
                    {
                        gridTemp.Add(mesh.Nodes[i + j * mesh.nu].Coordinate);
                    }
                    gridGroupU.Add(gridTemp);
                }

                for (int j = 0; j < mesh.nv; j++) // v grid
                {
                    List<Point3d> gridTemp = new List<Point3d>();
                    for (int i = 0; i < mesh.nu; i++)
                    {
                        gridTemp.Add(mesh.Nodes[i + j * mesh.nu].Coordinate);
                    }
                    gridGroupV.Add(gridTemp);
                }

                mesh.GridInformation = new List<List<List<Point3d>>>() { gridGroupU, gridGroupV };
            }


            // Update gridInformation
            if (gridInformation.Count == 0)
            {
                for (int i = 0; i < mesh.GridInformation.Count; i++)
                {
                    List<List<Point3d>> newGroup = new List<List<Point3d>>();
                    for (int j = 0; j < mesh.GridInformation[i].Count; j++)
                    {
                        List<Point3d> newGrid = new List<Point3d>();
                        for (int k = 0; k < mesh.GridInformation[i][j].Count; k++)
                        {
                            newGrid.Add(mesh.GridInformation[i][j][k]);
                        }
                        newGroup.Add(newGrid);
                    }
                    gridInformation.Add(newGroup);
                }
            }
            else
            {
                // Merge grids

                List<List<Point3d>> newGridGroup = new List<List<Point3d>>();
                List<List<Point3d>> oldGridGroup = new List<List<Point3d>>();

                // Get new grid group to merge
                for (int p = 0; p < mesh.GridInformation.Count; p++) // loop new gridinformation
                {
                    for (int n = 0; n < mesh.GridInformation[p].Count; n++) // loop grids 
                    {
                        List<Point3d> grid = mesh.GridInformation[p][n];

                        // Check if new grid is edge to merge
                        int edgeCheck = 0;
                        foreach (Node mergeNode in nodesOnEdge2)
                        {
                            foreach (Point3d gridNode in grid)
                            {
                                if ((gridNode - mergeNode.Coordinate).Length < 0.001)
                                {
                                    edgeCheck++;
                                }
                            }
                        }

                        if (edgeCheck != 1)
                        {
                            gridInformation.Add(mesh.GridInformation[p]);
                            break;
                        }
                        else
                        {
                            newGridGroup = mesh.GridInformation[p];
                            break;
                        }
                    }
                }


                // Get existing grid group to merge
                int numGridGroups = gridInformation.Count;
                for (int k = 0; k < numGridGroups; k++)
                {
                    int gridCheck = 0;
                    oldGridGroup = gridInformation[k];
                    for (int j = 0; j < oldGridGroup.Count; j++)
                    {
                        List<Point3d> oldGrid = oldGridGroup[j];

                        // Check if grid is edge to merge
                        int edgeCheck = 0;
                        foreach (Node mergeNode in nodesOnEdge1)
                        {
                            foreach (Point3d gridNode in oldGrid)
                            {
                                if ((gridNode - mergeNode.Coordinate).Length < 0.001)
                                {
                                    edgeCheck++;
                                }
                            }
                        }

                        if (edgeCheck != 1)
                        {
                            continue;
                        }
                        gridCheck++;
                    }

                    if (gridCheck == newGridGroup.Count)
                    {
                        break;
                    }
                }

                foreach (List<Point3d> newGrid in newGridGroup)
                {
                    for (int counterOldGrids = 0; counterOldGrids < oldGridGroup.Count; counterOldGrids++)
                    {
                        int oldGridCount = oldGridGroup[counterOldGrids].Count;
                        for (int counterOldNodes = 0; counterOldNodes < oldGridCount; counterOldNodes++)
                        {
                            Point3d oldPoint = oldGridGroup[counterOldGrids][counterOldNodes];
                            if ((oldPoint - newGrid[0]).Length < 0.001 | (oldPoint - newGrid.Last()).Length < 0.001)
                            {
                                foreach (Point3d pointToAdd in newGrid)
                                {
                                    if ((pointToAdd - oldPoint).Length < 0.001) { continue; }
                                    oldGridGroup[counterOldGrids].Add(pointToAdd);
                                }
                                oldGridCount = oldGridGroup.Count;
                                counterOldNodes = oldGridCount;
                            }
                        }

                    }
                }
            }
        }
        private void UpdateGeometry(SmartMesh mesh1, SmartMesh mesh2, SmartMesh mergedSmartMesh)
        {
            List<BrepEdge> brepEdges = new List<BrepEdge>();
            var edgesToRemove = FindEdgeToMerge(mesh1, mesh2);
            brepEdges.AddRange(mesh1.Geometry.Edges);
            brepEdges.AddRange(mesh2.Geometry.Edges);
            brepEdges.Remove(edgesToRemove.Item1);
            brepEdges.Remove(edgesToRemove.Item2);

            List<BrepFace> faces = new List<BrepFace>();
            faces.AddRange(mesh1.Geometry.Faces);
            faces.AddRange(mesh2.Geometry.Faces);

            Brep totalBrep = null;
            totalBrep = mesh1.Geometry.Brep;
            totalBrep.Join(mesh2.Geometry.Brep, 0.001, true);

            Geometry mergedGeometry = new Geometry(totalBrep, faces, brepEdges, totalBrep.Vertices.ToList());
            mergedSmartMesh.Geometry = mergedGeometry;
        }
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Icon_MergeMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("4eede3a5-ff72-42e4-bac3-5c35cbeff345"); }
        }
    }
}
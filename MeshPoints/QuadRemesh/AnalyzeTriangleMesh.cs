using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using Rhino.Geometry.Collections;

namespace MeshPoints.QuadRemesh
{
    public class AnalyzeTriangleMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the AnalyzeTriangleMesh class.
        /// </summary>
        public AnalyzeTriangleMesh()
          : base("Analyze Triangle Mesh", "atm",
              "Analyse elements and edges in a triangle mesh",
              "MyPlugIn", "QuadRemesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Triangle mesh", "trimesh", "Input a trinagle mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("# elements to remesh", "number", "Input an integer", GH_ParamAccess.item, 1);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Front edges", "element", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("List 1-1", "11", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("List 1-0", "10", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("List 0-1", "01", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("List 0-0", "00", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("qElements", "element", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("E_front", "element", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("E_k_left", "element", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("E_k_right ", "element", "", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // variables
            Mesh mesh = new Mesh();
            double numberElementsToRemesh = 0;
            List<qElement> elementList = new List<qElement>();
            List<qEdge> edgeList = new List<qEdge>();
            qElement element = new qElement();
            List<qEdge> frontEdges = new List<qEdge>();
            List<double> angleList = new List<double>();

            // input
            DA.GetData(0, ref mesh);
            DA.GetData(1, ref numberElementsToRemesh);

            // NOTE: Topological vertices and vertices are not sorted the same way. Make use to topological vertices in these codes.
            edgeList = CreateEdges(mesh);
            elementList = CreateElements(mesh, edgeList);
            SetNeighborElements(mesh, elementList, edgeList);
            frontEdges = GetFrontEdges(mesh, edgeList);
            SetNeighorFrontEdges(frontEdges);
            var edgeStates = CreateEdgeStateList(mesh, frontEdges); // todo: check bug. Angles.
            var list11 = edgeStates.Item1;
            var list10 = edgeStates.Item2;
            var list01 = edgeStates.Item3;
            var list00 = edgeStates.Item4;

            qEdge E_k_left = new qEdge(); // left side edge of quad
            qEdge E_k_right = new qEdge(); // right side edge of quad
            qEdge E_front = new qEdge(); // bottom of quad


            List<qElement> sortedConnectedElements = new List<qElement>();
            List<qEdge> sortedConnectedEdges = new List<qEdge>();

            qEdge E_i = new qEdge();
            qEdge E_k = new qEdge();
            qEdge E_k1 = new qEdge();
            List<qEdge> intersectingS = new List<qEdge>();
            List<qElement> connectedElements = new List<qElement>();

            for (int n = 0; n < numberElementsToRemesh; n++)
            {
                // get E_front
                // todo: extand the hieraki
                int edgeState = 0;
                if (list11.Count != 0) { E_front = list11[0]; edgeState = 11; list11.RemoveAt(0); }
                else if (list01.Count != 0) { E_front = list01[0]; edgeState = 01; list01.RemoveAt(0); }
                else if (list10.Count != 0) { E_front = list10[0]; edgeState = 10; list10.RemoveAt(0); }
                else { E_front = list00[0]; edgeState = 00; list00.RemoveAt(0); }

                // get E_k
                if (edgeState == 11)
                {
                    E_k_left = E_front.LeftFrontNeighbor;
                    E_k_right = E_front.RightFrontNeighbor;
                }
                else if (edgeState == 01 | edgeState == 00) 
                {
                    int nodeToEvaluate = 0; // evaluate left node
                    E_k_left = GetSideEdge(mesh, elementList, edgeList, nodeToEvaluate, E_front);
                    E_k_right = E_front.RightFrontNeighbor;
                }
                else
                {
                    int nodeToEvaluate = 1; // evaluate right node
                    E_k_right = GetSideEdge(mesh, elementList, edgeList, nodeToEvaluate, E_front);
                    E_k_left = E_front.LeftFrontNeighbor;
                }

                // todo: update edgestate list after top recovery


                #region top recovery:
                


                // todo: do not work with the swap and split functions because not updated edges and surrounding objects
                qNode leftNode = GetNodeOfFrontEdge(0, E_front);
                qNode rightNode = GetNodeOfFrontEdge(1, E_front);
                qNode N_d = new qNode();
                qNode N_c = new qNode();


                // Get N_d and N_c
                if (E_k_left.StartNode == leftNode) { N_d = E_k_left.EndNode; }
                else { N_d = E_k_left.StartNode; }

                if (E_k_right.StartNode == rightNode) { N_c = E_k_right.EndNode; }
                else { N_c = E_k_right.StartNode; }

                Vector3d V_s = N_d.Coordinate - N_c.Coordinate;
                

                // Get connected elements
                connectedElements = new List<qElement>();
                List<qEdge> connectedEdges = GetConnectedEdges(N_c, edgeList);
                foreach (qEdge connectedEdge in connectedEdges)
                {
                    if (!connectedElements.Contains(connectedEdge.Element1)) 
                    {
                        connectedElements.Add(connectedEdge.Element1);
                    }
                    if (connectedEdge.Element2 != null)
                    {
                        if (!connectedElements.Contains(connectedEdge.Element2))
                        {
                            connectedElements.Add(connectedEdge.Element2);
                        }
                    }
                }
                
                #region Sort connected elements ccw
                // create vectors from N_c to center of elements
                List<Vector3d> vectorToCenterList = new List<Vector3d>();
                foreach (qElement connectedElement in connectedElements)
                {
                    Point3d elementCenter = GetElementCenter(connectedElement); 
                    vectorToCenterList.Add(elementCenter - N_c.Coordinate);      
                }
                
                // calculate angle between a vector and the first center-vector
                Vector3d intialVector = vectorToCenterList[0];
                List<double> vecAngleList = new List<double>() { 0 };
                for (int i = 1; i < vectorToCenterList.Count; i++)
                {
                    double angle = Vector3d.VectorAngle(intialVector, vectorToCenterList[i], Vector3d.ZAxis); // todo: change to normal of intialface
                    vecAngleList.Add(angle);
                }

                // sort connected elements ccw from the first element in connected element list
                sortedConnectedElements = new List<qElement>();
                for (int j = 0; j < vectorToCenterList.Count; j++)
                {
                    double minAngle = vecAngleList[0];
                    int minAngleIndex = 0;
                    for (int i = 1; i < vecAngleList.Count; ++i)
                    {
                        if (vecAngleList[i] < minAngle)
                        {
                            minAngle = vecAngleList[i];
                            minAngleIndex = i;
                        }
                    }
                    sortedConnectedElements.Add(connectedElements[minAngleIndex]);
                    connectedElements.RemoveAt(minAngleIndex);
                    vecAngleList.RemoveAt(minAngleIndex);
                }
                #endregion


                #region Find intersected edges
                qElement T_k = new qElement();
                qElement T_k1 = new qElement();

                foreach (qElement connectedElement in sortedConnectedElements)
                {
                    T_k = connectedElement;
                    if (connectedElement.IsQuad) 
                    {
                        // todo: Abort the front edge, try with a new front edge. 
                        { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Top edge recovery fails because of intersection of quads."); }
                    }

                    // get edges
                    E_i = new qEdge();
                    E_k = new qEdge();
                    E_k1 = new qEdge();
                    Vector3d V_k = Vector3d.Zero;
                    Vector3d V_k1 = Vector3d.Zero;

                    List<qEdge> E_kCandidates = new List<qEdge>();
                    foreach (qEdge elementEdge in connectedElement.EdgeList)
                    {
                        if (elementEdge.StartNode == N_c | elementEdge.EndNode == N_c) { E_kCandidates.Add(elementEdge); }
                        else { E_i = elementEdge; }
                    }

                    var vectors = CalculateVectorsFromSharedNode(E_kCandidates[0], E_kCandidates[1]);
                    Vector3d vec1 = vectors.Item1;
                    Vector3d vec2 = vectors.Item2;
                    Vector3d cros = Vector3d.CrossProduct(vec1, vec2);
                    if (cros.Z > 0) // todo: change cros.Z > 0 to make more general
                    { 
                        E_k = E_kCandidates[0]; 
                        E_k1 = E_kCandidates[1];
                        V_k = vec1;
                        V_k1 = vec2;
                    } 
                    else
                    {
                        E_k = E_kCandidates[1];
                        E_k1 = E_kCandidates[0];
                        V_k = vec2;
                        V_k1 = vec1;
                    }

                    // check if intersection
                    // todo: change if 3d
                    if (Vector3d.Multiply(V_s, V_k) > 0 & Vector3d.Multiply(V_s, V_k1) > 0) // todo: check if ok with different last criterior than kok.
                    {
                        break;
                    }
                    else if (Math.Abs(V_s.IsParallelTo(V_k)) == 1 | Math.Abs(V_s.IsParallelTo(V_k1)) == 1)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "E_k or E_k1 in top edge recovery failed because V_s || V_k or V_k1.");
                        break;
                    }
                    
                }


                if (!frontEdges.Contains(E_i))
                {
                    intersectingS.Add(E_i);
                }
                else 
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Top edge recovery failed.");
                    break;
                }
                
                bool done = false;
                while (!done)
                {
                    if (E_i.Element1 != T_k)
                    {
                        T_k1 = E_i.Element1;
                    }
                    else if (E_i.Element2 != T_k)
                    {
                        T_k1 = E_i.Element2;
                    }

                    List<qNode> elementNodes = GetNodesOfElement(T_k1);
                    if (elementNodes.Contains(N_d))
                    {
                        done = true;
                        continue;
                    }

                    T_k = T_k1;

                    qNode N_i = new qNode();
                    foreach (qNode node in elementNodes)
                    {
                        if ((node != E_i.StartNode) & (node != E_i.EndNode))
                        {
                            N_i = node;
                        }
                    }

                    Vector3d V_i = N_i.Coordinate - N_c.Coordinate;

                    List<qEdge> E_nCandidates = new List<qEdge>();
                    foreach (qEdge edge in T_k.EdgeList)
                    {
                        if (edge.StartNode == N_i | edge.EndNode == N_i) { E_nCandidates.Add(edge); }
                    }

                    var vectors = CalculateVectorsFromSharedNode(E_nCandidates[0], E_nCandidates[1]);
                    Vector3d vec1 = vectors.Item1;
                    Vector3d vec2 = vectors.Item2;
                    qEdge E_n = new qEdge();
                    qEdge E_n1 = new qEdge();

                    Vector3d cros = Vector3d.CrossProduct(vec1, vec2);
                    if (cros.Z > 0) // todo: change cros.Z > 0 to make more general
                    {
                        E_n = E_nCandidates[1];
                        E_n1 = E_nCandidates[0];
                    }
                    else
                    {
                        E_n = E_nCandidates[0];
                        E_n1 = E_nCandidates[1];
                    }

                    if (Vector3d.Multiply(V_s, V_i) < 0) // todo: change if 3d
                    {
                        E_i = E_n;
                    }
                    else 
                    {
                        E_i = E_n1;
                    }
                    if (!frontEdges.Contains(E_i)) { intersectingS.Add(E_i); }
                }
                #endregion

                // edge recovery process
           
                while (intersectingS.Count > 0)
                {
                    qEdge edge = intersectingS[0];
                    SwapEdge(edge); // swap the edge
                   
                    //List<qElement> Tswaped = new List<qElement>() { edge.Element1, edge.Element2 }; // todo: need to update mesh first

                    // todo: check area of adjacent elements
                    bool areaOK = true;
                    if (areaOK)
                    {
                        intersectingS.RemoveAt(0);
                    }
                    else
                    {
                        intersectingS.Add(edge);
                        intersectingS.RemoveAt(0);
                    }
                }
                #endregion

            }

            // output
            DA.SetDataList(0, frontEdges);
            DA.SetDataList(1, list11);
            DA.SetDataList(2, list10);
            DA.SetDataList(3, list01);
            DA.SetDataList(4, list00);
            DA.SetDataList(5, intersectingS);
            DA.SetData(6, E_front); 
            DA.SetData(7, E_k_right);
            DA.SetData(8, E_i);
        }
        #region Methods
        private List<qNode> CreateNodes(Mesh mesh)
        {
            // summary: create nodes from initial mesh
            List<qNode> nodeList = new List<qNode>();
            MeshTopologyVertexList topologyVertexList = mesh.TopologyVertices;
            for (int i = 0; i < topologyVertexList.Count; i++)
            {
                qNode node = new qNode();
                node.TopologyVertexIndex = i; // todo: delete
                node.MeshVertexIndex = topologyVertexList.MeshVertexIndices(i)[0]; // todo: delete
                node.Coordinate = mesh.Vertices.Point3dAt(node.MeshVertexIndex);
                nodeList.Add(node);
            }
            return nodeList;
        }
        private List<qEdge> GetConnectedEdges(qNode node, List<qEdge> edgeList)
        {
            List<qEdge> connectedEdges = new List<qEdge>();

            foreach (qEdge edge in edgeList)
            {
                if ((edge.StartNode == node) | (edge.EndNode == node))
                {
                    connectedEdges.Add(edge);
                }
            }
            return connectedEdges;
        }
        private List<qElement> GetConnectedElements(qEdge edge)
        {
            List<qElement> connectedElements = new List<qElement>();
            connectedElements.Add(edge.Element1);
            if (edge.Element2 != null)
            {
                connectedElements.Add(edge.Element2);
            }
            return connectedElements;
        }
        private List<qEdge> CreateEdges(Mesh mesh)
        {
            // summary: create edges from intial mesh

            List<qEdge> edgeList = new List<qEdge>();
            MeshTopologyEdgeList topologyEdgeList = mesh.TopologyEdges;
            List<qNode> nodeList = CreateNodes(mesh);

            for (int i = 0; i < topologyEdgeList.Count; i++)
            {
                Rhino.IndexPair edgeTopoVerticesIndex = topologyEdgeList.GetTopologyVertices(i);
                qNode startNode = nodeList[edgeTopoVerticesIndex[0]];
                qNode endNode = nodeList[edgeTopoVerticesIndex[1]];

                qEdge edge = new qEdge(startNode, endNode); //sjekk om i er nødvendig
                edgeList.Add(edge);
            }
            return edgeList;
        }
        private List<qElement> CreateElements(Mesh mesh, List<qEdge> edgeList)
        {
            qElement element = new qElement();
            List<qElement> elementList = new List<qElement>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                List<qEdge> edgeListElement = new List<qEdge>();

                int[] elementEdgesIndex = mesh.TopologyEdges.GetEdgesForFace(i);
                foreach (int n in elementEdgesIndex)
                {
                    edgeListElement.Add(edgeList[n]);
                }
                element = new qElement(i, edgeListElement);
                elementList.Add(element);
            }
            return elementList;
        }
        private void SetNeighborElements(Mesh mesh, List<qElement> elementList, List<qEdge> edgeList)
        {
            for (int i = 0; i < mesh.TopologyEdges.Count; i++)
            {
                int[] connectedElementsIndex = mesh.TopologyEdges.GetConnectedFaces(i);
                edgeList[i].Element1 = elementList[connectedElementsIndex[0]];

                if (connectedElementsIndex.Length == 2)
                {
                    edgeList[i].Element2 = elementList[connectedElementsIndex[1]];
                }
            }
            return;
        }
        private List<qEdge> GetFrontEdges(Mesh mesh, List<qEdge> edgeList)
        {
            List<qEdge> frontEdges = new List<qEdge>();
            for (int i = 0; i < edgeList.Count; i++)
            {
                int[] connectedElementsIndex = mesh.TopologyEdges.GetConnectedFaces(i);
                if ((connectedElementsIndex.Length == 1) & !edgeList[i].Element1.IsQuad)
                {
                    frontEdges.Add(edgeList[i]);
                }
            }
            return frontEdges;
        }
        private void SetNeighorFrontEdges(List<qEdge> frontEdges)
        {
            for (int i = 0; i < frontEdges.Count; i++)
            {
                qEdge edge = frontEdges[i];
                qEdge edg1 = new qEdge();
                qEdge edg2 = new qEdge();
                for (int j = 0; j < frontEdges.Count; j++)
                {
                    qEdge testEdge = frontEdges[j];
                    if (i != j)
                    {
                        if ((edge.StartNode == testEdge.StartNode) | (edge.StartNode == testEdge.EndNode)) // check edge start node if connected
                        {
                            edg1 = frontEdges[j];
                        }
                        else if ((edge.EndNode == testEdge.StartNode) | (edge.EndNode == testEdge.EndNode)) // check edge end point if connected
                        {
                            edg2 = frontEdges[j];
                        }
                    }
                }

                Point3d midPointEdg = (edge.StartNode.Coordinate + edge.EndNode.Coordinate) / 2.0; // mid point of edge

                Point3d center1 = GetElementCenter(edg1.Element1); // center of element to neighbor edge 1

                Point3d center2 = GetElementCenter(edg2.Element1); // center of elemenet to neighbor edge 2

                Vector3d vec1 = center1 - midPointEdg;
                Vector3d vec2 = center2 - midPointEdg;

                // check what is left and right neighbor for edge
                if (Vector3d.CrossProduct(vec1, vec2).Z > 0) // todo: modify to more general relative to face normals
                {
                    edge.RightFrontNeighbor = edg1;
                    edge.LeftFrontNeighbor = edg2;
                }
                else
                {
                    edge.RightFrontNeighbor = edg2;
                    edge.LeftFrontNeighbor = edg1;
                }
            }
            return;
        }
        Tuple<List<qEdge>, List<qEdge>, List<qEdge>, List<qEdge>> CreateEdgeStateList(Mesh mesh, List<qEdge> frontEdges)
        {
            // Return the EdgeStates in the format of the list: list11, list01, list10, andlist00.
            List<qEdge> list11 = new List<qEdge>();
            List<qEdge> list01 = new List<qEdge>();
            List<qEdge> list10 = new List<qEdge>();
            List<qEdge> list00 = new List<qEdge>();

            double angleTolerance = 0.75 * Math.PI;
            double leftAngle = 0;
            double rightAngle = 0;
            int edgeState = 0;
            
            for (int i = 0; i < frontEdges.Count; i++)
            {
                for (int nodeToCalculate = 0; nodeToCalculate < 2; nodeToCalculate++) // nodeToCalculate = 0 is left neighbor, otherwise right neighbor
                {
                    if (nodeToCalculate == 0)
                    {
                        leftAngle = CalculateAngleOfAdjecentFrontEdges(nodeToCalculate, frontEdges[i]);
                    }
                    else
                    {
                        rightAngle = CalculateAngleOfAdjecentFrontEdges(nodeToCalculate, frontEdges[i]);
                    }

                    if (leftAngle < angleTolerance & rightAngle < angleTolerance) { edgeState = 11; }
                    else if (leftAngle >= angleTolerance & rightAngle < angleTolerance) { edgeState = 01; }
                    else if (leftAngle < angleTolerance & rightAngle >= angleTolerance) { edgeState = 10; }
                    else if (leftAngle >= angleTolerance & rightAngle >= angleTolerance) { edgeState = 00; }
                }

                switch (edgeState)
                {
                    case 11:
                        list11.Add(frontEdges[i]);
                        break;
                    case 10:
                        list10.Add(frontEdges[i]);
                        break;
                    case 01:
                        list01.Add(frontEdges[i]);
                        break;
                    case 00:
                        list00.Add(frontEdges[i]);
                        break;
                }
            }
            return Tuple.Create(list11, list10, list01, list00);
        }
        private double CalculateAngleOfAdjecentFrontEdges(int nodeToCalculate, qEdge edge)
        {
            // Calculate the angle between front edges with a shared point, i.e. neighbor front edges.

            List<qElement> connectedElements = GetConnectedElements(edge);
            if ( connectedElements.Count > 1 | edge.Element1.IsQuad )
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The edge is not a front edge.");}

            qEdge edgeNeighbor = new qEdge();
            double angle = 0;

            // Get neighbor nodes
            if (nodeToCalculate == 0) { edgeNeighbor = edge.LeftFrontNeighbor; }
            else { edgeNeighbor = edge.RightFrontNeighbor; }

            // Create vectors from shared node
            var vectors = CalculateVectorsFromSharedNode(edge, edgeNeighbor);
            Vector3d vec1 = vectors.Item1;
            Vector3d vec2 = vectors.Item2;
            qNode sharedNode= vectors.Item3;
            Vector3d V_k = vec1.Length * vec2 + vec2.Length * vec1; // angle bisector
            if (V_k == Vector3d.Zero) // create V_k if zero vector
            {
                if (nodeToCalculate == 0) { V_k = vec1; }
                else { V_k = vec2; }

                Vector3d rotationAxis = Vector3d.CrossProduct(vec1, vec2);
                V_k.Rotate(0.5 * Math.PI, rotationAxis);
            }
            V_k = V_k / V_k.Length; // normalize

            // check with domain
            Point3d endPointV_k = Point3d.Add(sharedNode.Coordinate, V_k); // endpoint of V_k from sharedNode
            Point3d endPointV_kNegative = Point3d.Add(sharedNode.Coordinate, -V_k); // endpoint of -V_k from sharedNode

            // distance to domain
            int edgeElementFaceId = edge.Element1.FaceIndex;
            int edgeNeighborElementFaceId = edgeNeighbor.Element1.FaceIndex;
            Point3d edgeFaceCenter = GetElementCenter(edge.Element1);
            Point3d edgeNeighborFaceCenter = GetElementCenter(edgeNeighbor.Element1);
            Point3d midFaceCenter = (edgeFaceCenter + edgeNeighborFaceCenter) / 2; // mid face center

            double distanceEndPoint = midFaceCenter.DistanceTo(endPointV_k); // todo: ok for all cases?
            double distanceEndPointNegative = midFaceCenter.DistanceTo(endPointV_kNegative); // todo: ok for all cases?
            double alpha = Vector3d.VectorAngle(vec1, vec2);

            if (distanceEndPoint < distanceEndPointNegative)
            {
                // V_k is inside domain
                angle = alpha;
            }
            else 
            {
                // V_k is outside domain
                angle = 2 * Math.PI - alpha;
            }
            return angle;
        }
        private qNode GetNodeOfFrontEdge(int nodeToEvaluate, qEdge edge)
        {
            // get a node of a front edge, 0=left, 1=right
            qNode sharedNode = new qNode();
            qEdge neighborFrontEdge = new qEdge();
            if (nodeToEvaluate == 0) { neighborFrontEdge = edge.LeftFrontNeighbor; }
            else { neighborFrontEdge = edge.RightFrontNeighbor; }

            if (edge.StartNode == neighborFrontEdge.StartNode) { sharedNode = edge.StartNode; }
            else if (edge.StartNode == neighborFrontEdge.EndNode) { sharedNode = edge.StartNode; }
            else if (edge.EndNode == neighborFrontEdge.StartNode) { sharedNode = edge.EndNode; }
            else if (edge.EndNode == neighborFrontEdge.EndNode) { sharedNode = edge.EndNode; }
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No shared nodes"); }
            return sharedNode;
        }
        private Tuple<Vector3d, Vector3d, qNode> CalculateVectorsFromSharedNode(qEdge edge1, qEdge edge2)
        {
            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;
            qNode sharedNode = new qNode();
            if (edge1.StartNode == edge2.StartNode)
            {
                sharedNode = edge1.StartNode;
                vec1 = edge1.EndNode.Coordinate - edge1.StartNode.Coordinate;
                vec2 = edge2.EndNode.Coordinate - edge2.StartNode.Coordinate;
            }
            else if (edge1.StartNode == edge2.EndNode)
            {
                sharedNode = edge1.StartNode;
                vec1 = edge1.EndNode.Coordinate - edge1.StartNode.Coordinate;
                vec2 = edge2.StartNode.Coordinate - edge2.EndNode.Coordinate;
            }
            else if (edge1.EndNode == edge2.StartNode)
            {
                sharedNode = edge1.EndNode;
                vec1 = edge1.StartNode.Coordinate - edge1.EndNode.Coordinate;
                vec2 = edge2.EndNode.Coordinate - edge2.StartNode.Coordinate;
            }
            else if (edge1.EndNode == edge2.EndNode)
            {
                sharedNode = edge1.EndNode;
                vec1 = edge1.StartNode.Coordinate - edge1.EndNode.Coordinate;
                vec2 = edge2.StartNode.Coordinate - edge2.EndNode.Coordinate;
            }
            return Tuple.Create(vec1, vec2, sharedNode);
        }
        private Vector3d GetVectorOfEdgeFromNode(qEdge edge, qNode node)
        {
            Vector3d vec = Vector3d.Zero;
            if (edge.StartNode == node) { vec = edge.EndNode.Coordinate - node.Coordinate; }
            else if (edge.EndNode == node) { vec = edge.StartNode.Coordinate - node.Coordinate; }
            return vec;
        }
        private qEdge FindEdge(List<qEdge> edgeList, qNode node1, qNode node2)
        {
            qEdge foundEdge = new qEdge();
            foreach (qEdge edge in edgeList)
            {
                if (edge.StartNode == node1 | edge.EndNode == node1)
                {
                    if (edge.StartNode == node2 | edge.EndNode == node2)
                    { foundEdge = edge; } // todo: break?
                }
            }
            if (foundEdge == null) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No edge contains the given nodes."); }
            return foundEdge;
        }
        private void SwapEdge(qEdge E_0)
        {
            // update edge
            List<qNode> swapedNodes = GetSwapedNodes(E_0);
            E_0.StartNode = swapedNodes[0];
            E_0.EndNode = swapedNodes[1];
            E_0.EdgeLine = E_0.VisualizeLine(E_0.StartNode, E_0.EndNode);
            E_0.Length = E_0.CalculateLength(E_0.StartNode, E_0.EndNode);
            // todo: update mesh and neighbor elements
        }
        private Tuple<qEdge, qEdge> SplitEdge(qEdge E_0, Vector3d V_k, qNode N_k, qNode N_m)
        {
            qEdge E_k = new qEdge();
            qEdge E_m = new qEdge();
            qNode N_n = new qNode();

            Line V_k_line = new Line(N_k.Coordinate, V_k);
            Line E_0_line = E_0.EdgeLine;
            NurbsCurve V_k_curve = V_k_line.ToNurbsCurve();
            NurbsCurve E_0_curve = E_0_line.ToNurbsCurve();
            double[] V_k_param = V_k_curve.DivideByCount(100, true);
            Point3d testPoint = new Point3d();
            Point3d V_k_point = new Point3d();
            double distanceToCurve = 100;
            double minDistance = 100;
            double minDistanceParam = 0;
            foreach (double n in V_k_param)
            {
                V_k_point = V_k_curve.PointAt(n);
                E_0_curve.ClosestPoint(V_k_point, out double pointParamOnCurve); // closest parameter on curve
                testPoint = E_0_curve.PointAt(pointParamOnCurve);  // make test point 
                distanceToCurve = testPoint.DistanceTo(V_k_point);

                if (distanceToCurve < minDistance)
                {
                    minDistance = distanceToCurve;
                    minDistanceParam = pointParamOnCurve;
                }
            }

            Point3d N_n_coordinate = E_0_curve.PointAt(minDistanceParam);
            N_n = new qNode() { Coordinate = N_n_coordinate }; // todo: indexing

            E_k = new qEdge(N_n, N_m); 
            E_m = new qEdge(N_k, N_n);
            
            // todo: update mesh, divide E_0 into two edges, create new elements..

            return Tuple.Create(E_k, E_m);
        }
        private qEdge GetSideEdge(Mesh mesh, List<qElement> elementList, List<qEdge> edgeList, double nodeToEvaluate, qEdge E_front)
        {
            // todo: update new elements, new nodes, new connections etc. Might be done later. Have not used N_n and E_m.

            double thetaTolerance = 0.16667 * Math.PI; // todo: an assumption
            qEdge E_neighborFront = new qEdge();
            if (nodeToEvaluate == 0) { E_neighborFront = E_front.LeftFrontNeighbor; ; }
            else { E_neighborFront = E_front.RightFrontNeighbor; }

            #region Get V_k
            var vectorsAndSharedNode = CalculateVectorsFromSharedNode(E_front, E_neighborFront);
            Vector3d vec1 = vectorsAndSharedNode.Item1; // get vector for E_front
            Vector3d vec2 = vectorsAndSharedNode.Item2; // get vector for E_neighborFront
            qNode N_k = vectorsAndSharedNode.Item3; // shared node
            Vector3d V_k = vec1.Length * vec2 + vec2.Length * vec1; // angle bisector
            if (Math.Round(V_k.Length, 2) == 0)
            {
                if (nodeToEvaluate == 0) { V_k = vec1; }
                else { V_k = vec2; }
                V_k.Rotate(0.5 * Math.PI, Vector3d.ZAxis); // todo: not for 3d surface, make axis for normal to plane (vec1, vec2)
            }
            #endregion

            #region Get E_i (connected edges to shared node)
            //int[] connectedEdgesIndex = N_k.ConnectedEdges; // get connected edges indices of shared node
            List<qEdge> connectedEdges = GetConnectedEdges(N_k, edgeList); // get connected edges of shared node
            // get E_i
            List<qEdge> E_i_list = new List<qEdge>();
            foreach (qEdge edge in connectedEdges) // loop connected edges
            {
                if ((edge == E_front) | (edge == E_neighborFront)) { continue; } // continue if connected edge is E_frontEdge or E_neighborFront 

                // check if edge is only connected with triangle element
                List<qElement> connectedElements = new List<qElement>();// get connected element of edge
                connectedElements.Add(edge.Element1);
                if (edge.Element2 == null) { connectedElements.Add(edge.Element2); }

                bool connectedToQuadElement = false;
                foreach (qElement connectedElement in connectedElements) // loop connected elements
                {
                    if (connectedElement.IsQuad)
                    {
                        connectedToQuadElement = true;
                    }
                }
                if (!connectedToQuadElement) { E_i_list.Add(edge); } // add edge if not connected with quad
            }
            #endregion

            #region Find smallest theta to edge
            // calcualte thetas 
            List<double> theta_i_list = new List<double>();
            foreach (qEdge E_i in E_i_list)
            {
                Vector3d E_i_vec = GetVectorOfEdgeFromNode(E_i, N_k);
                double theta_i = Vector3d.VectorAngle(V_k, E_i_vec);
                theta_i_list.Add(theta_i);
            }

            // find smallest theta
            double minTheta = theta_i_list[0];
            int minThetaIndex = 0;
            for (int i = 1; i < theta_i_list.Count; ++i)
            {
                if (theta_i_list[i] < minTheta)
                {
                    minTheta = theta_i_list[i];
                    minThetaIndex = i;
                }
            }
            #endregion

            #region Get E_k
            qEdge E_k = new qEdge();
            qEdge E_0 = new qEdge();
            qEdge E_m = new qEdge();

            if (minTheta < thetaTolerance) // use existing edge
            {
                E_k = E_i_list[minThetaIndex];
            }
            else { E_k = E_i_list[minThetaIndex]; } // todo: change back to swap/split

            /*
            else // swap or split
            {
                if (E_i_list.Count != 2) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Assumption of E_i_list.Count = 2 is wrong."); }

                #region Find E_0
                // find shared node and connected egdes
                qNode E_1_NotSharedNode = new qNode();
                qNode E_2_NotSharedNode = new qNode();
                if (E_i_list[0].StartNode == N_k) { E_1_NotSharedNode = E_i_list[0].EndNode; }
                else { E_1_NotSharedNode = E_i_list[0].StartNode; }
                if (E_i_list[1].StartNode == N_k) { E_2_NotSharedNode = E_i_list[1].EndNode; }
                else { E_2_NotSharedNode = E_i_list[1].StartNode; }

                E_0 = FindEdge(edgeList, E_1_NotSharedNode, E_2_NotSharedNode); // find edge
                #endregion

                // Get N_m
                List<qNode> swapedNodes = GetSwapedNodes(E_0);
                qNode N_m = new qNode();
                if (swapedNodes[0] == N_k) { N_m = swapedNodes[1]; }
                else { N_m = swapedNodes[0]; }

                // Swap or split
                double lengthN_kN_m = N_k.Coordinate.DistanceTo(N_m.Coordinate);
                double beta = Vector3d.VectorAngle(V_k, N_m.Coordinate - N_k.Coordinate);
                if ((lengthN_kN_m < Math.Sqrt(3) * (E_front.Length + E_neighborFront.Length) * 0.5) & beta < thetaTolerance)
                {
                    SwapEdge(E_0);
                    E_k = E_0;
                }
                else
                {
                    var E_kAndE_m = SplitEdge(E_0, V_k, N_k, N_m);
                    E_k = E_kAndE_m.Item1;
                    E_m = E_kAndE_m.Item2;
                }
            }*/
            #endregion
            return E_k;
        }
        private Point3d GetElementCenter(qElement element) // todo: if not use faces..
        {
            double sx = 0;
            double sy = 0;
            double sz = 0;

            List<qEdge> edgeList = element.EdgeList;
            foreach (qEdge edge in edgeList)
            {
                Point3d startPoint = edge.StartNode.Coordinate;
                Point3d endPoint = edge.EndNode.Coordinate;
                List<Point3d> pts = new List<Point3d>() { startPoint, endPoint };
                foreach (Point3d pt in pts)
                {
                    sx = sx + pt.X;
                    sy = sy + pt.Y;
                    sz = sz + pt.Z;
                }
            }
            int n = edgeList.Count*2;
            Point3d centerPt = new Point3d(sx / n, sy / n, sz / n);
            return centerPt;
        }
        private List<qNode> GetNodesOfElement(qElement element)
        {
            List<qNode> nodeList = new List<qNode>();
            foreach (qEdge edge in element.EdgeList)
            {
                if (!nodeList.Contains(edge.StartNode))
                {
                    nodeList.Add(edge.StartNode);
                }
                else if (!nodeList.Contains(edge.EndNode))
                {
                    nodeList.Add(edge.EndNode);
                }
            }
            return nodeList;
        }
        private List<qNode> GetSwapedNodes(qEdge E_0)
        {
            // summary: get the new nodes of E_0 if swaped. Assume E_0 is swapable, i.e. two adjacent triangles.

            List<qNode> swapedNodes = new List<qNode>(); // list of new nodes if swaping
            List<qNode> nodeCandidates = new List<qNode>(); // list of all adjacent element nodes
            List<qElement> adjecentElements = new List<qElement>() { E_0.Element1, E_0.Element2 };

            // find candidate nodes
            foreach (qElement element in adjecentElements)
            {
                List<qNode> elementNodes = GetNodesOfElement(element);
                foreach (qNode node in elementNodes)
                {
                    nodeCandidates.Add(node);
                }
            }

            // find swaped nodes
            foreach (qNode node in nodeCandidates)
            { 
                if ((node != E_0.StartNode) & (node != E_0.EndNode))
                {
                    swapedNodes.Add(node);
                }
            }
            return swapedNodes;
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
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("71512089-a4f6-45fd-8b73-c1c42c18e59e"); }
        }
    }
}
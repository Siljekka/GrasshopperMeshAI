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
            pManager.AddGenericParameter("testList", "tL", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("testItem1", "tI", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("testItem2", "tI", "", GH_ParamAccess.item);
            pManager.AddNumberParameter("number iteration before stop", "temp", "Temoprary: number of iteration to perform", GH_ParamAccess.item, 2);

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
            pManager.AddGenericParameter("testList ", "tL", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("testItem ", "tI", "", GH_ParamAccess.item);

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
            double iterationsToPerformBeforeStop = 0;

            // input
            DA.GetData(0, ref mesh);
            DA.GetData(1, ref numberElementsToRemesh);
            //DA.GetData(2, testList); //to test 
            //DA.GetData(3, ref testItem1); // to test
            //DA.GetData(4, ref testItem2); // to test

            DA.GetData(5, ref iterationsToPerformBeforeStop);
            #region Code

            // Get initial edges and elements using mesh topology properties
            var initialEdgeAndElementList = GetInitialEdgesAndElements(mesh);
            List<qEdge> globalEdgeList = initialEdgeAndElementList.Item1;
            List<qElement> globalElementList = initialEdgeAndElementList.Item2;

            // Get intial front edges
            List<qEdge> frontEdges = GetFrontEdges(globalEdgeList);

            // temporary variables
            qEdge E_k_left = new qEdge(); // left side edge of quad
            qEdge E_k_right = new qEdge(); // right side edge of quad
            qEdge E_front = new qEdge(); // bottom of quad
            qEdge E_top = new qEdge(); // top if quad
            qEdge E_frontFail = null;
            List<qEdge> listE_frontFailed = new List<qEdge>();
            int iterationCounter = 0;

            for (int n = 0; n < numberElementsToRemesh; n++)
            {
                // Temporary stop
                if (iterationCounter == iterationsToPerformBeforeStop) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Iteration stop"); break; } 
                iterationCounter++;

                // back up if selected E_front is not to useable
                List<qElement> globalElementListBackUp = globalElementList;
                List<qEdge> globalEdgeListBackUp = globalEdgeList;
                
                // check for special cases
                // * seams
                // * transision

                // select next front edge
                var E_frontAndEdgeState = SelectNextFrontEdge(frontEdges);
                E_front = E_frontAndEdgeState.Item1;
                var edgeState = E_frontAndEdgeState.Item2;

                // to do: temporay solution for E_frontFail
                /*
                if (E_front == E_frontFail) 
                { 
                    E_front = E_frontFail.RightFrontNeighbor;
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "E_front == E_frontFail, switch to rightneighbor.");
                    edgeState[0] = 0;
                    edgeState[1] = 0;
                }*/

                // get left edge
                switch (edgeState[0])
                {
                    case 0:
                        E_k_left = GetSideEdge(globalElementList, globalEdgeList, 0, E_front); break;
                    case 1:
                        E_k_left = E_front.LeftFrontNeighbor; break;
                }

                // get right edge
                switch (edgeState[1])
                {
                    case 0:
                        E_k_right = GetSideEdge(globalElementList, globalEdgeList, 1, E_front); break;
                    case 1:
                        E_k_right = E_front.RightFrontNeighbor; break;
                }

                // get top edge
                var topEdgeValue = GetTopEdge(E_front, E_k_left, E_k_right, globalEdgeList, globalElementList, frontEdges);
                E_top = topEdgeValue.Item1;
                bool E_top_performCheck = topEdgeValue.Item2;

                // if not possible to perform top recovery skip the chosen E_front and choose new
                if (!E_top_performCheck)
                {
                    listE_frontFailed.Add(E_front);
                    E_frontFail = E_front;
                    globalElementList = globalElementListBackUp; // reset changes made in the iteration
                    globalEdgeList = globalEdgeListBackUp; // reset changes made in the iteration
                    E_frontFail.Level = 100 ; // to do: temporary

                    n--;
                    continue;
                }

                E_frontFail = null;

                // quadrilateral formation
                List<qEdge> quadEdges = new List<qEdge>() { E_front, E_k_right, E_k_left, E_top };
                CreateQuadElement(quadEdges, globalEdgeList, globalElementList);

                // Mesh modification
                frontEdges = GetFrontEdges(globalEdgeList);

                // local smoothing
            }
            #endregion End Code

            List<qEdge> test = new List<qEdge>() { E_front, E_k_left, E_k_right, E_top };

            // todo: when new Level: change qEdge.IsQuadSideEdge = false;

            // output
            DA.SetDataList(0, globalEdgeList);
            DA.SetDataList(1, frontEdges); 
            DA.SetDataList(2, globalElementList); //list10
            DA.SetDataList(3, globalElementList); //list01
            DA.SetDataList(4, globalElementList);
            DA.SetDataList(5, test);
            DA.SetData(6, E_front);
            DA.SetData(7, E_k_left);
            DA.SetData(8, E_k_right);
            
            //DA.SetDataList(9, );
            //DA.SetData(10, );
        }



        #region Methods
        // _____________________________________ for ininital mesh _________________________________________
        private Tuple<List<qEdge>, List<qElement>> GetInitialEdgesAndElements(Mesh mesh)
        {
            // summary: create edges and elements from the initial mesh topology
            List<qEdge> edgeList = CreateInitialEdges(mesh);
            List<qElement> elementList = CreateInitialElements(mesh, edgeList);
            SetNeighborElements(mesh, elementList, edgeList);

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

            bool[] meshVertexBool = mesh.GetNakedEdgePointStatus();
            for (int i = 0; i < nodeList.Count; i++)
            {
                if (meshVertexBool[i] == true) { nodeList[i].BoundaryNode = true; }
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

                qEdge edge = new qEdge(startNode, endNode); //sjekk om i er nødvendig
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
        } // todo: if needed for more cases than initial mesh, make more general
        
        // ________________________________ for front definiton and classification ___________________________________
        private List<qEdge> GetFrontEdges(List<qEdge> globalEdgeList)
        {
            // summary: get all edges connected to a single triangle element
            List<qEdge> frontEdges = new List<qEdge>();
            foreach (qEdge edge in globalEdgeList)
            {
                List<qElement> connectedElements = GetConnectedElements(edge);
                if (connectedElements.Count == 1)
                {
                    if (!edge.Element1.IsQuad)
                    {
                        frontEdges.Add(edge);
                    }
                }
                else if (connectedElements.Count != 1)
                {
                    if (!edge.Element1.IsQuad & edge.Element2.IsQuad) { frontEdges.Add(edge); }
                    else if (edge.Element1.IsQuad & !edge.Element2.IsQuad) { frontEdges.Add(edge); }
                }
            }
            SetNeighorFrontEdges(frontEdges);
            return frontEdges;
        }
        private void SetNeighorFrontEdges(List<qEdge> frontEdges)
        {
            // summary: assign neighbor front edges to each front edge
            for (int i = 0; i < frontEdges.Count; i++)
            {
                qEdge edge = frontEdges[i];
                qEdge neigborEdgeToStartNode = new qEdge();
                qEdge neigborEdgeToEndNode = new qEdge();
                for (int j = 0; j < frontEdges.Count; j++)
                {
                    qEdge testEdge = frontEdges[j];
                    if (i != j)
                    {
                        if ((edge.StartNode == testEdge.StartNode) | (edge.StartNode == testEdge.EndNode)) // check edge start node if connected
                        {
                            neigborEdgeToStartNode = frontEdges[j];
                        }
                        else if ((edge.EndNode == testEdge.StartNode) | (edge.EndNode == testEdge.EndNode)) // check edge end point if connected
                        {
                            neigborEdgeToEndNode = frontEdges[j];
                        }
                    }
                }

                Point3d midPointEdg = (edge.StartNode.Coordinate + edge.EndNode.Coordinate) / 2.0; // mid point of edge

                Point3d centerPoint = GetElementCenter(GetFrontElement(edge));

                Vector3d centerToMidVector = midPointEdg - centerPoint;

                Vector3d centerToEndNodeVector = edge.EndNode.Coordinate - centerPoint;

                Vector3d centerToStartNodeVector = edge.StartNode.Coordinate - centerPoint;

                double startAngle = Vector3d.VectorAngle(centerToMidVector, centerToStartNodeVector, Vector3d.ZAxis); // todo: make normal more general

                double endAngle = Vector3d.VectorAngle(centerToMidVector, centerToEndNodeVector, Vector3d.ZAxis); // todo: make normal more general


                if (endAngle < startAngle)
                {
                    edge.LeftFrontNeighbor = neigborEdgeToStartNode;
                    edge.RightFrontNeighbor = neigborEdgeToEndNode;
                }
                else
                {
                    edge.LeftFrontNeighbor = neigborEdgeToEndNode;
                    edge.RightFrontNeighbor = neigborEdgeToStartNode;
                }
            }

            return;
        }
        Tuple<List<qEdge>, List<qEdge>, List<qEdge>, List<qEdge>> CreateEdgeStateList(List<qEdge> frontEdges)
        {
            // summary: create edge state lists in the format: list11, list01, list10, and list00.
            List<qEdge> list11 = new List<qEdge>();
            List<qEdge> list01 = new List<qEdge>();
            List<qEdge> list10 = new List<qEdge>();
            List<qEdge> list00 = new List<qEdge>();

            double angleTolerance = 0.75 * Math.PI;
            double leftAngle = 0;
            double rightAngle = 0;

            for (int i = 0; i < frontEdges.Count; i++)
            {
                leftAngle = CalculateAngleOfNeighborFrontEdges(0, frontEdges[i]);
                rightAngle = CalculateAngleOfNeighborFrontEdges(1, frontEdges[i]);

                if (leftAngle < angleTolerance & rightAngle < angleTolerance) { list11.Add(frontEdges[i]); }
                else if (leftAngle >= angleTolerance & rightAngle < angleTolerance) { list01.Add(frontEdges[i]); }
                else if (leftAngle < angleTolerance & rightAngle >= angleTolerance) { list10.Add(frontEdges[i]); }
                else if (leftAngle >= angleTolerance & rightAngle >= angleTolerance) { list00.Add(frontEdges[i]); }
            }
            return Tuple.Create(list11, list10, list01, list00);
        }
        private double CalculateAngleOfNeighborFrontEdges(int nodeToCalculate, qEdge edge)
        {
            // summary: calculate the angle between front edges with a shared point, i.e. neighbor front edges

            if (!IsFrontEdge(edge))
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The edge is not a front edge."); }

            qEdge edgeNeighbor = new qEdge();
            double angle = 0;

            // Get neighbor nodes
            if (nodeToCalculate == 0) { edgeNeighbor = edge.LeftFrontNeighbor; }
            else { edgeNeighbor = edge.RightFrontNeighbor; }

            // Create vectors from shared node
            var vectors = CalculateVectorsFromSharedNode(edge, edgeNeighbor);
            Vector3d vec1 = vectors.Item1;
            Vector3d vec2 = vectors.Item2;

            if (nodeToCalculate == 0)
            {
                angle = Vector3d.VectorAngle(vec1, vec2, Vector3d.ZAxis); // to do: make normal more general
            }
            else
            {
                angle = Vector3d.VectorAngle(vec2, vec1, Vector3d.ZAxis); // to do: make normal more general
            }

            /* to do: CHECK IF NEEDED
            qNode sharedNode = vectors.Item3;
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
            Point3d edgeFaceCenter = GetElementCenter(GetFrontElement(edge));
            Point3d edgeNeighborFaceCenter = GetElementCenter(GetFrontElement(edgeNeighbor));
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
            }*/
            return angle;
        }
        private bool IsFrontEdge(qEdge edge)
        {
            bool check = false;
            List<qElement> connectedElements = GetConnectedElements(edge);
            if (connectedElements.Count == 1)
            {
                if (!edge.Element1.IsQuad)
                {
                    check = true;
                }
            }
            else if (connectedElements.Count != 1)
            {
                if (!edge.Element1.IsQuad & edge.Element2.IsQuad) { check = true; }
                else if (edge.Element1.IsQuad & !edge.Element2.IsQuad) { check = true; }
            }
            return check;
        }

        private qEdge SelectFrontEdgeFromList(List<qEdge> edgeStateList)
        {
            // summary: select an edge from edge state list input. Choose lowest level, if multiple choose shortest of them.

            // get edges of lowest level
            List<qEdge> edgesOfLowestLevel = new List<qEdge>();
            int lowestLevel = edgeStateList[0].Level;
            foreach (qEdge edge in edgeStateList)
            {
                if (edge.Level == lowestLevel)
                {
                    edgesOfLowestLevel.Add(edge);
                }
                else if (edge.Level < lowestLevel)
                {
                    edgesOfLowestLevel.Clear();
                    lowestLevel = edge.Level;
                    edgesOfLowestLevel.Add(edge);
                }
            }

            // get shortes edge
            qEdge E_front = new qEdge();
            double shortestLength = edgesOfLowestLevel[0].Length;
            foreach (qEdge edge in edgesOfLowestLevel)
            {
                if (edge.Length <= shortestLength)
                {
                    E_front = edge;
                }
            }
            return E_front;
        }

        private Tuple<qEdge, int[] > SelectNextFrontEdge(List<qEdge> frontEdges)
        {
            // summary: select next front edge with the hierarchy: list -> lowest level -> shortest length and switch to neighbor if large transision
            double transitionTolerance = 2.5;
            qEdge E_front = new qEdge();

            var edgeStates = CreateEdgeStateList(frontEdges);
            var list11 = edgeStates.Item1;
            var list10 = edgeStates.Item2;
            var list01 = edgeStates.Item3;
            var list00 = edgeStates.Item4;
            int[] edgeState = { 0, 0 };

            // select
            if (list11.Count != 0)
            {
                E_front = SelectFrontEdgeFromList(list11);
                edgeState[0] = 1; edgeState[1] = 1;

            }
            else if (list01.Count != 0 | list10.Count != 0)
            {
                // merge list01 and 10
                List<qEdge> list0110 = new List<qEdge>();
                foreach (qEdge edge in list01)
                {
                    list0110.Add(edge);
                }
                foreach (qEdge edge in list10)
                {
                    list0110.Add(edge);
                }

                E_front = SelectFrontEdgeFromList(list0110);

                // set edge state
                if (list01.Contains(E_front))
                {
                    edgeState[0] = 0; edgeState[1] = 1;
                }
                else
                {
                    edgeState[0] = 1; edgeState[1] = 0;
                }
            }
            else if (list00.Count == 0)
            {
                E_front = SelectFrontEdgeFromList(list00);
                edgeState[0] = 0; edgeState[1] = 0;
            }
            else
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "SelectNextFrontEdge: No more edges to select."); }

            // switch E_front to neighbor front if transision to one of the front neighbors is large and neighbor is selectable
            if (E_front.Length / E_front.LeftFrontNeighbor.Length > transitionTolerance)
            {
                if (!list11.Contains(E_front.LeftFrontNeighbor))
                {
                    E_front = E_front.LeftFrontNeighbor;
                    { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "SelectNextFrontEdge: Large transision to left neighor. Switched."); }
                }
            }
            else if (E_front.Length / E_front.RightFrontNeighbor.Length > transitionTolerance)
            {
                if (!list11.Contains(E_front.RightFrontNeighbor))
                {
                    E_front = E_front.RightFrontNeighbor;
                    { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "SelectNextFrontEdge: Large transision to right neighor. Switched."); }
                }
            }

            return Tuple.Create(E_front, edgeState);
        }

        // _________________________________________ for topology  _______________________________________________
        private List<qEdge> GetConnectedEdges(qNode node, List<qEdge> globalEdgeList)
        {
            // summary: get connected edges to a node
            List<qEdge> connectedEdges = new List<qEdge>();
            foreach (qEdge edge in globalEdgeList)
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
            // summary: get conneccted elements to an edge. Assume edge has updated elements element 1 and/or element 2.
            List<qElement> connectedElements = new List<qElement>();
            connectedElements.Add(edge.Element1);
            if (edge.Element2 != null)
            {
                connectedElements.Add(edge.Element2);
            }
            return connectedElements;
        }
        private qNode GetNodeOfFrontEdge(int nodeToEvaluate, qEdge edge)
        {
            // summary: get a node of a front edge, 0=left, 1=right
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
            // summary: calculate vectors for two edges with a shared node with direction away from the shared node.
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
        private Point3d GetElementCenter(qElement element)
        {
            // summary: get center of an element
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
            int n = edgeList.Count * 2;
            Point3d centerPt = new Point3d(sx / n, sy / n, sz / n);
            return centerPt;
        }
        private List<qNode> GetNodesOfElement(qElement element)
        {
            // summary: get nodes of an element
            List<qNode> nodeList = new List<qNode>();
            if (!element.IsQuad)
            {

                foreach (qEdge edge in element.EdgeList)
                {
                    if (!nodeList.Contains(edge.StartNode))
                    {
                        nodeList.Add(edge.StartNode);
                    }
                   
                    if (!nodeList.Contains(edge.EndNode))
                    {
                        nodeList.Add(edge.EndNode);
                    }
                }
            }
            else if (element.IsQuad) //todo: test function
            {
                List<qEdge> quadEdges = element.EdgeList;

                qEdge baseEdge = quadEdges[0];
                qEdge rightEdge = quadEdges[1];
                qEdge leftEdge = quadEdges[2];
                qEdge topEdge = quadEdges[3];

                qNode node1 = new qNode();
                qNode node2 = new qNode();
                qNode node3 = new qNode();
                qNode node4 = new qNode();

                if (baseEdge.StartNode == leftEdge.StartNode || baseEdge.StartNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.StartNode;
                    node2 = baseEdge.EndNode;
                }
                else if (baseEdge.EndNode == leftEdge.StartNode || baseEdge.EndNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.EndNode;
                    node2 = baseEdge.StartNode;
                }

                if (topEdge.StartNode == leftEdge.StartNode || topEdge.StartNode == leftEdge.EndNode)
                {
                    node3 = topEdge.EndNode;
                    node4 = topEdge.StartNode;
                }
                else if (topEdge.EndNode == leftEdge.StartNode || topEdge.EndNode == leftEdge.EndNode)
                {
                    node3 = topEdge.StartNode;
                    node4 = topEdge.EndNode;
                }

                nodeList = new List<qNode> { node1, node2, node3, node4 }; // n1: bottom left, n2: bottom right, n2: top right, n3: top left
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
        private Vector3d GetVectorOfEdgeFromNode(qEdge edge, qNode node)
        {
            Vector3d vec = Vector3d.Zero;
            if (edge.StartNode == node) { vec = edge.EndNode.Coordinate - node.Coordinate; }
            else if (edge.EndNode == node) { vec = edge.StartNode.Coordinate - node.Coordinate; }
            return vec;
        }
        private qEdge FindEdge(List<qEdge> globalEdgeList, qNode node1, qNode node2)
        {
            qEdge foundEdge = new qEdge();
            foreach (qEdge edge in globalEdgeList)
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
        private List<qElement> GetNeighborElements(qElement element)
        {
            List<qElement> neighborElements = new List<qElement>();

            foreach (qEdge edge in element.EdgeList)
            {
                List<qElement> connectedElements = GetConnectedElements(edge);
                foreach (qElement elementCandidate in connectedElements)
                {
                    if (elementCandidate != element)
                    {
                        neighborElements.Add(elementCandidate);
                    }
                }

            }
            return neighborElements;
        }
        private qNode GetOppositeNode(qNode node, qEdge edge)
        {
            qNode oppositeNode = new qNode();
            if (node == edge.StartNode) { oppositeNode = edge.EndNode; }
            else if (node == edge.EndNode) { oppositeNode = edge.StartNode; }
            return oppositeNode;
        }
        private qElement GetFrontElement(qEdge edge)
        {
            // summary: get the connected triangle element of a front edge
            qElement triangleElement = new qElement();
            if (GetConnectedElements(edge).Count == 1) { triangleElement = edge.Element1; }
            else if (!edge.Element1.IsQuad) { triangleElement = edge.Element1; }
            else { triangleElement = edge.Element2; }

            return triangleElement;
        }

        // _______________________________________ for mesh modification __________________________________________________
        private qEdge GetSideEdge(List<qElement> globalElementList, List<qEdge> globalEdgeList, double nodeToEvaluate, qEdge E_front)
        {
            // summary: get side edge of a new quad; nodeToEvaluate: 0 = left, 1 = right;
            double thetaTolerance = 0.16667 * Math.PI; // todo: an assumption
            qEdge E_neighborFront = new qEdge();
            if (nodeToEvaluate == 0) { E_neighborFront = E_front.LeftFrontNeighbor; ; }
            else { E_neighborFront = E_front.RightFrontNeighbor; }

            #region Get V_k
            var vectorsAndSharedNode = CalculateVectorsFromSharedNode(E_front, E_neighborFront);
            Vector3d vec1 = vectorsAndSharedNode.Item1; // get vector for E_front
            Vector3d vec2 = vectorsAndSharedNode.Item2; // get vector for E_neighborFront
            qNode N_k = vectorsAndSharedNode.Item3; // shared node
            Vector3d V_k = Vector3d.Zero;

            if (nodeToEvaluate == 0)
            {
                double angle = Vector3d.VectorAngle(vec1, vec2, Vector3d.ZAxis); // to do: make normal mor general
                V_k = vec1;
                V_k.Rotate(0.5 * angle, Vector3d.ZAxis); // todo: not for 3d surface, make axis for normal to plane (vec1, vec2)
            }
            else
            {
                double angle = Vector3d.VectorAngle(vec2, vec1, Vector3d.ZAxis); // to do: make normal mor general
                V_k = vec2;
                V_k.Rotate(0.5 * angle, Vector3d.ZAxis); // todo: not for 3d surface, make axis for normal to plane (vec1, vec2)
            }
            #endregion

            #region Get E_i (connected edges to shared node)
            List<qEdge> connectedEdges = GetConnectedEdges(N_k, globalEdgeList); // get connected edges of shared node
            // get E_i
            List<qEdge> E_i_list = new List<qEdge>();
            foreach (qEdge edge in connectedEdges) // loop connected edges
            {
                if ((edge == E_front) | (edge == E_neighborFront)) { continue; } // continue if connected edge is E_frontEdge or E_neighborFront 

                // check if edge is only connected with triangle element
                List<qElement> connectedElements = GetConnectedElements(edge);

                bool connectedToQuadElement = false;
                foreach (qElement connectedElement in connectedElements) // loop element 1 and (if it exists) element 2 of a conncted edge
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
                double theta_i = Vector3d.VectorAngle(V_k, E_i_vec); // to do: make more general, assume there are edges close to V_K
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
            else // swap or split
            {
                if (E_i_list.Count != 2) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "GetSideEdge: Assumption of E_i_list.Count = 2 is wrong."); }

                #region Find E_0
                // find shared node and connected egdes
                qNode E_1_NotSharedNode = new qNode();
                qNode E_2_NotSharedNode = new qNode();
                if (E_i_list[0].StartNode == N_k) { E_1_NotSharedNode = E_i_list[0].EndNode; }
                else { E_1_NotSharedNode = E_i_list[0].StartNode; }
                if (E_i_list[1].StartNode == N_k) { E_2_NotSharedNode = E_i_list[1].EndNode; }
                else { E_2_NotSharedNode = E_i_list[1].StartNode; }

                E_0 = FindEdge(globalEdgeList, E_1_NotSharedNode, E_2_NotSharedNode); // find edge
                #endregion

                // Get N_m
                qNode N_m = new qNode();
                List<qNode> swapedNodes = GetSwapedNodes(E_0);
                if (swapedNodes[0] == N_k) { N_m = swapedNodes[1]; }
                else { N_m = swapedNodes[0]; }

                // Swap or split
                double lengthN_kN_m = N_k.Coordinate.DistanceTo(N_m.Coordinate);
                double beta = Vector3d.VectorAngle(V_k, N_m.Coordinate - N_k.Coordinate);
                
                if ((lengthN_kN_m < Math.Sqrt(3) * (E_front.Length + E_neighborFront.Length) * 0.5) & beta < thetaTolerance)
                {
                    SwapEdge(E_0, globalElementList);
                    E_k = E_0;
                }
                else
                {
                    E_k = SplitEdge(E_0, V_k, N_k, globalEdgeList, globalElementList);
                }
            }
            #endregion
            return E_k;
        }
        private void SwapEdge(qEdge E_0, List<qElement> globalElementList)
        {
            // is swapable
            if (E_0.Element1.IsQuad | E_0.Element2.IsQuad)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Edge is not swapable.");
            }

            // get nodes of edge to swap
            qNode node1 = E_0.StartNode;
            qNode node2 = E_0.EndNode;
            List<qNode> swapedNodes = GetSwapedNodes(E_0);

            // get edges of adjacent elements
            List<qEdge> element1edges = E_0.Element1.EdgeList;
            List<qEdge> element2edges = E_0.Element2.EdgeList;

            // get list of all edges for new elements
            List<qEdge> edgeForNewElements = new List<qEdge>();
            element1edges.RemoveAt(element1edges.IndexOf(E_0));
            element2edges.RemoveAt(element2edges.IndexOf(E_0));
            for (int i = 0; i < 2; i++)
            {
                edgeForNewElements.Add(element1edges[i]);
                edgeForNewElements.Add(element2edges[i]);
            }

            qElement oldElement1 = E_0.Element1;
            qElement oldElement2 = E_0.Element2;

            // update edge
            E_0.StartNode = swapedNodes[0];
            E_0.EndNode = swapedNodes[1];
            E_0.EdgeLine = E_0.VisualizeLine(E_0.StartNode, E_0.EndNode);
            E_0.Length = E_0.CalculateLength(E_0.StartNode, E_0.EndNode);
            E_0.Element1 = null;
            E_0.Element2 = null;

            // create new elements
            List<qEdge> element1edgesNew = new List<qEdge>();
            List<qEdge> element2edgesNew = new List<qEdge>();
            element1edgesNew.Add(E_0);
            element2edgesNew.Add(E_0);

            foreach (qEdge edge in edgeForNewElements)
            {
                if (edge.StartNode == node1 | edge.EndNode == node1)
                {
                    element1edgesNew.Add(edge);
                }
                else if (edge.StartNode == node2 | edge.EndNode == node2)
                {
                    element2edgesNew.Add(edge);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Swaping failed.");
                }
            }

            qElement newElement1 = new qElement(element1edgesNew);
            qElement newElement2 = new qElement(element2edgesNew);
            globalElementList.Add(newElement1);
            globalElementList.Add(newElement2);

            // update edge elements
            E_0.Element1 = newElement1;
            E_0.Element2 = newElement2;

            // element update for element edges
            foreach (qEdge edge in element1edgesNew)
            {
                if ((edge.Element1 == oldElement1) | (edge.Element1 == oldElement2))
                {
                    edge.Element1 = newElement1;
                }
                else if ((edge.Element2 == oldElement1) | (edge.Element2 == oldElement2))
                {
                    edge.Element2 = newElement1;

                }
            }

            foreach (qEdge edge in element2edgesNew)
            {
                if ((edge.Element1 == oldElement1) | (edge.Element1 == oldElement2))
                {
                    edge.Element1 = newElement2;
                }
                else if ((edge.Element2 == oldElement1) | (edge.Element2 == oldElement2))
                {
                    edge.Element2 = newElement2;

                }
            }

            // remove old elements from list of elements
            globalElementList.RemoveAt(globalElementList.IndexOf(oldElement1));
            globalElementList.RemoveAt(globalElementList.IndexOf(oldElement2));
        }
        private qEdge SplitEdge(qEdge E_0, Vector3d V_k, qNode N_k, List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            // summary: split element

            // get N_m
            qNode N_m = new qNode();
            List<qNode> swapedNodes = GetSwapedNodes(E_0);
            if (swapedNodes[0] == N_k) { N_m = swapedNodes[1]; }
            else { N_m = swapedNodes[0]; }

            // get N_n
            Line V_k_line = new Line(N_k.Coordinate, V_k);
            Line E_0_line = E_0.EdgeLine;
            NurbsCurve V_k_curve = V_k_line.ToNurbsCurve();
            NurbsCurve E_0_curve = E_0_line.ToNurbsCurve();
            double[] V_k_param = V_k_curve.DivideByCount(100, true);
            double minDistance = 100;
            double minDistanceParam = 0;
            foreach (double n in V_k_param)
            {
                Point3d V_k_point = V_k_curve.PointAt(n);
                E_0_curve.ClosestPoint(V_k_point, out double pointParamOnCurve); // closest parameter on curve
                Point3d testPoint = E_0_curve.PointAt(pointParamOnCurve);  // make test point 
                double distanceToCurve = testPoint.DistanceTo(V_k_point);

                if (distanceToCurve < minDistance)
                {
                    minDistance = distanceToCurve;
                    minDistanceParam = pointParamOnCurve;
                }
            }
            Point3d N_n_coordinate = E_0_curve.PointAt(minDistanceParam);
            qNode N_n = new qNode(N_n_coordinate, false);

            // create new edges
            qEdge E_k = new qEdge(N_n, N_k); 
            qEdge E_m = new qEdge(N_m, N_n);
            qEdge E_n_1 = new qEdge(E_0.EndNode, N_n);
            qEdge E_n_2 = new qEdge(E_0.StartNode, N_n);

            List<qEdge> newEdges = new List<qEdge>() { E_k, E_m, E_n_1, E_n_2 };


            // get kept edges from old elements
            List<qEdge> element1edges = E_0.Element1.EdgeList;
            List<qEdge> element2edges = E_0.Element2.EdgeList;
            element1edges.RemoveAt(element1edges.IndexOf(E_0));
            element2edges.RemoveAt(element2edges.IndexOf(E_0));

            List<qEdge> edgesFromOldToNewElements = new List<qEdge>();
            for (int i = 0; i < 2; i++)
            {
                edgesFromOldToNewElements.Add(element1edges[i]);
                edgesFromOldToNewElements.Add(element2edges[i]);
            }


            // create edge lists for new elements
            List<qEdge> newElement1Edges = new List<qEdge>() { E_k, E_n_1 };
            List<qEdge> newElement2Edges = new List<qEdge>() { E_k, E_n_2 };
            List<qEdge> newElement3Edges = new List<qEdge>() { E_m, E_n_2 };
            List<qEdge> newElement4Edges = new List<qEdge>() { E_m, E_n_1 };

            foreach (qEdge edge in edgesFromOldToNewElements)
            {
                if (edge.StartNode == N_k | edge.EndNode == N_k)
                {
                    if (edge.StartNode == GetOppositeNode(N_n, E_n_1) | edge.EndNode == GetOppositeNode(N_n, E_n_1))
                    {
                        newElement1Edges.Add(edge); // element 1
                    }
                    else
                    {
                        newElement2Edges.Add(edge); // element 2 
                    }
                }
                else if (edge.StartNode == N_m | edge.EndNode == N_m)
                {
                    if (edge.StartNode == GetOppositeNode(N_n, E_n_1) | edge.EndNode == GetOppositeNode(N_n, E_n_1))
                    {
                        newElement4Edges.Add(edge); // element 4
                    }
                    else
                    {
                        newElement3Edges.Add(edge); // element 3
                    }
                }
            
            }


            // create new elements
            qElement newElement1 = new qElement(newElement1Edges);
            qElement newElement2 = new qElement(newElement2Edges);
            qElement newElement3 = new qElement(newElement3Edges);
            qElement newElement4 = new qElement(newElement4Edges);

            List<qElement> newElements = new List<qElement>() { newElement1, newElement2, newElement3, newElement4 };


            // locate old elements
            qElement oldElement1 = E_0.Element1;
            qElement oldElement2 = E_0.Element2;

            List<qNode> nodesOfOldElement1 = new List<qNode>();
            List<qNode> nodesOfOldElement2 = new List<qNode>();

            foreach (qEdge edge in oldElement1.EdgeList)
            {
                nodesOfOldElement1.Add(edge.EndNode);
                nodesOfOldElement1.Add(edge.StartNode);
            }

            foreach (qEdge edge in oldElement2.EdgeList)
            {
                nodesOfOldElement2.Add(edge.EndNode);
                nodesOfOldElement2.Add(edge.StartNode);
            }

            qElement elementWithN_k = new qElement();
            qElement elementWithN_m = new qElement();

            if (nodesOfOldElement1.Contains(N_k))
            { 
                elementWithN_k = oldElement1;
                elementWithN_m = oldElement2;
            }
            else if (nodesOfOldElement2.Contains(N_k))
            { 
                elementWithN_k = oldElement2;
                elementWithN_m = oldElement1;
            }

            // update elements to kept edges 
            for (int i = 0; i < 4; i++)
            {
                qElement newElement = newElements[i];
                qEdge keptEdge = newElement.EdgeList[2]; // get the edge kept from old to new elements
                if (i < 2) // if element 1 or 2: replace oldElement 
                {
                    if (keptEdge.Element1 == elementWithN_k)
                    {
                        keptEdge.Element1 = newElement;
                    }
                    else if (keptEdge.Element2 == elementWithN_k)
                    {
                        keptEdge.Element2 = newElement;
                    }
                    else
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Split: update old edges failed with elementWithN_k");
                    }
                }
                else
                {
                    if (keptEdge.Element1 == elementWithN_m)
                    {
                        keptEdge.Element1 = newElement;
                    }
                    else if (keptEdge.Element2 == elementWithN_m)
                    {
                        keptEdge.Element2 = newElement;
                    }
                    else
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Split: update old edges failed with elementWithN_m");
                    }
                }
            }

            // update elements for new edges
            E_n_1.Element1 = newElement1;
            E_n_1.Element2 = newElement4;

            E_n_2.Element1 = newElement2;
            E_n_2.Element2 = newElement3;

            E_k.Element1 = newElement1;
            E_k.Element2 = newElement2;

            E_m.Element1 = newElement3;
            E_m.Element2 = newElement4;


            // update edgeList
            globalEdgeList.RemoveAt(globalEdgeList.IndexOf(E_0));
            foreach (qEdge newEdge in newEdges)
            {
                globalEdgeList.Add(newEdge);
            }

            // update elementList
            globalElementList.RemoveAt(globalElementList.IndexOf(E_0.Element1));
            globalElementList.RemoveAt(globalElementList.IndexOf(E_0.Element2));
            foreach (qElement newElement in newElements)
            {
                globalElementList.Add(newElement);
            }

            return E_k;
        }
        private Tuple<qEdge, bool> GetTopEdge(qEdge E_front, qEdge E_k_left, qEdge E_k_right, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: top edge recovery
            bool performed = true;
            qEdge E_top = new qEdge();

            // Get N_d and N_c
            qNode N_d = new qNode();
            qNode N_c = new qNode();

            qNode leftNode = GetNodeOfFrontEdge(0, E_front);
            qNode rightNode = GetNodeOfFrontEdge(1, E_front);

            if (E_k_left.StartNode == leftNode) { N_d = E_k_left.EndNode; }
            else { N_d = E_k_left.StartNode; }

            if (E_k_right.StartNode == rightNode) { N_c = E_k_right.EndNode; }
            else { N_c = E_k_right.StartNode; }

            
            // to do: only temporary:If left and right edge share node, need to find a solution
            if (N_c == N_d)
            {
                performed = false;
                { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Top edge recovery not performed because N_C and N_d is equal."); }
                return Tuple.Create(E_top, performed);
            }

            var E_topAndPerfomed = EdgeRecoveryProcess( N_c, N_d, globalEdgeList, globalElementList, frontEdges);
            E_top = E_topAndPerfomed.Item1;
            performed = E_topAndPerfomed.Item2;
            return Tuple.Create(E_top, performed);
        }
        private Tuple<qEdge, bool> EdgeRecoveryProcess(qNode N_c, qNode N_d, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            qEdge E_recovered = new qEdge();
            bool performed = true;

            #region Get connected elements to node _N_c
            List<qElement> connectedElements = new List<qElement>();
            List<qEdge> connectedEdges = GetConnectedEdges(N_c, globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                List<qElement> edgeElements = GetConnectedElements(edge);
                foreach (qElement element in edgeElements)
                {
                    if (!connectedElements.Contains(element))
                    {
                        connectedElements.Add(element);
                    }
                }
            }
            #endregion

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
            List<qElement> sortedConnectedElements = new List<qElement>();
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

            #region Find intersected edges of line S from N_c to N_d
            qElement T_k = new qElement();
            qElement T_k1 = new qElement();
            Vector3d V_k = Vector3d.Zero;
            Vector3d V_k1 = Vector3d.Zero;
            Vector3d V_s = N_d.Coordinate - N_c.Coordinate;

            qEdge E_i = new qEdge();
            qEdge E_k = new qEdge();
            qEdge E_k1 = new qEdge();
            foreach (qElement connectedElement in sortedConnectedElements)
            {
                T_k = connectedElement;
                if (connectedElement.IsQuad)
                {
                    continue;
                }

                // get edges
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

                if (V_k.IsParallelTo(V_s) == 1)
                {
                    E_recovered = E_k;
                    return Tuple.Create(E_recovered, performed);
                }
                else if (V_k1.IsParallelTo(V_s) == 1)
                {
                    E_recovered = E_k1;
                    return Tuple.Create(E_recovered, performed);
                }
                else if (Vector3d.Multiply(V_s, V_k) > 0 & Vector3d.Multiply(V_s, V_k1) > 0) // todo: riktig ift Owen??
                {
                    break;
                }
            }

            // add E_i to intersection of S if not a front edge
            List<qEdge> intersectingS = new List<qEdge>();
            if (!frontEdges.Contains(E_i))
            {
                intersectingS.Add(E_i);
            }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Top edge recovery failed.");
            }

            bool done = false;
            while (!done) // todo: bug. E_i may only have 1 elements..
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

                if (E_nCandidates.Count < 2)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Top edge: Fails because E_nCandidates is zero");
                    performed = false;
                    return Tuple.Create(E_recovered, performed);
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
                if (!IsFrontEdge(E_i)) { intersectingS.Add(E_i); }
            }
            #endregion

            #region Edge recovery process
            if (intersectingS.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Top edge: No intersection of S");
                performed = false;
                return Tuple.Create(E_recovered, performed);
            }

            foreach (qEdge edge in intersectingS)
            {
                if (edge.Element1.IsQuad | edge.Element2.IsQuad)
                {
                    performed = false;
                    { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Top edge recovery not performed because of intersection of quads."); }
                    return Tuple.Create(E_recovered, performed);
                }

            }

            while (intersectingS.Count > 0)
            {
                E_recovered = intersectingS[0];
                SwapEdge(E_recovered, globalElementList); // swap the edge

                List<qElement> Tswaped = new List<qElement>() { E_recovered.Element1, E_recovered.Element2 };

                // todo: check area of Tswaped
                // check if inverted...
                bool areaOK = true;
                if (areaOK)
                {
                    intersectingS.RemoveAt(0);
                }
                else
                {
                    intersectingS.Add(E_recovered);
                    intersectingS.RemoveAt(0);
                }
            }
            #endregion

            return Tuple.Create(E_recovered, performed);
        }
        private void Seam(qEdge rightEdge, qEdge leftEdge, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: seam right edge and left edge together to a new seam edge with the adjacent quad elements

            // get needed nodes for seaming
            var vectorsAndShareNode = CalculateVectorsFromSharedNode(rightEdge, leftEdge);
            qNode N_k = vectorsAndShareNode.Item3;
            qNode N_k_left = GetOppositeNode(N_k, leftEdge);
            qNode N_k_right = GetOppositeNode(N_k, rightEdge);

            // edge recovery
            var edgeAndPerformed = EdgeRecoveryProcess(N_k_right, N_k_left, globalEdgeList, globalElementList, frontEdges);
            qEdge E_0 = edgeAndPerformed.Item1;

            // check if edge recovery failed
            bool performed = edgeAndPerformed.Item2;
            if (!performed) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Seam: Edge recovery failed."); return; }

            #region Merge elements
            // get N_new
            Point3d newNodeCoordinate = 0.5 * (E_0.EndNode.Coordinate + E_0.StartNode.Coordinate);
            qNode N_new = new qNode(newNodeCoordinate, false);

            // get N_t
            qNode N_t = new qNode();
            List<qNode> swapedNodes = GetSwapedNodes(E_0);
            if (N_k == swapedNodes[0]) { N_t = swapedNodes[1]; }
            else if (N_k == swapedNodes[1]) { N_t = swapedNodes[0]; }

            // create new seam edges
            qEdge seamEdge = new qEdge(N_k, N_new);
            qEdge connectedSeamEdge = new qEdge(N_new, N_t);

            // get elements to E_0 
            qElement elementToN_k = new qElement(); 
            qElement elementToN_t = new qElement();
            List<qElement> connectedElementsToE_0 = GetConnectedElements(E_0);
            if (GetNodesOfElement(connectedElementsToE_0[0]).Contains(N_k))
            {
                elementToN_k = connectedElementsToE_0[0]; elementToN_t = connectedElementsToE_0[1];
            }
            else if (GetNodesOfElement(connectedElementsToE_0[0]).Contains(N_t))
            {
                elementToN_t = connectedElementsToE_0[0]; elementToN_k = connectedElementsToE_0[1];
            }
            else  { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Seam: Problems with elements connected to E_0."); return; }

            #region Merge elements
            List<qEdge> edgesToSeam = new List<qEdge>();

            // fix elements to N_k
            int indexToUpdate = 0;
            leftEdge = FindEdge(globalEdgeList, N_k, N_k_left); // in case any changes to edge
            rightEdge = FindEdge(globalEdgeList, N_k, N_k_right); // in case any changes to edge
            foreach (qEdge edge in elementToN_k.EdgeList)
            {
                if (edge == E_0) { continue; }

                if (edge == leftEdge)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(leftEdge);
                        edge.Element1.EdgeList[indexToUpdate] = seamEdge;
                        edgesToSeam.Add(edge);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(leftEdge);
                        edge.Element2.EdgeList[indexToUpdate] = seamEdge;
                        edgesToSeam.Add(edge);
                    }
                }
                
                if (edge == rightEdge)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(rightEdge);
                        edge.Element1.EdgeList[indexToUpdate] = seamEdge;
                        edgesToSeam.Add(edge);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(rightEdge);
                        edge.Element2.EdgeList[indexToUpdate] = seamEdge;
                        edgesToSeam.Add(edge);
                    }
                }
            }
            foreach (qEdge edge in edgesToSeam)
            {
                // set element 1 and element 2 to edgesToSeam
            }

            // fix elements to N_t
            indexToUpdate = 0;
            edgesToSeam.Clear();
            leftEdge = FindEdge(globalEdgeList, N_t, N_k_left); // in case any changes to edge
            rightEdge = FindEdge(globalEdgeList, N_t, N_k_right); // in case any changes to edge
            foreach (qEdge edge in elementToN_t.EdgeList)
            {
                if (edge == E_0) { continue; }

                if (edge == leftEdge)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(leftEdge);
                        edge.Element1.EdgeList[indexToUpdate] = connectedSeamEdge;
                        edgesToSeam.Add(edge);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(leftEdge);
                        edge.Element2.EdgeList[indexToUpdate] = connectedSeamEdge;
                        edgesToSeam.Add(edge);
                    }
                }

                if (edge == rightEdge)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(rightEdge);
                        edge.Element1.EdgeList[indexToUpdate] = connectedSeamEdge;
                        edgesToSeam.Add(edge);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(rightEdge);
                        edge.Element2.EdgeList[indexToUpdate] = connectedSeamEdge;
                        edgesToSeam.Add(edge);
                    }
                }
            }

            foreach (qEdge edge in edgesToSeam)
            {
                // set element 1 and element 2 to edgesToSeam
            }
            #endregion 


            // update global element list
            int indexToRemove = globalElementList.IndexOf(connectedElementsToE_0[0]);
            globalElementList.RemoveAt(indexToRemove);

            indexToRemove = globalElementList.IndexOf(connectedElementsToE_0[1]);
            globalElementList.RemoveAt(indexToRemove);

            // update global edge list
            indexToRemove = globalEdgeList.IndexOf(E_0);
            globalEdgeList.RemoveAt(indexToRemove);
            #endregion
        }
        private void CreateQuadElement(List<qEdge> quadEdge, List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            // get inside elements
            qEdge E_front = quadEdge[0];
            List<qElement> elementInside = new List<qElement>() { GetFrontElement(E_front) };
            qElement startElement = GetFrontElement(E_front); // start with front Edge 
            bool done = false;
            int count = 0;
            while (!done)
            {
                foreach (qEdge edge in startElement.EdgeList)
                {
                    if (quadEdge.Contains(edge)) { continue; }

                    List<qElement> connectedElements = GetConnectedElements(edge);
                    foreach (qElement element in connectedElements)
                    {
                        if (!elementInside.Contains(element)) { elementInside.Add(element); }
                    }
                }
                qElement nextElement = elementInside[elementInside.Count - 1];

                if (nextElement == startElement) { done = true; }
                else { startElement = nextElement; }

                if (count > 20) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Infinity loop"); break; }
                count++;
            }
            // to do: find a better way to solve "hided" elements. Ask Silje.
            startElement = elementInside[1]; // start with element that has maybe not been used in loop above
            done = false;
            count = 0;
            while (!done)
            {
                foreach (qEdge edge in startElement.EdgeList)
                {
                    if (quadEdge.Contains(edge)) { continue; }

                    List<qElement> connectedElements = GetConnectedElements(edge);
                    foreach (qElement element in connectedElements)
                    {
                        if (!elementInside.Contains(element)) { elementInside.Add(element); }
                    }
                }
                qElement nextElement = elementInside[elementInside.Count - 1];

                if (nextElement == startElement) { done = true; }
                else { startElement = nextElement; }

                if (count > 20) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Infinity loop"); break; }
                count++;
            }

            // remove inside elements
            foreach (qElement element in elementInside)
            {
                List<qEdge> elementEdgeCody = element.EdgeList;
                foreach (qEdge edge in elementEdgeCody)
                {

                    if (edge.Element1 == element) { edge.Element1 = null; }
                    else { edge.Element2 = null; }

                    if (!quadEdge.Contains(edge) & globalEdgeList.Contains(edge))
                    {
                        globalEdgeList.RemoveAt(globalEdgeList.IndexOf(edge));
                    }
                }
                globalElementList.RemoveAt(globalElementList.IndexOf(element));
            }

            // create new element
            qElement newQuadElement = new qElement(quadEdge);
            globalElementList.Add(newQuadElement);

            // update connected edges
            foreach (qEdge edge in quadEdge)
            {
                if (edge.Element1 == null)
                {
                    edge.Element1 = newQuadElement;
                }
                else if (edge.Element2 == null)
                {
                    edge.Element2 = newQuadElement;
                }
                else { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Quadelement not assigned to edges"); }
            }
            quadEdge[1].IsQuadSideEdge = true;
            quadEdge[2].IsQuadSideEdge = true; // make sure not pick as new front edge, but is in frontEdges
            quadEdge[3].Level++; // update level for top edge of quad
        }

        // __________________________________________ Local smoothing ______________________________________________________

        private void DoLocalSmoothing(qElement QuadElement, List<qEdge> globalEdgeList, List<qEdge> frontEdges)
        {
            qNode smoothNode = new qNode();
            List<qNode> adjacentNodes = new List<qNode>();
            List<qEdge> connectedEdges = new List<qEdge>();

            List<qEdge> QuadEdges = QuadElement.EdgeList;
            List<qNode> QuadNodes = GetNodesOfElement(QuadElement);


            qEdge baseEdge = QuadEdges[0];
            qEdge rightEdge = QuadEdges[1];
            qEdge leftEdge = QuadEdges[2];
            qEdge topEdge = QuadEdges[3];

            qNode node1 = QuadNodes[0];
            qNode node2 = QuadNodes[1];
            qNode node3 = QuadNodes[2];
            qNode node4 = QuadNodes[3];

            // Smooth nodes on Quad
            foreach (qNode n in QuadNodes)
            {
                bool isFrontNode = IsFrontNode(n, frontEdges);
                if (isFrontNode & !n.BoundaryNode)
                {
                    // smooth front node
                    smoothNode = FrontNodeSmooth(n, globalEdgeList);
                }
                else if (!isFrontNode & !n.BoundaryNode)
                {
                    // do laplacian smooth
                    smoothNode = ModifiedLengthWeightedLaplacianSmooth(n, globalEdgeList);
                }
                else
                {
                    smoothNode = n;
                }
            }

            adjacentNodes = GetNeighborNodesToElement(QuadElement, globalEdgeList);

            foreach (qNode n in adjacentNodes)
            {
                connectedEdges = GetConnectedEdges(n, globalEdgeList);
                bool isFrontNode = IsFrontNode(n, frontEdges);
                if (isFrontNode & !n.BoundaryNode)
                {
                    // smooth front node
                    smoothNode = FrontNodeSmooth(n, connectedEdges);
                }
                else if (!isFrontNode & !n.BoundaryNode)
                {
                    // do laplacian smooth
                    smoothNode = ModifiedLengthWeightedLaplacianSmooth(n, globalEdgeList);
                }
                else
                {
                    smoothNode = n;
                }
            }



        }

        private List<qNode> GetNeighborNodesToElement(qElement element, List<qEdge> globalEdgeList) // todo: check if this is OK.
        {
            List<qNode> adjacentNodes = new List<qNode>();

            List<qNode> elementNodes = GetNodesOfElement(element);
            List<qEdge> elementEdges = element.EdgeList;

            for (int i = 0; i < elementNodes.Count; i++)
            {
                qNode node = elementNodes[i];
                List<qEdge> connectedEdges = GetConnectedEdges(node, globalEdgeList);

                for (int j = 0; j < connectedEdges.Count; j++)
                {
                    bool IsElementEdge = false;
                    qEdge edge = connectedEdges[j];
                    for (int k = 0; k < elementEdges.Count; k++)
                    {
                        if (edge == elementEdges[k]) { IsElementEdge = true; }
                    }

                    if (!IsElementEdge)
                    {
                        qNode adjacentNode = GetOppositeNode(node, edge);
                        adjacentNodes.Add(adjacentNode);
                    }
                }
            }
            return adjacentNodes;
        }
        private bool IsFrontNode(qNode node, List<qEdge> frontEdges) // todo: test if this works. Check each node of quad.
        {
            bool isFrontNode = false;
            foreach (qEdge front in frontEdges)
            {
                if (node == front.StartNode || node == front.EndNode)
                {
                    isFrontNode = true;
                }
            }
            return isFrontNode;
        }
        private qNode ModifiedLengthWeightedLaplacianSmooth(qNode Ni, List<qEdge> globalEdgeList) // todo: check if this works. Check if each method osv do what I want
        {
            Vector3d lengthCjVectorCj = Vector3d.Zero;
            qNode smoothNode = new qNode();

            List<qEdge> connectedEdgesToNi = GetConnectedEdges(Ni, globalEdgeList); // todo: check if it finds connected edges
            double lengthCj = 0;

            for (int i = 0; i < connectedEdgesToNi.Count; i++)
            {
                qEdge edge = connectedEdgesToNi[i];
                qNode Nj = GetOppositeNode(Ni, edge);
                Vector3d Vi = Nj.Coordinate - Ni.Coordinate;
                Vector3d Cj = Vector3d.Zero;

                if (Nj.BoundaryNode)
                {
                    //Vector3d deltaCj = GetAngularSmoothness(Nj); //todo: fix angular smoothness
                    Cj = Vi; //+ deltaCj;
                }
                else
                {
                    Cj = Vi;
                }

                lengthCj = Cj.Length + lengthCj;
                lengthCjVectorCj = Cj.Length * Cj;
            }

            Vector3d delta = lengthCjVectorCj / lengthCj;
            smoothNode.Coordinate = new Point3d(Ni.Coordinate.X + delta.X, Ni.Coordinate.Y + delta.Y, Ni.Coordinate.Z + delta.Z);

            return smoothNode;
        }
        private qEdge GetSharedEdge(List<qElement> quadElements)
        {
            qEdge sharedEdge = new qEdge();
            qElement quadElement1 = quadElements[0];
            qElement quadElement2 = quadElements[0];

            foreach (qEdge edge in quadElement1.EdgeList)
            {
                if (quadElement2.EdgeList.Contains(edge))
                {
                    sharedEdge = edge;
                }
            }
            return sharedEdge;
        }
        private qNode FrontNodeSmooth(qNode Ni, List<qEdge> globalEdgeList)
        {
            // Find front edges connected to front node: todo: test if this works as I want.
            List<qEdge> nodeFrontEdges = new List<qEdge>();
            List<qEdge> connectedEdges = GetConnectedEdges(Ni, globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                bool isFrontEdge = IsFrontEdge(edge);
                if (isFrontEdge)
                {
                    nodeFrontEdges.Add(edge);
                }
            }

            if (nodeFrontEdges.Count != 2) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, String.Format("Node is connected to {0} front edges", nodeFrontEdges.Count)); }
            else
            {
                qEdge front1 = nodeFrontEdges[0];
                qEdge front2 = nodeFrontEdges[1];

                if (front1.RightFrontNeighbor.StartNode == Ni || front1.RightFrontNeighbor.EndNode == Ni)
                {
                    qEdge nodeLeftFront = front1;
                    qEdge nodeRightFront = front2;
                }
                else
                {
                    qEdge nodeLeftFront = front2;
                    qEdge nodeRightFront = front1;
                }
            }


            // Get quads that are connected to node. Todo: check if this do what I want.
            List<qElement> quadElements = GetQuadsConnectedToNode(Ni, globalEdgeList);
            int numberOfConnectedQuads = quadElements.Count;


            // True isoparametric smooth:
            if (numberOfConnectedQuads != 0)
            {
                Vector3d Vmj = Vector3d.Zero;
                Vector3d Vmk = Vector3d.Zero;
                Vector3d Vml = Vector3d.Zero;
                Vector3d vectorSum = Vector3d.Zero;
                Point3d origo = new Point3d(0, 0, 0);

                foreach (qElement quadElement in quadElements)
                {
                    List<qNode> quadNodes = GetNodesOfElement(quadElement);
                    qNode node1 = quadNodes[0];
                    qNode node2 = quadNodes[1];
                    qNode node3 = quadNodes[2];
                    qNode node4 = quadNodes[3];

                    if (node1 == Ni)
                    {
                        Vmj = node2.Coordinate - origo;
                        Vmk = node3.Coordinate - origo;
                        Vml = node4.Coordinate - origo;
                    }
                    if (node2 == Ni)
                    {
                        Vmj = node3.Coordinate - origo;
                        Vmk = node4.Coordinate - origo;
                        Vml = node1.Coordinate - origo;
                    }
                    if (node3 == Ni)
                    {
                        Vmj = node4.Coordinate - origo;
                        Vmk = node1.Coordinate - origo;
                        Vml = node2.Coordinate - origo;
                    }
                    if (node4 == Ni)
                    {
                        Vmj = node1.Coordinate - origo;
                        Vmk = node2.Coordinate - origo;
                        Vml = node3.Coordinate - origo;
                    }

                    vectorSum = (Vmj + Vml - Vmk) + vectorSum;
                }

                Vector3d Vi = Ni.Coordinate - origo;
                Vector3d Vi_mark = (1 / numberOfConnectedQuads) * vectorSum;

                Vector3d deltaA = Vi_mark - Vi;


                // If numberOfConnectedQuads = 2, do utvidet smooth 
                Vector3d deltaI = Vector3d.Zero;
                Vector3d deltaB = Vector3d.Zero;
                Vector3d deltaC = Vector3d.Zero;
                if (numberOfConnectedQuads == 2)
                {
                    // Length modification:
                    qEdge sharedEdge = null;
                    sharedEdge = GetSharedEdge(quadElements);

                    qNode Nj = GetOppositeNode(Ni, sharedEdge);
                    Vector3d Vj = Nj.Coordinate - origo;

                    double la = (Vi_mark - Vj).Length;

                    // Find ld
                    double ld = 0;
                    double tr = 2.5; //todo: fix tr!!
                    double sumSurroundingEdges = 0;
                    double n = 0;
                    if (tr <= 2.5)
                    {
                        ld = (Ni.Coordinate - Nj.Coordinate).Length; // todo: check assumption
                    }
                    else if (tr > 2.5 & tr <= 20)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            qEdge edge1 = quadElements[0].EdgeList[i];
                            qEdge edge2 = quadElements[1].EdgeList[i];

                            // length of edges of element1
                            if (!IsFrontEdge(edge1) || edge1 != sharedEdge)
                            {
                                sumSurroundingEdges = edge1.Length + sumSurroundingEdges;
                            }

                            // length of edges of element2
                            if (!IsFrontEdge(edge2) || edge2 != sharedEdge)
                            {
                                sumSurroundingEdges = edge2.Length + sumSurroundingEdges;
                            }
                        }

                        // length of edges not connected to elements or is frontedge
                        foreach (qEdge edge in connectedEdges)
                        {
                            if (!IsFrontEdge(edge) || edge != sharedEdge)
                            {
                                sumSurroundingEdges = edge.Length + sumSurroundingEdges;
                                n++;
                            }
                        }

                        ld = sumSurroundingEdges / (4 + n);
                    }
                    else
                    {
                        foreach (qEdge edge in connectedEdges)
                        {
                            sumSurroundingEdges = edge.Length + sumSurroundingEdges;
                            n++;
                        }

                        ld = sumSurroundingEdges / n;
                    }

                    if (IsFrontEdge(sharedEdge)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Shared edge in Front Node smooth (two quads) is a front edge.."); }

                    deltaB = Vj - Vi + (deltaA + Vi - Vj) * ld / la;


                    // Angle modification:




                }
                else { deltaI = deltaA; }


            }
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "numberOfConnectedQuads is zero. Blacker smooth can not be performed!!"); }
            return null;
        }


        private Vector3d GetAngularSmoothness(qNode Ni, qNode Nj, List<qEdge> globalEdgeList) // todo: check if this works
        {
            List<qElement> quadElements = GetQuadsConnectedToNode(Ni, globalEdgeList);
            List<qEdge> connectedEdges = GetConnectedEdges(Ni, globalEdgeList);
            
            //Hvis node til quadElements er på BC - gjør slik.
            if (!IsFrontNode(Ni, globalEdgeList))
            {

            }

            if (quadElements.Count == 2)
            {
                qEdge sharedEdge = GetSharedEdge(quadElements);

                qElement quadElement1 = quadElements[0];
                qEdge topEdge1 = quadElement1.EdgeList[3];
                List<qNode> nodes1 = GetNodesOfElement(quadElement1);

                qElement quadElement2 = quadElements[1];
                qEdge topEdge2 = quadElement2.EdgeList[3];

                Vector3d vectorLeft = Vector3d.Zero; // dummy vector
                Vector3d vectorRight = Vector3d.Zero; // dummy vector

                if (GetOppositeNode(Ni, topEdge1) == nodes1[4])
                {
                    vectorLeft = GetOppositeNode(Ni, topEdge1).Coordinate - Nj.Coordinate;
                    vectorRight = GetOppositeNode(Ni, topEdge2).Coordinate - Nj.Coordinate;
                }
                else
                {
                    vectorLeft = GetOppositeNode(Ni, topEdge2).Coordinate - Nj.Coordinate;
                    vectorRight = GetOppositeNode(Ni, topEdge1).Coordinate - Nj.Coordinate;
                }

                Vector3d P_B1 = GetBisectingVector(vectorRight, vectorLeft);

                Vector3d Pi = Ni.Coordinate - Nj.Coordinate;
                Vector3d P_B2 = P_B1.Length * Pi + Pi.Length * P_B1; // Assume angle always less than 180 degree.
                if (Vector3d.Multiply(P_B1, Pi) < 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "P_B1 and Pi have angle larger than 90 degree. Do not know it is above 180."); }

                Line line1 = new Line(Nj.Coordinate, P_B2, 100);
                Line line2 = new Line(GetOppositeNode(Ni, topEdge1).Coordinate, GetOppositeNode(Ni, topEdge2).Coordinate);


            }


            return Vector3d.Zero;
        }

        private Vector3d GetBisectingVector(Vector3d VectorRight, Vector3d VectorLeft)
        {
            double angle = Vector3d.VectorAngle(VectorRight, VectorLeft, Vector3d.ZAxis); // todo: make normal mor general
            Vector3d bisectVector = VectorRight;
            bisectVector.Rotate(0.5 * angle, Vector3d.ZAxis); // todo: not for 3d surface, make axis for normal to plane (vec1, vec2)
            bisectVector.Unitize();
            bisectVector = bisectVector * (VectorRight.Length * VectorLeft + VectorLeft.Length * VectorRight).Length;
            return bisectVector;
        }
        private List<qEdge> GetFrontEdgesConnectedToNode(qNode node, List<qEdge> globalEdgeList)
        {
            List<qEdge> nodeFrontEdges = new List<qEdge>();
            qEdge nodeLeftFront = new qEdge();
            qEdge nodeRightFront = new qEdge();

            List<qEdge> connectedEdges = GetConnectedEdges(node, globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                bool isFrontEdge = IsFrontEdge(edge);
                if (isFrontEdge)
                {
                    nodeFrontEdges.Add(edge);
                }
            }

            if (nodeFrontEdges.Count != 2) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, String.Format("Node is connected to {0} front edges", nodeFrontEdges.Count)); }
            else
            {
                qEdge front1 = nodeFrontEdges[0];
                qEdge front2 = nodeFrontEdges[1];

                if (front1.RightFrontNeighbor.StartNode == node || front1.RightFrontNeighbor.EndNode == node)
                {
                    nodeLeftFront = front1;
                    nodeRightFront = front2;
                }
                else
                {
                    nodeLeftFront = front2;
                    nodeRightFront = front1;
                }
            }
            nodeFrontEdges.Clear();
            nodeFrontEdges.Add(nodeLeftFront);
            nodeFrontEdges.Add(nodeRightFront);

            return nodeFrontEdges;

        }

        private List<qElement> GetQuadsConnectedToNode(qNode node, List<qEdge> globalEdgeList)
        {
            // Get how quads that are connected to node. Todo: check if this do what I want.
            List<qElement> quadElements = new List<qElement>();

            var connectedEdges = GetConnectedEdges(node, globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                if (edge.Element1.IsQuad) { quadElements.Add(edge.Element1); }
                if (edge.Element2.IsQuad) { quadElements.Add(edge.Element2); }
            }

            if (quadElements.Count > 0)
            {
                foreach (qElement element in quadElements)
                {
                    if (quadElements.Contains(element)) { quadElements.Remove(element); }
                }
            }
            if (quadElements.Count == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "number of connected quad-elements to frontNode is zero."); }

            return quadElements;
        }

        private List<qNode> GetQuadNodesOfElement(qElement element)
        {
            List<qEdge> quadEdges = element.EdgeList;
            
            qEdge baseEdge = quadEdges[0];
            qEdge rightEdge = quadEdges[1];
            qEdge leftEdge = quadEdges[2];
            qEdge topEdge = quadEdges[3];

            qNode node1 = new qNode();
            qNode node2 = new qNode();
            qNode node3 = new qNode();
            qNode node4 = new qNode();

            if (baseEdge.StartNode == leftEdge.StartNode || baseEdge.StartNode == leftEdge.EndNode)
            {
                node1 = baseEdge.StartNode;
                node2 = GetOppositeNode(node1, baseEdge);
            }
            else if (baseEdge.EndNode == leftEdge.StartNode || baseEdge.EndNode == leftEdge.EndNode)
            {
                node1 = baseEdge.EndNode;
                node2 = GetOppositeNode(node1, baseEdge);
            }

            if (topEdge.StartNode == leftEdge.StartNode || topEdge.StartNode == leftEdge.EndNode)
            {
                node3 = baseEdge.StartNode;
                node4 = GetOppositeNode(node3, topEdge);
            }
            else if (topEdge.EndNode == leftEdge.StartNode || topEdge.EndNode == leftEdge.EndNode)
            {
                node3 = baseEdge.EndNode;
                node4 = GetOppositeNode(node3, topEdge);
            }

            List<qNode> quadNodes = new List<qNode> { node1, node2, node3, node4 };

            return quadNodes;
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
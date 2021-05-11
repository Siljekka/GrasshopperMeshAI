using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using Rhino.Geometry.Collections;
using System.Linq;
using Rhino.Geometry.Intersect;
using System.Drawing;

namespace MeshPoints.QuadRemesh
{
    public class QuadRemesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the AnalyzeTriangleMesh class.
        /// </summary>
        public QuadRemesh()
          : base("Quad Remesh", "atm",
              "Remesh from triangle mesh to quadrilateral mesh.",
              "SmartMesh", "QuadRemesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Triangle mesh", "trimesh", "Input a trinagle mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("# elements to remesh", "number", "Input an integer", GH_ParamAccess.item, 1);
            pManager.AddGenericParameter("Local smoothing", "tI", "", GH_ParamAccess.item);
            pManager.AddNumberParameter("number iteration before stop", "temp", "Temoprary: number of iteration to perform", GH_ParamAccess.item, 2);
            pManager[0].Optional = true;
            pManager[1].Optional = true;
            pManager[2].Optional = true;
            pManager[3].Optional = true;

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Front edges", "element", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Edge list", "11", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("Element list", "10", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("qElements", "element", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("E_front", "element", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("E_k_left", "element", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("E_k_right ", "element", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("avg quality ", "q", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("badest quality ", "bq", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("color mesh", "bq", "m", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Mesh mesh = new Mesh();
            double numberElementsToRemesh = 0;
            double iterationsToPerformBeforeStop = 0;
            bool performeLocalSmoothing = false;
            
            DA.GetData(0, ref mesh);
            DA.GetData(1, ref numberElementsToRemesh);
            DA.GetData(2, ref performeLocalSmoothing);
            DA.GetData(3, ref iterationsToPerformBeforeStop);

            if (!DA.GetData(0, ref mesh)) return; 
            if (!DA.GetData(1, ref numberElementsToRemesh)) return; 
            if (!DA.GetData(2, ref performeLocalSmoothing)) return; 
            if (!DA.GetData(3, ref iterationsToPerformBeforeStop)) return;

            #region Code


            // Get initial edges and elements using mesh topology properties
            var initialEdgeAndElementList = GetInitialEdgesAndElements(mesh);
            List<qEdge> globalEdgeList = initialEdgeAndElementList.Item1;
            List<qElement> globalElementList = initialEdgeAndElementList.Item2;

            // Get intial front edges
            List<qEdge> frontEdges = GetFrontEdges(globalEdgeList);

            // check if even number of boundary nodes
            if (!IsFrontLoopsEven(frontEdges, null, globalEdgeList).Item1) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Need the initial mesh to have an even number of boundary nodes to make it an all quad mesh"); return; }
            if (numberElementsToRemesh == 0 | iterationsToPerformBeforeStop == 0) return;

            // temporary variables
            qEdge E_k_left = new qEdge(); // left side edge of quad
            qEdge E_k_right = new qEdge(); // right side edge of quad
            qEdge E_front = new qEdge(); // bottom of quad
            qEdge E_top = new qEdge(); // top if quad
            int iterationCounter = 0;
            qElement quadElement = new qElement();
            qEdge unselectedEdge = new qEdge();

            for (int n = 0; n < numberElementsToRemesh; n++)
            {
                // Temporary stop
                if (iterationCounter == iterationsToPerformBeforeStop) { AddRuntimeMessage(GH_RuntimeMessageLevel.Blank, "Iteration stop"); break; }
                iterationCounter++;
                ClearRuntimeMessages();


                // Logical values for succesful selection of edges
                bool E_k_left_performed = true;
                bool E_k_right_performed = true;
                bool E_top_performed = true;

                // Back up if selected E_front is not to selectable
                List<qElement> globalElementListBackUp = globalElementList;
                List<qEdge> globalEdgeListBackUp = globalEdgeList;

                // Stop if quad remesh is finished
                if (frontEdges.Count == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Quad Remesh succeeded!");
                    break;
                }

                //________________ select next front edge________________
                var E_frontAndEdgeState = SelectNextFrontEdge(frontEdges);
                E_front = E_frontAndEdgeState.Item1;
                var edgeState = E_frontAndEdgeState.Item2;
                unselectedEdge.Unselectable = false;

                bool loopOf4Edges = LoopControll(frontEdges, globalEdgeList, globalElementList); // to do: legg inn krit når kun 8 edges igjen
                CleanUpChevorns(globalEdgeList, globalElementList, frontEdges);

                if (loopOf4Edges)
                {
                    continue;
                }
                /*
                 if (loop.Count == 8)
                {
                    foreach (qEdge edge in loop)
                    {
                        qElement frontElement = edge.GetFrontElement();
                        foreach (qEdge elementEdge in frontElement.EdgeList)
                        {
                            if (IsFrontLoopsEven(frontEdges, elementEdge, globalEdgeList).Item1)
                            { 
                                E_front = 
                            }
                        }
                    }              
                }
                 
                 */
                //________________ check special case________________
                var specialCaseValues = CheckSpecialCase(E_front, globalEdgeList, globalElementList, frontEdges);
                bool seamAnglePerformed = specialCaseValues.Item1;
                bool isSpecialCase = specialCaseValues.Item2;

                if (seamAnglePerformed)
                {
                    continue;
                }

                //________________ get side edges ________________
                if (isSpecialCase) // if special case
                {
                    E_front = specialCaseValues.Item3;
                    E_k_right = specialCaseValues.Item4;
                    E_k_left = specialCaseValues.Item5;
                }
                else if  (!isSpecialCase) // if not special case
                {
                    // get left edge
                    switch (edgeState[0])
                    {
                        case 0:

                            var leftSideEdgeValues = GetSideEdge(globalElementList, globalEdgeList, 0, E_front, frontEdges, false);
                            E_k_left = leftSideEdgeValues.Item1;
                            E_k_left_performed = leftSideEdgeValues.Item2;  // false only when split/swap is performed on an closing edge, i.e. E_0 is a front edge, and swap edge fails
                            break;
                        case 1:
                            qNode sharedNode = E_front.GetSharedNode(E_front.LeftFrontNeighbor);
                            if (IsEdgeClosingFront(sharedNode, E_front.LeftFrontNeighbor, frontEdges, globalEdgeList, E_front))
                            {
                                E_k_left = CloseFront(frontEdges, globalEdgeList, globalElementList, E_front.LeftFrontNeighbor, sharedNode);
                                break;
                            }
                            else
                            {
                                E_k_left = E_front.LeftFrontNeighbor; 
                                break;
                            }
                    }

                    // get right edge
                    switch (edgeState[1])
                    {
                        case 0:
                            var rightSideEdgeValues = GetSideEdge(globalElementList, globalEdgeList, 1, E_front, frontEdges, false);
                            E_k_right = rightSideEdgeValues.Item1;
                            E_k_right_performed = rightSideEdgeValues.Item2; // false only when split/swap is performed on an closing edge, i.e. E_0 is a front edge, and swap edge fails
                            break;
                        case 1:
                            qNode sharedNode = E_front.GetSharedNode(E_front.RightFrontNeighbor);
                            if (IsEdgeClosingFront(sharedNode, E_front.RightFrontNeighbor, frontEdges, globalEdgeList, E_front))
                            {
                                E_k_right = CloseFront(frontEdges, globalEdgeList, globalElementList, E_front.RightFrontNeighbor, sharedNode);
                                break;
                            }
                            else
                            {
                                E_k_right = E_front.RightFrontNeighbor;
                                break;
                            }
                    }
                }

                //________________get top edge________________
                if (E_k_left_performed & E_k_right_performed)
                {
                    var topEdgeValue = GetTopEdge(E_front, E_k_left, E_k_right, globalEdgeList, globalElementList, frontEdges);
                    E_top = topEdgeValue.Item1;
                    E_k_left = topEdgeValue.Item2; // in case triangle
                    E_k_right = topEdgeValue.Item3; // in case traingle
                    E_top_performed = topEdgeValue.Item4;
                }
                else { E_top_performed = false; }

                // If not possible to perform top recovery, skip the selected E_front and select a new front edge
                if (!E_top_performed)
                {
                    globalElementList = globalElementListBackUp; // reset changes made in the iteration
                    globalEdgeList = globalEdgeListBackUp; // reset changes made in the iteration
                    E_front.Unselectable = true;
                    unselectedEdge = E_front;
                    n--;
                    continue;
                }


                //________________ quadrilateral formation________________
                List<qEdge> quadEdges = new List<qEdge>() { E_front, E_k_right, E_k_left, E_top };
                quadElement = CreateQuadElement(quadEdges, globalEdgeList, globalElementList, frontEdges);


                // ________________Local smoothing________________
                if (performeLocalSmoothing)
                { DoLocalSmoothing(quadElement, globalEdgeList, frontEdges, globalElementList); }

                // to do: fix....
                if (!IsFrontLoopsEven(frontEdges, null, globalEdgeList).Item1)
                {
                    globalElementList = globalElementListBackUp; // reset changes made in the iteration
                    globalEdgeList = globalEdgeListBackUp; // reset changes made in the iteration
                    E_front.Unselectable = true;
                    unselectedEdge = E_front;
                    n--;
                    continue;
                }

                // to do: what if closing front og special case?
                // to do: apply local smoothing for seamAngle
                // to do: temporay solution for E_frontFail
            }

            DoGlobalSmoothing(globalEdgeList, globalElementList);
            

            // ___________Transform to main mesh classes_______________
            var meshProperties = ConvertToMainMeshClasses(globalElementList);
            List<Node> nodes = meshProperties.Item1;
            List<Element> elements = meshProperties.Item2;



            List<qEdge> test = new List<qEdge>() { E_front, E_k_right, E_k_left, E_top };
            var meshValues = CalculateQuality(globalElementList);
            double avgQuality = meshValues.Item1;
            double badestQuality = meshValues.Item2;
            Mesh colorMesh = meshValues.Item3;

            // Assign properties to surfaceMesh:
            SmartMesh surfaceMesh = new SmartMesh(nodes, elements, colorMesh, "Surface");

            #endregion End Code



           
            DA.SetDataList(0, frontEdges);
            DA.SetDataList(1, globalEdgeList);
            DA.SetDataList(2, globalElementList);
            DA.SetData(3, quadElement);
            DA.SetData(4, E_front);
            //DA.SetData(5, vec);
            /*DA.SetDataList(6, nodesTest);
            DA.SetDataList(7, v);
            DA.SetDataList(8, p);
            DA.SetData(9, e);*/

        }

        #region Methods
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
            for (int  i = 0; i < globalqElementList.Count; i++)
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

        // _____________________________________ for ininital mesh _________________________________________
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

        // ________________________________ for front definiton and classification ___________________________________
        private List<qEdge> GetFrontEdges(List<qEdge> globalEdgeList)
        {
            // summary: get all edges connected to a single triangle element
            List<qEdge> frontEdges = new List<qEdge>();
            foreach (qEdge edge in globalEdgeList)
            {
                List<qElement> connectedElements = edge.GetConnectedElements();
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
                
                var edgesToOrientate = edge.OrientateNeigborEdges(neigborEdgeToStartNode, neigborEdgeToEndNode);
                edge.LeftFrontNeighbor = edgesToOrientate.Item1;
                edge.RightFrontNeighbor = edgesToOrientate.Item2;
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

            double angleTolerance = 0.75 * Math.PI; // constant

            for (int i = 0; i < frontEdges.Count; i++)
            {
                double leftAngle = frontEdges[i].CalculateAngleOfNeighborFrontEdges("left");
                double rightAngle = frontEdges[i].CalculateAngleOfNeighborFrontEdges("right");

                if (leftAngle < angleTolerance & rightAngle < angleTolerance) { list11.Add(frontEdges[i]); }
                else if (leftAngle >= angleTolerance & rightAngle < angleTolerance) { list01.Add(frontEdges[i]); }
                else if (leftAngle < angleTolerance & rightAngle >= angleTolerance) { list10.Add(frontEdges[i]); }
                else if (leftAngle >= angleTolerance & rightAngle >= angleTolerance) { list00.Add(frontEdges[i]); }
            }
            return Tuple.Create(list11, list10, list01, list00);
        }
        private int[] GetEdgeState(qEdge E_front, List<qEdge> frontEdges)
        {
            int[] edgeState = { 0, 0 };

            var edgeStates = CreateEdgeStateList(frontEdges);
            var list11 = edgeStates.Item1;
            var list10 = edgeStates.Item2;
            var list01 = edgeStates.Item3;
            var list00 = edgeStates.Item4;

            if (list11.Contains(E_front))
            {edgeState[0] = 1; edgeState[1] = 1; }
            else if (list01.Contains(E_front))
            { edgeState[0] = 0; edgeState[1] = 1; }
            else if (list10.Contains(E_front))
            { edgeState[0] = 1; edgeState[1] = 0; }
            else if (list00.Contains(E_front))
            { edgeState[0] = 0; edgeState[1] = 0; }

            return edgeState;
        }
        private Tuple<qEdge, int[]> SelectNextFrontEdge(List<qEdge> frontEdges)
        {
            // summary: select next front edge with the hierarchy: list -> lowest level -> shortest length and switch to neighbor if large transision
            qConstants constants = new qConstants();
            double transitionTolerance = constants.GetTransitionTolerance();
            qEdge E_front = new qEdge();

            var edgeStates = CreateEdgeStateList(frontEdges);
            var list11 = edgeStates.Item1;
            var list10 = edgeStates.Item2;
            var list01 = edgeStates.Item3;
            var list00 = edgeStates.Item4;

            // to do: make more efficient
            // get lowest global level
            int lowestGlobalLevel = 1000000;
            foreach (qEdge edge in frontEdges)
            { if (edge.Level < lowestGlobalLevel & !edge.Unselectable) { lowestGlobalLevel = edge.Level; } }

            // get lowest list level
            int lowestLevelList11 = 100000;
            foreach (qEdge edge in list11)
            { if (edge.Level < lowestLevelList11 & !edge.Unselectable) { lowestLevelList11 = edge.Level; } }

            int lowestLevelList01 = 100000;
            foreach (qEdge edge in list01)
            { if (edge.Level < lowestLevelList01 & !edge.Unselectable) { lowestLevelList01 = edge.Level; } }

            int lowestLevelList10 = 100000;
            foreach (qEdge edge in list10)
            { if (edge.Level < lowestLevelList10 & !edge.Unselectable) { lowestLevelList10 = edge.Level; } }

            int lowestLevelList00 = 100000;
            foreach (qEdge edge in list00)
            { if (edge.Level < lowestLevelList00 & !edge.Unselectable) { lowestLevelList00 = edge.Level; } }


            // select
            if (list11.Count != 0 & lowestLevelList11 == lowestGlobalLevel)
            {
                 E_front = SelectFrontEdgeFromList(list11, lowestGlobalLevel);
            }
            else if ((list01.Count != 0 | list10.Count != 0) & (lowestLevelList01 == lowestGlobalLevel | lowestLevelList10 == lowestGlobalLevel))
            {
                // merge list01 and 10
                List<qEdge> list0110 = new List<qEdge>(list01);
                list0110.AddRange(list10);
                E_front = SelectFrontEdgeFromList(list0110, lowestGlobalLevel);
            }
            else if (list00.Count != 0)
            {
                E_front = SelectFrontEdgeFromList(list00, lowestGlobalLevel);
            }
            else
            {
                // to do: temporary end
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "SelectNextFrontEdge: No more edges to select.");
                E_front = null;
                int[] tempList = { 0, 0 }; 
                return Tuple.Create(E_front, tempList);
            }

            // switch E_front to neighbor front if transision to one of the front neighbors is large and neighbor is selectable
            
            if (E_front.Length / E_front.LeftFrontNeighbor.Length > transitionTolerance & !E_front.LeftFrontNeighbor.Unselectable)
            {
                if (!list11.Contains(E_front))
                {
                    E_front = E_front.LeftFrontNeighbor;
                    { AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "SelectNextFrontEdge: Large transision to left neighor. Switched."); }
                }
            }
            else if (E_front.Length / E_front.RightFrontNeighbor.Length > transitionTolerance & !E_front.RightFrontNeighbor.Unselectable)
            {
                if (!list11.Contains(E_front))
                {
                    E_front = E_front.RightFrontNeighbor;
                    { AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "SelectNextFrontEdge: Large transision to right neighor. Switched."); }
                }
            }
            

            int[] edgeState = GetEdgeState(E_front, frontEdges);

            return Tuple.Create(E_front, edgeState);
        }
        private qEdge SelectFrontEdgeFromList(List<qEdge> edgeStateList, int lowestGlobalLevel)
        {
            // summary: select an edge from edge state list input. Choose lowest level, if multiple choose shortest of them.

            // ensure side edge is not selected
            List<qEdge> selectableEdges = new List<qEdge>(edgeStateList);

            // get edges of lowest level
            List<qEdge> edgesOfLowestLevel = new List<qEdge>();
            foreach (qEdge edge in selectableEdges)
            {
                if (edge.Unselectable) { continue; } // remove unselectabel edges
                if (edge.Level == lowestGlobalLevel) { edgesOfLowestLevel.Add(edge); }
            }

            // get shortes edge
            qEdge E_front = new qEdge();
            double shortestLength = edgesOfLowestLevel[0].Length;
            foreach (qEdge edge in edgesOfLowestLevel)
            {
                if (edge.Length <= shortestLength) { E_front = edge; shortestLength = edge.Length; }
            }
            return E_front;
        }

        // _________________________________________ for topology  _______________________________________________
        private qNode GetNodeOfFrontEdge(int nodeToEvaluate, qEdge edge) // to do: flytt til edge
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
        private List<qNode> GetSwapedNodes(qEdge E_0)
        {
            // summary: get the new nodes of E_0 if swaped. Assume E_0 is swapable, i.e. two adjacent triangles.

            List<qNode> swapedNodes = new List<qNode>(); // list of new nodes if swaping
            List<qNode> nodeCandidates = new List<qNode>(); // list of all adjacent element nodes
            List<qElement> adjecentElements = new List<qElement>() { E_0.Element1, E_0.Element2 };

            // find candidate nodes
            foreach (qElement element in adjecentElements)
            {
                List<qNode> elementNodes = element.GetNodesOfElement();
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

        private Tuple<bool, bool, qEdge, qEdge, qEdge> CheckSpecialCase(qEdge E_front, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: check if selected front edge is a special case and performe the needed operations

            double epsilon_1 = 0.04 * Math.PI; // constant: can be changed
            double epsilon_2 = 0.09 * Math.PI; // constant: can be changed

            bool seamAnglePerformed = false;
            bool specialCase = false;
            qEdge E_k_left = new qEdge();
            qEdge E_k_right = new qEdge();

            // to do: check if correct to say
            if (E_front.StartNode.BoundaryNode | E_front.EndNode.BoundaryNode) { return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_left, E_k_right); }

            #region Check left side

            qNode leftNode = GetNodeOfFrontEdge(0, E_front);
            double leftAngle = E_front.CalculateAngleOfNeighborFrontEdges("left");
            int numberConnectedQuad = leftNode.GetQuadsConnectedToNode(globalEdgeList).Count;

            // ensure length ratio is > 1
            double lengthRatio = 0;
            if (E_front.Length < E_front.LeftFrontNeighbor.Length) { lengthRatio = E_front.LeftFrontNeighbor.Length / E_front.Length; }
            else { lengthRatio = E_front.Length / E_front.LeftFrontNeighbor.Length; }


            if ((leftAngle < epsilon_2 & numberConnectedQuad <= 4) | ((leftAngle < epsilon_1) & numberConnectedQuad > 4))
            {
                if (lengthRatio > 2.5)
                {
                    // perform seam transition
                    var edgesOfQuad = TransitionSeam(E_front, E_front.LeftFrontNeighbor, globalEdgeList, globalElementList, frontEdges);
                    E_front = edgesOfQuad.Item1;
                    E_k_right = edgesOfQuad.Item2;
                    E_k_left = edgesOfQuad.Item3;

                    if (E_k_left != null) 
                    { 
                        specialCase = true;
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "CheckSpecialCase: performed a transition seam of left edge.");
                        return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_right, E_k_left);
                    }
                    else { specialCase = false; }
                }
                else
                {
                    // perform angle seam
                    Seam(E_front.LeftFrontNeighbor, E_front, globalEdgeList, globalElementList, frontEdges);
                    specialCase = true;
                    seamAnglePerformed = true;
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "CheckSpecialCase: performed a seam of left edge.");
                    return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_right, E_k_left);
                }

            }
            else if (lengthRatio > 2.5)
            {
                // perform transition split
                var edgesOfQuad = TransitionSplit(E_front, E_front.LeftFrontNeighbor, globalEdgeList, globalElementList, frontEdges);
                E_front = edgesOfQuad.Item1;
                E_k_right = edgesOfQuad.Item2;
                E_k_left = edgesOfQuad.Item3;

                if (E_k_left != null)
                { 
                    specialCase = true;
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "CheckSpecialCase: performed a split of left edge.");
                    return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_right, E_k_left);
                }
                else { specialCase = false; }
            }
            #endregion check left side


            #region Check right side

            qNode rightNode = GetNodeOfFrontEdge(1, E_front);
            double rightAngle = E_front.CalculateAngleOfNeighborFrontEdges("right");
            numberConnectedQuad = rightNode.GetQuadsConnectedToNode(globalEdgeList).Count;

            // ensure length ratio is > 1
            lengthRatio = 0;
            if (E_front.Length < E_front.RightFrontNeighbor.Length) { lengthRatio = E_front.RightFrontNeighbor.Length / E_front.Length; }
            else { lengthRatio = E_front.Length / E_front.RightFrontNeighbor.Length; }

            if ((rightAngle < epsilon_2 & numberConnectedQuad <= 4) | ((rightAngle < epsilon_1) & numberConnectedQuad > 4))
            {
                if (lengthRatio > 2.5)
                {
                    // perform seam transition
                    var edgesOfQuad = TransitionSeam(E_front, E_front.RightFrontNeighbor, globalEdgeList, globalElementList, frontEdges);
                    E_front = edgesOfQuad.Item1;
                    E_k_right = edgesOfQuad.Item2;
                    E_k_left = edgesOfQuad.Item3;

                    if (E_k_right != null)
                    {
                        specialCase = true;
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "CheckSpecialCase: performed a transition seam of  right edge.");
                        return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_right, E_k_left);
                    }
                    else { specialCase = false; }
                }
                else
                {
                    // perform angle seam
                    Seam(E_front, E_front.RightFrontNeighbor, globalEdgeList, globalElementList, frontEdges);
                    specialCase = true;
                    seamAnglePerformed = true;
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "CheckSpecialCase: performed a seam of  right edge.");
                    return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_right, E_k_left);
                }

            }
            else if (lengthRatio > 2.5)
            {
                // perform split
                var edgesOfQuad = TransitionSplit(E_front, E_front.RightFrontNeighbor, globalEdgeList, globalElementList, frontEdges);
                E_front = edgesOfQuad.Item1;
                E_k_right = edgesOfQuad.Item2;
                E_k_left = edgesOfQuad.Item3;

                if (E_k_right != null)
                {
                    specialCase = true;
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "CheckSpecialCase: performed a split of right edge.");
                    return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_right, E_k_left);
                }
                else { specialCase = false; }
            }

            #endregion Check right side

            return Tuple.Create(seamAnglePerformed, specialCase, E_front, E_k_right, E_k_left);

        }
        
        private Tuple<bool, List<List<qEdge>>> IsFrontLoopsEven(List<qEdge> frontEdges, qEdge checkSideEdge, List<qEdge> globalEdgeList)
        {
            // summary: check if front loops are comprised of an even number of edges. If checkSideEdge != null, check if new loops to be formed are even loops.
            bool evenEdgesInLoops = true;
            List<int> loopCount = new List<int>();
            List<List<qEdge>> frontLoopList = new List<List<qEdge>>();

            #region Get fronLoopList
            if (checkSideEdge == null)
            {
                List<qEdge> remainingFrontEdges = new List<qEdge>(frontEdges);
                List<qNode> frontNodes = new List<qNode>();

                // check nodes of front
                foreach (qEdge edge in frontEdges)
                {
                    if (!frontNodes.Contains(edge.StartNode))
                    {
                        frontNodes.Add(edge.StartNode);
                    }
                    if (!frontNodes.Contains(edge.EndNode))
                    {
                        frontNodes.Add(edge.EndNode);
                    }
                }

                List<List<qEdge>> loopNodes = new List<List<qEdge>>();
                foreach (qNode node in frontNodes)
                {
                    List<qEdge> connectedEdges = node.GetConnectedEdges(globalEdgeList);
                    List<qEdge> connectedFrontEdges = new List<qEdge>();
                    foreach (qEdge edge in connectedEdges)
                    {
                        if (edge.IsFrontEdge())
                        {
                            connectedFrontEdges.Add(edge);
                        }
                    }

                    if (connectedFrontEdges.Count > 2)
                    {
                        loopNodes.Add(connectedFrontEdges);
                    }
                }

                if (loopNodes.Count == 0)
                {
                    loopCount.Add(frontEdges.Count);
                }
                else 
                {
                    foreach (List<qEdge> frontMeet in loopNodes) // loop meeting points 
                    {
                        qNode sharedNode = frontMeet[0].GetSharedNode(frontMeet[1]);

                        // find loops
                        for (int i = 0; i < 2; i++) // assume only two loops can be connected in one node
                        {
                            qNode startNode = sharedNode;
                            qEdge startEdge = frontMeet[0];
                            qNode nextNode = new qNode();
                            qEdge nextEdge = new qEdge();

                            List<qEdge> currentLoop = new List<qEdge>();
                            bool done = false;
                            while (!done) // loop a edge loop
                            {

                                currentLoop.Add(startEdge);
                                // Find next node
                                nextNode = startEdge.GetOppositeNode(startNode);

                                // Find next edge
                                remainingFrontEdges.Remove(startEdge);

                                foreach (qEdge edge in remainingFrontEdges)
                                {
                                    if (edge.StartNode == nextNode | edge.EndNode == nextNode)
                                    {
                                        nextEdge = edge; 
                                        break;
                                    }
                                }

                                if (nextNode == sharedNode)
                                {
                                    done = true; continue;
                                }

                                startNode = nextNode;
                                startEdge = nextEdge;
                            }

                            frontMeet.RemoveAt(0);
                            frontMeet.Remove(startEdge); // assume be connected to node
                            frontLoopList.Add(currentLoop);
                            loopCount.Add(frontLoopList.Count);
                        }
                    }

                }
            }
            else if (checkSideEdge != null)
            {
                // get startNode and endNode, and connected front edges to startNode
                qNode startNode = checkSideEdge.StartNode;
                qNode endNode = checkSideEdge.EndNode;
                List<qEdge> startEdges = GetFrontEdgesConnectedToNode(startNode, globalEdgeList);

                // find left and right front edge connected to startNode
                qEdge startEdgeLeft = new qEdge();
                qEdge startEdgeRight = new qEdge();
                if (startEdges[0] == startEdges[1].LeftFrontNeighbor)
                {
                    startEdgeLeft = startEdges[0];
                    startEdgeRight = startEdges[1];
                }
                else if (startEdges[0] == startEdges[1].RightFrontNeighbor)
                {
                    startEdgeLeft = startEdges[1];
                    startEdgeRight = startEdges[0];
                }

                // get left loop
                List<qEdge> edgesOfLeftLoop = new List<qEdge>() { checkSideEdge, startEdgeLeft };
                bool leftLoopDone = false;
                int counter = 0;
                while (!leftLoopDone)
                {
                    edgesOfLeftLoop.Add(startEdgeLeft.LeftFrontNeighbor);
                    if (endNode != startEdgeLeft.LeftFrontNeighbor.StartNode & endNode != startEdgeLeft.LeftFrontNeighbor.EndNode)
                    {
                        startEdgeLeft = startEdgeLeft.LeftFrontNeighbor;
                    }
                    else { leftLoopDone = true; }

                    if (counter > 1000) { leftLoopDone = true; AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "IsFrontLoopEven: failed left loop"); }
                    counter++;
                }

                // get right loop
                List<qEdge> edgesOfRightLoop = new List<qEdge>() { checkSideEdge, startEdgeRight };
                bool rightLoopDone = false;
                counter = 0;
                while (!rightLoopDone)
                {
                    edgesOfRightLoop.Add(startEdgeRight.RightFrontNeighbor);
                    if (endNode != startEdgeRight.RightFrontNeighbor.StartNode & endNode != startEdgeRight.RightFrontNeighbor.EndNode)
                    {
                        startEdgeRight = startEdgeRight.RightFrontNeighbor;
                    }
                    else { rightLoopDone = true; }
                    if (counter > 1000) { rightLoopDone = true; AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "IsFrontLoopEven: failed right loop"); }
                    counter++;
                }
                loopCount.Add(edgesOfLeftLoop.Count);
                loopCount.Add(edgesOfRightLoop.Count);
            }
            #endregion Get frontLoopList

            foreach (int numEdges in loopCount)
            {
                if (numEdges % 2 != 0) { evenEdgesInLoops = false; break; }
            }

            return Tuple.Create(evenEdgesInLoops, frontLoopList);
        }
        private bool IsEdgeClosingFront(qNode N_k, qEdge edgeToCheck, List<qEdge> frontEdges, List<qEdge> globalEdgeList, qEdge E_front)
        {
            bool existingEdgeIsClosingFront = false;
            qNode N_m = edgeToCheck.GetOppositeNode(N_k);
            if (N_m.IsFrontNode(frontEdges))
            {
                List<qEdge> edgesToIgnore = new List<qEdge>() { E_front, E_front.LeftFrontNeighbor, E_front.RightFrontNeighbor };
                foreach (qEdge edge in edgesToIgnore)
                {
                    if (edge.StartNode == N_m | edge.EndNode == N_m) { existingEdgeIsClosingFront = false; break; }
                    else { existingEdgeIsClosingFront = true; }
                }
            }


            // to do: temporary, need to find a correct solution            
            // to do: delete
            /*
            List<qNode> notNodeCandidates = new List<qNode>();
                foreach (qEdge edge in GetFrontEdgesConnectedToNode(N_k, globalEdgeList))
                {
                    foreach (qEdge edgeConncected in GetConnectedEdges(GetOppositeNode(N_k, edge), globalEdgeList))
                    {
                        notNodeCandidates.Add(edgeConncected.StartNode);
                        notNodeCandidates.Add(edgeConncected.EndNode);
                    }
                }

                if (theta_i_list_sorted[0] < thetaToleranceForClosing)
                {
                    N_m = GetOppositeNode(N_k, E_i_candidates_sorted[0]);
                    if (IsFrontNode(N_m, frontEdges) & !notNodeCandidates.Contains(N_m)) // to do: change this to another criterion.
                    { 
                        existingEdgeIsClosingFront = IsFrontNode(N_m, frontEdges);
                    }
                }*/
                return existingEdgeIsClosingFront;
        }
        private Tuple<double, double, Mesh> CalculateQuality(List<qElement> globalElementList) // to do: more qualities?
        {
            // summary: calculate quality aspect ratio and create color mesh
            double avgQuality = 0;
            double badestQuality = 0;
            List<double> ARList = new List<double>();

            foreach (qElement element in globalElementList)
            {
                double maxEdgeLength = 0;
                double minEdgeLength = 10000;
                foreach(qEdge edge in element.EdgeList)
                {
                    if (edge.Length < minEdgeLength) { minEdgeLength = edge.Length;}
                    if (edge.Length > maxEdgeLength) { maxEdgeLength = edge.Length; }
                }
                double AR = minEdgeLength / maxEdgeLength;
                ARList.Add(AR);
                avgQuality = avgQuality + AR;
            }

            badestQuality = ARList.Min();
            avgQuality = Math.Round(avgQuality / globalElementList.Count, 3);

            Mesh colorMesh = new Mesh();
            for (int i = 0; i < globalElementList.Count; i++)
            {
                Mesh mesh = new Mesh();
                List<qNode> nodes = globalElementList[i].GetNodesOfElement();
                foreach (qNode node in nodes) {mesh.Vertices.Add(node.Coordinate);}

                if (globalElementList[i].IsQuad) { mesh.Faces.AddFace(0, 1, 2, 3);  }
                else { mesh.Faces.AddFace(0, 1, 2); }
     
                if (ARList[i] > 0.9) { mesh.VertexColors.CreateMonotoneMesh(Color.Green); }
                else if (ARList[i] > 0.7) { mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);}
                else if (ARList[i] > 0.6){mesh.VertexColors.CreateMonotoneMesh(Color.Orange);}
                else if (ARList[i] > 0){mesh.VertexColors.CreateMonotoneMesh(Color.Red);}
                colorMesh.Append(mesh);
            }
            return Tuple.Create(avgQuality, badestQuality, colorMesh);
        }
        private bool LoopControll(List<qEdge> frontEdges, List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            // summary: check if E_k_left and E_k_right is in same loop
            qEdge E_front = new qEdge();
            qEdge E_k_left = new qEdge();
            qEdge E_k_right = new qEdge();
            qEdge E_top = new qEdge();

            List<List<qEdge>> frontLoops = IsFrontLoopsEven(frontEdges, null, globalEdgeList).Item2;
            foreach (List<qEdge> loop in frontLoops)
            {
                if (loop.Count == 4) // prioritize
                {
                    E_front = loop[0];

                    qEdge connectedToStartNode = new qEdge();
                    qEdge connectedToEndNode = new qEdge();

                    for (int i = 1; i < 4; i++)
                    {
                        if (E_front.StartNode == loop[i].StartNode | E_front.StartNode == loop[i].EndNode)
                        {
                            connectedToStartNode = loop[i];
                        }
                        else if (E_front.EndNode == loop[i].StartNode | E_front.EndNode == loop[i].EndNode)
                        {
                            connectedToEndNode = loop[i];
                        }
                        else
                        {
                            E_top = loop[i];
                        }
                    }
                    var sideEdges = E_front.OrientateNeigborEdges(connectedToStartNode, connectedToEndNode);
                    E_k_left = sideEdges.Item1;
                    E_k_right = sideEdges.Item2;

                    List<qEdge> quadEdgeList = new List<qEdge>() { E_front, E_k_right, E_k_left, E_top };
                    qElement quadElement = CreateQuadElement(quadEdgeList, globalEdgeList, globalElementList, frontEdges);
                    DoLocalSmoothing(quadElement, globalEdgeList, frontEdges, globalElementList);
                    return true;
                }
            }
            return false;
        }

        // _______________________________________ for mesh modification __________________________________________________
        private Tuple<qEdge, bool> GetSideEdge(List<qElement> globalElementList, List<qEdge> globalEdgeList, double nodeToEvaluate, qEdge E_front, List<qEdge> frontEdges, bool triangleSidesCase)
        {
            // summary: get side edge of a new quad; nodeToEvaluate: 0 = left, 1 = right;
            qConstants constants = new qConstants();
            qEdge E_k = new qEdge();
            bool performed = true;
            double thetaTolerance = constants.GetThetaTolerance(); // constant 
            
            qEdge E_neighborFront = new qEdge();
            if (nodeToEvaluate == 0) { E_neighborFront = E_front.LeftFrontNeighbor; }
            else { E_neighborFront = E_front.RightFrontNeighbor; }

            #region Get V_k
            qNode N_k = E_front.GetSharedNode(E_neighborFront); // shared node
            var vectors = E_front.CalculateVectorsFromSharedNode(E_neighborFront);
            Vector3d vec1 = vectors.Item1; // get vector for E_front
            Vector3d vec2 = vectors.Item2; // get vector for E_neighborFront
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
            List<qEdge> connectedEdges = N_k.GetConnectedEdges(globalEdgeList); // get connected edges of shared node
            // get E_i
            List<qEdge> E_i_candidates = new List<qEdge>();
            foreach (qEdge edge in connectedEdges) // loop connected edges
            {
                // E_i
                if (!edge.IsFrontEdge()) 
                {
                     E_i_candidates.Add(edge);
                }
            }
            #endregion

            // to do: check if this is correct
            if (N_k.BoundaryNode & E_i_candidates.Count == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "SideEdge: N_k is boundary node and no E_i. No edge to select for side edge."); }

            #region Find smallest theta to edge
            // calcualte thetas 
            List<double> theta_i_list = new List<double>();
            foreach (qEdge E_i in E_i_candidates)
            {
                Vector3d E_i_vec = GetVectorOfEdgeFromNode(E_i, N_k);
                double theta_i = Vector3d.VectorAngle(V_k, E_i_vec); // to do: make more general, assume there are edges close to V_K
                theta_i_list.Add(theta_i);
            }

            // sort edges wrt theta
            double minTheta = 0;
            int minThetaIndex = 0;
            List<qEdge> E_i_candidates_sorted = new List<qEdge>(); // sorted from smallest theta angle to largest
            List<double> theta_i_list_sorted = new List<double>();
            int numEdges = theta_i_list.Count;
            for (int j = 0; j < numEdges; j++)
            {
                minTheta = theta_i_list[0];
                minThetaIndex = 0;
                for (int i = 1; i < theta_i_list.Count; ++i)
                {
                    if (theta_i_list[i] < minTheta)
                    {
                        minTheta = theta_i_list[i];
                        minThetaIndex = i;
                    }
                }
                E_i_candidates_sorted.Add(E_i_candidates[minThetaIndex]);
                theta_i_list_sorted.Add(minTheta);
                E_i_candidates.RemoveAt(minThetaIndex);
                theta_i_list.RemoveAt(minThetaIndex);
            }

            #endregion

            #region Get E_k
            qEdge E_0 = new qEdge();
            qEdge E_m = new qEdge();

            bool existingEdgeIsClosingFront = false;
            if (theta_i_list_sorted[0] < constants.GetThetaToleranceForClosing())
            {
                 existingEdgeIsClosingFront = IsEdgeClosingFront(N_k, E_i_candidates_sorted[0], frontEdges, globalEdgeList, E_front);
            }

            // if trangle case, assumption
            bool forceNextBestExisting = false;
            bool forceSwapOrSplit = false;
            if (triangleSidesCase)
            {
                if (E_i_candidates_sorted.Count > 1)
                {
                    if (theta_i_list_sorted[1] < thetaTolerance) { forceNextBestExisting = true; }
                    else {forceSwapOrSplit = true; }
                }
                else { forceSwapOrSplit = true; }
            }

            if ( ((theta_i_list_sorted[0] < thetaTolerance) | existingEdgeIsClosingFront) & !forceSwapOrSplit) // use existing edge
            {
                if (forceNextBestExisting & !existingEdgeIsClosingFront) // assume triangle case is not a closing front
                {
                    E_k = E_i_candidates_sorted[1]; // use next best existing edge b/c triangle case
                }
                else if (existingEdgeIsClosingFront) // if closing edge
                {

                    E_k = CloseFront(frontEdges, globalEdgeList, globalElementList, E_i_candidates_sorted[0], N_k);
                }
                else 
                { 
                    E_k = E_i_candidates_sorted[0]; // use existing edge as side edge
                }
            }
            else // swap or split
            {

                // sort edges wrt angle from V_k
                List<qEdge> connectedFrontEdges = GetFrontEdgesConnectedToNode(N_k, globalEdgeList);
                E_i_candidates_sorted.AddRange(connectedFrontEdges);
                
                List<double> angleFromV_kToE_i_candidates = new List<double>();
                foreach (qEdge E_i in E_i_candidates_sorted)
                {
                    Vector3d E_i_vec = GetVectorOfEdgeFromNode(E_i, N_k);
                    double theta_i = Vector3d.VectorAngle(V_k, E_i_vec, Vector3d.ZAxis); // to do: make more general
                    angleFromV_kToE_i_candidates.Add(theta_i);
                }

                double minAngle = 0;
                int minAngleIndex = 0;
                List<qEdge> E_i_candidates_sortedBasedOnAngleFromV_k = new List<qEdge>(); // sorted from smallest to largest angle from V_k ccw
                numEdges = E_i_candidates_sorted.Count;
                for (int j = 0; j < numEdges; j++)
                {
                    minAngle = angleFromV_kToE_i_candidates[0];
                    minAngleIndex = 0;
                    for (int i = 1; i < angleFromV_kToE_i_candidates.Count; ++i)
                    {
                        if (angleFromV_kToE_i_candidates[i] < minAngle)
                        {
                            minAngle = angleFromV_kToE_i_candidates[i];
                            minAngleIndex = i;
                        }
                    }
                    E_i_candidates_sortedBasedOnAngleFromV_k.Add(E_i_candidates_sorted[minAngleIndex]);
                    E_i_candidates_sorted.RemoveAt(minAngleIndex);
                    angleFromV_kToE_i_candidates.RemoveAt(minAngleIndex);
                }
                
                qEdge E_1 = E_i_candidates_sortedBasedOnAngleFromV_k[0];
                qEdge E_2 = E_i_candidates_sortedBasedOnAngleFromV_k[E_i_candidates_sortedBasedOnAngleFromV_k.Count - 1];

                Vector3d E_1_vec = GetVectorOfEdgeFromNode(E_1, N_k);
                Vector3d E_2_vec = GetVectorOfEdgeFromNode(E_2, N_k);
                Vector3d cross1 = Vector3d.CrossProduct(V_k, E_1_vec);
                Vector3d cross2 = Vector3d.CrossProduct(V_k, E_2_vec);

                if ((cross1.Z * cross2.Z) > 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "SideEdge: Split or swap not performed because V_k does not intersect E_0.");
                }
            

                // get edge between the two edges closest to V_k
                qNode E_1_NotSharedNode = E_1.GetOppositeNode(N_k);
                qNode E_2_NotSharedNode = E_2.GetOppositeNode(N_k);
                E_0 = FindEdge(globalEdgeList, E_1_NotSharedNode, E_2_NotSharedNode);

                
                if (E_0.IsFrontEdge())
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "GetSideEdge: E_0 is a front edge. Closing front is performed");
                    //E_k = CloseFront(frontEdges, globalEdgeList, globalElementList, E_i_candidates_sorted[0], N_k); // to do: assumption... må sjekke ut dette
                    performed = false;
                    return Tuple.Create(E_k, performed);
                }
                /*
       
                // check if edge is closing front
                if (E_0_list[E_0_list.Count-1].IsFrontEdge())
                {
                    if (E_0_list.Count == 1)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "GetSideEdge: E_0 is a front edge. Side edge is aborted");
                        performed = false;
                        return Tuple.Create(E_k, performed);
                    }
                    for (int i = 0; i < E_0_list.Count - 1; i++)
                    {
                        if (!E_0_list[i].IsFrontEdge())
                        {
                            E_0 = E_0_list[i];
                        }
                    }
                }
                */

                // Get N_m
                qNode N_m = new qNode();
                List<qNode> swapedNodes = GetSwapedNodes(E_0);
                if (swapedNodes[0] == N_k) { N_m = swapedNodes[1]; }
                else { N_m = swapedNodes[0]; }

                //heipådeg
                // Swap or split
                double lengthN_kN_m = N_k.Coordinate.DistanceTo(N_m.Coordinate);
                double beta = Vector3d.VectorAngle(V_k, N_m.Coordinate - N_k.Coordinate);

                if ((lengthN_kN_m < Math.Sqrt(3) * (E_front.Length + E_neighborFront.Length) * 0.5) & beta < thetaTolerance)
                {
                    performed = SwapEdge(E_0, globalElementList);
                    E_k = E_0;
                }
                else
                {
                    E_k = SplitEdge(E_0, V_k, N_k, globalEdgeList, globalElementList);
                }
                //performed = !IsFrontLoopsEven(frontEdges, E_k, globalEdgeList); // to do: slett?
            }
            #endregion

            return Tuple.Create(E_k, performed);
        }
        private bool SwapEdge(qEdge E_0, List<qElement> globalElementList)
        {
            // summary: swap edge if swapable

            bool performed = true ;
            // is swapable
            if (E_0.Element1.IsQuad | E_0.Element2.IsQuad)
            {
                performed = false;
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "SwapEdge: Edge is not swapable.");
                return performed;
            }

            // get nodes of edge to swap
            qNode node1 = E_0.StartNode;
            qNode node2 = E_0.EndNode;
            List<qNode> swapedNodes = GetSwapedNodes(E_0);

            // is swapable
            Vector3d vector1 = swapedNodes[0].Coordinate - node1.Coordinate;
            Vector3d vector2 = swapedNodes[1].Coordinate - node1.Coordinate;

            Vector3d vector3 = swapedNodes[0].Coordinate - node2.Coordinate;
            Vector3d vector4 = swapedNodes[1].Coordinate - node2.Coordinate;
            double check12 = Vector3d.VectorAngle(vector1, vector2, Vector3d.ZAxis); // to do: make more general
            double check34 = Vector3d.VectorAngle(vector3, vector4, Vector3d.ZAxis);  // to do: make more general

            if (check12 >= Math.PI & check34 >= Math.PI)
            {
                performed = false;
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "SwapEdge: Edge is not swapable b/c of elements will be inverted.");
                return performed;
            }
            else if (check12 <= Math.PI & check34 <= Math.PI)
            {
                performed = false;
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "SwapEdge: Edge is not swapable b/c of elements will be inverted.");
                return performed;
            }

            // get edges of adjacent elements
            List<qEdge> element1edges = E_0.Element1.EdgeList;
            List<qEdge> element2edges = E_0.Element2.EdgeList;

            // get list of all edges for new elements
            List<qEdge> edgeForNewElements = new List<qEdge>();
            element1edges.Remove(E_0);
            element2edges.Remove(E_0);
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
            globalElementList.Remove(oldElement1);
            globalElementList.Remove(oldElement2);
            
            return performed;
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
            element1edges.Remove(E_0);
            element2edges.Remove(E_0);

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
                    if (edge.StartNode == E_n_1.GetOppositeNode(N_n) | edge.EndNode == E_n_1.GetOppositeNode(N_n))
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
                    if (edge.StartNode == E_n_1.GetOppositeNode(N_n) | edge.EndNode == E_n_1.GetOppositeNode(N_n))
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
            List<qEdge> keptEdges = new List<qEdge>() { newElement1Edges[2], newElement2Edges[2], newElement3Edges[2], newElement4Edges[2] };
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
                qEdge keptEdge = keptEdges[i];
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
            globalEdgeList.Remove(E_0);
            foreach (qEdge newEdge in newEdges)
            {
                globalEdgeList.Add(newEdge);
            }

            // update elementList
            globalElementList.Remove(E_0.Element1);
            globalElementList.Remove(E_0.Element2);
            foreach (qElement newElement in newElements)
            {
                globalElementList.Add(newElement);
            }
            return E_k;
        }
        private Tuple<qEdge, qEdge, qEdge, bool> GetTopEdge(qEdge E_front, qEdge E_k_left, qEdge E_k_right, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: top edge recovery
            bool performed = true;
            bool E_k_left_updated = false;
            bool E_k_right_updated = false;

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

            if (N_c == N_d) // triangle case
            {
                // to do: legg inn prioritering på left/right edge ut fra antall edges connected to node..

                // find new side edge
                if (GetEdgeState(E_front, frontEdges)[0] == 0)
                {
                    // select E_k_left based on triangle case
                    var leftSideEdgeValues = GetSideEdge(globalElementList, globalEdgeList, 0, E_front, frontEdges, true);
                    E_k_left = leftSideEdgeValues.Item1;
                    E_k_left_updated = leftSideEdgeValues.Item2; // false only when split/swap is performed on an closing edge, i.e. E_0 is a front edge

                    if (E_k_left.StartNode == leftNode) { N_d = E_k_left.EndNode; }
                    else { N_d = E_k_left.StartNode; }
                }
                               
                if (GetEdgeState(E_front, frontEdges)[1] == 0 & !E_k_left_updated)
                {
                    // select E_k_right based on triangle case
                    var rightSideEdgeValues = GetSideEdge(globalElementList, globalEdgeList, 1, E_front, frontEdges, true);
                    E_k_right = rightSideEdgeValues.Item1;
                    E_k_right_updated = rightSideEdgeValues.Item2; // false only when split/swap is performed on an closing edge, i.e. E_0 is a front edge

                    if (E_k_right.StartNode == rightNode) { N_c = E_k_right.EndNode; }
                    else { N_c = E_k_right.StartNode; }
                }

                if (!E_k_left_updated & !E_k_right_updated) // triangle case that could not be fixed
                {
                    performed = false;
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Top edge recovery not performed because N_C and N_d is equal.");
                    return Tuple.Create(E_top, E_k_left, E_k_right, performed);
                }
            }

            var E_topAndPerfomed = EdgeRecoveryProcess( N_c, N_d, globalEdgeList, globalElementList, frontEdges);
            E_top = E_topAndPerfomed.Item1;
            performed = E_topAndPerfomed.Item2;
            return Tuple.Create(E_top, E_k_left, E_k_right, performed);
        }
        private Tuple<qEdge, bool> EdgeRecoveryProcess(qNode N_c, qNode N_d, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            qEdge E_recovered = new qEdge();
            bool performed = true;

            #region Get connected elements to node _N_c
            List<qElement> connectedElements = new List<qElement>();
            List<qEdge> connectedEdges = N_c.GetConnectedEdges(globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                List<qElement> edgeElements = edge.GetConnectedElements();
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
                Point3d elementCenter = connectedElement.GetElementCenter();
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

                var vectors = E_kCandidates[0].CalculateVectorsFromSharedNode(E_kCandidates[1]);
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

                
                if (V_k.IsParallelTo(V_s) == 1 & V_k.Length == V_s.Length)  // the recovered edge can consist of more than one parallel edge 
                {
                    E_recovered = E_k;
                    return Tuple.Create(E_recovered, performed);
                }
                else if (V_k1.IsParallelTo(V_s) == 1 & V_k1.Length == V_s.Length)  // the recovered edge can consist of more than one parallel edge 
                {
                    E_recovered = E_k1;
                    return Tuple.Create(E_recovered, performed);
                }
                else if (Vector3d.CrossProduct(V_s, V_k).Z <= -0.0001 & Vector3d.CrossProduct(V_s, V_k1).Z >= -0.0001) // todo: correct compared to Owen??
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
                
                List<qNode> elementNodes = T_k1.GetNodesOfElement();
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

                var vectors = E_nCandidates[0].CalculateVectorsFromSharedNode(E_nCandidates[1]);
                Vector3d vec1 = vectors.Item1;
                Vector3d vec2 = vectors.Item2;
                qEdge E_n = new qEdge();
                qEdge E_n1 = new qEdge();

                if (Vector3d.CrossProduct(vec1, vec2).Z > 0) // to do: make more general for 3d?
                {
                    E_n = E_nCandidates[1];
                    E_n1 = E_nCandidates[0];
                }
                else
                {
                    E_n = E_nCandidates[0];
                    E_n1 = E_nCandidates[1];
                }

                if (Vector3d.CrossProduct(V_i, V_s).Z <= 0) // OBS: changes from Owen's article, to do: make more general for 3d?
                {
                    E_i = E_n;
                }
                else
                {
                    E_i = E_n1;
                }
                if (!E_i.IsFrontEdge()) { intersectingS.Add(E_i); }
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
            int n = 0;
            while (intersectingS.Count > 0 & n < 100)
            {
                E_recovered = intersectingS[0];
                bool swapPerformed = SwapEdge(E_recovered, globalElementList); // swap the edge
                
                // check if inverted
                if (swapPerformed)
                {
                    if (IntersectionCheck(E_recovered, N_c, N_d))
                    {
                        intersectingS.Add(E_recovered);
                        intersectingS.RemoveAt(0);
                    }
                    else 
                    {
                        intersectingS.RemoveAt(0);
                    }

                }
                else
                {
                    intersectingS.Add(E_recovered);
                    intersectingS.RemoveAt(0);
                }
                n++;
            }
            #endregion

            return Tuple.Create(E_recovered, performed);
        }
        private bool IntersectionCheck(qEdge E_i, qNode N_c, qNode N_d)
        {
            Vector3d V_s = N_d.Coordinate - N_c.Coordinate;
            Vector3d V_endNode = E_i.EndNode.Coordinate - N_c.Coordinate;
            Vector3d V_startNode = E_i.StartNode.Coordinate - N_c.Coordinate;

            if (Vector3d.CrossProduct(V_s, V_endNode).Z*Vector3d.CrossProduct(V_s, V_startNode).Z <  0 ) 
            {
                return true;
            }
            return false;
        }
        private void Seam(qEdge rightEdgeToSeam, qEdge leftEdgeToSeam, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: seam right edge and left edge together to a new seam edge with the adjacent quad elements

            // get needed nodes for seaming
            qNode N_k = rightEdgeToSeam.GetSharedNode(leftEdgeToSeam);
            qNode N_k_left = leftEdgeToSeam.GetOppositeNode(N_k);
            qNode N_k_right = rightEdgeToSeam.GetOppositeNode(N_k);

            if( N_k_left.BoundaryNode | N_k_right.BoundaryNode)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Seam: failed because of boundary nodes.");
                return;
            }

            // edge recovery
            var edgeAndPerformed = EdgeRecoveryProcess(N_k_right, N_k_left, globalEdgeList, globalElementList, frontEdges);
            qEdge E_0 = edgeAndPerformed.Item1;

            // check if edge recovery failed
            bool performed = edgeAndPerformed.Item2;
            if (!performed) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Seam: Edge recovery failed."); return; }

            // get N_new
            Point3d newNodeCoordinate = 0.5 * (E_0.EndNode.Coordinate + E_0.StartNode.Coordinate);
            qNode N_new = new qNode(newNodeCoordinate, false);

            // get N_t
            qNode N_t = new qNode();
            List<qNode> swapedNodes = GetSwapedNodes(E_0);
            if (N_k == swapedNodes[0]) { N_t = swapedNodes[1]; }
            else if (N_k == swapedNodes[1]) { N_t = swapedNodes[0]; }

            // create new seam edges
            qEdge seamEdgeN_K = new qEdge(N_k, N_new);
            qEdge seamEdgeN_t = new qEdge(N_new, N_t);

            // get elements to E_0 
            qElement elementToN_k = new qElement(); 
            qElement elementToN_t = new qElement();
            List<qElement> elementsToE_0 = E_0.GetConnectedElements();
            if (elementsToE_0[0].GetNodesOfElement().Contains(N_k))
            {
                elementToN_k = elementsToE_0[0]; elementToN_t = elementsToE_0[1];
            }
            else if (elementsToE_0[0].GetNodesOfElement().Contains(N_t))
            {
                elementToN_t = elementsToE_0[0]; elementToN_k = elementsToE_0[1];
            }
            else  { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Seam: Problems with elements connected to E_0."); return; }

            #region Merge elements
            List<qElement> elementsToSeamAtN_k = new List<qElement>();

            // fix elements to N_k
            int indexToUpdate = 0;
            qEdge leftEdgeToN_k = FindEdge(globalEdgeList, N_k, N_k_left); // in case any changes to edge
            qEdge rightEdgeToN_k = FindEdge(globalEdgeList, N_k, N_k_right); // in case any changes to edge
            foreach (qEdge edge in elementToN_k.EdgeList)
            {
                if (edge == E_0) { continue; }

                if (edge == leftEdgeToN_k)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(leftEdgeToN_k);
                        edge.Element1.EdgeList[indexToUpdate] = seamEdgeN_K;
                        elementsToSeamAtN_k.Add(edge.Element1);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(leftEdgeToN_k);
                        edge.Element2.EdgeList[indexToUpdate] = seamEdgeN_K;
                        elementsToSeamAtN_k.Add(edge.Element2);
                    }
                }
                
                if (edge == rightEdgeToN_k)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(rightEdgeToN_k);
                        edge.Element1.EdgeList[indexToUpdate] = seamEdgeN_K;
                        elementsToSeamAtN_k.Add(edge.Element1);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(rightEdgeToN_k);
                        edge.Element2.EdgeList[indexToUpdate] = seamEdgeN_K;
                        elementsToSeamAtN_k.Add(edge.Element2);
                    }
                }
            }
            seamEdgeN_K.Element1 = elementsToSeamAtN_k[0];
            seamEdgeN_K.Element2 = elementsToSeamAtN_k[1];

            // fix elements to N_t
            indexToUpdate = 0;
            List<qElement> elementsToSeamAtN_t = new List<qElement>();
            qEdge leftEdgeToN_t = FindEdge(globalEdgeList, N_t, N_k_left); // in case any changes to edge
            qEdge rightEdgeToN_t = FindEdge(globalEdgeList, N_t, N_k_right); // in case any changes to edge
            foreach (qEdge edge in elementToN_t.EdgeList)
            {
                if (edge == E_0) { continue; }

                if (edge == leftEdgeToN_t)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(leftEdgeToN_t);
                        edge.Element1.EdgeList[indexToUpdate] = seamEdgeN_t;
                        elementsToSeamAtN_t.Add(edge.Element1);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(leftEdgeToN_t);
                        edge.Element2.EdgeList[indexToUpdate] = seamEdgeN_t;
                        elementsToSeamAtN_t.Add(edge.Element2);
                    }
                }

                if (edge == rightEdgeToN_t)
                {
                    if (edge.Element1.IsQuad)
                    {
                        indexToUpdate = edge.Element1.EdgeList.IndexOf(rightEdgeToN_t);
                        edge.Element1.EdgeList[indexToUpdate] = seamEdgeN_t;
                        elementsToSeamAtN_t.Add(edge.Element1);
                    }
                    else
                    {
                        indexToUpdate = edge.Element2.EdgeList.IndexOf(rightEdgeToN_t);
                        edge.Element2.EdgeList[indexToUpdate] = seamEdgeN_t;
                        elementsToSeamAtN_t.Add(edge.Element2);
                    }
                }
            }
            seamEdgeN_t.Element1 = elementsToSeamAtN_t[0];
            seamEdgeN_t.Element2 = elementsToSeamAtN_t[1];
            #endregion 

           

            // update global element list
            globalElementList.Remove(elementsToE_0[0]);
            globalElementList.Remove(elementsToE_0[1]);

            // update global edge list
            globalEdgeList.Remove(E_0);
            globalEdgeList.Remove(leftEdgeToN_k);
            globalEdgeList.Remove(rightEdgeToN_k);
            globalEdgeList.Remove(leftEdgeToN_t);
            globalEdgeList.Remove(rightEdgeToN_t);

            // update frontEdges
            frontEdges.Remove(leftEdgeToN_k);
            frontEdges.Remove(rightEdgeToN_k);
            List<qEdge> edgesToCheck =  N_new.GetConnectedEdges(globalEdgeList);

            foreach (qEdge edge in edgesToCheck)
            {
                if (edge.IsFrontEdge())
                {
                    frontEdges.Add(edge);
                    indexToUpdate = globalEdgeList.IndexOf(edge);
                    globalEdgeList[indexToUpdate].Level++;
                }
            }
            SetNeighorFrontEdges(frontEdges);

            // local smoothing of new quads
            DoLocalSmoothing(elementsToSeamAtN_k[0], globalEdgeList, frontEdges, globalElementList);
            DoLocalSmoothing(elementsToSeamAtN_k[1], globalEdgeList, frontEdges, globalElementList);
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Seam: is performed.");
            return;
        }
        private Tuple<qEdge, qEdge, qEdge> TransitionSeam(qEdge edge1, qEdge edge2, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: if small angles between connected front edges and large transitions, perform transition seam. Return E_front, E_k_left, E_k_right.

            qEdge E_front = new qEdge();
            qEdge E_k_left = new qEdge();
            qEdge E_k_right = new qEdge();

            // find shortest and longest edge
            qEdge E_long = new qEdge();
            qEdge E_short = new qEdge();
            if (edge1.Length < edge2.Length) { E_short = edge1; E_long = edge2; }
            else { E_short = edge2; E_long = edge1; }

            // nodes to use
            qNode N_k = edge1.GetSharedNode(edge2);
            qNode N_k_short = E_short.GetOppositeNode(N_k);
            qNode N_k_long = E_long.GetOppositeNode(N_k);

            // split E_long at mid point
            Point3d midPointE_long = 0.5 * (E_long.EndNode.Coordinate + E_long.StartNode.Coordinate);
            qNode N_new = new qNode(midPointE_long, false);

            if (E_long.EndNode.BoundaryNode & E_long.StartNode.BoundaryNode)
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Transition seam: Wrong to assumed no boundary nodes of E_long."); }

            // get E_front and elements connected to E_long 
            qElement quadElementOfE_long = new qElement();
            qElement triElementOfE_long = new qElement();

            List<qElement> connectedElementsToE_long = E_long.GetConnectedElements(); // connected elements to E_long
            List<qEdge> connectedEdgesToN_k = N_k.GetConnectedEdges(globalEdgeList); // connected edges to N_k

            foreach (qElement element in connectedElementsToE_long)
            {
                if (element.IsQuad)
                {
                    quadElementOfE_long = element; 
                    foreach (qEdge edge in element.EdgeList)
                    {
                        if (connectedEdgesToN_k.Contains(edge) & edge != E_long)
                        { 
                            E_front = edge;
                        }
                    }
                }
                else
                {
                    triElementOfE_long = element;
                }             
            }

            if (quadElementOfE_long.EdgeList == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition seam: Cannot perform because E_long has no connected quad element.");
                return Tuple.Create(E_front, E_k_right, E_k_left);
            }

            // create new edges

            qNode N_triangle_tip = new qNode();
            var nodesOfTriElementOfE_long = triElementOfE_long.GetNodesOfElement();
            foreach (qNode node in nodesOfTriElementOfE_long)
            {
                if (node != N_k & node != N_k_long)
                {
                    N_triangle_tip = node;
                }
            }

            if (nodesOfTriElementOfE_long.Count > 1)
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "TransitionSeam: Nodes of triangle element of E_long is not correctly found.");}

            //1.
            qEdge newEdgeInTriangle = new qEdge(N_new, nodesOfTriElementOfE_long[0]);
            
            //2.
            qEdge newEdgeInQuad = new qEdge(N_new, E_front.GetOppositeNode(N_k));

            // 3. 4.
            qEdge E_long_part1 = new qEdge(N_k, N_new);
            qEdge E_long_part2 = new qEdge(N_new, N_k_long);

            // new elements
            List<qEdge> newElementFromQuadEdges = new List<qEdge>() { E_front, E_long_part1, newEdgeInQuad };
            List<qEdge> newElementFromTriEdges_part1 = new List<qEdge>() { E_long_part1, newEdgeInTriangle, FindEdge(globalEdgeList, N_k, N_triangle_tip) };
            List<qEdge> newElementFromTriEdges_part2 = new List<qEdge>() { E_long_part2, newEdgeInTriangle, FindEdge(globalEdgeList, N_k_long, N_triangle_tip) };

            qElement newElementFromQuad = new qElement(newElementFromQuadEdges);
            qElement newElementFromTri_part1 = new qElement(newElementFromTriEdges_part1);
            qElement newElementFromTri_part2 = new qElement(newElementFromTriEdges_part2);
        
            List<qEdge> modifiedEdgeListFromRemainingQuad = new List<qEdge>(quadElementOfE_long.EdgeList);
            int indexToUpdate = modifiedEdgeListFromRemainingQuad.IndexOf(E_front);
            if (indexToUpdate == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Transition seam: did not find E_front in edgelist. E_front is only a copy of original edge??");
            }
            modifiedEdgeListFromRemainingQuad[indexToUpdate] = newEdgeInQuad;
            indexToUpdate = modifiedEdgeListFromRemainingQuad.IndexOf(E_long);
            modifiedEdgeListFromRemainingQuad[indexToUpdate] = E_long_part1;
            quadElementOfE_long.EdgeList = modifiedEdgeListFromRemainingQuad;

            // fix elements
            // fix newEdgeInQuad
            newEdgeInQuad.Element1 = quadElementOfE_long;
            newEdgeInQuad.Element2 = newElementFromQuad;

            // fix E_front
            if (E_front.Element1 == quadElementOfE_long)
            { E_front.Element1 = newElementFromQuad; }
            else if (E_front.Element2 == quadElementOfE_long)
            { E_front.Element2 = newElementFromQuad; }
            else
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition seam: problems with assigning neighbors to E_front"); }

            // fix E_long_part1
            E_long_part1.Element1 = newElementFromQuad;
            E_long_part1.Element2 = newElementFromTri_part1;

            // fix E_long_part2
            E_long_part2.Element1 = quadElementOfE_long;
            E_long_part2.Element2 = newElementFromTri_part2;

            // fix newEdgeInTriangle
            newEdgeInTriangle.Element1 = newElementFromTri_part1;
            newEdgeInTriangle.Element2 = newElementFromTri_part2;

            // fix old edges from triangle
            qEdge oldEdgeTriangle_part1 = FindEdge(globalEdgeList, N_k, N_triangle_tip);
            qEdge oldEdgeTriangle_part2 = FindEdge(globalEdgeList, N_k_long, N_triangle_tip);

            if (oldEdgeTriangle_part1.Element1 == triElementOfE_long)
            { oldEdgeTriangle_part1.Element1 = newElementFromTri_part1; }
            else if (oldEdgeTriangle_part1.Element2 == triElementOfE_long)
            { oldEdgeTriangle_part1.Element2 = newElementFromTri_part1; }
            else
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition seam: problems with assigning neighbors to oldEdgeTriangle_part1"); }


            if (oldEdgeTriangle_part2.Element1 == triElementOfE_long)
            { oldEdgeTriangle_part2.Element1 = newElementFromTri_part2; }
            else if (oldEdgeTriangle_part2.Element2 == triElementOfE_long)
            { oldEdgeTriangle_part2.Element2 = newElementFromTri_part2; }
            else
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition seam: problems with assigning neighbors to oldEdgeTriangle_part2"); }

            // check
            if (!(E_front.Element1 == null & E_front.Element2 == null))
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition seam: problems with assigning neighbors. E_front is not ok."); }

            /* to do: fixslett
            // fix edge order
            newElementFromQuad.FixElementEdgeAndAngle();
            newElementFromTri_part1.FixElementEdgeAndAngle();
            newElementFromTri_part2.FixElementEdgeAndAngle();
            */

            // update global lists
            globalElementList.Remove(triElementOfE_long);
            globalElementList.Add(newElementFromQuad);
            globalElementList.Add(newElementFromTri_part1);
            globalElementList.Add(newElementFromTri_part2);

            globalEdgeList.Remove(E_long);
            globalEdgeList.Add(newEdgeInQuad);
            globalEdgeList.Add(newEdgeInTriangle);
            globalEdgeList.Add(E_long_part1);
            globalEdgeList.Add(E_long_part2);

            // update frontEdges
            frontEdges.Add(E_long_part2);
            frontEdges.Add(newEdgeInQuad);
            frontEdges.Add(E_front);
            frontEdges.Remove(E_long);
            SetNeighorFrontEdges(frontEdges);

            // local smoothing
            //DoLocalSmoothing(quadElementOfE_long, globalEdgeList, frontEdges, globalElementList);

            // add a controll for all edges in elements are found in global lists

            qEdge neigborEdgeToStartNode = new qEdge();
            qEdge neigborEdgeToEndNode = new qEdge();

            List<qEdge> sideEdgeCandidates = new List<qEdge>() { E_short, newEdgeInQuad };
            foreach (qEdge sideEdgeCandidate in sideEdgeCandidates)
            {
                if ((E_front.StartNode == sideEdgeCandidate.StartNode) | (E_front.StartNode == sideEdgeCandidate.EndNode)) // check edge start node if connected
                {
                    neigborEdgeToStartNode = sideEdgeCandidate;
                }
                else if ((E_front.EndNode == sideEdgeCandidate.StartNode) | (E_front.EndNode == sideEdgeCandidate.EndNode)) // check edge end point if connected
                {
                    neigborEdgeToEndNode = sideEdgeCandidate;
                }
            }
            
            var sideEdges= E_front.OrientateNeigborEdges(neigborEdgeToStartNode, neigborEdgeToEndNode);
            E_k_left = sideEdges.Item1;
            E_k_right = sideEdges.Item2;

            return Tuple.Create(E_front, E_k_right, E_k_left);
        }

        private Tuple<qEdge, qEdge, qEdge> TransitionSplit(qEdge edge1, qEdge edge2, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: transition split if large transition between adjacent edges, with no criterior for angle

            qEdge E_front = new qEdge();
            qEdge E_k_left = new qEdge();
            qEdge E_k_right = new qEdge();

            // find shortest and longest edge
            qEdge E_long = new qEdge();
            qEdge E_short = new qEdge();
            if (edge1.Length < edge2.Length) { E_short = edge1; E_long = edge2; }
            else { E_short = edge2; E_long = edge1; }

            // nodes to use
            qNode N_k = edge1.GetSharedNode(edge2);
            qNode N_k_short = E_short.GetOppositeNode(N_k);
            qNode N_k_long = E_long.GetOppositeNode(N_k);

            if (E_long.EndNode.BoundaryNode & E_long.StartNode.BoundaryNode)
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Transition split: Wrong to assumed no boundary nodes of E_long."); }

            // get elements connected to E_long 
            qElement quadElementOfE_long = new qElement();
            qElement triElementOfE_long = new qElement();

            List<qElement> connectedElementsToE_long = E_long.GetConnectedElements(); // connected elements to E_long
            List<qEdge> connectedEdgesToN_k = N_k.GetConnectedEdges(globalEdgeList); // connected edges to N_k

            foreach (qElement element in connectedElementsToE_long)
            {
                if (element.IsQuad)
                {
                    quadElementOfE_long = element;
                }
                else
                {
                    triElementOfE_long = element;
                }
            }

            if (quadElementOfE_long.EdgeList == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition split: Cannot perform because E_long has no connected quad element.");
                return Tuple.Create(E_front, E_k_right, E_k_left);
            }

            Point3d midPointE_long = 0.5 * (E_long.EndNode.Coordinate + E_long.StartNode.Coordinate);
            qNode N_new_midLong = new qNode(midPointE_long, false);

            qNode N_triangle_tip = new qNode();
            var nodesOfTriElementOfE_long = triElementOfE_long.GetNodesOfElement();
            foreach (qNode node in nodesOfTriElementOfE_long)
            { 
                if (node != N_k  & node != N_k_long)
                {
                    N_triangle_tip = node;
                }
            }

            qNode N_new_quadCorner = new qNode();
            foreach (qEdge edge in quadElementOfE_long.EdgeList)
            {
                List<qNode> nodes = new List<qNode>() { edge.EndNode, edge.StartNode };
                if (nodes.Contains(N_k_long) & !nodes.Contains(N_k))
                {
                    N_new_quadCorner = edge.GetOppositeNode(N_k_long);
                }
            }

            Point3d N_new_quadMid_coord = 0.5 * (N_k.Coordinate + N_new_quadCorner.Coordinate);
            qNode N_new_quadMid = new qNode(N_new_quadMid_coord, false);

            // create new edges
            qEdge newEdgeInTriangle = new qEdge(N_new_midLong, N_triangle_tip);
            qEdge newEdgeInQuad_part1 = new qEdge(N_k, N_new_quadMid);
            qEdge newEdgeInQuad_part2 = new qEdge(N_new_quadMid, N_new_quadCorner);
            qEdge newEdgeFromQuadToTriangle = new qEdge(N_new_quadMid, N_new_midLong);
            qEdge newEdge_long_part1 = new qEdge(N_k, N_new_midLong);
            qEdge newEdge_long_part2 = new qEdge(N_new_midLong, N_k_long);

            // new elements
            // get edges to elements

            List<qEdge> newElementFromQuadEdges_part2 = new List<qEdge>() {newEdgeFromQuadToTriangle, newEdgeInQuad_part1, newEdge_long_part1 }; // done
            List<qEdge> newElementFromQuadEdges_part3 = new List<qEdge>() { newEdgeFromQuadToTriangle, newEdgeInQuad_part2, newEdge_long_part2 }; // not done
            List<qEdge> newElementFromQuadEdges_part1 = new List<qEdge>() { newEdgeInQuad_part1, newEdgeInQuad_part2 }; // not done
            foreach (qEdge edge in quadElementOfE_long.EdgeList)
            {
                if (edge == E_long) { continue; }
                List<qNode> nodes = new List<qNode>() { edge.EndNode, edge.StartNode };
                if (!nodes.Contains(N_k_long))
                {
                    newElementFromQuadEdges_part1.Add(edge);
                }
                else 
                {
                    newElementFromQuadEdges_part3.Add(edge);
                }
            }

            List<qEdge> newElementFromTriEdges_part1 = new List<qEdge>() { newEdge_long_part1, newEdgeInTriangle, FindEdge(globalEdgeList, N_k, N_triangle_tip) };
            List<qEdge> newElementFromTriEdges_part2 = new List<qEdge>() { newEdge_long_part2, newEdgeInTriangle, FindEdge(globalEdgeList, N_k_long, N_triangle_tip) };

            qElement newElementFromQuad_part1 = new qElement(newElementFromQuadEdges_part1);
            qElement newElementFromQuad_part2 = new qElement(newElementFromQuadEdges_part2);
            qElement newElementFromQuad_part3 = new qElement(newElementFromQuadEdges_part3);
            qElement newElementFromTri_part1 = new qElement(newElementFromTriEdges_part1);
            qElement newElementFromTri_part2 = new qElement(newElementFromTriEdges_part2);

            // fix elements to new edges
            // fix newEdgeInQuad_part1
            newEdgeInQuad_part1.Element1 = newElementFromQuad_part1;
            newEdgeInQuad_part1.Element2 = newElementFromQuad_part2;

            // fix newEdgeInQuad_part2
            newEdgeInQuad_part2.Element1 = newElementFromQuad_part1;
            newEdgeInQuad_part2.Element2 = newElementFromQuad_part3;

            // fix newEdgeFromQuadToTri
            newEdgeFromQuadToTriangle.Element1 = newElementFromQuad_part2;
            newEdgeFromQuadToTriangle.Element2 = newElementFromQuad_part3;

            // fix newEdge_long_part1
            newEdge_long_part1.Element1 = newElementFromQuad_part2;
            newEdge_long_part1.Element2 = newElementFromTri_part1;

            // fix newEdge_long_part2
            newEdge_long_part2.Element1 = newElementFromQuad_part3;
            newEdge_long_part2.Element2 = newElementFromTri_part2;

            // fix newEdgeInTri
            newEdgeInTriangle.Element1 = newElementFromTri_part1;
            newEdgeInTriangle.Element2 = newElementFromTri_part2;

            // fix old edges from quad
            foreach (qEdge edge in quadElementOfE_long.EdgeList)
            {
                List<qNode> nodes = new List<qNode>() { edge.EndNode, edge.StartNode };
                if (edge.Element1 == quadElementOfE_long & !nodes.Contains(N_k_long))
                {
                    edge.Element1 = newElementFromQuad_part1;
                }
                else if (edge.Element2 == quadElementOfE_long & !nodes.Contains(N_k_long))
                {
                    edge.Element2 = newElementFromQuad_part1;
                }
                else if (edge.Element1 == quadElementOfE_long & nodes.Contains(N_k_long)) // note that the edge of quad connected to both part 2 and 3 will be deleted
                {
                    edge.Element1 = newElementFromQuad_part3;
                }
                else if (edge.Element2 == quadElementOfE_long & nodes.Contains(N_k_long))
                {
                    edge.Element2 = newElementFromQuad_part3;
                }
            }

            // fix old edges from triangle
            qEdge oldEdgeTriangle_part1 = FindEdge(globalEdgeList, N_k, N_triangle_tip);
            qEdge oldEdgeTriangle_part2 = FindEdge(globalEdgeList, N_k_long, N_triangle_tip);

            if (oldEdgeTriangle_part1.Element1 == triElementOfE_long)
            { oldEdgeTriangle_part1.Element1 = newElementFromTri_part1; }
            else if (oldEdgeTriangle_part1.Element2 == triElementOfE_long)
            { oldEdgeTriangle_part1.Element2 = newElementFromTri_part1; }
            else
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition seam: problems with assigning neighbors to oldEdgeTriangle_part1"); }

            if (oldEdgeTriangle_part2.Element1 == triElementOfE_long)
            { oldEdgeTriangle_part2.Element1 = newElementFromTri_part2; }
            else if (oldEdgeTriangle_part2.Element2 == triElementOfE_long)
            { oldEdgeTriangle_part2.Element2 = newElementFromTri_part2; }
            else
            { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Transition seam: problems with assigning neighbors to oldEdgeTriangle_part2"); }

            /* to do: fixslett
            // fix edge order
            newElementFromQuad_part1.FixElementEdgeAndAngle();
            newElementFromQuad_part3.FixElementEdgeAndAngle();
            newElementFromQuad_part2.FixElementEdgeAndAngle();
            newElementFromTri_part1.FixElementEdgeAndAngle();
            newElementFromTri_part2.FixElementEdgeAndAngle();
            */

            // update global lists
            globalElementList.Remove(quadElementOfE_long);
            globalElementList.Remove(triElementOfE_long);
            globalElementList.Add(newElementFromQuad_part1);
            globalElementList.Add(newElementFromQuad_part2);
            globalElementList.Add(newElementFromQuad_part3);
            globalElementList.Add(newElementFromTri_part1);
            globalElementList.Add(newElementFromTri_part2);

            globalEdgeList.Remove(E_long);
            globalEdgeList.Add(newEdgeInQuad_part1);
            globalEdgeList.Add(newEdgeInQuad_part2);
            globalEdgeList.Add(newEdge_long_part1);
            globalEdgeList.Add(newEdge_long_part2);
            globalEdgeList.Add(newEdgeFromQuadToTriangle);
            globalEdgeList.Add(newEdgeInTriangle);

            // update frontEdges
            frontEdges.Add(newEdgeFromQuadToTriangle);
            frontEdges.Add(newEdgeInQuad_part1);
            frontEdges.Add(newEdge_long_part2);
            frontEdges.Remove(E_long);
            SetNeighorFrontEdges(frontEdges);

            // local smoothing
            //DoLocalSmoothing(newElementFromQuad_part1, globalEdgeList, frontEdges, globalElementList);
            //DoLocalSmoothing(newElementFromQuad_part3, globalEdgeList, frontEdges, globalElementList);

            // to do: add a controll for all edges in elements are found in global lists

            #region Find left and right side edge
            E_front = newEdgeInQuad_part1;
            qEdge neigborEdgeToStartNode = new qEdge();
            qEdge neigborEdgeToEndNode = new qEdge();

            List<qEdge> sideEdgeCandidates = new List<qEdge>() { E_short, newEdgeFromQuadToTriangle };
            foreach (qEdge sideEdgeCandidate in sideEdgeCandidates)
            {
                if ((E_front.StartNode == sideEdgeCandidate.StartNode) | (E_front.StartNode == sideEdgeCandidate.EndNode)) // check edge start node if connected
                {
                    neigborEdgeToStartNode = sideEdgeCandidate;
                }
                else if ((E_front.EndNode == sideEdgeCandidate.StartNode) | (E_front.EndNode == sideEdgeCandidate.EndNode)) // check edge end point if connected
                {
                    neigborEdgeToEndNode = sideEdgeCandidate;
                }
            }

            Point3d midPointEdg = 0.5 * (E_front.StartNode.Coordinate + E_front.EndNode.Coordinate); // mid point of edge

            Point3d centerPoint = E_front.GetFrontElement().GetElementCenter();

            Vector3d centerToMidVector = midPointEdg - centerPoint;

            Vector3d centerToEndNodeVector = E_front.EndNode.Coordinate - centerPoint;

            Vector3d centerToStartNodeVector = E_front.StartNode.Coordinate - centerPoint;

            double startAngle = Vector3d.VectorAngle(centerToMidVector, centerToStartNodeVector, Vector3d.ZAxis); // todo: make normal more general

            double endAngle = Vector3d.VectorAngle(centerToMidVector, centerToEndNodeVector, Vector3d.ZAxis); // todo: make normal more general

            if (endAngle < startAngle)
            {
                E_k_left = neigborEdgeToStartNode;
                E_k_right = neigborEdgeToEndNode;
            }
            else
            {
                E_k_left = neigborEdgeToEndNode;
                E_k_right = neigborEdgeToStartNode;
            }
            #endregion Find left and right side edge


            return Tuple.Create(E_front, E_k_right, E_k_left);
        }
        private qElement CreateQuadElement(List<qEdge> quadEdge, List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // get inside elements
            qEdge E_front = quadEdge[0];
            List<qElement> elementInside = new List<qElement>() { E_front.GetFrontElement() };
            qElement startElement = E_front.GetFrontElement(); // start with front Edge 
            bool done = false;
            int count = 0;
            while (!done)
            {
                foreach (qEdge edge in startElement.EdgeList)
                {
                    if (quadEdge.Contains(edge)) { continue; }

                    List<qElement> connectedElements = edge.GetConnectedElements();
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
            startElement = elementInside[1]; // start with element that has perhaps not been used in loop above
            done = false;
            count = 0;
            while (!done)
            {
                foreach (qEdge edge in startElement.EdgeList)
                {
                    if (quadEdge.Contains(edge)) { continue; }

                    List<qElement> connectedElements = edge.GetConnectedElements();
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
                        globalEdgeList.Remove(edge);
                    }
                }
                globalElementList.Remove(element);
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

            int newLevel = quadEdge[0].Level + 1;
            //foreach (qEdge edge in quadEdge) { if (edge.Level > maxLevel) { maxLevel = edge.Level; } } // to do: Silje, slett?       

            // update frontEdges
            foreach (qEdge edge in quadEdge)
            {
                if (edge.IsFrontEdge() & !frontEdges.Contains(edge))
                {
                    frontEdges.Add(edge);
                    int indexToUpdate = globalEdgeList.IndexOf(edge);
                    globalEdgeList[indexToUpdate].Level = newLevel;
                }
                else if (frontEdges.Contains(edge))
                {
                    frontEdges.Remove(edge); // remove edge from front edge if not a front edge, but is in the frontEdges list
                }
                else
                { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "CreateQuad: problems with front edge criterion or is closing front."); } // to do: temporary check 
            }
            SetNeighorFrontEdges(frontEdges);

            return newQuadElement;
        }
        private void CleanUpChevorns(List<qEdge> globalEdgeList, List<qElement> globalElementList, List<qEdge> frontEdges)
        {
            // summary: clean up chevrons if there are any. Assume only quads
            foreach (qElement element in globalElementList)
            {

                if (element.IsChevron())
                {
                    DoLocalSmoothing(element, globalEdgeList, frontEdges, globalElementList);
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "CleanUpChevron: fixed chevron");
                }
          

                /*
                double maxAngle = 0;
                int id = 0;
                double chevronConstant = (double)200 / (double)180 * Math.PI;
                List<qEdge> edgeList = element.EdgeList;
                List<double> angleList = element.AngleList;
                for (int i = 0; i < angleList.Count ; i++)
                {
                    if (angleList[i] > maxAngle) { maxAngle = angleList[i] ; id = i; }
                }

                //
                if (maxAngle > chevronConstant)
                {
                    qNode concavNode = element.GetNodesOfElement()[id]; // assume same order 
                    bool laplaceSmoothFailed = false;
                    if (!concavNode.BoundaryNode)
                    {
                        //DoLocalSmoothing()

                        // check if ok
                        foreach (double angle in element.AngleList)
                        {
                            if (angle > chevronConstant) { laplaceSmoothFailed = true; }
                        }

                    }
                 

                    if (concavNode.BoundaryNode | laplaceSmoothFailed)
                    {
                        // fill
                        List<qEdge> conncectedEdgesOfconcaveNode = concavNode.GetConnectedEdges(globalEdgeList);
                        List<qEdge> edgeCandidates = new List<qEdge>();
                        List<qElement> elementCandidates = new List<qElement>();
                        foreach (qEdge edge in element.EdgeList)
                        {
                            if (!conncectedEdgesOfconcaveNode.Contains(edge)) 
                            {
                                edgeCandidates.Add(edge);
                                if (edge.Element1 == element) { elementCandidates.Add(edge.Element2); }
                                else { elementCandidates.Add(edge.Element1); }
                            }
                        }

                        // check node valence
                        // temporary: select first element
                        qEdge edgeToUse = edgeCandidates[0];
                        qElement elementToUse = elementCandidates[0];
                        
                        // find n1
                        qNode n1 = new qNode();
                        if (edgeToUse.StartNode == edgeCandidates[1].StartNode | edgeToUse.StartNode == edgeToUse.EndNode)
                        { n1 = edgeToUse.EndNode; }
                        else { n1 = edgeToUse.StartNode; }

                        // check if fill_4
                        bool fill_3 = true;
                        foreach (qEdge edge2 in elementToUse.EdgeList)
                        {
                            if ((edge2.StartNode == concavNode | edge2.EndNode == concavNode) & edge2 != edgeToUse)
                            {
                                qEdge edge1 = FindEdge(globalEdgeList, n1, concavNode);
                                var vectors = edge1.CalculateVectorsFromSharedNode(edge2);
                                Vector3d vec1 = vectors.Item1;
                                Vector3d vec2 = vectors.Item2;
                                double angleToCheck = 0;
                                if (Vector3d.CrossProduct(vec1, vec2).Z >= 0)
                                {
                                    angleToCheck =Vector3d.VectorAngle(vec1, vec2, Vector3d.ZAxis); // to do: make more general
                                }
                                else
                                {
                                    angleToCheck = Vector3d.VectorAngle(vec2, vec1, Vector3d.ZAxis); // to do: make more general
                                }
                                if (angleToCheck > chevronConstant) { fill_3 = false; }
                                break;
                            }
                        }

                        if (fill_3)
                        {
                            qNode n5 = new qNode();
                            qNode n4 = new qNode();
                            qNode n3 = edgeToUse.GetOppositeNode(n1);
                            qEdge n1n5 = new qEdge();
                            qEdge n3n4 = new qEdge();


                            foreach (qEdge edge in elementToUse.EdgeList)
                            {
                                if ((edge.StartNode == n3 | edge.EndNode == n3) & edge != edgeToUse)
                                {
                                    n3n4 = edge;
                                    n4 = edge.GetOppositeNode(n3);
                                }
                                if ((edge.StartNode == n1 | edge.EndNode == n1) & edge != edgeToUse)
                                {
                                    n1n5 = edge;
                                    n5 = edge.GetOppositeNode(n1);
                                }

                            }

                            // new node
                            qNode n6 = new qNode(0.5 * (n4.Coordinate + n1.Coordinate), false); // assume new node is given by this

                            // edges to keep
                            List<qEdge> edgesToKeep = new List<qEdge>(elementToUse.EdgeList);
                            edgesToKeep.AddRange(element.EdgeList);
                            edgesToKeep.Remove(edgeToUse);

                            // new edges
                            qEdge n6c = new qEdge(n6, concavNode);
                            qEdge n6n5 = new qEdge(n6, n5);
                            qEdge n6n3 = new qEdge(n6, n3);
                            
                            // new elements
                           

                        }

                    }

                }*/
            }
            return;
        }// to do: finish this
        private qEdge CloseFront(List<qEdge> frontEdges, List<qEdge> globalEdgeList, List<qElement> globalElementList, qEdge edgeToClose, qNode N_k)
        {
            qEdge E_k = new qEdge();

            if (!IsFrontLoopsEven(frontEdges, edgeToClose, globalEdgeList).Item1) // potential loops are not even loops, need to split
            {
                // split the edge because not an even number edges in loops
                qEdge edgeToSplit = edgeToClose;
                List<qNode> nodes = edgeToClose.Element1.GetNodesOfElement();
                qNode newNode = new qNode(0.5 * (edgeToSplit.StartNode.Coordinate + edgeToSplit.EndNode.Coordinate), false);

                qNode nodeToSpitFrom = new qNode();
                foreach (qNode node in nodes)
                {
                    if (node != edgeToSplit.StartNode & node != edgeToSplit.EndNode) { nodeToSpitFrom = node; break; }
                }
                qEdge E_temp = SplitEdge(edgeToSplit, newNode.Coordinate - nodeToSpitFrom.Coordinate, nodeToSpitFrom, globalEdgeList, globalElementList);
                E_k = FindEdge(globalEdgeList, N_k, E_temp.GetOppositeNode(nodeToSpitFrom));
            }
            else
            {
                E_k = edgeToClose; // use existing edge as side edge in closing
            }
            return E_k;
        }
        
        // __________________________________________ Local smoothing ______________________________________________________
        private void DoLocalSmoothing(qElement quadElement, List<qEdge> globalEdgeList, List<qEdge> frontEdges, List<qElement> globalElementList)
        {
            
            Point3d newCoordinate = new Point3d();

            // Smooth nodes of Quad
            List<qNode> quadNodes = quadElement.GetNodesOfElement();
            List<qEdge> globalEdgeListOld = new List<qEdge>(globalEdgeList);

            foreach (qNode node in quadNodes) //assume order is irrelevant
            {
                qNode oldNode = new qNode(node.Coordinate, node.BoundaryNode);
                bool isFrontNode = node.IsFrontNode(frontEdges);
                
                if (isFrontNode & !node.BoundaryNode)
                {
                    // smooth front node
                    newCoordinate = FrontNodeSmooth(node, globalEdgeListOld);
                }
                else if (!isFrontNode & !node.BoundaryNode)
                {
                    // do laplacian smooth
                    newCoordinate = ModifiedLengthWeightedLaplacianSmooth(node, globalEdgeListOld);
                }
                else
                {
                    // assume that node is not moved
                    newCoordinate = node.Coordinate;
                    continue;
                }
                
                var changedEdges = UpdateGlobalEdgeList_NodePosition(node, newCoordinate, globalEdgeList);
                List<qEdge> oldEdges = changedEdges.Item1;
                List<qEdge> newEdges = changedEdges.Item2;
                UpdateGlobalElementList_ChangedEdges(newEdges, oldEdges, globalElementList);

                List<qElement> connectedTrinagles = node.GetTrianglesConnectedToNode(globalEdgeList);
                if (connectedTrinagles.Count == 0) { continue; }
                
                Vector3d movingVector = oldNode.Coordinate - newCoordinate;
                var newCoordinateInverted = InvertedElementsCleanUp(node, connectedTrinagles, movingVector);
                if (newCoordinateInverted != newCoordinate)
                {
                    changedEdges = UpdateGlobalEdgeList_NodePosition(node, newCoordinateInverted, globalEdgeList);
                    oldEdges = changedEdges.Item1;
                    newEdges = changedEdges.Item2;
                    UpdateGlobalElementList_ChangedEdges(newEdges, oldEdges, globalElementList);
                }
            }

            // Smooth adjacent nodes
            List<qNode> adjacentNodes = GetNeighborNodesToElement(quadElement, globalEdgeList);
            List<qEdge> globalEdgeListOldUpdated = new List<qEdge>(globalEdgeList);
            
            foreach (qNode adjNode in adjacentNodes)
            {
                qNode oldAdjNode = new qNode(adjNode.Coordinate, adjNode.BoundaryNode);
                bool isFrontNode = adjNode.IsFrontNode(frontEdges);
                if (isFrontNode & !adjNode.BoundaryNode)
                {
                    // smooth front node
                    newCoordinate = FrontNodeSmooth(adjNode, globalEdgeListOldUpdated);
                }
                else if (!isFrontNode & !adjNode.BoundaryNode)
                {
                    // do laplacian smooth
                    newCoordinate = ModifiedLengthWeightedLaplacianSmooth(adjNode, globalEdgeListOldUpdated);
                }
                else
                {
                    // assume that node is not moved
                    newCoordinate = adjNode.Coordinate;
                    continue;
                }
                var changedEdges = UpdateGlobalEdgeList_NodePosition(adjNode, newCoordinate, globalEdgeList);
                List<qEdge> oldEdges = changedEdges.Item1;
                List<qEdge> newEdges = changedEdges.Item2;
                UpdateGlobalElementList_ChangedEdges(newEdges, oldEdges, globalElementList);

                List<qElement> connectedTrinagles = adjNode.GetTrianglesConnectedToNode(globalEdgeList);
                if (connectedTrinagles.Count == 0) { continue; }

                Vector3d movingVector = oldAdjNode.Coordinate - adjNode.Coordinate;
                Point3d newCoordinateInverted = InvertedElementsCleanUp(adjNode, connectedTrinagles, movingVector);
                if (newCoordinateInverted != newCoordinate)
                {
                    changedEdges = UpdateGlobalEdgeList_NodePosition(adjNode, newCoordinateInverted, globalEdgeList);
                    oldEdges = changedEdges.Item1;
                    newEdges = changedEdges.Item2;
                    UpdateGlobalElementList_ChangedEdges(newEdges, oldEdges, globalElementList);
                }
            }
        } 

        private Point3d InvertedElementsCleanUp(qNode smoothNode, List<qElement> connectedTriangles, Vector3d movingVector)
        {
            Point3d newCoordinate = new Point3d(smoothNode.Coordinate);
            for (int i = 0; i < connectedTriangles.Count; i++)
            {
                // to do: dummy solution
                qElement triangle = new qElement();
                triangle.EdgeList = connectedTriangles[i].EdgeList;
                triangle.CalculateAngles();
                triangle.IsQuad = false;
                triangle.GetContourOfElement();

                int counter = 0;
                while (triangle.IsInverted() & counter < 1000)
                {
                    newCoordinate = new Point3d(newCoordinate.X + movingVector.X * 0.25, newCoordinate.Y + movingVector.Y * 0.25, newCoordinate.Z + movingVector.Z * 0.25);
                    for (int j = 0; j < triangle.EdgeList.Count; j++)
                    {
                        qEdge edge = triangle.EdgeList[j];
                        if (edge.StartNode == smoothNode) { triangle.EdgeList[j].StartNode.Coordinate = newCoordinate; }
                        else if (edge.EndNode == smoothNode) { triangle.EdgeList[j].EndNode.Coordinate = newCoordinate; }
                    }
                    smoothNode.Coordinate = newCoordinate;
                    counter++;
                }
                if (counter >= 1000) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "InvertedElementsCleanUp: failed"); }
            }
            return newCoordinate;
        }
        private Tuple<List<qEdge>,List<qEdge>> UpdateGlobalEdgeList_NodePosition(qNode smoothNode, Point3d newCoordinate, List<qEdge> globalEdgeList)
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
                    int edgeId2 = globalEdgeListCopy[id].Element2.EdgeList.IndexOf(edge);
                    globalEdgeList[id].Element1.EdgeList[edgeId1] = globalEdgeList[id];
                    
                    globalEdgeList[id].Element1.GetContourOfElement();
                    globalEdgeList[id].Element1.CalculateAngles(); // to do:isquad?


                    globalEdgeList[id].Element2.EdgeList[edgeId2] = globalEdgeList[id];
                    globalEdgeList[id].Element2.GetContourOfElement();
                    globalEdgeList[id].Element2.CalculateAngles(); // to do: isquad?


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
        private List<qNode> GetNeighborNodesToElement(qElement element, List<qEdge> globalEdgeList) 
        {
            // summary: get neighbor nodes of element
            List<qNode> adjacentNodes = new List<qNode>();

            List<qNode> elementNodes = element.GetNodesOfElement();
            List<qEdge> elementEdges = element.EdgeList;

            // silje comment:
            /*
             foreach(qNode ndoe in elementNodes)
            {
                 List<qEdge> connectedEdges = GetConnectedEdges(node, globalEdgeList);
                 foreach (qEdge edge in connectedEdges)
                 {
                     if (!elementEdges.Contains(edge))
                     {
                         qNode adjacentNode = GetOppositeNode(node, edge);
                         adjacentNodes.Add(adjacentNode);
                     }
                 }
            }
            */
            for (int i = 0; i < elementNodes.Count; i++)
            {
                qNode node = elementNodes[i];
                List<qEdge> connectedEdges = node.GetConnectedEdges(globalEdgeList);

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
                        qNode adjacentNode = edge.GetOppositeNode(node);
                        adjacentNodes.Add(adjacentNode);
                    }
                }
            }
            return adjacentNodes;
        }  
        private Point3d ModifiedLengthWeightedLaplacianSmooth(qNode Ni, List<qEdge> globalEdgeList) // todo: check if this works. Check if each method osv do what I want
        {
            Vector3d lengthCjVectorCj = Vector3d.Zero;
            double lengthCj = 0;

            List<qEdge> connectedEdgesToNi = Ni.GetConnectedEdges(globalEdgeList); // todo: check if it finds connected edges

            for (int i = 0; i < connectedEdgesToNi.Count; i++)
            {
                qEdge edge = connectedEdgesToNi[i];
                qNode Nj = edge.GetOppositeNode(Ni);
                Vector3d Vj = Nj.Coordinate - Ni.Coordinate;
                Vector3d Cj = Vector3d.Zero;

                if (Nj.BoundaryNode)
                {
                    Vector3d deltaCj = GetAngularSmoothness(Ni, Nj, (Ni.Coordinate-Nj.Coordinate).Length, false, globalEdgeList); //to do: fix angular smoothness
                    Cj = Vj + deltaCj;
                }
                else
                {
                    Cj = Vj;
                }

                lengthCj = Cj.Length + lengthCj;
                lengthCjVectorCj = Cj.Length * Cj + lengthCjVectorCj; // todo: check if last is OK
            }

            Vector3d delta = lengthCjVectorCj / lengthCj;
            Point3d smoothNode = new Point3d(Ni.Coordinate.X + delta.X, Ni.Coordinate.Y + delta.Y, Ni.Coordinate.Z + delta.Z);

            return smoothNode;
        }
        private qEdge GetSharedEdge(List<qElement> quadElements)
        {
            // summary: get shared edge from adjacent elements
            // silje comment: endre navn til GetSharedEdgeFromAdjacentQuads?
            qEdge sharedEdge = new qEdge();
            if (quadElements.Count == 2)
            {
                qElement quadElement1 = quadElements[0];
                qElement quadElement2 = quadElements[1];

                foreach (qEdge edge in quadElement1.EdgeList)
                {
                    if (quadElement2.EdgeList.Contains(edge))
                    {
                        sharedEdge = edge;
                    }
                }
            }
            if (sharedEdge.Length == 0) { sharedEdge = null; }

            return sharedEdge;
        } 
        private Point3d FrontNodeSmooth(qNode Ni, List<qEdge> globalEdgeList)
        {
            // summary: Get new coordinate of front node to be smooth
            Point3d origo = new Point3d(0, 0, 0);

            // Find front edges connected to front node: todo: test if this works as I want.
            //List<qEdge> nodeFrontEdges = GetFrontEdgesConnectedToNode(Ni, globalEdgeList); // to do: hilde slett?
            List<qEdge> connectedEdges = Ni.GetConnectedEdges(globalEdgeList);

            // Get quads that are connected to node. Todo: check if this do what I want.
            List<qElement> quadElements = Ni.GetQuadsConnectedToNode(globalEdgeList);
            int numberOfConnectedQuads = quadElements.Count;
            Vector3d deltaI = Vector3d.Zero;

            if (numberOfConnectedQuads != 0)
            {
                // True isoparametric smooth:
                var vectors = TrueIsoparametricSmooth(Ni, quadElements);
                Vector3d Vi_mark = vectors.Item1;
                Vector3d deltaA = vectors.Item2;
                
                Vector3d deltaB = Vector3d.Zero;
                Vector3d deltaC = Vector3d.Zero;
                qEdge sharedEdge = GetSharedEdge(quadElements);

                if (numberOfConnectedQuads == 2) 
                {
                    if (sharedEdge == null) { deltaI = deltaA; AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "FrontNodeSmooth: No shared edge between the two quads."); }
                    else if (sharedEdge.IsFrontEdge()) { deltaI = deltaA; AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "FrontNodeSmooth: Shared edge in Front Node smooth (two quads) is a front edge.."); }
                    else
                    {
                        // Length modification:
                        qNode Nj = sharedEdge.GetOppositeNode(Ni);
                        Vector3d Vi = Ni.Coordinate - origo;
                        Vector3d Vj = Nj.Coordinate - origo;

                        double la = (Vi_mark - Vj).Length;

                        // Find ld
                        double ld = 0;
                        double tr = 3; //todo: fix tr!!
                        double sumSurroundingEdges = 0;
                        double n = 0;
                        if (tr <= 2.5)
                        {
                            ld = (Ni.Coordinate - Nj.Coordinate).Length; // todo: check assumption
                        }
                        else if (tr > 2.5 & tr <= 20)
                        {
                            // silje comment:  Gjør koden mindre: Ta ut edge-lister fra elementene og connectedEdges. Slå sammen til en list. Loop foreach edge add dersom edge ikke er sharedEdge eller front?.
                            /*
                                List<qEdge> edgesToLoop = new List<qEdge>();
                                edgesToLoop.AddRange(quadElements[0].EdgeList);
                                edgesToLoop.AddRange(quadElements[1].EdgeList);
                                edgesToLoop.AddRange(connectedEdges);

                            foreach(qEdge edge in edgesToLoop)
                            {
                                (!IsFrontEdge(edge) & edge != sharedEdge) { sumSurroundingEdges = edge1.Length + sumSurroundingEdges; }
                                if (connectedEdges.Contain(edge)) {n++;}
                            }
                            */
                            for (int i = 0; i < 4; i++)
                            {
                                qEdge edge1 = quadElements[0].EdgeList[i];
                                qEdge edge2 = quadElements[1].EdgeList[i];

                                // length of edges of element1
                                if (!edge1.IsFrontEdge() & edge1 != sharedEdge) { sumSurroundingEdges = edge1.Length + sumSurroundingEdges; }

                                // length of edges of element2
                                if (!edge2.IsFrontEdge() & edge2 != sharedEdge) { sumSurroundingEdges = edge2.Length + sumSurroundingEdges; }
                            }

                            // length of edges not connected to elements or is frontedge
                            foreach (qEdge edge in connectedEdges)
                            {
                                if (!edge.IsFrontEdge() & edge != sharedEdge)
                                {
                                    sumSurroundingEdges = edge.Length + sumSurroundingEdges;
                                    n++;
                                }
                            }

                            ld = sumSurroundingEdges / ((double)4 + n);
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

                        // Angle modification:
                        deltaB = Vj - Vi + (deltaA + Vi - Vj) * ld / la;
                        deltaC = GetAngularSmoothness(Ni, Nj, ld, true, globalEdgeList);
                        deltaI = (deltaB + deltaC) / (double)2;
                    }
                }
                else { deltaI = deltaA; }
            }
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "FrontNodeSmooth: numberOfConnectedQuads is zero. Blacker smooth can not be performed!!"); }
            
            Point3d smoothNode = new Point3d(Ni.Coordinate.X + deltaI.X, Ni.Coordinate.Y + deltaI.Y, Ni.Coordinate.Z + deltaI.Z);
            return smoothNode;
        } // todo: check if this is OK
        private Tuple<Vector3d,Vector3d> TrueIsoparametricSmooth(qNode Ni, List<qElement> quadElements)
        {
            // summary: get translations vectors to perfome true isoparametric smoothing
            Vector3d Vi = Vector3d.Zero;
            Vector3d Vi_mark = Vector3d.Zero;
            Vector3d deltaA = Vector3d.Zero;
            Vector3d vectorSum = Vector3d.Zero;
            Point3d origo = new Point3d(0, 0, 0);
            int numberOfConnectedQuads = quadElements.Count;

            foreach (qElement quadElement in quadElements)
            {
                List<qNode> quadNodes = quadElement.GetNodesOfElement();
                qNode node1 = quadNodes[0];
                qNode node2 = quadNodes[1];
                qNode node3 = quadNodes[2];
                qNode node4 = quadNodes[3];
                Vector3d Vmj = Vector3d.Zero;
                Vector3d Vmk = Vector3d.Zero;
                Vector3d Vml = Vector3d.Zero;

                if (node1 == Ni)
                {
                    Vmj = node2.Coordinate - origo;
                    Vmk = node3.Coordinate - origo;
                    Vml = node4.Coordinate - origo;
                }
                else if (node2 == Ni)
                {
                    Vmj = node3.Coordinate - origo;
                    Vmk = node4.Coordinate - origo;
                    Vml = node1.Coordinate - origo;
                }
                else if (node3 == Ni)
                {
                    Vmj = node4.Coordinate - origo;
                    Vmk = node1.Coordinate - origo;
                    Vml = node2.Coordinate - origo;
                }
                else if (node4 == Ni)
                {
                    Vmj = node1.Coordinate - origo;
                    Vmk = node2.Coordinate - origo;
                    Vml = node3.Coordinate - origo;
                }

                vectorSum = (Vmj + Vml - Vmk) + vectorSum;
            }
            Vi = Ni.Coordinate - origo;
            Vi_mark = ((double)1 / (double)numberOfConnectedQuads) * vectorSum;
            deltaA =  Vi_mark - Vi;
            return Tuple.Create(Vi_mark, deltaA);
        } // todo: check if this is OK
        private Vector3d GetAngularSmoothness(qNode Ni, qNode Nj, double ld, bool Ni_IsFrontNode, List<qEdge> globalEdgeList) 
        {
            // summary: get a translation vector to performe angular smoothing
            Vector3d P_B1 = Vector3d.Zero;
            Vector3d P_B2 = Vector3d.Zero;
            Vector3d Pi = Ni.Coordinate - Nj.Coordinate;
            Point3d pointQ = new Point3d();
            Vector3d deltaC = Vector3d.Zero;

            List<qElement> quadElements = new List<qElement>();
            
            if (!Ni_IsFrontNode)
            {
                quadElements = Nj.GetQuadsConnectedToNode(globalEdgeList);
            }
            else
            {
                quadElements = Ni.GetQuadsConnectedToNode(globalEdgeList);
            }

            if (quadElements.Count == 2)
            {
                qEdge topEdge1 = new qEdge();
                qEdge topEdge2 = new qEdge();
                qEdge sharedEdge = GetSharedEdge(quadElements);

                if (sharedEdge == null) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "AngularSmoothing: shared edge = null"); }
                else if (!sharedEdge.IsFrontEdge())
                {
                    for (int i = 0; i < 4; i++)
                    {
                        qEdge edge1 = quadElements[0].EdgeList[i];
                        qEdge edge2 = quadElements[1].EdgeList[i];
                        if (edge1 != sharedEdge & (edge1.StartNode == Ni | edge1.EndNode == Ni)) { topEdge1 = edge1; }
                        if (edge2 != sharedEdge & (edge2.StartNode == Ni | edge2.EndNode == Ni)) { topEdge2 = edge2; }
                    }

                    /*
                    List<qNode> node1 = GetNodesOfElement(quadElements[0]);
                    
                    List<qNode> testNodes = GetNodesOfElement(quadElements[0]);
                    int testNodeIndex = testNodes.IndexOf(Nj) + 2;
                    if (testNodeIndex > 3) { testNodeIndex = testNodeIndex - 4; }
                    

                    if (GetOppositeNode(Ni, topEdge1) == node1[3]) //testNodes[testNodeIndex]
                    {
                        Vector3d vectorLeft = GetOppositeNode(Ni, topEdge1).Coordinate - Nj.Coordinate;
                        Vector3d vectorRight = GetOppositeNode(Ni, topEdge2).Coordinate - Nj.Coordinate;
                        P_B1 = GetBisectingVector(vectorRight, vectorLeft);
                    }
                    else
                    {
                        Vector3d vectorLeft = GetOppositeNode(Ni, topEdge2).Coordinate - Nj.Coordinate;
                        Vector3d vectorRight = GetOppositeNode(Ni, topEdge1).Coordinate - Nj.Coordinate;
                        P_B1 = GetBisectingVector(vectorRight, vectorLeft);
                    }
                    */
                    Vector3d vectorLeft = topEdge1.GetOppositeNode(Ni).Coordinate - Nj.Coordinate;
                    Vector3d vectorRight = topEdge2.GetOppositeNode(Ni).Coordinate - Nj.Coordinate;
                    if (Vector3d.VectorAngle(vectorLeft, vectorRight, Vector3d.ZAxis) < Math.PI) // todo: make for not planar 
                    {
                        vectorRight = topEdge1.GetOppositeNode(Ni).Coordinate - Nj.Coordinate;
                        vectorLeft = topEdge2.GetOppositeNode(Ni).Coordinate - Nj.Coordinate;
                    }
                    P_B1 = GetBisectingVector(vectorRight, vectorLeft);
                    P_B2 = (double)P_B1.Length * Pi + (double)Pi.Length * P_B1; // Assume angle always less than 180 degree.
                    P_B2.Unitize();
                    if (Vector3d.Multiply(P_B1, Pi) < 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "P_B1 and Pi have angle larger than 90 degree. Do not know it is above 180."); }

                    NurbsCurve line1 = new Line(Nj.Coordinate, P_B2, 100).ToNurbsCurve();
                    NurbsCurve line2 = new Line(topEdge1.GetOppositeNode(Ni).Coordinate, topEdge2.GetOppositeNode(Ni).Coordinate).ToNurbsCurve();
                    var placesWithIntersection = Intersection.CurveCurve(line1, line2, 0.0001, 0.0001);
                    bool isIntersecting = true;
                    if (placesWithIntersection.Count == 0)
                    {
                        isIntersecting = false;
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "AngularSmooth: Vector P_B2 do not intersect with line Ni-1 to Ni+1.");
                    }
                    else
                    {
                        pointQ = placesWithIntersection[0].PointA;
                    }

                    // make P_B2:
                    double lq = (pointQ - Nj.Coordinate).Length;
                    if (ld > lq & isIntersecting)
                    {
                        double lengthP_B2 = (lq + ld) / (double)2;
                        P_B2 = P_B2 * lengthP_B2;
                    }
                    else
                    {
                        double lengthP_B2 = ld;
                        P_B2 = P_B2 * lengthP_B2;
                    }

                    deltaC = P_B2 - Pi;
                }
            }
            return deltaC;
        } 
        private Vector3d GetBisectingVector(Vector3d VectorRight, Vector3d VectorLeft)
            {
                // summary: get bisecting vector, assuming the source vectors are not parallel 
                double angle = Vector3d.VectorAngle(VectorRight, VectorLeft, Vector3d.ZAxis); // todo: make normal mor general
                Vector3d bisectVector = VectorRight;
                bisectVector.Rotate(0.5 * angle, Vector3d.ZAxis); // todo: not for 3d surface, make axis for normal to plane (vec1, vec2)
                bisectVector.Unitize();
                return bisectVector;
            } // todo: test if this is OK
        private List<qEdge> GetFrontEdgesConnectedToNode(qNode node, List<qEdge> globalEdgeList)
            {
                // summary: get front edges connected to a node
                List<qEdge> nodeFrontEdges = new List<qEdge>();
                qEdge nodeLeftFront = new qEdge();
                qEdge nodeRightFront = new qEdge();

                List<qEdge> connectedEdges = node.GetConnectedEdges(globalEdgeList);

                foreach (qEdge edge in connectedEdges)
                {
                    if (edge.IsFrontEdge()) { nodeFrontEdges.Add(edge); }
                }

                if (nodeFrontEdges.Count != 2) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, String.Format("GetFrontEdgesConnectedToNode: Node is connected to {0} front edges", nodeFrontEdges.Count)); } // to do: fjerne denne ? closing front er et tilfelle der dette er aktuelt
                else
                {
                    qEdge front1 = nodeFrontEdges[0];
                    qEdge front2 = nodeFrontEdges[1];

                    if (front1.RightFrontNeighbor.StartNode == node | front1.RightFrontNeighbor.EndNode == node)
                    {
                        nodeLeftFront = front1; nodeRightFront = front2;
                    }
                    else
                    {
                        nodeLeftFront = front2; nodeRightFront = front1;
                    }
                }
                nodeFrontEdges.Clear();
                nodeFrontEdges.Add(nodeLeftFront);
                nodeFrontEdges.Add(nodeRightFront);

                return nodeFrontEdges;
        } //todo: test if this is OK

        // __________________________________________ Global smoothing ______________________________________________________
        
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

        // After paper by Cannan, with modifications by Qmorp-fyren for the quads.
        private void DoGlobalSmoothing(List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            List<qNode> globalNodeList = GetGlobalNodeList(globalEdgeList);
            List<qNode> smoothableNodes = new List<qNode>(globalNodeList);
            double moveTolerance = 0.01; // constant: fix: todo.
            double maxDistanceMoved = 0; // todo: fix this so that it is bigger than 1.75*moveTol
            double maxModelDimension = 0;

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
            while (continueSmooth)
            {
                bool nodeIsMoved = false;
                maxDistanceMoved = 0;
                List<qNode> smoothableNodesCopy = new List<qNode>(smoothableNodes);
                foreach (qNode node in smoothableNodesCopy)
                {
                    if (node.BoundaryNode | !smoothableNodesCopy.Contains(node)) { smoothableNodes.Remove(node); continue; } // todo: forlag: ta BC-edge og flytt vektor.
                    if (iterations == 1) { node.OBS = false; }

                    // 4.1. Perform Constrained Laplacian Smooth
                    if (!node.OBS)
                    {
                        // 4.1.1. Move node
                        qNode movedNode = ConstrainedLaplacianSmooth(node, globalEdgeList, globalElementList);
                        double distanceMoved = (movedNode.Coordinate - node.Coordinate).Length;

                        // 4.1.2. Check if move shall be allowed;
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
                    }
                }
                iterations++;
                if (!nodeIsMoved | maxDistanceMoved <= 1.75 * moveTolerance) { continueSmooth = false; }
                n++;
            }
        }
        private double CalculateDistortionMetric(qElement element) //OK  move to qElement
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
        } // OK
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
        } // OK
        private qNode ConstrainedLaplacianSmooth(qNode node, List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            /*// 0. Defining variables used to check if move is acceptable:
            int posN = 0;
            int negN = 0;
            int upN = 0;
            int downN = 0;
            int invN = 0;
            double theta = 0;
            double thetaMax = 200 * Math.PI / 180; // constant chosen
            double deltaDistMetric = 0;*/

            // 1. Move node with Laplacian smooth
            Vector3d laplacianVector = LaplacianSmooth(node, globalEdgeList);
            Point3d newPoint = new Point3d(node.Coordinate.X + laplacianVector.X, node.Coordinate.Y + laplacianVector.Y, node.Coordinate.Z + laplacianVector.Z);
            qNode newNode = new qNode(newPoint, node.BoundaryNode);
            qNode oldNode = new qNode(node.Coordinate, node.BoundaryNode);

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
                //var update = UpdateGlobalEdgeList_NodePosition(node, newNode.Coordinate, globalEdgeList);
                //UpdateGlobalElementList_ChangedEdges(update.Item1, update.Item2, globalElementList);
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
                    //if (i == 19) { newNode = node; }
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
                Vector3d vector = edge.GetOppositeNode(node).Coordinate - node.Coordinate;
                vectorSum = vectorSum + vector;
            }
            return vectorSum;
        } // OK
        private qNode OptimizationBasedSmoothing(qNode node, double maxModelDimension, List<qEdge> globalEdgeList)
        {
            qNode newNode = new qNode(node.Coordinate, node.BoundaryNode);
            double myMin = 100; // dummy-value
            Vector3d g = Vector3d.Zero;
            List<Vector3d> gi = new List<Vector3d>();
            List<qElement> connectedElements = node.GetConnectedElements(globalEdgeList);

            // 1. Estimate gradient vector for each element connected to node
            double delta = Math.Pow(10, -5) * maxModelDimension; // constant form paper

            foreach (qElement element in connectedElements)
            {
                double giX = CalculateGradient(element, node, delta, "x");
                double giY = CalculateGradient(element, node, delta, "y");
                double giZ = CalculateGradient(element, node, delta, "z");
                if ((new Vector3d(giX, giY, giZ)).Length > 0.00001 & element.DistortionMetric < myMin)
                {
                    myMin = element.DistortionMetric;
                    g = new Vector3d(giX, giY, giZ);
                }
                gi.Add(new Vector3d(giX, giY, giZ));
            }

            // 2. Calculate gamma used to move node
            double gamma = 100;
            bool gammaLimited = false;
            for (int i = 0; i < connectedElements.Count; i++)
            {
                if (Vector3d.Multiply(g, gi[i]) < 0) 
                {
                    gammaLimited = true;
                    double gamma_i = (connectedElements[i].DistortionMetric - myMin) / (Vector3d.Multiply(g, g) - Vector3d.Multiply(g, gi[i]));
                    if (gamma_i < gamma) { gamma = gamma_i; }
                }
            }
            if (!gammaLimited) { gamma = 0.8; } // Constant as in qmorph-fyr

            // 3. Move node:
            for (int i = 0; i <= 4; i++) // Constant "4" as proposed in paper
            {
                Point3d newPoint = node.Coordinate + gamma * g;
                double myMinNew = 100;

                foreach (qElement element in connectedElements)
                {
                    List<qNode> elementNodes = element.GetNodesOfElement();
                    elementNodes[elementNodes.IndexOf(node)].Coordinate = newPoint;
                    qElement newElement = CreateElementFromNodes(elementNodes);
                    double my = CalculateDistortionMetric(newElement);
                    if (my < myMinNew) { myMinNew = my; }
                }

                if (myMinNew >= myMin + 0.0001) // Constant as proposed in paper
                {
                    newNode.Coordinate = newPoint;
                    newNode.OBS = true;
                    break;
                }
                else { gamma = gamma / 2; newNode.OBS = false; }
            }
            return newNode;
        }
        private double CalculateGradient(qElement element, qNode node, double delta, string direction)
        {
            if (direction == "x")
            {
                List<qNode> elementNodes = element.GetNodesOfElement();
                Point3d point = new Point3d(node.Coordinate.X + delta, node.Coordinate.Y, node.Coordinate.Z);
                elementNodes[elementNodes.IndexOf(node)].Coordinate = point;
                qElement newElement = CreateElementFromNodes(elementNodes);
                double myPertubed = CalculateDistortionMetric(newElement); // pertubed distortion metric. 
                double gi = (myPertubed - element.DistortionMetric) / delta;
                return gi;
            }
            else if (direction == "y")
            {
                List<qNode> elementNodes = element.GetNodesOfElement();
                Point3d point = new Point3d(node.Coordinate.X, node.Coordinate.Y + delta, node.Coordinate.Z);
                elementNodes[elementNodes.IndexOf(node)].Coordinate = point;
                qElement newElement = CreateElementFromNodes(elementNodes);
                double myPertubed = CalculateDistortionMetric(newElement); // pertubed distortion metric. 
                double gi = (myPertubed - element.DistortionMetric) / delta;
                return gi;
            }
            else if (direction == "z")
            {
                List<qNode> elementNodes = element.GetNodesOfElement();
                Point3d point = new Point3d(node.Coordinate.X, node.Coordinate.Y, node.Coordinate.Z + delta);
                elementNodes[elementNodes.IndexOf(node)].Coordinate = point;
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
                qEdge edge2 = new qEdge(nodes[1], nodes[2]);
                qEdge edge3 = new qEdge(nodes[2], nodes[0]);
                List<qEdge> edgeList = new List<qEdge>() { edge1, edge2, edge3 };
                qElement element = new qElement(edgeList);
                return element;
            }
        } // OK
        #endregion
        /*
         * private qNode ConstrainedLaplacianSmooth(qNode node, List<qEdge> globalEdgeList, List<qElement> globalElementList)
        {
            // 0. Defining variables used to check if move is acceptable:
            int posN = 0;
            int negN = 0;
            int upN = 0;
            int downN = 0;
            int invN = 0;
            double theta = 0;
            double thetaMax = 200 * Math.PI / 180; // constant chosen
            double deltaDistMetric = 0;

            // 1. Move node with Laplacian smooth
            Vector3d laplacianVector = LaplacianSmooth(node, globalEdgeList);
            Point3d newPoint = new Point3d(node.Coordinate.X + laplacianVector.X, node.Coordinate.Y + laplacianVector.Y, node.Coordinate.Z + laplacianVector.Z);
            qNode newNode = new qNode(newPoint, node.BoundaryNode);
            qNode oldNode = new qNode(node.Coordinate, node.BoundaryNode);

            // 2. Loop to find final node position
            for (int i = 0; i < 20; i++) // constant proposed in paper: 20
            {
                // 2.1. Update edge and element list
                var update = UpdateGlobalEdgeList_NodePosition(node, newNode.Coordinate, globalEdgeList);
                UpdateGlobalElementList_ChangedEdges(update.Item1, update.Item2, globalElementList);
                List<qElement> connectedElements = oldNode.GetConnectedElements(globalEdgeList);
                int N = connectedElements.Count;

                // 2.2. Calculate Acceptance Criteria for each element 
                foreach (qElement element in connectedElements)
                {
                    double newDisMetric = CalculateDistortionMetric(element);
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
                    if (element.IsInverted())
                    {
                        invN++;
                    }
                    double testAngle = element.AngleList.Max();
                    if (testAngle > theta) { theta = testAngle; }

                    deltaDistMetric = (newDisMetric - oldDisMetric) + deltaDistMetric;
                }
                
                // 2.4. Check node location
                deltaDistMetric = deltaDistMetric / N;
                if (negN == N | invN > 0 | downN > posN | deltaDistMetric < -0.05 | theta > thetaMax) // constant as proposed in paper, 
                {
                    laplacianVector = laplacianVector * 0.5;
                    newPoint = new Point3d(oldNode.Coordinate.X + laplacianVector.X, oldNode.Coordinate.Y * laplacianVector.Y, oldNode.Coordinate.Z * laplacianVector.Z);
                    //node = new qNode(newPoint, newNode.BoundaryNode);
                }
                else { break; }
            }
            return newNode;
        }
         */
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Icon_QuadRemesh;
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
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
            pManager.AddGenericParameter("item", "element", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // variables
            Mesh mesh = new Mesh();
            List < qElement > elementList = new List<qElement>();
            List<qEdge> edgeList = new List<qEdge>();
            qElement element = new qElement();
            List<qEdge> frontEdges = new List<qEdge>();
            List<double> angleList = new List<double>();

            // input
            DA.GetData(0, ref mesh);

            // NOTE: Topological vertices and vertices are not sorted the same way. Make use to topological vertices in these codes.
            edgeList = CreateEdges(mesh);
            elementList = CreateElements(mesh, edgeList);
            SetNeighborElements(mesh, elementList, edgeList);
            frontEdges = GetFrontEdges(mesh, edgeList);
            SetNeighorFrontEdges(mesh, frontEdges);
            var edgeStates = CreateEdgeStateList(frontEdges); // todo: check bug. Angles...
            var list11 = edgeStates.Item1;
            var list10 = edgeStates.Item2;
            var list01 = edgeStates.Item3;
            var list00 = edgeStates.Item4;
            /*
            
            //-----side edge definiton-------
            qEdge choosenEdge = new qEdge();
            int edgeState = 0;
            
            // get edge
            if (list11.Count != 0) { choosenEdge = list11[0]; edgeState = 11; }
            else if (list10.Count != 0) { choosenEdge = list10[0]; edgeState = 10; }
            else if (list01.Count != 0) { choosenEdge = list01[0]; edgeState = 01; }
            else { choosenEdge = list00[0]; edgeState = 00; }

            // delete: temporary
            choosenEdge = list01[0]; 
            edgeState = 01;

            // get node, if not edgestate == 11
            qNode rightNode = GetRightNode(choosenEdge);
            qNode leftNode = GetLeftNode(choosenEdge);

            // fix sides
            // if edgestate 10: fix rightSide
            // if edgestate 01: fix leftSide
            // if edgestate 00: fix first leftSiden and then rightSide

            //__fix leffSide:
            qEdge E_fl = choosenEdge.LeftFrontNeighbor;
            qEdge E_fr = choosenEdge.RightFrontNeighbor;
            qEdge E_f = choosenEdge;
            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;

            // leftSide
            if (leftNode == E_fl.StartNode) // for left side
            {
                vec1 = rightNode.Coordinate - leftNode.Coordinate;
                vec2 = E_fl.EndNode.Coordinate - E_fl.StartNode.Coordinate;
            }
            else if (leftNode == E_fl.EndNode) // for left side
            {
                vec1 = rightNode.Coordinate - leftNode.Coordinate;
                vec2 = E_fl.StartNode.Coordinate - E_fl.EndNode.Coordinate;
            }
            qNode N_k = leftNode;
            Vector3d V_k = vec1.Length * vec2 + vec2.Length * vec1; // angle bisector
            qEdge E_1 = new qEdge();
            qEdge E_2 = new qEdge();

            int[] adjecentEdges = N_k.AdjacentEdges;

            foreach (int i in adjecentEdges) // loop adjacent edges to N_k
            {
                int[] adjacentElements = mesh.TopologyEdges.GetConnectedFaces(i); 
                foreach (int j in adjacentElements) // loop adjacent elements to adjacent edges
                {
                    // check if element !IsQuad
                    // check if element is equal to E_f element 1 or 2
                    // check if element is equal to E_fl element 1 or 2

                }

                List<int> sharedEdges = new List<int>();
                if (!edgeList[i].Element1.IsQuad) // if triangles elements on both sides (i.e. not front edges or quads)
                {
                    for (int j = 0; j < choosenEdge.Element1.EdgeList.Count; j++) // loop element 1
                    {
                        if (edgeList[i] == E_f.Element1.EdgeList[j]) { E_1 = edgeList[i]; }
                        else if (edgeList[i] == E_fl.Element1.EdgeList[j]) { E_2 = edgeList[i]; }
                    }

                }
            }
            */
            // output
            DA.SetDataList(0, edgeList);
            DA.SetDataList(1, list11);
            DA.SetDataList(2, list10);
            DA.SetDataList(3, list01);
            DA.SetDataList(4, list00);
            DA.SetDataList(5, elementList);
            DA.SetData(6,1);
        }
        #region Methods
        private List<qNode> CreateNodes(Mesh mesh)
        {
            List<qNode> nodeList = new List<qNode>();
            MeshTopologyVertexList topologyVertexList = mesh.TopologyVertices;
            for (int i = 0; i < topologyVertexList.Count; i++)
            {
                qNode node = new qNode();
                node.TopologyVertexIndex = i;
                node.MeshVertexIndex = topologyVertexList.MeshVertexIndices(i)[0]; // todo: check if ok
                node.Coordinate = mesh.Vertices.Point3dAt(node.MeshVertexIndex);
                node.AdjacentEdges = mesh.TopologyVertices.ConnectedEdges(node.TopologyVertexIndex);
                nodeList.Add(node);
            }
            return nodeList;
        }
        private List<qEdge> CreateEdges(Mesh mesh)
        {
            List<qEdge> edgeList = new List<qEdge>();
            MeshTopologyEdgeList topologyEdgeList = mesh.TopologyEdges;
            MeshTopologyVertexList topologyVertexList = mesh.TopologyVertices;
            List<qNode> nodeList = CreateNodes(mesh);

            for (int i = 0; i < topologyEdgeList.Count; i++)
            {
                Rhino.IndexPair edgeTopoVerticesIndex = topologyEdgeList.GetTopologyVertices(i);
                qNode startNode = nodeList[edgeTopoVerticesIndex[0]];
                qNode endNode = nodeList[edgeTopoVerticesIndex[1]];

                qEdge edge = new qEdge(i, startNode, endNode); //sjekk om i er nødvendig
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
        private void SetNeighorFrontEdges(Mesh mesh, List<qEdge> frontEdges)
        {
            for (int i = 0; i < frontEdges.Count; i++)
            {
                qEdge edg = frontEdges[i];
                qEdge edg1 = new qEdge();
                qEdge edg2 = new qEdge();
                qNode[] edgeNodes = new qNode[] { edg.StartNode, edg.EndNode };
                for (int j = 0; j < frontEdges.Count; j++)
                {
                    if (i != j)
                    {
                        qNode[] testNodes = new qNode[] { frontEdges[j].StartNode, frontEdges[j].EndNode };
                        if ((edgeNodes[0] == testNodes[0]) | (edgeNodes[0] == testNodes[1]))
                        {
                             edg1 = frontEdges[j];
                        }
                        else if ((edgeNodes[1] == testNodes[0]) | (edgeNodes[1] == testNodes[1]))
                        {
                             edg2 = frontEdges[j];
                        }
                        Point3d midPointEdg = (edgeNodes[0].Coordinate + edgeNodes[1].Coordinate) / 2.0;

                        int[] connectedElementsIndex1 = mesh.TopologyEdges.GetConnectedFaces(edg1.Index); // get the face index to edg1
                        Point3d center1 = mesh.Faces.GetFaceCenter(connectedElementsIndex1[0]);           // get the facecenter of edg1

                        int[] connectedElementsIndex2 = mesh.TopologyEdges.GetConnectedFaces(edg2.Index); // get the face index to edg1
                        Point3d center2 = mesh.Faces.GetFaceCenter(connectedElementsIndex2[0]);           // get the facecenter of edg1

                        Vector3d vec1 = center1 - midPointEdg;
                        Vector3d vec2 = center2 - midPointEdg;
                        if (Vector3d.CrossProduct(vec1, vec2).Z > 0) // todo: modify to more general relative to face normals
                        {
                            edg.RightFrontNeighbor = edg1;
                            edg.LeftFrontNeighbor = edg2;
                        }
                        else
                        {
                            edg.RightFrontNeighbor = edg2;
                            edg.LeftFrontNeighbor = edg1;
                        }
                    }
                }
            }
            return;
        }
        Tuple<List<qEdge>, List<qEdge>, List<qEdge>, List<qEdge>> CreateEdgeStateList(List<qEdge> frontEdges)
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
                        leftAngle = CalculateAngleOfAdjecentEdges(nodeToCalculate, frontEdges[i]);
                    }
                    else
                    {
                        rightAngle = CalculateAngleOfAdjecentEdges(nodeToCalculate, frontEdges[i]);
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
        private double CalculateAngleOfAdjecentEdges(int nodeToCalculate, qEdge edge)
        {
            // Calculate the angle between vector shared point --> other edgPoint and shared point --> other neigbor point.
            // Angles less than zero are to the left. Angles greater than
            // zero are to the right.
            qNode edgeStartNode = new qNode();
            qNode edgeEndNode = new qNode();
            qNode neighborStartNode = new qNode();
            qNode neighborEndNode = new qNode();

            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;
            double angle = 0;

            edgeStartNode = edge.StartNode;
            edgeEndNode = edge.EndNode;

            // Get neighbor nodes
            if (nodeToCalculate == 0)
            {
                neighborStartNode = edge.LeftFrontNeighbor.StartNode;
                neighborEndNode = edge.LeftFrontNeighbor.EndNode;
            }
            else
            {
                neighborStartNode = edge.RightFrontNeighbor.StartNode;
                neighborEndNode = edge.RightFrontNeighbor.EndNode;
            }

            // Create vectors from shared node
            if (neighborStartNode == edgeStartNode)
            {
                vec1 = edgeEndNode.Coordinate - edgeStartNode.Coordinate;
                vec2 = neighborEndNode.Coordinate - neighborStartNode.Coordinate;
            }
            else if (neighborStartNode == edgeEndNode)
            {
                vec1 = edgeStartNode.Coordinate - edgeEndNode.Coordinate;
                vec2 = neighborEndNode.Coordinate - neighborStartNode.Coordinate;
            }
            else if (neighborEndNode == edgeStartNode)
            {
                vec1 = edgeEndNode.Coordinate - edgeStartNode.Coordinate;
                vec2 = neighborStartNode.Coordinate - neighborEndNode.Coordinate;
            }
            else
            {
                vec1 = edgeStartNode.Coordinate - edgeEndNode.Coordinate;
                vec2 = neighborStartNode.Coordinate - neighborEndNode.Coordinate;
            }

            // Change vectors to 2d
            Vector2d ve1 = new Vector2d(vec1.X, vec1.Y); // todo 3d: calculate angles differently in 3d
            Vector2d ve2 = new Vector2d(vec2.X, vec2.Y);

            // Calculate lengths
            double len1 = Math.Sqrt(ve1.X * ve1.X + ve1.Y * ve1.Y);
            double len2 = Math.Sqrt(ve2.X * ve2.X + ve2.Y * ve2.Y);

            // Use the dot product to get the cosine.
            double dot_product = ve1.X * ve2.X + ve1.Y * ve2.Y;
            double cos = dot_product / len1 / len2;

            // Use the cross product to get the sine.
            double cross_product = ve1.X * ve2.Y - ve1.Y * ve2.X;
            double sin = cross_product / len1 / len2;

            // Find the angle.
            angle = Math.Acos(cos);
            if (sin < 0 & nodeToCalculate == 0)
            {
                angle = -angle;
            }
            return angle;
        }
        private qNode GetRightNode(qEdge edge)
        {
            // get the right node of a front edge
            qNode sharedNode = new qNode();
            if (edge.StartNode == edge.RightFrontNeighbor.StartNode) { sharedNode = edge.StartNode; }
            else if (edge.StartNode == edge.RightFrontNeighbor.EndNode) { sharedNode = edge.StartNode; }
            else if (edge.EndNode == edge.RightFrontNeighbor.StartNode) { sharedNode = edge.EndNode; }
            else if (edge.EndNode == edge.RightFrontNeighbor.EndNode) { sharedNode = edge.EndNode; }
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No shared nodes"); }
            return sharedNode;
        }
        private qNode GetLeftNode(qEdge edge)
        {
            // get the right node of a front edge
            qNode sharedNode = new qNode();
            if (edge.StartNode == edge.LeftFrontNeighbor.StartNode) { sharedNode = edge.StartNode; }
            else if (edge.StartNode == edge.LeftFrontNeighbor.EndNode) { sharedNode = edge.StartNode; }
            else if (edge.EndNode == edge.LeftFrontNeighbor.StartNode) { sharedNode = edge.EndNode; }
            else if (edge.EndNode == edge.LeftFrontNeighbor.EndNode) { sharedNode = edge.EndNode; }
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No shared nodes"); }
            return sharedNode;
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
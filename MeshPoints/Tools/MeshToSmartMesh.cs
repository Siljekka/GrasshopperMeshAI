using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using Rhino.Geometry.Collections;
using System.Linq;

namespace MeshPoints.Tools
{
    public class MeshToSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MeshToSmarMesh class.
        /// </summary>
        public MeshToSmartMesh()
          : base("TrinagleMesh To SmarMesh", "ToSM",
              "Transforms the class Mesh to class SmartMesh. The mesh must be composed of triangle elements.",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Triangle Mesh", "mesh", "Surface mesh with triangle elements.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "mesh", "SmartMesh for surfaces.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh mesh = new Mesh();
            DA.GetData(0, ref mesh);

            // 1. Get initial edges and elements of mesh using mesh topology properties
            var initialEdgeAndElementList = GetInitialEdgesAndElements(mesh);
            List<qElement> globalElementList = initialEdgeAndElementList.Item2;
            foreach (qElement element in globalElementList)
            {
                element.FixEdgeOrder();
            }

            // 2. Convert mesh to SmartMesh Class and create final mesh
            var meshProperties = ConvertToMainMeshClasses(globalElementList);
            List<Node> nodes = meshProperties.Item1;
            List<Element> elements = meshProperties.Item2;
            SmartMesh surfaceMesh = new SmartMesh(nodes, elements, "Surface");

            // Output
            DA.SetData(0, surfaceMesh);

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
                return Properties.Resources.Icon_MeshToSM;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("c303644f-24eb-4ad9-a32a-bbe2adb74a9a"); }
        }
    }
}
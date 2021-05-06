using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.Tools
{
    public class MergeSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MergeSmartMesh class.
        /// </summary>
        public MergeSmartMesh()
          : base("Merge SmartMesh", "merge",
              "Merge multiple SmartMesh to one SmartMesh",
              "MyPlugIn", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("List of SmartMesh", "smartMeshes", "Input a list of SmartMesh to merge", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "smartMesh", "Merged SmartMesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Assumptions: merge only surface edges, edges to merge must share start/edge point, must have equal number of nodes. For solid, sweep the merged SmartMesh.

            // Input
            List<SmartMesh> smartMeshList = new List<SmartMesh>(); 
            DA.GetDataList(0,  smartMeshList);

            // Code
            SmartMesh mergedSmartMesh = new SmartMesh();

            for (int i = 0; i < smartMeshList.Count - 1; i++) // loop the meshes to merge
            {
                SmartMesh mesh1 = smartMeshList[i];
                SmartMesh mesh2 = smartMeshList[i + 1];

                BrepEdge edge1 = null;
                BrepEdge edge2 = null;

                // Find edges to merge
                for (int j = 0; j < 4; j++)
                {
                    BrepEdge edgeToTest1 = mesh1.Geometry.Edges[j];
                    for (int k = 0; k < 4; k++)
                    {
                        BrepEdge edgeToTest2 = mesh2.Geometry.Edges[k];
                        List<Point3d> edgePoints1 = new List<Point3d>() { edgeToTest1.PointAtStart, edgeToTest1.PointAtEnd };
                        List<Point3d> edgePoints2 = new List<Point3d>() { edgeToTest2.PointAtStart, edgeToTest2.PointAtEnd };

                        if (edgePoints1.Contains(edgePoints2[0]) & edgePoints1.Contains(edgePoints2[1])) // if edges share start and end points
                        {
                            edge1 = edgeToTest1;
                            edge2 = edgeToTest2;
                            i = 4; // break
                            j = 4; // break
                        }
                    }
                }
                if (edge1 == null | edge2 == null) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No edges to merge."); return; }
                
                // Find nodes to merge
                List<Node> nodesOnEdge1 = new List<Node>();
                List<Node> nodesOnEdge2 = new List<Node>();

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

                if (nodesOnEdge1.Count != nodesOnEdge2.Count) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Number of nodes on edges to merge do not match."); return; }

                // Pare nodes to merge
                List<List<Node>> nodesToMerge = new List<List<Node>>();
                int numNodesToMerge = nodesOnEdge1.Count; // to do: sjekk at denne ikke oppdateres
                for (int j = 0; j < numNodesToMerge; j++)
                {
                    Node node1 = nodesOnEdge1[j];
                    Node node2 = FindClosestNode(node1, nodesOnEdge2);

                    //nodesToMerge[j] = new List<Node>(){node1, node2};

                    // sjekk ut dette:
                    Point3d newNodeCoordinate = (node1.Coordinate + node2.Coordinate) / (double)2;
                    Node newNode = new Node(1000000 + j, newNodeCoordinate); // dummy solution for global id: 1000000 +j

                    int node1Index = mesh1.Nodes.IndexOf(node1);
                    int node2Index = mesh2.Nodes.IndexOf(node2);

                    // to do: assign BC

                    // Replace old nodes with new nodes
                    mesh1.Nodes[node1Index] = newNode;
                    mesh2.Nodes[node2Index] = newNode;


                    // to do: check if elements are corretly changed
                    // if not: find elements to change and update the nodes and connectivity

                }            

                // Update global id of nodes
                int numNodesMesh1 = mesh1.Nodes.Count;
                foreach (Node node2 in mesh2.Nodes)
                {
                    if (node2.GlobalId < 1000000)
                    {
                        node2.GlobalId += numNodesMesh1;
                    }
                }

                // Merge nodes
                List<Node> nodes = new List<Node>();
                nodes.AddRange(mesh1.Nodes);
                nodes.AddRange(mesh2.Nodes);

                // to do: edge nodes are dublicated...

                // Fix elementes adjacent to edge

                // Merge elements
                List<Element> elements = new List<Element>();
                elements.AddRange(mesh1.Elements);
                elements.AddRange(mesh2.Elements);

                // Merge mesh
                Mesh mergedGlobalMesh = null;
                mergedGlobalMesh.Append(mesh1.Mesh);
                mergedGlobalMesh.Append(mesh1.Mesh);

                // Fix geometry
            }

            DA.SetData(0, mergedSmartMesh);

        }


        // Methods
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
            get { return new Guid("4eede3a5-ff72-42e4-bac3-5c35cbeff345"); }
        }
    }
}
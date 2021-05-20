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
              "Merge multiple SmartMesh to one SmartMesh",
              "SmartMesh", "Tools")
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
            // Assumption: U and V direction match

            // Input
            List<SmartMesh> smartMeshList = new List<SmartMesh>(); 
            DA.GetDataList(0,  smartMeshList);

            // Code
            SmartMesh mergedSmartMesh = new SmartMesh();
            SmartMesh mesh1 = smartMeshList[0];


            for (int i = 1; i < smartMeshList.Count; i++) // loop the meshes to merge
            {
                SmartMesh mesh2 = smartMeshList[i];

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
                            j = 4; // break
                            k = 4; // break
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

                // Merge nodes
                int numNodesToMerge = nodesOnEdge1.Count; 
                List<Node> newNodes = new List<Node>();
                newNodes.AddRange(mesh1.Nodes);
                newNodes.AddRange(mesh2.Nodes);

                int startIdNewNodes = mesh1.Nodes.Count + mesh2.Nodes.Count - nodesOnEdge1.Count - nodesOnEdge2.Count; 
                for (int j = 0; j < numNodesToMerge; j++)
                {
                    Node node1 = nodesOnEdge1[j];
                    Node node2 = FindClosestNode(node1, nodesOnEdge2);

                    Point3d newNodeCoordinate = (node1.Coordinate + node2.Coordinate) / (double)2;
                    Node newNode = new Node(1000000 + j, newNodeCoordinate); // dummy solution for global id: 1000000 +j

                    if (node1.Type == "Corner")
                    {
                        newNode.BC_U = true;
                        newNode.BC_V = true;
                    }

                    int node1Index = mesh1.Nodes.IndexOf(node1);
                    int node2Index = mesh2.Nodes.IndexOf(node2);

                    // Find elements to change and update the nodes and connectivity
                    // to do: fix bug for element in middel for multiple SmartMesh
                    foreach (Element element in mesh1.Elements)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            if (element.Nodes[k] == mesh1.Nodes[node1Index])
                            {
                                element.Nodes[k] = newNode;
                                element.Connectivity[k] = startIdNewNodes + j;
                            }
                        }
                        UpdateElementMesh(element);
                    }

                    foreach (Element element in mesh2.Elements)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            if (element.Nodes[k] == mesh2.Nodes[node2Index])
                            {
                                element.Nodes[k] = newNode;
                                element.Connectivity[k] = startIdNewNodes + j;
                            }
                        }
                        UpdateElementMesh(element);
                    }

                    newNodes.Add(newNode);
                    newNodes.Remove(node1);
                    newNodes.Remove(node2);
                }
 
                // Merge elements
                List<Element> newElements = new List<Element>();
                newElements.AddRange(mesh1.Elements);
                newElements.AddRange(mesh2.Elements);

                // Update global id of nodes in mesh2
                for (int j = mesh1.Nodes.Count - numNodesToMerge; j < (mesh1.Nodes.Count + mesh2.Nodes.Count - nodesOnEdge1.Count); j++ ) // loop nodes from mesh2
                {
                    newNodes[j].GlobalId = j;
                }

                // Update elements in mesh2
                for (int j = mesh1.Elements.Count; j < mesh1.Elements.Count + mesh2.Elements.Count; j++) // loop elements from mesh2
                {
                    for (int k = 0; k < 4; k++)
                    {
                        if (newElements[j].Nodes[k].GlobalId < startIdNewNodes)
                        {
                            newElements[j].Connectivity[k] += mesh1.Nodes.Count - nodesOnEdge1.Count;
                            newElements[j].Nodes[k] = newNodes[newElements[j].Connectivity[k]];
                        }
                    }
                }

                // Merge mesh, with weld
                Mesh mergedGlobalMesh = new Mesh();
                
                foreach (Element element in newElements)
                {
                    mergedGlobalMesh.Append(element.Mesh);
                }
                mergedGlobalMesh.Weld(0.001);

                // Create merged SmartMesh
                mergedSmartMesh = new SmartMesh(newNodes, newElements, mergedGlobalMesh, "Surface");



                //  Fix Geometry

                List<BrepEdge> edges = new List<BrepEdge>();
                edges.AddRange(mesh1.Geometry.Edges);
                edges.AddRange(mesh2.Geometry.Edges);
                edges.Remove(edge1);
                edges.Remove(edge2);

                List<BrepFace> faces = new List<BrepFace>();
                faces.AddRange(mesh1.Geometry.Faces);
                faces.AddRange(mesh2.Geometry.Faces);

                Brep totalBrep = mesh1.Geometry.Brep;
                totalBrep.Join(mesh2.Geometry.Brep,0.001, true);

                Geometry mergedGeometry = new Geometry(totalBrep, faces, edges, totalBrep.Vertices.ToList());
                mergedSmartMesh.Geometry = mergedGeometry;

                mesh1 = mergedSmartMesh;
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
        


        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Icon_MergeNodes;
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
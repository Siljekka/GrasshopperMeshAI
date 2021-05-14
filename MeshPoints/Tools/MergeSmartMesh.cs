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

                // Create new node list
                List<Node> newNodes = new List<Node>();
                List<Node> oldNodes = new List<Node>();
                oldNodes.AddRange(mesh1.Nodes);
                oldNodes.AddRange(mesh2.Nodes);
                foreach (Node node in oldNodes)
                {
                    Node newNode = new Node(node.GlobalId, node.Coordinate, node.BC_U, node.BC_V);
                    newNodes.Add(newNode);
                }

                // Create new element list
                List<Element> newElements = new List<Element>();
                List<Element> oldElements = new List<Element>();
                oldElements.AddRange(mesh1.Elements);
                oldElements.AddRange(mesh2.Elements);
                foreach (Element element in oldElements)
                {
                    Element newElement = new Element(element.Id, element.Nodes, element.Connectivity);
                    newElement.Mesh = element.Mesh;
                    newElements.Add(newElement);
                }

                // Find nodes to merge
                var edges = FindEdgeToMerge(mesh1, mesh2);
                BrepEdge edge1 = edges.Item1;
                BrepEdge edge2 = edges.Item2;

                var nodesToMerge = FindNodesToMerge(mesh1, mesh2, edge1, edge2);
                List<Node> nodesOnEdge1 = nodesToMerge.Item1;
                List<Node> nodesOnEdge2 = nodesToMerge.Item2;

                if (nodesOnEdge1[0] == null) { return; }

                // Merge nodes
                int numNodesToMerge = nodesOnEdge1.Count; 
                int startIdNewNodes = mesh1.Nodes.Count + mesh2.Nodes.Count - numNodesToMerge * 2; 

                for (int j = 0; j < numNodesToMerge; j++)
                {
                    // Find node from mesh2 that is equal to node from mesh1, and delete it. Update connected elements.
                    Node node1 = nodesOnEdge1[j];
                    Node node2 = FindClosestNode(node1, nodesOnEdge2);

                    Point3d newNodeCoordinate = (node1.Coordinate + node2.Coordinate) / (double)2;
                    Node newNode = new Node(startIdNewNodes + j, newNodeCoordinate);

                    if (node1.Type == "Corner") // keep BC if corner node
                    {
                        newNode.BC_U = true;
                        newNode.BC_V = true;
                    }

                    // Find elements to change and update the nodes and connectivity
                    foreach (Element element in newElements)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            if (element.Nodes[k] == node1)
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
                            else if (element.Nodes[k] == node2)
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
                        //UpdateElementMesh(element);
                    }
                    newNodes.Add(newNode);
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
                Mesh mergedGlobalMesh = new Mesh();
                int counter1 = 0;
                int counter2 = 0;

                for (int j = 0; j < newElements.Count; j++) // loop elements from mesh1
                {
                    for (int k = 0; k < 4; k++)
                    {
                        // loop through - old to new id of nodes to get the updated id
                        if (j < mesh1.Elements.Count)
                        {
                            for (int n = 0; n < mesh1.Nodes.Count - numNodesToMerge; n++)
                            {
                                List<int> oldToNew = oldToNewId[n];
                                if (oldToNew[0] == mesh1.Elements[counter1].Connectivity[k])
                                {
                                    newElements[j].Nodes[k] = newNodes[oldToNew[1]];
                                    newElements[j].Connectivity[k] = oldToNew[1];
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (int n = mesh1.Nodes.Count - numNodesToMerge; n < mesh1.Nodes.Count + mesh2.Nodes.Count - numNodesToMerge*2; n++)
                            {
                                List<int> oldToNew = oldToNewId[n];
                                if (oldToNew[0] == mesh2.Elements[counter2].Connectivity[k])
                                {
                                    newElements[j].Nodes[k] = newNodes[oldToNew[1]];
                                    newElements[j].Connectivity[k] = oldToNew[1];
                                    break;
                                }
                            }

                        } 
                    }
                    if (j < mesh1.Elements.Count) { counter1++; }
                    else { counter2++; }
                    mergedGlobalMesh.Append(newElements[j].Mesh); // add element mesh to global mesh
                }
                mergedGlobalMesh.Weld(0.001); // delete overlapping vertices and faces

                // Create merged SmartMesh
                mergedSmartMesh = new SmartMesh(newNodes, newElements, mergedGlobalMesh, "Surface");

                //  Fix Geometry
                List<BrepEdge> brepEdges = new List<BrepEdge>();
                brepEdges.AddRange(mesh1.Geometry.Edges);
                brepEdges.AddRange(mesh2.Geometry.Edges);
                brepEdges.Remove(edge1);
                brepEdges.Remove(edge2);

                List<BrepFace> faces = new List<BrepFace>();
                faces.AddRange(mesh1.Geometry.Faces);
                faces.AddRange(mesh2.Geometry.Faces);

                Brep totalBrep = null;
                totalBrep =  mesh1.Geometry.Brep;
                totalBrep.Join(mesh2.Geometry.Brep,0.001, true);

                Geometry mergedGeometry = new Geometry(totalBrep, faces, brepEdges, totalBrep.Vertices.ToList());
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
        private Tuple<BrepEdge, BrepEdge> FindEdgeToMerge(SmartMesh mesh1, SmartMesh mesh2)
        {
            BrepEdge edge1 = null;
            BrepEdge edge2 = null;
            List<Node> nodesOnEdge1 = new List<Node>();
            List<Node> nodesOnEdge2 = new List<Node>();

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
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;


namespace MeshPoints.Classes
{
    class Element
    {
        public List<Node> Nodes { get; set; }
        public List<int> Connectivity { get; set; }
        public string Type { get; set; }
        public Quality MeshQuality { get; set; }
        public int Id { get; set; }
        public Mesh Mesh { get; set; }

        //Constructer
        public Element()
        {
            //Empty constructor
        }
        public Element(int _id, List<Node> _nodes, List<int> _connectivity)
        {
            Id = _id;
            Nodes = _nodes;
            Connectivity = _connectivity;
            GetElementType();
            GetElementMesh();
        }

        // Methods
        public void GetElementType()
        {
            string type = "null";
            switch (this.Nodes.Count)
            {
                case 3:
                    type = "Triangle";
                    break;
                case 4:
                    type = "Quad";
                    break;
                case 6:
                    type = "Tet";
                    break;
                case 8:
                    type = "Hex";
                    break;
            }
            this.Type = type;
        }

        public List<List<Node>> GetFaces()
        {
            List<List<Node>> Faces = new List<List<Node>>();
            List<Node> nodes = this.Nodes;

            if (this.Type == "Quad" | this.Nodes.Count == 3) // surface
            {
                Faces.Add(nodes);
            }
            else // solid
            {
                List<List<int>> nodeIndex = new List<List<int>>
                {
                      new List<int> {0, 1, 5, 4}, new List<int> {1, 2, 6, 5}, new List<int> {2, 3, 7, 6},
                      new List<int> { 3, 0, 4, 7 }, new List<int> { 0, 3, 2, 1 },  new List<int> { 4, 5, 6, 7 }
                };
  
                foreach (List<int> indices in nodeIndex)
                {
                    List<Node> nodesOfFace = new List<Node>();
                    foreach (int i in indices)
                    {
                        nodesOfFace.Add(nodes[i]);
                    }
                    Faces.Add(nodesOfFace);
                }
            }
            return Faces;
        }

        public void GetElementMesh()
        {
            Mesh mesh = new Mesh();
            foreach (Node node in this.Nodes)
            {
                mesh.Vertices.Add(node.Coordinate);
            }
            if (this.Type == "Quad")
            {
                mesh.Faces.AddFace(0, 1, 2, 3);
            }
            else if (this.Type == "Hex")
            {
                mesh.Faces.AddFace(0, 1, 5, 4);
                mesh.Faces.AddFace(1, 2, 6, 5);
                mesh.Faces.AddFace(2, 3, 7, 6);
                mesh.Faces.AddFace(3, 0, 4, 7);
                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.Faces.AddFace(4, 5, 6, 7);
            }
            else if (this.Type == "Triangle")
            {
                mesh.Faces.AddFace(0, 1, 2);
            }
            else if (this.Type == "Tet")
            {
                mesh.Faces.AddFace(0, 1, 2);
                mesh.Faces.AddFace(3, 4, 5);
                mesh.Faces.AddFace(0, 1, 4, 3);
                mesh.Faces.AddFace(1, 2, 5, 4);
                mesh.Faces.AddFace(2, 0, 3, 5);
            }

            mesh.Compact(); //to ensure that it calculate
            mesh.FaceNormals.ComputeFaceNormals();
            mesh.UnifyNormals();
            this.Mesh = mesh;
        }
    }
}

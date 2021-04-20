using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;


namespace MeshPoints.Classes
{
    class Element
    {
        public Node Node1 { get; set; } // delete
        public Node Node2 { get; set; } // delete
        public Node Node3 { get; set; } // delete
        public Node Node4 { get; set; } // delete
        public Node Node5 { get; set; } // delete
        public Node Node6 { get; set; } // delete
        public Node Node7 { get; set; } // delete
        public Node Node8 { get; set; } // delete
        public List<Node> Nodes { get; set; }
        public List<int> Connectivity { get; set; }

        public string Type { get; set; }

        public Quality MeshQuality { get; set; }
        public int Id { get; set; }
        public Mesh mesh { get; set; }

        //Constructer

        //_empty
        public Element()
        {
            //Empty constructor
        }

        public Element(int _id, List<Node> _nodes, List<int> _connectivity) // new
        {
            Id = _id;
            Nodes = _nodes;
            Connectivity = _connectivity;
            GetType(Nodes);
        }


        //_for 2D - triangle
        public Element(int _id, Node _node1, Node _node2, Node _node3, Mesh _mesh)
        {
            Id = _id;
            Node1 = _node1;
            Node2 = _node2;
            Node3 = _node3;
            mesh = _mesh;
        }
        //_for 2D - quad
        public Element(int _id, Node _node1, Node _node2, Node _node3, Node _node4, Mesh _mesh)
        {
            Id = _id;
            Node1 = _node1;
            Node2 = _node2;
            Node3 = _node3;
            Node4 = _node4;
            mesh = _mesh;
        }

        //_for 3D
        public Element(int _id, Node _node1, Node _node2, Node _node3, Node _node4, Node _node5, Node _node6, Node _node7, Node _node8)
        {
            Id = _id;
            Node1 = _node1;
            Node2 = _node2;
            Node3 = _node3;
            Node4 = _node4;
            Node5 = _node5;
            Node6 = _node6;
            Node7 = _node7;
            Node8 = _node8;
        }


        // Methods
        private void GetType(List<Node> Nodes)
        {
            string type = "null";
            switch (Nodes.Count)
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

            if (this.Type == "Quad") // surface
            {
                Faces.Add(nodes);
            }
            else // solid
            {
                List<List<int>> nodeIndex = new List<List<int>>
                {
                      new List<int> {0, 1, 5, 4}, new List<int> {1, 2, 6, 5}, new List<int> {2, 3, 7, 6},
                      new List<int> { 3, 0, 4, 7 }, new List<int> { 0, 1, 2, 3 },  new List<int> { 4, 5, 6, 7 }
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

    }
}

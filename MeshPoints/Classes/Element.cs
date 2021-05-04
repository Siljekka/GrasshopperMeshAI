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
        public List<Node> Nodes { get; set; }
        public List<int> Connectivity { get; set; }
        public string Type { get; set; }
        public Quality MeshQuality { get; set; }
        public int Id { get; set; }
        public Mesh Mesh { get; set; }

        //Constructer

        //_empty
        public Element()
        {
            //Empty constructor
        }

        public Element(int _id, List<Node> _nodes, List<int> _connectivity)
        {
            Id = _id;
            Nodes = _nodes;
            Connectivity = _connectivity;
            GetType(Nodes);
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

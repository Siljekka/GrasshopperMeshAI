using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeshPoints.Classes
{
    class Mesh2D
    {
        public List<Element> Elements{ get; set; } //list of nodes
        public List<Node> Nodes { get; set; } //list of nodes
        public int Nx{ get; set; } //number of nodes in x-dir
        public int Ny { get; set; } //number of nodes in y-dir

        public Mesh2D()
        {
            //Empty constructor
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeshPoints.Classes
{
    class Mesh3D
    {
        public List<Element> Elements { get; set; } //list of elements
        public List<Node> Nodes { get; set; } //list of nodes
        public int Nx { get; set; } //number of nodes in x-dir
        public int Ny { get; set; } //number of nodes in y-dir
        public int Nz { get; set; } //number of nodes in z-dir

        public Mesh3D()
        {
            //Empty constructor
        }

    }
}

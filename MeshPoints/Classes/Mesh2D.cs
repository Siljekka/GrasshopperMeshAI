using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    class Mesh2D
    {
        public List<Element> Elements{ get; set; } //list of nodes
        public List<Node> Nodes { get; set; } //list of nodes
        public Mesh mesh { get; set; } //list of nodes
        public int nu{ get; set; } //number of nodes in u-dir
        public int nv { get; set; } //number of nodes in v-dir

        public Mesh2D()
        {
            //Empty constructor
        }
        public Mesh2D(int _nu, int _nv, List<Node> _nodes, List<Element> _elements, Mesh _mesh)
        {
            nu = _nu;
            nv = _nv;
            Nodes = _nodes;
            Elements = _elements;
            mesh = _mesh;
        }
    }
}

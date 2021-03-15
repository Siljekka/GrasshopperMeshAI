using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    class Mesh3D
    {
        public List<Element> Elements { get; set; } //list of elements
        public List<Node> Nodes { get; set; } //list of nodes
        public Mesh mesh { get; set; } //mesh
        public int nu { get; set; } //number of nodes in x-dir
        public int nv { get; set; } //number of nodes in y-dir
        public int nw { get; set; } //number of nodes in z-dir
        public bool inp { get; set; }
        public Mesh3D()
        {
            //Empty constructor
        }

        public Mesh3D(int _nu, int _nv, int _nw, List<Node> _nodes, List<Element> _elements, Mesh _mesh)
        {
            nu = _nu;
            nv = _nv;
            nw = _nw;
            Nodes = _nodes;
            Elements = _elements;
            mesh = _mesh;
        }
    }
}

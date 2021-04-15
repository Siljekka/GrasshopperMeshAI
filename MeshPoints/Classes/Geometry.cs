using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    class Geometry
    {
        public Brep Brep { get; set; }
        public List<BrepFace> Faces { get; set; }
        public List<BrepEdge> Edges { get; set; }
        public List<BrepVertex> Vertices { get; set; }

        // Constructors
        public Geometry()
        { 
            // empty constructor
        }
        public Geometry(Brep _brep, List<BrepFace> _faces, List<BrepEdge> _edges, List<BrepVertex> _vertices)
        {
            Brep = _brep;
            Faces = _faces;
            Edges = _edges;
            Vertices = _vertices;
        }
    }
}

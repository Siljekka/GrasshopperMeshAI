using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    class GeometryInfo
    {
        public List<BrepFace> Faces { get; set; }
        public List<BrepEdge> Edges { get; set; }
        public List<BrepVertex> Vertices { get; set; }


        public GeometryInfo()
        { 
            // empty constructor
        }
    }
}

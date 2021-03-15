using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    class qNode
    {
        public Point3d Coordinate { get; set; }
        public int TopologyVertexIndex { get; set; } // vertex index in topology
        public int MeshVertexIndex { get; set; } // vertex index in mesh
        public int[] ConnectedEdges { get; set; } // index of connected edges

        public qNode()
        {
            // empty constructor
        }
    }
}

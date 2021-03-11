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
        public int[] AdjacentEdges { get; set; } // index of adjacent edges

        public qNode()
        {
            // empty constructor
        }

        // contains only a coordinate. Use Point3d class.
    }
}

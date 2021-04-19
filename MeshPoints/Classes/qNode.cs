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
        public bool BoundaryNode { get; set; }

        public qNode()
        {
            // empty constructor
        }

        public qNode(Point3d _coordinate, bool _boundaryNode)
        {
            Coordinate = _coordinate;
            BoundaryNode = _boundaryNode;
        }
    }
}

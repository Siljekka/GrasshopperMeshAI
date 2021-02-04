using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using Grasshopper.Kernel;

namespace MeshPoints.Classes
{
    class Node
    {
        public int LocalId { get; set; }
        public int GlobalId { get; set; }
        public Point3d Coordinate { get; set; }

        public bool BC_X { get; set; }
        public bool BC_Y { get; set; }
        public bool CornerNode { get; set; }


        //Constructor
        public Node()
        {
            //Empty constructor
        }

        public Node(int _globalId, Point3d _coord)
        {
            GlobalId = _globalId;
            Coordinate = _coord;
        }
        public Node(int _globalId, Point3d _coord, bool _BC_X, bool _BC_Y)
        {
            GlobalId = _globalId;
            Coordinate = _coord;
            BC_X = _BC_X;
            BC_Y = _BC_Y;
        }

        //Methods

    }
}

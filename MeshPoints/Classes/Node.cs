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

        public bool BC_U { get; set; }
        public bool BC_V { get; set; }
        
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
        public Node(int _globalId, Point3d _coord, bool _BC_U, bool _BC_V)
        {
            GlobalId = _globalId;
            Coordinate = _coord;
            BC_U = _BC_U;
            BC_V = _BC_V;
        }

        public Node(int _locald, int _globalId, Point3d _coord, bool _BC_U, bool _BC_V)
        {
            LocalId = _locald;
            GlobalId = _globalId;
            Coordinate = _coord;
            BC_U = _BC_U;
            BC_V = _BC_V;
        }

        //Methods

    }
}

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
        public int LocalId { get; set; } // delete
        public int GlobalId { get; set; }
        public Point3d Coordinate { get; set; }

        public bool BC_U { get; set; } // to do: change if time
        public bool BC_V { get; set; } // to do: change if time
        public bool BC_W { get; set; } // to do: change if time

        public string Type { get; set; }

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
            BC_W = true;
            this.SetType();
        }

        public Node(int _locald, int _globalId, Point3d _coord, bool _BC_U, bool _BC_V) // to do: slett
        {
            LocalId = _locald;
            GlobalId = _globalId;
            Coordinate = _coord;
            BC_U = _BC_U;
            BC_V = _BC_V;
            BC_W = true; // to do: sjekk om ok
            this.SetType();
        }

        public Node(int _globalId, Point3d _coord, bool _BC_U, bool _BC_V, bool  _BC_W)
        {
            GlobalId = _globalId;
            Coordinate = _coord;
            BC_U = _BC_U;
            BC_V = _BC_V;
            BC_W = _BC_W;
            this.SetType();
        }

        public Node(int _locald, int _globalId, Point3d _coord, bool _BC_U, bool _BC_V, bool _BC_W) // to do: slett
        {
            LocalId = _locald;
            GlobalId = _globalId;
            Coordinate = _coord;
            BC_U = _BC_U;
            BC_V = _BC_V;
            BC_W = _BC_W;
            this.SetType();
        }

        //Methods
        public bool IsOnFace(BrepFace face)
        {
            Point3d point = this.Coordinate;
            bool isOnFace = false;

            face.ClosestPoint(point, out double PointOnCurveU, out double PointOnCurveV);
            Point3d testPoint = face.PointAt(PointOnCurveU, PointOnCurveV);  // make test point 
            double distanceToFace = (testPoint - point).Length; // calculate distance between testPoint and node
            if (distanceToFace <= 0.0001 & distanceToFace >= -0.0001) // if distance = 0: node is on edge
            {
                isOnFace = true;
            }
            return isOnFace;
        }

        public bool IsOnEdge(BrepEdge edge)
        {
            Point3d point = this.Coordinate;
            bool isOnEdge = false;

            edge.ClosestPoint(point, out double PointOnCurve);
            Point3d testPoint = edge.PointAt(PointOnCurve);  // make test point 
            double distanceToEdge = (testPoint - point).Length; // calculate distance between testPoint and node
            if (distanceToEdge <= 0.0001 & distanceToEdge >= -0.0001) // if distance = 0: node is on edge
            {
                isOnEdge = true;               
            }
            return isOnEdge;
        }

        public void SetType()
        {
            int counter = 0;
            if (this.BC_U) { counter++; }
            if (this.BC_V) { counter++; }
            if (this.BC_W) { counter++; }

            switch (counter)
            {
                case 3: // corner node
                    this.Type = "Corner"; break;
                case 2: // edge node
                    this.Type = "Edge"; break;
                case 1: // midle node
                    this.Type = "Face"; break;
                case 0:
                    this.Type = "Interior"; break;
            }
        }
    }

}

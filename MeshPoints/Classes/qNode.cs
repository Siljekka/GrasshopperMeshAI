﻿using System;
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
        public bool boundaryNode { get; set; }

        public qNode()
        {
            // empty constructor
        }
    }
}

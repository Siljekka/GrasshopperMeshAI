using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using Grasshopper.Kernel;
using MeshPoints.Classes;

namespace MeshPoints.Classes
{
    class Material
    {
        public double YoungModulus { get; set; }
        public double PossionRatio { get; set; }
        public double YieldingStress { get; set; }

        // constructer
        public Material()
        { 
             // empty constructer
        }

        public Material(double _youngModulus, double _possionRatio, double _yieldingStress)
        {
            YoungModulus = _youngModulus;
            PossionRatio = _possionRatio;
            YieldingStress = _yieldingStress;
        }

    }
}

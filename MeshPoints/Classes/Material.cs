using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using Grasshopper.Kernel;
using MeshPoints.Classes;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;


namespace MeshPoints.Classes
{
    class Material
    {
        public double YoungModulus { get; set; }
        public double PossionRatio { get; set; }
        public double ShellThickness { get; set; }

        // constructer
        public Material()
        { 
             // empty constructer
        }

        public Material(double _youngModulus, double _possionRatio)
        {
            YoungModulus = _youngModulus;
            PossionRatio = _possionRatio;
        }

    }
}

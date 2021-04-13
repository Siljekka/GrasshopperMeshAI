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
        public double YoungModulus { get; set; } // delete
        public double PossionRatio { get; set; }

        public Matrix<double> MaterialConstant {get; set;}

        // constructer
        public Material()
        { 
             // empty constructer
        }

        public Material(double _youngModulus, double _possionRatio)
        {
            YoungModulus = _youngModulus;
            PossionRatio = _possionRatio;
            MaterialConstant = GetMaterialConstant(_youngModulus, _possionRatio);
        }


        // method
        private Matrix<double> GetMaterialConstant(double _youngModulus, double _possionRatio)
        {
            Matrix<double> C = DenseMatrix.OfArray(new double[,]
            {
                    {1-_possionRatio, _possionRatio, _possionRatio, 0, 0, 0},
                    {_possionRatio, 1-_possionRatio, _possionRatio, 0, 0, 0},
                    {_possionRatio, _possionRatio, 1- _possionRatio, 0, 0, 0},
                    {0, 0, 0, (1-2*_possionRatio)/2, 0, 0},
                    {0, 0, 0, 0, (1-2*_possionRatio)/2, 0},
                    {0, 0, 0, 0, 0, (1-2*_possionRatio)/2},
            });
            C.Multiply((double) _youngModulus / ((1 + _possionRatio) * (1 - 2 * _possionRatio)));

            return C;
        }
    }
}

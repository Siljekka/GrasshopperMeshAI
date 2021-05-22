using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;


namespace MeshPoints.Classes
{
    class FEM
    {
        public Vector<double> GetShapeFunctions(double r, double s, double t, int nodeDOFS)
        {
            // Shapefunctions on matrix form
            if (nodeDOFS == 2)
            {
                Vector<double> N = DenseVector.OfArray(new double[]
                {
                    (1-r)*(1-s),
                    (1+r)*(1-s),
                    (1+r)*(1+s),
                    (1-r)*(1+s),
                });
                N = N.Multiply(0.25);
                return N;
            }
            else
            {
                Vector<double> N = DenseVector.OfArray(new double[]
                {
                    (1-r)*(1-s)*(1-t),
                    (1+r)*(1-s)*(1-t),
                    (1+r)*(1+s)*(1-t),
                    (1-r)*(1+s)*(1-t),
                    (1-r)*(1-s)*(1+t),
                    (1+r)*(1-s)*(1+t),
                    (1+r)*(1+s)*(1+t),
                    (1-r)*(1+s)*(1+t)
                });
                N = N.Multiply(0.125);
                return N;
            }

        }
        public Matrix<double> DerivateWithNatrualCoordinates(double r, double s, double t, int nodeDOFS)
        {
            if (nodeDOFS == 2)
            {
                Matrix<double> shapeFunctionsDerivatedNatural = DenseMatrix.OfArray(new double[,]
                {
                    {-(1-s), (1-s), (1+s), -(1+s)},
                    {-(1-r), -(1+r), (1+r), (1-r)}
                });
                shapeFunctionsDerivatedNatural = shapeFunctionsDerivatedNatural.Multiply(0.25);
                return shapeFunctionsDerivatedNatural;

            }
            else
            {
                Matrix<double> shapeFunctionsDerivatedNatural = DenseMatrix.OfArray(new double[,]
                {
                        {-(1-s)*(1-t), (1-s)*(1-t), (1+s)*(1-t),-(1+s)*(1-t),-(1-s)*(1+t),(1-s)*(1+t),(1+s)*(1+t),-(1+s)*(1+t)},
                        {-(1-r)*(1-t), -(1+r)*(1-t), (1+r)*(1-t),(1-r)*(1-t),-(1-r)*(1+t),-(1+r)*(1+t),(1+r)*(1+t),(1-r)*(1+t)},
                        {-(1-r)*(1-s), -(1+r)*(1-s), -(1+r)*(1+s),-(1-r)*(1+s),(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)},
                });
                shapeFunctionsDerivatedNatural = shapeFunctionsDerivatedNatural.Multiply(0.125);
                return shapeFunctionsDerivatedNatural;
            }

        }
        public Matrix<double> GetGaussPoints(double scaleFactor, int nodeDOFS)
        {
            double gp = scaleFactor;

            if (nodeDOFS == 2)
            {
                Matrix<double> gaussNodes = DenseMatrix.OfArray(new double[,]
                {
                    {-gp,-gp},
                    {gp,-gp},
                    {gp, gp},
                    {-gp, gp},
                });
                return gaussNodes;
            }
            else
            {
                Matrix<double> gaussNodes = DenseMatrix.OfArray(new double[,]
                {
                    {-gp,-gp,-gp},
                    {gp,-gp,-gp},
                    {gp, gp,-gp},
                    {-gp, gp,-gp},
                    {-gp,-gp, gp},
                    {gp,-gp, gp},
                    {gp, gp, gp},
                    {-gp, gp,gp},
                });
                return gaussNodes;
            }
        }
    }
}
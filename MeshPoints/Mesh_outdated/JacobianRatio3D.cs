using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MathNet.Symbolics;
using MathNet.Numerics.LinearAlgebra;
using Expr = MathNet.Symbolics.SymbolicExpression;
using MeshPoints.Classes;
using System.Linq;

namespace MeshPoints.MeshQuality
{
    public class JacobianRatio3D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the JacobianRatio3D class.
        /// </summary>
        public JacobianRatio3D()
          : base("JacobianRatio3D", "jr3d",
              "Jacobian Ratio for Mesh3D element",
              "MyPlugIn", "Quality")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh3D", "m", "Insert Mesh3D class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Avg. Jacobian", "jb", "Average Jacobian ratio of all elements", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input 
            Mesh3D mesh3D = new Mesh3D();
            DA.GetData(0, ref mesh3D);

            List<double> jacobianRatios = CalculateJacobianRatioOfmesh3D(mesh3D);

            // Output
            DA.SetDataList(0, jacobianRatios);

        }

        private List<double> CalculateJacobianRatioOfmesh3D(Mesh3D mesh)
        {
            // Corner nodes of the isoparametric element
            var naturalNodes = new List<List<Double>>
            {
                new List<double> { -1, -1, -1 }, new List<double> { 1, -1, -1}, new List<double> { 1, 1, -1 }, new List<double> { -1, 1, -1 },
                new List<double> { -1, -1, 1 }, new List<double> { 1, -1, 1 }, new List<double> { 1, 1, 1 }, new List<double> { -1, 1, 1 }
            };

            var jacobianRatios = new List<double>();
            foreach (Element e in mesh.Elements)
            {
                // Global X, Y, and Z-coordinates of the corner nodes of the actual element
                List<double> gX = new List<double>()
                {
                    e.Node1.Coordinate.X, e.Node2.Coordinate.X, e.Node3.Coordinate.X, e.Node4.Coordinate.X,
                    e.Node5.Coordinate.X, e.Node6.Coordinate.X, e.Node7.Coordinate.X, e.Node8.Coordinate.X
                };
                List<double> gY = new List<double>()
                {
                    e.Node1.Coordinate.Y, e.Node2.Coordinate.Y, e.Node3.Coordinate.Y, e.Node4.Coordinate.Y,
                    e.Node5.Coordinate.Y, e.Node6.Coordinate.Y, e.Node7.Coordinate.Y, e.Node8.Coordinate.Y
                };
                List<double> gZ = new List<double>()
                {
                    e.Node1.Coordinate.Z, e.Node2.Coordinate.Z, e.Node3.Coordinate.Z, e.Node4.Coordinate.Z,
                    e.Node5.Coordinate.Z, e.Node6.Coordinate.Z, e.Node7.Coordinate.Z, e.Node8.Coordinate.Z
                };

                List<double> jacobiansOfElement = new List<double>();
                foreach (List<Double> node in naturalNodes)
                {
                    // Substitute the natural coordinates into the symbolic expression
                    var r = node[0];
                    var s = node[1];
                    var t = node[2];

                    // Partial derivatives of the shape functions
                    var N1Dr = -0.125 * (s - 1) * (t - 1);
                    var N1Ds = -0.125 * (r - 1) * (t - 1);
                    var N1Dt = -0.125 * (r - 1) * (s - 1);
                    var N2Dr = 0.125 * (s - 1) * (t - 1);
                    var N2Ds = 0.125 * (r + 1) * (t - 1);
                    var N2Dt = 0.125 * (r + 1) * (s - 1);
                    var N3Dr = -0.125 * (s + 1) * (t - 1);
                    var N3Ds = -0.125 * (r + 1) * (t - 1);
                    var N3Dt = -0.125 * (r + 1) * (s + 1);
                    var N4Dr = 0.125 * (s + 1) * (t - 1);
                    var N4Ds = 0.125 * (r - 1) * (t - 1);
                    var N4Dt = 0.125 * (r - 1) * (s + 1);
                    var N5Dr = 0.125 * (s - 1) * (t + 1);
                    var N5Ds = 0.125 * (r - 1) * (t + 1);
                    var N5Dt = 0.125 * (r - 1) * (s - 1);
                    var N6Dr = -0.125 * (s - 1) * (t + 1);
                    var N6Ds = -0.125 * (r + 1) * (t + 1);
                    var N6Dt = -0.125 * (r + 1) * (s - 1);
                    var N7Dr = 0.125 * (s + 1) * (t + 1);
                    var N7Ds = 0.125 * (r + 1) * (t + 1);
                    var N7Dt = 0.125 * (r + 1) * (s + 1);
                    var N8Dr = -0.125 * (s + 1) * (t + 1);
                    var N8Ds = -0.125 * (r - 1) * (t + 1);
                    var N8Dt = -0.125 * (r - 1) * (s + 1);

                    var sfDr = new List<double>
                    {
                        N1Dr, N2Dr, N3Dr, N4Dr, N5Dr, N6Dr, N7Dr, N8Dr
                    };
                    var sfDs = new List<double>
                    {
                        N1Ds, N2Ds, N3Ds, N4Ds, N5Ds, N6Ds, N7Ds, N8Ds
                    };
                    var sfDt = new List<double>
                    {
                        N1Dt, N2Dt, N3Dt, N4Dt, N5Dt, N6Dt, N7Dt, N8Dt
                    };

                    // Evaluates each partial derivative in the isoparametric node
                    var calcDerivs = new List<Double>
                    {
                        MultiplyLists(gX, sfDr),
                        MultiplyLists(gX, sfDs),
                        MultiplyLists(gX, sfDt),

                        MultiplyLists(gY, sfDr),
                        MultiplyLists(gY, sfDs),
                        MultiplyLists(gY, sfDt),

                        MultiplyLists(gZ, sfDr),
                        MultiplyLists(gZ, sfDs),
                        MultiplyLists(gZ, sfDt)
                    };

                    // Helper function to piecewise multiply elements of two lists of length 8
                    double MultiplyLists(List<double> a, List<double> b)
                    {
                        double sum = 0.0;
                        for (int i = 0; i<8; i++)
                        {
                            sum += (a[i] * b[i]);
                        }
                        return sum;
                    }

                    // Structure data in the form of a Jacobian matrix
                    Matrix<double> jacobianMatrix = DenseMatrixModule.ofArray2(new double[,]
                    {
                        {calcDerivs[0], calcDerivs[3], calcDerivs[6] },
                        {calcDerivs[1], calcDerivs[4], calcDerivs[7] },
                        {calcDerivs[2], calcDerivs[5], calcDerivs[8] },
                    });

                    var jacobianDeterminant = jacobianMatrix.Determinant();
                    jacobiansOfElement.Add(Math.Abs(jacobianDeterminant));

                    // This might indicate something wrong with the ordering of nodes in elements (???)
                    if (jacobianDeterminant < 0) {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "The Jacobian is negative.");
                    }
                }

                // A value of 1 denotes a cuboid element.
                double jacobianRatio = jacobiansOfElement.Min() / jacobiansOfElement.Max();

                jacobianRatios.Add(jacobianRatio);
            }
            return jacobianRatios;
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("45e31423-e9b0-492e-82bb-63a7cb0f74d8"); }
        }
    }
}
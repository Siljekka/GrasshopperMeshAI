using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MathNet.Symbolics;
using Expr = MathNet.Symbolics.SymbolicExpression;
using MeshPoints.Classes;


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
            pManager.AddGenericParameter("Avg. Jacobian", "jb", "Average Jacobian ratio of all elements", GH_ParamAccess.item);
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

            double jacobianRatio = CalculateJacobianRatioOf3dElement(mesh3D);

            // Output
            DA.SetData(0, jacobianRatio);

        }

        private double CalculateJacobianRatioOf3dElement(Mesh3D mesh3D)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Function 'CalculateJacobianRatioOf3dElement' not implemented.");

            // Corner nodes of the isoparametric element
            var naturalNodes = new List<List<Double>>
            {
                new List<double> { -1, -1, -1 }, new List<double> { 1, -1, -1}, new List<double> { 1, 1, -1 }, new List<double> { -1, 1, -1 },
                new List<double> { -1, -1, 1 }, new List<double> { 1, -1, 1 }, new List<double> { 1, 1, 1 }, new List<double> { -1, 1, 1 }
            };

            foreach (Element e in mesh3D.Elements)
            {
                // Collect the global X, Y, and Z-coordinates of the corner nodes of the actual element
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

                // Declare the three variables of the natural coordinate system
                var r = Expr.Variable("r"); // xi
                var s = Expr.Variable("s"); // nu
                var t = Expr.Variable("t"); // zeta

                // Shape functions
                Expr N1 = 0.125 * (1 - r) * (1 - s) * (1 - t);
                Expr N2 = 0.125 * (1 + r) * (1 - s) * (1 - t);
                Expr N3 = 0.125 * (1 + r) * (1 + s) * (1 - t);
                Expr N4 = 0.125 * (1 - r) * (1 + s) * (1 - t);
                Expr N5 = 0.125 * (1 - r) * (1 - s) * (1 + t);
                Expr N6 = 0.125 * (1 + r) * (1 - s) * (1 + t);
                Expr N7 = 0.125 * (1 + r) * (1 + s) * (1 + t);
                Expr N8 = 0.125 * (1 - r) * (1 + s) * (1 + t);

                // The function f(x, y, z) deconstructed into its x, y, and z-terms
                Expr x = N1 * gX[0] + N2 * gX[1] + N3 * gX[2] + N4 * gX[3] + N5 * gX[4] + N6 * gX[5] + N7 * gX[6] + N8 * gX[7];
                Expr y = N1 * gY[0] + N2 * gY[1] + N3 * gY[2] + N4 * gY[3] + N5 * gY[4] + N6 * gY[5] + N7 * gY[6] + N8 * gY[7];
                Expr z = N1 * gZ[0] + N2 * gZ[1] + N3 * gZ[2] + N4 * gZ[3] + N5 * gZ[4] + N6 * gZ[5] + N7 * gZ[6] + N8 * gZ[7];

                // The Jacobian matrix of f is the matrix of all its first-order partial derivatives
                Expr x_dr = x.Differentiate(r);
                Expr x_ds = x.Differentiate(s);
                Expr x_dt = x.Differentiate(t);
                     
                Expr y_dr = y.Differentiate(r);
                Expr y_ds = y.Differentiate(s);
                Expr y_dt = y.Differentiate(t);
                     
                Expr z_dr = z.Differentiate(r);
                Expr z_ds = z.Differentiate(s);
                Expr z_dt = z.Differentiate(t);

                List<Expr> derivatives = new List<Expr>
                {
                    x_dr, x_ds, x_dt, y_dr, y_ds, y_dt, z_ds, z_dr, z_dt
                };

                foreach (List<Double> node in naturalNodes)
                {
                    // Substitute the natural coordinates of a given node
                    var actual_values = new Dictionary<string, FloatingPoint>
                    {
                        { "r", node[0] },
                        { "s", node[1] },
                        { "t", node[2] },
                    };

                    var calculated_derivatives = new List<Double>();
                    foreach (var deriv in derivatives)
                    {
                        calculated_derivatives.Add(deriv.Evaluate(actual_values).RealValue);
                    }
                    // Structure data in completed Jacobian matrix
                }
                

            }

            return 0; // placeholder
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
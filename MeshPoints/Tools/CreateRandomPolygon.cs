using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MeshPoints.Tools
{
    public class CreateRandomPolygon : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateRandomNgon class.
        /// </summary>
        public CreateRandomPolygon()
          : base("CreateRandomNgon", "Nickname",
              "Creates a random N-gon roughly defined in the range x, y := [-1, 1].",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Edge Count", "ed", "The number of edges (integer) of the N-gon. Non-integer parameters are rounded down.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "sf", "A polygon surface with a given number of edges.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double edgeCountInput = 0;
            DA.GetData(0, ref edgeCountInput);

            if (edgeCountInput < 3) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Edge count must be atleast 3."); return; }
            int edgeCount = (int)Math.Floor(edgeCountInput);

            Brep polygon = CreateRandomNgon(edgeCount);

            DA.SetData(0, polygon);
        }

        /// <summary>
        /// Creates a random, simple n-gon with x and y between -1 and 1, where n is the given number of sides.
        /// The polygon is created using a "circle method" where the points of the polygon are placed randomly inside
        /// one of n quantiles.
        /// </summary>
        /// <param name="edgeCount"> An <see cref="int"/> denoting the number of sides of the polygon.</param>
        /// <returns>A randomly generated <see cref="Brep"/> polygon with the given number of sides.</returns>
        public Brep CreateRandomNgon(int edgeCount)
        {
            List<Point3d> contourPoints = new List<Point3d>();
            var random = new Random();

            double exclusion = 0.2;
            double quantile = 1 / edgeCount;

            for (int i = 0; i < edgeCount; i++)
            {
                double r = random.NextDouble() * (1 - exclusion) + exclusion;
                double theta = 2 * Math.PI * i / edgeCount + random.NextDouble() * quantile;
                contourPoints.Add(new Point3d(r * Math.Cos(theta), r * Math.Sin(theta), 0));
            }
            contourPoints.Add(contourPoints[0]);

            Curve curveLoop = Curve.CreateControlPointCurve(contourPoints, 1);

            Brep[] surfaces = Brep.CreatePlanarBreps(curveLoop, RhinoMath.ZeroTolerance);
            Brep surface = surfaces[0];

            return surface;
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
            get { return new Guid("0e84588d-8c79-48f3-ba73-99d47ce9cdc0"); }
        }
    }
}
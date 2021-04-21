using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Linq;
using System.Collections.Generic;

namespace MeshPoints.Tools
{
    public class NormalizeSurface : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the NormalizeSurface class.
        /// </summary>
        public NormalizeSurface()
          : base("NormalizeSurface", "Nickname",
              "Translates a Surface to origin and scales to [-1,1]",    
              "MyPlugIn", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "sf", "Surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "sf", "Translated and scaled Surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Surface inputSurface = null;
            DA.GetData(0, ref inputSurface);

            NurbsSurface surface = inputSurface.ToNurbsSurface();
            // Define parametric surface
            surface.SetDomain(0, new Interval(0, 1));
            surface.SetDomain(1, new Interval(0, 1));

            // Translate surface mid-point to origin.
            surface.Translate(Point3d.Origin - surface.PointAt(0.5, 0.5));

            // Scale surface to fit roughly between (-1, -1) and (1, 1)
            BoundingBox boundingBox = surface.GetBoundingBox(false);
            double[] boundingBoxPoints = { boundingBox.Max.MaximumCoordinate, boundingBox.Min.MaximumCoordinate };
            double maxValue = boundingBoxPoints.Max();
            surface.Scale(1 / maxValue);

            DA.SetData(0, surface);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.Icon_NormalizeSurface;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("c738214a-5e14-46c9-a13a-126ff4ae912f"); }
        }
    }
}
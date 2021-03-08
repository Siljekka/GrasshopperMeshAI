using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

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
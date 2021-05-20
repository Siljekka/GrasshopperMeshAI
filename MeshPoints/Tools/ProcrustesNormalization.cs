using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MeshPoints.Tools
{
    public class ProcrustesNormalization : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ProcrustesNormalization class.
        /// </summary>
        public ProcrustesNormalization()
          : base("ProcrustesNormalization", "Procrustes",
              "Transform a plane n-gon to best fit a unit n-gon",
              "SmartMesh", "Tools")
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
            pManager.AddGenericParameter("Surface", "sf", "Transformed Surface", GH_ParamAccess.item);
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

            NurbsSurface transformedSurface = ProcrustesSuperimposition(surface);

            DA.SetData(0, transformedSurface);
        }

        public NurbsSurface ProcrustesSuperimposition(NurbsSurface surface)
        {
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
            get { return new Guid("19464952-82db-4633-a673-7bdbdde86664"); }
        }
    }
}
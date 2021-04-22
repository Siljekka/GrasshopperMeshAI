using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.Tools
{
    public class MapPoints : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MapPoints class.
        /// </summary>
        public MapPoints()
          : base("MapPoints", "map",
              "Map points of a surface to a perfect quad/cube.",
              "MyPlugIn", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "sm", "SmartMesh Class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Normalized points", "np", "List with normalized points.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SmartMesh mesh = new SmartMesh();
            DA.GetData(0, ref mesh);

            Brep brep = mesh.Geometry.Brep;


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
            get { return new Guid("14c464b4-b62c-4dec-9420-8bf9826578f0"); }
        }
    }
}
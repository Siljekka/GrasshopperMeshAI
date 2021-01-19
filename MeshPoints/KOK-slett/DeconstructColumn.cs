using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace FirstPlugin
{
    public class DeconstructColumn : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the _5Plugincomponent class.
        /// </summary>
        public DeconstructColumn()
          : base("DeconstructColumn", "dc",
              "deconstruct properties of column class",
              "PC2021", "advanced")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("column", "c", "our own column class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("name", "n", "", GH_ParamAccess.item);
            pManager.AddLineParameter("axis", "a", "", GH_ParamAccess.item);
            pManager.AddBrepParameter("brep", "b", "", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //variables
            Column c = new Column();


            //input
            DA.GetData(0, ref c);

            //code

            //output
            DA.SetData(0,c.name);
            DA.SetData(1, c.axis);
            //DA.SetData(2, c.brep);

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
            get { return new Guid("9543b342-26bb-4636-9e90-73d78835785a"); }
        }
    }
}
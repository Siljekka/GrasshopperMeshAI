using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace FirstPlugin
{
    public class SecondPluginComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the SecondPluginComponent class.
        /// </summary>
        public SecondPluginComponent()
          : base("SecondPlugin", "sp",
              "second code ever",
              "PC2021", "simple")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("count", "c", "number of elements", GH_ParamAccess.item);
            pManager.AddNumberParameter("span", "s", "divison between number", GH_ParamAccess.item);
            pManager.AddNumberParameter("start", "st", "first number in list", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("series", "s", "series of numbers", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //variables
            double start = 0;
            double span = 10;
            int count = 5;

            //input
            DA.GetData(0, ref count); //first input
            DA.GetData(1, ref span);
            DA.GetData(2, ref start);

            //code
            List<double> series = new List<double>();

            for (int i = 0; i < count; i++)
            {
                double element = start + i * span;
                series.Add(element);
            }

            //output
            DA.SetDataList(0, series);
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
            get { return new Guid("485f81ce-4d57-406e-9bb4-853ad0b9a67e"); }
        }
    }
}
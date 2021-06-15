using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MeshPoints.Tools
{
    public class KillSwitch : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the KillSwitch class.
        /// </summary>
        public KillSwitch()
          : base("KillSwitch", "Engager",
              "Takes an input and repeatedly outputs it",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("input", "i", "input to be repeated", GH_ParamAccess.item);
            pManager.AddIntegerParameter("delay", "d", "delay in ms", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("output", "o", "output to be repeated", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double input = 0;
            int delay = 1000;
            DA.GetData(0, ref input);
            DA.GetData(1, ref delay);

            GH_Document a = this.OnPingDocument();
            a.ScheduleSolution(delay, ReturnNumber);

            var output = input;

            DA.SetData(0, output);

        }
        void ReturnNumber(GH_Document doc)
        {
            this.ExpireSolution(true);
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
            get { return new Guid("4108c327-53fe-4a4c-8113-7cf94752df25"); }
        }
    }
}
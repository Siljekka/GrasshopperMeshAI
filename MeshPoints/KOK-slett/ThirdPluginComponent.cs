using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace FirstPlugin
{
    public class ThirdPluginComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public ThirdPluginComponent()
          : base("ThirdPlugin", "tp",
              "third code ever",
              "PC2021", "simple")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("numbers", "n", "list of double numbers", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("points", "pts", "list of points", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
                 
            //variables
            List<double> numbers = new List<double>();



            //inputs
            DA.GetDataList(0, numbers);


            //code
            double limit = 3;
            List<Point3d> pts = new List<Point3d>();
            foreach (var n in numbers)
            {
                bool ifItIsSmaller = MyMethod1(n, limit);
                if (ifItIsSmaller)
                {
                    pts.Add(new Point3d(n, 0, 0));
                }
            }

            //outputs
            DA.SetDataList(0, pts);


            //my method
            bool MyMethod1(double number1, double number2)
            {

                bool result = false;
                if (number1 > number2)
                {
                    result = true;
                }
                return result;
            }
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
            get { return new Guid("a472a8cc-dd7c-4d74-b99b-25c6146e6be7"); }
        }
    }
}
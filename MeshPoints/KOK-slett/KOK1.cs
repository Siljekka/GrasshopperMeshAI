using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace BuildCreator
{
    public class BuildingParameters : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the BuildingParameters class.
        /// </summary>
        public BuildingParameters()
          : base("BuildingParameters", "Nickname",
              "Description",
              "PC2021", "YourBuilding")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Floor Height [cm]", "h", "", GH_ParamAccess.item, 300); // 0
            pManager.AddNumberParameter("Column diameter [cm]", "diameter", "", GH_ParamAccess.item, 30); // 1
            pManager.AddNumberParameter("Slab thickness [cm]", "t", "", GH_ParamAccess.item, 30); // 2
            pManager.AddNumberParameter("Column Spacing X [cm]", "h", "", GH_ParamAccess.item, 500); // 3
            pManager.AddNumberParameter("Column Spacing Y [cm]", "h", "", GH_ParamAccess.item, 500); // 4
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Parameters", "p", "List of parameters", GH_ParamAccess.list); // 0
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Assign input
            double floorH = 0.0;
            double colDiameter = 0.0;
            double slabThickness = 0.0;
            double spacingX = 0.0;
            double spacingY = 0.0;

            DA.GetData(0, ref floorH); // 0
            DA.GetData(1, ref colDiameter); // 1
            DA.GetData(2, ref slabThickness); // 2
            DA.GetData(3, ref spacingX); // 3
            DA.GetData(4, ref spacingY); // 4
            #endregion

            List<double> data = new List<double>();
            data.Add(floorH);
            data.Add(colDiameter);
            data.Add(slabThickness);
            data.Add(spacingX);
            data.Add(spacingY);

            DA.SetDataList(0, data);
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
            get { return new Guid("ff5a0848-6e5e-4879-b27e-8f974a801262"); }
        }
    }
}
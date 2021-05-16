using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.Tools
{
    public class PreviewGridInformation : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the PreviewGridInformation class.
        /// </summary>
        public PreviewGridInformation()
          : base("Preview grids", "grids",
              "Preview grid information.",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Grid Information", "grid info", "Input grid information.", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Grid group", "grid group", "Input grid group index.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Grid", "grid", "Input grid index in grid group.", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Grid groups", "smartMesh", "Merged SmartMesh", GH_ParamAccess.list);
            //pManager.AddGenericParameter("Grids", "smartMesh", "Merged SmartMesh", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            List<List<List<Node>>> gridInfo = new List<List<List<Node>>>();
            int gridGroupNum = 0;
            int gridNum = 0;

            DA.GetDataList(0, gridInfo);
            DA.GetData(1, ref gridGroupNum);
            DA.GetData(2, ref gridNum);

            // Code
            List<Node> gridGroups = new List<Node>();
            foreach (List<List<Node>> gridGroup in gridInfo)
            {
                foreach (List<Node> grid in gridGroup)
                {
                    gridGroups.AddRange(grid);
                }
            }

            // Output
            if (gridGroupNum >= gridInfo.Count) { return; }
            if (gridNum >= gridInfo[gridGroupNum].Count) { return; }

            DA.SetDataList(0, gridInfo[gridGroupNum][gridNum]);
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
            get { return new Guid("87546636-f76d-44f7-85fe-55da010e8bbe"); }
        }
    }
}
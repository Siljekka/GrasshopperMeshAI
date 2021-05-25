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
            pManager.AddGenericParameter("Grid groups", "group", "Index of grid group to return.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Grid in groups", "grid in group", "Index of grid in grid group to return.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Grid", "grid", "Index of grid to return, independent of grid groups.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Grid in groups", "grid in group", "Grid in grid group.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Grid", "grid", "Grid independent of grid groups.", GH_ParamAccess.item);
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
            int gridNumWOgroup = 0;

            DA.GetDataList(0, gridInfo);
            DA.GetData(1, ref gridGroupNum);
            DA.GetData(2, ref gridNum);

            // Code
            List<int> gridGroupsCount = new List<int>();
            List<Node> gridGroups = new List<Node>();
            List<Node> grid = new List<Node>();
            int counter = 0;
            foreach (List<List<Node>> gridGroup in gridInfo)
            {
                if (gridNumWOgroup == counter)
                { 
                
                }
                gridGroupsCount.Add(gridGroup.Count);
                foreach (List<Node> gridInGroup in gridGroup)
                {
                    gridGroups.AddRange(gridInGroup);
                    counter++;
                }
            }

            // Output
            if (gridGroupNum >= gridInfo.Count) { return; }
            if (gridNum >= gridInfo[gridGroupNum].Count) { return; }

            DA.SetDataList(0, gridInfo[gridGroupNum][gridNum]); 
            DA.SetData(1, grid);

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
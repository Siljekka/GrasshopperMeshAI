using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Text;
using System.IO;

namespace MeshPoints.Tools
{
    public class ExcelExport : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GH_DataToExcel class.
        /// </summary>
        public ExcelExport()
          : base("Excel Export", "excel",
              "Create data file that can be imported to excel.", //todo: fix description
              "SmartMesh", "Tools")
        {
        }
        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FilePath", "fp", "File path to where data are saved.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Save", "save", "True: data is written to file, False: data is not written to file.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Col1", "c1", "Column 1 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col2", "c2", "Column 2 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col3", "c3", "Column 3 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col4", "c4", "Column 4 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col5", "c5", "Column 5 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col6", "c6", "Column 6 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col7", "c7", "Column 7 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col8", "c8", "Column 8 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col9", "c9", "Column 9 in excel.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Col10", "c10", "Column 10 in excel.", GH_ParamAccess.list);

            pManager[2].Optional = true;
            pManager[3].Optional = true;
            pManager[4].Optional = true;
            pManager[5].Optional = true;
            pManager[6].Optional = true;
            pManager[7].Optional = true;
            pManager[8].Optional = true;
            pManager[9].Optional = true;
            pManager[10].Optional = true;
            pManager[11].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string filePath = "empty";
            bool writeData = false;
            List<double> col1 = new List<double>();
            List<double> col2 = new List<double>();
            List<double> col3 = new List<double>();
            List<double> col4 = new List<double>();
            List<double> col5 = new List<double>();
            List<double> col6 = new List<double>();
            List<double> col7 = new List<double>();
            List<double> col8 = new List<double>();
            List<double> col9 = new List<double>();
            List<double> col10 = new List<double>();

            DA.GetData(0, ref filePath);
            DA.GetData(1, ref writeData);
            DA.GetDataList(2, col1);
            DA.GetDataList(3, col2);
            DA.GetDataList(4, col3);
            DA.GetDataList(5, col4);
            DA.GetDataList(6, col5);
            DA.GetDataList(7, col6);
            DA.GetDataList(8, col7);
            DA.GetDataList(9, col8);
            DA.GetDataList(10, col9);
            DA.GetDataList(11, col10);


            if (!DA.GetData(0, ref filePath)) return;
            if (!DA.GetData(1, ref writeData)) return;
            if (!writeData) return;

            StringBuilder stringBuilder = new StringBuilder();
            List<string> text = new List<string>();
            int counter = 0;
            if (counter < col1.Count) { counter = col1.Count; }
            if (counter < col2.Count) { counter = col2.Count; }
            if (counter < col3.Count) { counter = col3.Count; }
            if (counter < col4.Count) { counter = col4.Count; }
            if (counter < col5.Count) { counter = col5.Count; }
            if (counter < col6.Count) { counter = col6.Count; }
            if (counter < col7.Count) { counter = col7.Count; }
            if (counter < col8.Count) { counter = col8.Count; }
            if (counter < col9.Count) { counter = col9.Count; }
            if (counter < col10.Count) { counter = col10.Count; }

            // 1. Add colums to string
            for (int i = 0; i < counter; i++)
            {
                if (DA.GetDataList(2, col1))
                    stringBuilder.Append(String.Format("{0}", col1[i]));
                if (DA.GetDataList(3, col2))
                    stringBuilder.Append(String.Format(",{0}", col2[i]));
                if (DA.GetDataList(4, col3))
                    stringBuilder.Append(String.Format(",{0}", col3[i]));
                if (DA.GetDataList(5, col4))
                    stringBuilder.Append(String.Format(",{0}", col4[i]));
                if (DA.GetDataList(6, col5))
                    stringBuilder.Append(String.Format(",{0}", col5[i]));
                if (DA.GetDataList(7, col6))
                    stringBuilder.Append(String.Format(",{0}", col6[i]));
                if (DA.GetDataList(8, col7))
                    stringBuilder.Append(String.Format(",{0}", col7[i]));
                if (DA.GetDataList(9, col8))
                    stringBuilder.Append(String.Format(",{0}", col8[i]));
                if (DA.GetDataList(10, col9))
                    stringBuilder.Append(String.Format(",{0}", col9[i]));
                if (DA.GetDataList(11, col10))
                    stringBuilder.Append(String.Format(",{0}", col10[i]));
                
                text.Add(stringBuilder.ToString());
                stringBuilder.Clear();
            }
            var file = @filePath;
            File.WriteAllLines(file, text.ToArray());
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
                return Properties.Resources.Icon_Excel;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("7218f9b3-5a2b-4f2c-abba-7c24ea98fbbf"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Text;
using System.IO;

namespace MeshPoints.Tools
{
    public class CSVExport : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateData class.
        /// </summary>
        public CSVExport()
          : base("CSV Export", "csv",
              "Export data to CSV file.",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh Class.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Structure", "cs", "1: Node coordinates (x,y) for internal nodes only + avg quality for AR and JR." +
                "2: All node coordinates (x,y) + avg quality for AR and JR." +
                "3: All node coordinates (x,y) except corner nodes + avg quality for AR and JR.", GH_ParamAccess.item);
            pManager.AddGenericParameter("FilePath", "fp", "File path to where data are saved.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Save", "w", "True: data is written to file, False: data is not written to file.", GH_ParamAccess.item);
            pManager.AddGenericParameter("AvgAR", "avgQ", "Average Aspect Ratio of mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("AvgJR", "minQ", "Average Jacobian Ratio of mesh.", GH_ParamAccess.item);;
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
            SmartMesh mesh = new SmartMesh();
            int structureType = 0;
            string filePath = "empty";
            bool writeData = false;
            double avgAR = 100;
            double avgJR = 100;

            DA.GetData(0, ref mesh);
            DA.GetData(1, ref structureType);
            DA.GetData(2, ref filePath);
            DA.GetData(3, ref writeData);
            DA.GetData(4, ref avgAR);
            DA.GetData(5, ref avgJR);

            // 0. Check input
            if (!DA.GetData(0, ref mesh)) return;
            if (!DA.GetData(1, ref structureType)) return;
            if (!DA.GetData(2, ref filePath)) return;
            if (!DA.GetData(4, ref avgAR)) return;
            if (!DA.GetData(5, ref avgJR)) return;

            StringBuilder stringBuilder = new StringBuilder();
            StringBuilder header = new StringBuilder();
            List<Node> nodes = mesh.Nodes;

            if (writeData)
            {
                if (structureType == 1)
                {
                    // 0. Add header
                    if (new FileInfo(@filePath).Length == 0)
                    {
                        header = AddHeader(mesh, filePath, 2);
                    }

                    // 1. Add quality measure to string.
                    stringBuilder.Append(String.Format("{0}", avgAR)); // AR
                    stringBuilder.Append(String.Format(",{0}", avgJR)); // JR

                    // 2. Feature for each x and y node coordinate
                    for (int i = 0; i < nodes.Count; i++)
                    {
                        if (nodes[i].BC_U | nodes[i].BC_V) { continue; }
                        double x = Math.Round(nodes[i].Coordinate.X, 2);
                        double y = Math.Round(nodes[i].Coordinate.Y, 2);
                        string text = String.Format(",{0},{1}", x, y); // temporary 0
                        stringBuilder.Append(text);
                    }
                }
                else if (structureType == 2)
                {
                    // 0. Add header
                    if (new FileInfo(@filePath).Length == 0)
                    {
                        header = AddHeader(mesh, filePath, 1);
                    }

                    // 1. Add quality measure to string.
                    stringBuilder.Append(String.Format("{0}", avgAR)); // AR
                    stringBuilder.Append(String.Format(",{0}", avgJR)); // JR

                    // 2. Feature for each x and y node coordinate
                    for (int i = 0; i < nodes.Count; i++)
                    {
                        //if (nodes[i].BC_U & nodes[i].BC_V) { continue; }
                        double x = Math.Round(nodes[i].Coordinate.X, 2);
                        double y = Math.Round(nodes[i].Coordinate.Y, 2);
                        string text = String.Format(",{0},{1}", x, y); // temporary 0
                        stringBuilder.Append(text);
                    }
                }
                else if (structureType == 3)
                {
                    // 0. Add header
                    if (new FileInfo(@filePath).Length == 0)
                    {
                        header = AddHeader(mesh, filePath, 1);
                    }

                    // 1. Add quality measure to string.
                    stringBuilder.Append(String.Format("{0}", avgAR)); // AR
                    stringBuilder.Append(String.Format(",{0}", avgJR)); // JR

                    // 2. Feature for each x and y node coordinate
                    for (int i = 0; i < nodes.Count; i++)
                    {
                        if (nodes[i].BC_U & nodes[i].BC_V) { continue; }
                        double x = Math.Round(nodes[i].Coordinate.X, 2);
                        double y = Math.Round(nodes[i].Coordinate.Y, 2);
                        string text = String.Format(",{0},{1}", x, y); // temporary 0
                        stringBuilder.Append(text);
                    }
                }

                // 4. Make CSV-file
                var data = Convert.ToString(stringBuilder);
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(@filePath, true))
                {
                    if (new FileInfo(@filePath).Length == 0) { var headerText = Convert.ToString(header); file.WriteLine(headerText); }
                    file.WriteLine(data);
                }

                //WriteTextFile(data, filePath);
            }
            else { return; }
        }
        private StringBuilder AddHeader(SmartMesh mesh, string filePath, int structureType)
        {
            StringBuilder header = new StringBuilder();

            if (structureType == 1)
            {
                header.Append("avgAR,avgJR");
                int nodes = 1;
                for (int i = 0; i < mesh.Nodes.Count; i++)
                {
                    if (mesh.Nodes[i].BC_U | mesh.Nodes[i].BC_V) { continue; }
                    header.Append(String.Format(",x{0},y{0}", nodes));
                    nodes++;
                }
            }
            else if (structureType == 2)
            {
                header.Append("avgAR,avgJR");
                int nodes = 1;
                for (int i = 0; i < mesh.Nodes.Count; i++)
                {
                    header.Append(String.Format(",x{0},y{0}", nodes));
                    nodes++;
                }
            }
            else if (structureType == 3)
            {
                header.Append("avgAR,avgJR");
                int nodes = 1;
                for (int i = 0; i < mesh.Nodes.Count; i++)
                {
                    if (mesh.Nodes[i].BC_U & mesh.Nodes[i].BC_V) { continue; }
                    header.Append(String.Format(",x{0},y{0}", nodes));
                    nodes++;
                }
            }
            return header;
        }
        // Referer til Youtube-video
        private void WriteTextFile(string variable, string filepath)
        {
            try
            {
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(@filepath, true))
                {
                    file.WriteLine(variable);
                }
            }
            catch (Exception exeption)
            {
                throw new ApplicationException("Something went wrong.", exeption);
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
                return Properties.Resources.Icon_CreateData;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6a4de1ca-ad9d-4b99-aedc-9d4be927d472"); }
        }
    }
}
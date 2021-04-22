using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Text;
using MathNet;
using System.Linq;

namespace MeshPoints
{
    public class CreateData : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateData class.
        /// </summary>
        public CreateData()
          : base("CreateCSV (Quality)", "data",
              "Export data to CSV file.",
              "MyPlugIn", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality", "q", "Average quality of mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("SmartMesh", "m", "Surface mesh.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Coordinate structure", "cs", "1: feautre for each x and y node coordinate. 2: one feature for all nodes.", GH_ParamAccess.item);
            pManager.AddGenericParameter("filePath", "fp", "File path to where data are saved", GH_ParamAccess.item);
            pManager.AddBooleanParameter("WriteData", "w", "True: data is written to file, False: data is not written to file.", GH_ParamAccess.item);

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
            double quality = 100;
            SmartMesh mesh = new SmartMesh();
            int structureType = 0;
            string filePath = "empty";
            bool writeData = false;


            DA.GetData(0, ref quality);
            DA.GetData(1, ref mesh);
            DA.GetData(2, ref structureType);
            DA.GetData(3, ref filePath);
            DA.GetData(4, ref writeData);

            // 0. Check input
            if (!DA.GetData(0, ref quality)) return;
            if (!DA.GetData(1, ref mesh)) return;
            if (!DA.GetData(2, ref structureType)) return;
            if (!DA.GetData(3, ref filePath)) return;


            StringBuilder stringBuilder = new StringBuilder();
            List<Node> nodes = mesh.Nodes;

            if (writeData)
            {
                // 1. Add quality measure to string.
                double qualityRound = Math.Round(quality, 3);
                stringBuilder.Append(String.Format("{0}", qualityRound));

                if (structureType == 1)
                {
                    // 2. Feature for each x and y node coordinate
                    for (int i = 0; i < nodes.Count; i++)
                    {
                        if (nodes[i].BC_U & nodes[i].BC_V) { continue; }
                        double x = Math.Round(nodes[i].Coordinate.X, 2);
                        double y = Math.Round(nodes[i].Coordinate.Y, 2);
                        string text = String.Format(",{0},{1},{2}", x, y, 0); // temporary 0
                        stringBuilder.Append(text);
                    }
                }
                else if (structureType == 2)
                {
                    // 2. One feature for all nodes
                    stringBuilder.Append(",");
                    for (int i = 0; i < nodes.Count; i++)
                    {
                        if (nodes[i].BC_U & nodes[i].BC_V) { continue; }
                        double x = Math.Round(nodes[i].Coordinate.X, 1) * (double) 10;
                        double y = Math.Round(nodes[i].Coordinate.Y, 1) * (double) 10;
                        double z = Math.Round(nodes[i].Coordinate.Z, 1) * (double) 10;
                        string.Format("{0:00}", x);
                        string.Format("{0:00}", y);
                        string.Format("{0:00}", z);
                        string xString = x.ToString("00");
                        string yString = y.ToString("00");
                        string zString = z.ToString("00");


                        string text = String.Format("{0}{1}{2}", xString, yString, zString);
                        stringBuilder.Append(text);

                        if (x > 99 | y > 99 | z > 99)
                        {
                           AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Check scaling. Coordinate value larger than 9.9.");
                        }
                    }


                }
                // 4. Make CSV-file
                var data = Convert.ToString(stringBuilder);
                WriteTextFile(data, filePath);
            }
            else { return; }
        }

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
                return null;
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
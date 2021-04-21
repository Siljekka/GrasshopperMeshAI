using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Text;

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
              "MyPlugIn", "Data")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality", "q", "Average quality of mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("SmartMesh", "m", "Surface mesh.", GH_ParamAccess.item);
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
            string filePath = "empty";
            bool writeData = false;

            DA.GetData(0, ref quality);
            DA.GetData(1, ref mesh);
            DA.GetData(2, ref filePath);
            DA.GetData(3, ref writeData);

            // 0. Check input
            if (!DA.GetData(0, ref quality)) return;
            if (!DA.GetData(1, ref mesh)) return;
            if (!DA.GetData(2, ref filePath)) return;


            List<Point3d> points = new List<Point3d>();
            StringBuilder stringBuilder = new StringBuilder();
            List<Node> nodes = mesh.Nodes;

            if (writeData)
            {
                // 1. Add quality measure to string.
                stringBuilder.Append(String.Format("{0}", quality));

                // 2. Normalize coordinates.
                NurbsSurface surface = mesh.Geometry.Brep.Faces[0].ToNurbsSurface();
                surface.SetDomain(0, new Interval(0, 1));
                surface.SetDomain(1, new Interval(0, 1));
                foreach (Node node in mesh.Nodes)
                {
                    surface.ClosestPoint(node.Coordinate, out double PointU, out double PointV);
                    Point3d newPoint = new Point3d(PointU, PointV, 0);
                    points.Add(newPoint);
                }

                // 3. Add normalized coordiantes to string.
                for (int i = 0; i < points.Count; i++)
                {
                    if (nodes[i].BC_U & nodes[i].BC_V) { continue; }
                    string text = String.Format(",{0},{1},{2}", points[i].X, points[i].Y, 0);
                    stringBuilder.Append(text);
                }

                // 4. Make CSV-file
                var data = Convert.ToString(stringBuilder);
                WriteTextFile(data, filePath);
                //CSV.addRecord(data, filePath);
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
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Text;
using MathNet;
using System.Linq;
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
            pManager.AddIntegerParameter("Structure", "cs", "1: feautre for each x and y node coordinate. 2: one feature for all nodes. 3: Input to NN Mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("FilePath", "fp", "File path to where data are saved.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Save", "w", "True: data is written to file, False: data is not written to file.", GH_ParamAccess.item);
            pManager.AddGenericParameter("AvgQuality", "avgQ", "Average quality of mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("MinQuality", "minQ", "Minimum quality of mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("TargetLength", "target", "Target length of mesh.", GH_ParamAccess.item);
            pManager[5].Optional = true;
            pManager[6].Optional = true;
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
            double avgQuality = 100;
            double minQuality = 100;
            double target = 0;

            DA.GetData(0, ref mesh);
            DA.GetData(1, ref structureType);
            DA.GetData(2, ref filePath);
            DA.GetData(3, ref writeData);
            DA.GetData(4, ref avgQuality);
            DA.GetData(5, ref minQuality);
            DA.GetData(6, ref target);

            // 0. Check input
            if (!DA.GetData(0, ref mesh)) return;
            if (!DA.GetData(1, ref structureType)) return;
            if (!DA.GetData(2, ref filePath)) return;
            if (!DA.GetData(4, ref avgQuality)) return;
            if (structureType == 3 & (!DA.GetData(5, ref minQuality) | !DA.GetData(6, ref target))) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Must have input minQuality and target when structureType is 3."); return; }

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
                        header = AddHeader(mesh, filePath, 1);
                    }

                    // 1. Add quality measure to string.
                    double qualityRound = Math.Round(avgQuality, 3);
                    double qualityRound1 = Math.Round(minQuality, 3);
                    stringBuilder.Append(String.Format("{0}", qualityRound)); // AR
                    stringBuilder.Append(String.Format(",{0}", qualityRound1)); // JR

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
                else if (structureType == 2)
                {
                    // 0. Add header
                    if (new FileInfo(@filePath).Length == 0)
                    {
                        header = AddHeader(mesh, filePath, 2);
                    }

                    // 1. Add quality measure to string.
                    stringBuilder.Append(String.Format("{0}", avgQuality)); // AR
                    stringBuilder.Append(String.Format(",{0}", minQuality)); // JR

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
                else if (structureType == 3)
                {
                    // 0. Add header
                    if (new FileInfo(@filePath).Length == 0)
                    {
                        header = AddHeader(mesh, filePath, 3);
                    }

                    // 1. Add corner nodes of surface to string.
                    List<Point3d> points = new List<Point3d>();
                    foreach (Node node in mesh.Nodes)
                    {
                        if (node.BC_U & node.BC_V)
                        {
                            points.Add(node.Coordinate);
                        }
                    }
                    stringBuilder.Append(String.Format("{0},{1}", points[0].X, points[0].Y));
                    stringBuilder.Append(String.Format(",{0},{1}", points[1].X, points[1].Y));
                    stringBuilder.Append(String.Format(",{0},{1}", points[3].X, points[3].Y));
                    stringBuilder.Append(String.Format(",{0},{1}", points[2].X, points[2].Y));

                    // 2. Add average quality
                    stringBuilder.Append(String.Format(",{0}", avgQuality));

                    // 3. Add minimum quality
                    stringBuilder.Append(String.Format(",{0}", Math.Round(minQuality, 3)));

                    // 4. Add target edge length
                    stringBuilder.Append(String.Format(",{0}", target));

                    // 5. Add number of nodes
                    stringBuilder.Append(String.Format(",{0}", mesh.Nodes.Count));

                    // 2. Add each x- and y- coordinate of nodes
                    for (int i = 0; i < nodes.Count; i++)
                    {
                        if (nodes[i].BC_U & nodes[i].BC_V) { continue; }
                        double x = Math.Round(nodes[i].Coordinate.X, 2);
                        double y = Math.Round(nodes[i].Coordinate.Y, 2);
                        string text = String.Format(",{0},{1}", x, y); // todo: fix Z?
                        stringBuilder.Append(text);
                    }
                }
                else if (structureType == 4)
                {
                    // 0. Add header
                    if (new FileInfo(@filePath).Length == 0)
                    {
                        header = AddHeader(mesh, filePath, 4);
                    }

                    // 1. Add corner nodes of surface to string.
                    List<Point3d> points = new List<Point3d>();
                    foreach (Node node in mesh.Nodes)
                    {
                        if (node.BC_U & node.BC_V)
                        {
                            points.Add(node.Coordinate);
                        }
                    }
                    stringBuilder.Append(String.Format("{0},{1}", points[0].X, points[0].Y));
                    stringBuilder.Append(String.Format(",{0},{1}", points[1].X, points[1].Y));
                    stringBuilder.Append(String.Format(",{0},{1}", points[3].X, points[3].Y));
                    stringBuilder.Append(String.Format(",{0},{1}", points[2].X, points[2].Y));

                    // 2. Add average quality
                    stringBuilder.Append(String.Format(",{0}", avgQuality));

                    // 3. Add minimum quality
                    stringBuilder.Append(String.Format(",{0}", Math.Round(minQuality, 3)));

                    // 4. Add number of nodes
                    stringBuilder.Append(String.Format(",{0}", mesh.Nodes.Count));

                    // 5. Add number of inner nodes
                    int counter = 0;
                    foreach (Node node in mesh.Nodes)
                    {
                        if (!node.BC_U | !node.BC_V)
                        {
                            counter++;
                        }
                    }
                    stringBuilder.Append(String.Format(",{0}", counter));

                    // 6. Add each x- and y- coordinate of nodes
                    for (int i = 0; i < nodes.Count; i++)
                    {
                        if (nodes[i].BC_U | nodes[i].BC_V) { continue; }
                        double x = Math.Round(nodes[i].Coordinate.X, 2);
                        double y = Math.Round(nodes[i].Coordinate.Y, 2);
                        string text = String.Format(",{0},{1}", x, y); // todo: fix Z?
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
                    //if (mesh.Nodes[i].BC_U & mesh.Nodes[i].BC_V) { continue; }
                    header.Append(String.Format(",x{0},y{0}", nodes));
                    nodes++;
                }
            }
            if (structureType == 2)
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
            else if (structureType == 3)
            {
                int corners = 1;
                foreach (Node node in mesh.Nodes)
                {
                    if (node.BC_U & node.BC_V)
                    {
                        if (new FileInfo(@filePath).Length == 0) { header.Append(String.Format("x{0},y{0},", corners)); }
                        corners++;
                    }
                }

                header.Append("avgQuality,minQuality");
                header.Append(",targetEdgeLength");
                header.Append(",totalNodes");
                int nodes = 1;
                for (int i = 0; i < mesh.Nodes.Count; i++)
                {
                    if (mesh.Nodes[i].BC_U & mesh.Nodes[i].BC_V) { continue; }
                    header.Append(String.Format(",x{0},y{0}", nodes));
                    nodes++;
                }
            }
            else if (structureType == 4)
            {
                // 1. For corners
                int corners = 1;
                foreach (Node node in mesh.Nodes)
                {
                    if (node.BC_U & node.BC_V)
                    {
                        if (new FileInfo(@filePath).Length == 0) { header.Append(String.Format("x{0},y{0},", corners)); }
                        corners++;
                    }
                }

                // 2. For quality and # nodes
                header.Append("avgQuality,minQuality");
                header.Append(",totalNodes");
                header.Append(",innerNodes");

                // 3. For inner nodes
                int innerNodes = 1;
                for (int i = 0; i < mesh.Nodes.Count; i++)
                {
                    if (mesh.Nodes[i].BC_U | mesh.Nodes[i].BC_V) { continue; }
                    header.Append(String.Format(",x{0},y{0}", innerNodes));
                    innerNodes++;
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
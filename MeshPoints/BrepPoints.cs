using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using Grasshopper.Kernel.Types;
using Rhino.Geometry.Intersect;
using MeshPoints.Classes;
using System.Linq;


namespace MeshPoints
{
    public class BrepPoints : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public BrepPoints()
          : base("BrepPoints", "B_pt",
              "Generates points on a Brep",
              "MyPlugIn", "Points")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "B", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("u", "u", "", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("v", "v", "", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("w", "w", "", GH_ParamAccess.item, 4);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "pt", "Points on Brep", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Variables
            Brep bp = new Brep();
            int u = 0;
            int v = 0;
            int w = 0;

            DA.GetData(0, ref bp);
            DA.GetData(1, ref u);
            DA.GetData(2, ref v);
            DA.GetData(3, ref w);


            //Code
            BrepVertexList bpVert = bp.Vertices;
            List<Point3d> pts = new List<Point3d>();

            var list = bpVert.OrderBy(f =>  f.Location.X).ToList();

            Point3d n1 = new Point3d(list[0].Location);
            Point3d n2 = new Point3d(list[2].Location);
            Point3d n3 = new Point3d(list[6].Location);
            Point3d n4 = new Point3d(list[4].Location);
            Point3d n5 = new Point3d(list[1].Location);
            Point3d n6 = new Point3d(list[3].Location);
            Point3d n7 = new Point3d(list[7].Location);
            Point3d n8 = new Point3d(list[5].Location);


            double spanW1 = (n5 - n1).Length / w;
            double spanW2 = (n6 - n2).Length / w;
            double spanW3 = (n7 - n3).Length / w;
            double spanW4 = (n8 - n4).Length / w;
            Vector3d vecW1 = new Vector3d((n5 - n1) / (n5 - n1).Length);
            Vector3d vecW2 = new Vector3d((n6 - n2) / (n6 - n2).Length);
            Vector3d vecW3 = new Vector3d((n7 - n3) / (n7 - n3).Length);
            Vector3d vecW4 = new Vector3d((n8 - n4) / (n8 - n4).Length);

            for (int k = 0; k < w + 1; k++)
            {
                Point3d p1 = new Point3d(n1.X + spanW1 * k * vecW1.X, n1.Y + spanW1 * k * vecW1.Y, n1.Z + spanW1 * k * vecW1.Z);
                Point3d p2 = new Point3d(n2.X + spanW2 * k * vecW2.X, n2.Y + spanW2 * k * vecW2.Y, n2.Z + spanW2 * k * vecW2.Z);
                Point3d p3 = new Point3d(n3.X + spanW3 * k * vecW3.X, n3.Y + spanW3 * k * vecW3.Y, n3.Z + spanW3 * k * vecW3.Z);
                Point3d p4 = new Point3d(n4.X + spanW4 * k * vecW4.X, n4.Y + spanW4 * k * vecW4.Y, n4.Z + spanW4 * k * vecW4.Z);

                double spanV1 = (p4 - p1).Length / v;
                double spanV2 = (p3 - p2).Length / v;
                Vector3d vecV1 = new Vector3d((p4 - p1) / (p4 - p1).Length);
                Vector3d vecV2 = new Vector3d((p3 - p2) / (p3 - p2).Length);

                for (int j = 0; j < v + 1; j++)
                {
                    p1 = new Point3d(p1.X + spanV1 * j * vecV1.X, p1.Y + spanV1 * j * vecV1.Y, p1.Z + spanV1 * j * vecV1.Z);
                    p2 = new Point3d(p2.X + spanV2 * j * vecV2.X, p2.Y + spanV2 * j * vecV2.Y, p2.Z + spanV2 * j * vecV2.Z);

                    double spanU = (p2 - p1).Length / u;
                    Vector3d vecU = new Vector3d((p2 - p1) / (p2 - p1).Length);

                    for (int i = 0; i < u + 1; i++)
                    {
                        Point3d pt = new Point3d(p1.X + spanU * i * vecU.X, p1.Y + spanU * i * vecU.Y, p1.Z + spanU * i * vecU.Z);
                        pts.Add(pt);
                    }
                    p1 = new Point3d(n1.X + spanW1 * k * vecW1.X, n1.Y + spanW1 * k * vecW1.Y, n1.Z + spanW1 * k * vecW1.Z);
                    p2 = new Point3d(n2.X + spanW2 * k * vecW2.X, n2.Y + spanW2 * k * vecW2.Y, n2.Z + spanW2 * k * vecW2.Z);
                }
            }

            //Output
            DA.SetDataList(0, pts);

        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.test;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("0c521427-0fca-4d1d-a5d4-fa002d459314"); }
        }
    }
}
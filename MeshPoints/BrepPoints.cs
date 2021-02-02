using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using Grasshopper.Kernel.Types;
using Rhino.Geometry.Intersect;

namespace MeshPoints
{
    public class BrepPoints : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public BrepPoints()
          : base("BrepPoints", "B_pt",
              "Makes points on a brep",
              "MyPlugIn", "Points")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "B", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("u", "u", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("v", "v", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("w", "w", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            //pManager.AddNumberParameter("Points", "pt", "Points on Brep", GH_ParamAccess.list);
            pManager.AddGenericParameter("test1", "t", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("test2", "", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Variables
            Brep bp = new Brep();
            int u = 5;
            int v = 5;
            int w = 5;
            // List<Brep> bpEdge = new List<Brep>();

            DA.GetData(0, ref bp);
            DA.GetData(1, ref u);
            DA.GetData(2, ref v);
            DA.GetData(3, ref w);


            //Code
            
            // Explode brep
            BrepSurfaceList bpSrf = bp.Surfaces;
            BrepEdgeList bpCrv = bp.Edges;
            BoundingBox bb = new BoundingBox();
            BrepVertexList bpVert = bp.Vertices;

            List<Point3d> pts11 = new List<Point3d>();
            List<Point3d> pts12 = new List<Point3d>();
            List<Point3d> pts21 = new List<Point3d>();
            List<Point3d> pts22 = new List<Point3d>();
            List<Point3d> pts31 = new List<Point3d>();
            List<Point3d> pts32 = new List<Point3d>();
            List<Point3d> pts41 = new List<Point3d>();
            List<Point3d> pts42 = new List<Point3d>();

            List<Point3d> pts1 = new List<Point3d>();
            List<Point3d> pts2 = new List<Point3d>();
            List<Point3d> pts3 = new List<Point3d>();
            List<Point3d> pts4 = new List<Point3d>();
            List<Point3d> pts5 = new List<Point3d>();
            List<Point3d> pts6 = new List<Point3d>();
            List<Point3d> pts7 = new List<Point3d>();
            List<Point3d> pts8 = new List<Point3d>();

            List<Line> lines_a = new List<Line>();
            List<Line> lines_b = new List<Line>();
            List<Line> lines = new List<Line>();


            for (int i = 0; i < u; i++)
            {
                pts11.Add(bpCrv[0].PointAtLength((bpCrv[0].GetLength() / u) * i));
                pts12.Add(bpCrv[2].PointAtLength((bpCrv[2].GetLength() / u) * i));
                Line line1 = new Line(pts11[pts11.Count-1], pts12[pts11.Count-1]);
                lines_a.Add(line1);

                pts1.Add(bpCrv[1].PointAtLength((bpCrv[1].GetLength() / u) * i));
                pts2.Add(bpCrv[3].PointAtLength((bpCrv[3].GetLength() / u) * i));
                Line line5 = new Line(pts1[pts1.Count - 1], pts2[pts2.Count - 1]);
                lines_b.Add(line5);

                var ua = Intersection.CurveCurve(line1.ToNurbsCurve(), line5.ToNurbsCurve(), 0.01, 0.01);

            }


            for (int i = 0; i < u; i++)
            {
                pts21.Add(bpCrv[4].PointAtLength((bpCrv[4].GetLength() / u) * i));
                pts22.Add(bpCrv[6].PointAtLength((bpCrv[6].GetLength() / u) * i));
                Line line2 = new Line(pts21[pts21.Count - 1], pts22[pts22.Count - 1]);
                lines.Add(line2);

                pts3.Add(bpCrv[5].PointAtLength((bpCrv[5].GetLength() / u) * i));
                pts4.Add(bpCrv[1].PointAtLength((bpCrv[1].GetLength() / u) * i));
                Line line6 = new Line(pts3[pts3.Count - 1], pts4[pts4.Count - 1]);
                lines.Add(line6);
            }

            for (int i = 0; i < u; i++)
            {
                pts31.Add(bpCrv[7].PointAtLength((bpCrv[7].GetLength() / u) * i));
                pts32.Add(bpCrv[9].PointAtLength((bpCrv[9].GetLength() / u) * i));
                Line line3 = new Line(pts31[pts31.Count - 1], pts32[pts32.Count - 1]);
                lines.Add(line3);

                pts5.Add(bpCrv[8].PointAtLength((bpCrv[8].GetLength() / u) * i));
                pts6.Add(bpCrv[5].PointAtLength((bpCrv[5].GetLength() / u) * i));
                Line line7 = new Line(pts5[pts5.Count - 1], pts6[pts6.Count - 1]);
                lines.Add(line7);
            }
            for (int i = 0; i < u; i++)
            {
                pts41.Add(bpCrv[10].PointAtLength((bpCrv[10].GetLength() / u) * i));
                pts42.Add(bpCrv[11].PointAtLength((bpCrv[11].GetLength() / u) * i));
                Line line4 = new Line(pts41[pts41.Count - 1], pts42[pts42.Count - 1]);
                lines.Add(line4);

                pts7.Add(bpCrv[8].PointAtLength((bpCrv[8].GetLength() / u) * i));
                pts8.Add(bpCrv[3].PointAtLength((bpCrv[3].GetLength() / u) * i));
                Line line8 = new Line(pts7[pts7.Count - 1], pts8[pts8.Count - 1]);
                lines.Add(line8);
            }

            /*
             bb = bpSrf[0].GetBoundingBox(true);
             Point3d pt = new Point3d(bb.Min);

             bpSrf[0].Domain(0);
             bpSrf[0].Domain(1);
             double stepU = 1 / ((double)u);
             double stepV = 1 / ((double)v);

             for (int i = 0; i < u; i++)
             {
                 for (int j = 0; j < v; j++)
                 {
                     pts.Add(bpSrf[0].PointAt(i/stepV, j/stepU));
                 }
             }

             */

            //Output
            DA.SetDataList(0, r);
            //DA.SetData(1, r);

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
            get { return new Guid("0c521427-0fca-4d1d-a5d4-fa002d459314"); }
        }
    }
}
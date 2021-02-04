using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Drawing;

namespace MeshPoints
{
    public class Mesh_Quality : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Mesh_Quality class.
        /// </summary>
        public Mesh_Quality()
          : base("Mesh Quality", "mq",
              "Mesh Quality",
              "MyPlugIn", "Quality")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Insert Mesh2D class", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality", "mq", "Mesh Quality for elements", GH_ParamAccess.list);
            pManager.AddGenericParameter("Aspect Ratio", "ar", "Sum apect ratio", GH_ParamAccess.item);
            pManager.AddGenericParameter("Skewness", "sk", "Sum skewness", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //variables
            Mesh2D m = new Mesh2D();
            Quality q = new Quality();
            List<Quality> qualityList = new List<Quality>();


            //_var for Quality Check
            List<double> dist = new List<double>(); //list distacens between vertices in a mesh face, following mesh edges CCW
            List<double> qualityValueList = new List<double>();
            double angleIdeal = 90; //ideal angle in degrees
            double angleRad = 0; //angle in radians
            List<double> angle = new List<double>(); //list of angles in a element
            int neigbourPt = 3; //variable used in skweness calcualtion

            double sumAR = 0;
            double sumSK = 0;

            //input
            DA.GetData(0, ref m);

           

            //code
            foreach (Element e in m.Elements)
            {
                List < Point3d > pts = new List<Point3d>()
                { 
                        e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                        e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                };

                
                for (int n = 0; n < pts.Count / 2; n++)
                {   
                    //Aspect Ratio
                    dist.Add(pts[n].DistanceTo(pts[n + 1])); //Add the distance between the points, following mesh edges CCW

                    Vector3d a = new Vector3d(pts[n].X - pts[n + 1].X, pts[n].Y - pts[n + 1].Y, pts[n].Z - pts[n + 1].Z); //creat a vector from a vertice to a neighbour vertice
                    Vector3d b = new Vector3d(pts[n].X - pts[n + neigbourPt].X, pts[n].Y - pts[n + neigbourPt].Y, pts[n].Z - pts[n + neigbourPt].Z); //creat a vector from a vertice to the other neighbour vertice
                    angleRad = Math.Abs(Math.Acos(Vector3d.Multiply(a, b) / (a.Length * b.Length))); //calc angles in radians between vectors
                    angle.Add(angleRad * 180 / Math.PI); //convert from rad to deg
                }
                
                dist.Sort();
                angle.Sort();

                q.element = e;
                q.AspectRatio = (dist[0] / dist[dist.Count - 1]);
                q.Skewness = 1 - Math.Max((angle[angle.Count - 1] - angleIdeal) / (180 - angleIdeal), (angleIdeal - angle[0]) / (angleIdeal));
                e.quality = q;

                sumAR += q.AspectRatio;
                sumSK += q.Skewness;
                qualityList.Add(q);
                
                q = new Quality();
                dist.Clear();
                angle.Clear();
                //pts.Clear();
            }
 
            //output
            DA.SetDataList(0, qualityList);
            DA.SetData(1, sumAR);
            DA.SetData(2, sumSK);

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
            get { return new Guid("31a76520-576a-4b37-a64f-8d51178e4be7"); }
        }
    }
}
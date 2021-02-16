using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Drawing;


namespace MeshPoints
{
    public class MeshQuality : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Mesh_Quality class.
        /// </summary>
        public MeshQuality()
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
            pManager.AddGenericParameter("Mesh2D", "m", "Insert Mesh2D class", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Quality metric", "q", "Aspect Ratio = 1, Skewness = 2, Jacobian = 3", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality", "mq", "Mesh Quality for elements", GH_ParamAccess.list);
            pManager.AddGenericParameter("Avg. Aspect Ratio", "ar", "Average aspect ratio of all elements.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Avg. Skewness", "sk", "Average skewness of all elements.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Avg. Jacobian", "jb", "Average Jacobian ratio of all elements", GH_ParamAccess.item);
            pManager.AddGenericParameter("Color Mesh", "cm", "Color map of quality check", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Variables
            //variables
            Mesh2D mesh = new Mesh2D();
            Mesh colorMesh = new Mesh();

            List<Quality> qualityList = new List<Quality>(); // list of Quality for each element in the mesh
            List<double> vertexDistance = new List<double>(); //list distances between vertices in a mesh face, following mesh edges CCW (counter-clockwise)
            List<double> elementAngles = new List<double>(); //list of angles in a element

            // determines which quality check type to color mesh with
            // 1 = aspect ratio, 2 = skewness, 3 = jacobian
            int qualityCheckType = 0; 

            double idealAngle = 90; //ideal angle in degrees
            int neighborPoint = 3; //variable used in skewness calculation

            double sumAspectRatio = 0;
            double sumSkewness = 0;
            double sumJacobian = 0; // todo: implement jacobian
            #endregion

            //input
            DA.GetData(0, ref mesh);
            DA.GetData(1, ref qualityCheckType);

            #region Calculate mesh quality
            foreach (Element e in mesh.Elements)
            {
                Quality elementQuality = new Quality();

                List < Point3d > pts = new List<Point3d>()
                { 
                        e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                        e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                };

                for (int n = 0; n < pts.Count / 2; n++)
                {   
                    //Aspect Ratio
                    vertexDistance.Add(pts[n].DistanceTo(pts[n + 1])); //Add the distance between the points, following mesh edges CCW

                    //Skewness
                    //create a vector from a vertex to a neighbouring vertex
                    Vector3d a = new Vector3d(pts[n].X - pts[n + 1].X, pts[n].Y - pts[n + 1].Y, pts[n].Z - pts[n + 1].Z); 
                    Vector3d b = new Vector3d(pts[n].X - pts[n + neighborPoint].X, pts[n].Y - pts[n + neighborPoint].Y, pts[n].Z - pts[n + neighborPoint].Z);

                    //calculate angles in radians between vectors
                    double angleRad = Math.Abs(Math.Acos(Vector3d.Multiply(a, b) / (a.Length * b.Length))); 
                    elementAngles.Add(angleRad * 180 / Math.PI); //convert from rad to deg
                }
                
                vertexDistance.Sort();
                elementAngles.Sort();

                elementQuality.AspectRatio = (vertexDistance[0] / vertexDistance[vertexDistance.Count - 1]);
                elementQuality.Skewness = 1 - Math.Max((elementAngles[elementAngles.Count - 1] - idealAngle) / (180 - idealAngle), (idealAngle - elementAngles[0]) / (idealAngle));
                
                elementQuality.element = e;
                e.quality = elementQuality;

                sumAspectRatio += elementQuality.AspectRatio;
                sumSkewness += elementQuality.Skewness;
                qualityList.Add(elementQuality);
        
                vertexDistance.Clear();
                elementAngles.Clear();
            }

            double avgAspectRatio = sumAspectRatio / mesh.Elements.Count;
            double avgSkewness = sumSkewness / mesh.Elements.Count;
            double avgJacobian = 123;
            #endregion

            #region Color the mesh based on quality type
            // 1 = aspect ratio
            if (qualityCheckType == 1) 
            {
                foreach (Quality q in qualityList)
                {
                    if (q.AspectRatio > 0.9)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                    }
                    else if (q.AspectRatio > 0.7)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                    }
                    else if (q.AspectRatio > 0.6)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                    }
                    else if (q.AspectRatio > 0)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                    }
                    colorMesh.Append(q.element.mesh);
                }
            }
            // 2 = skewness
            else if (qualityCheckType == 2)
            {
                foreach (Quality q in qualityList)
                {
                    if (q.Skewness > 0.9)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                    }
                    else if (q.Skewness > 0.7)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                    }
                    else if (q.Skewness > 0.6)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                    }
                    else if (q.Skewness > 0)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                    }
                    colorMesh.Append(q.element.mesh);
                }
            }
            #endregion

            #region Outputs
            DA.SetDataList(0, qualityList);
            DA.SetData(1, avgAspectRatio);
            DA.SetData(2, avgSkewness);
            DA.SetData(3, avgJacobian);
            DA.SetData(4, colorMesh);
            #endregion
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
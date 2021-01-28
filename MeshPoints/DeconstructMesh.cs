using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;
using Rhino.Geometry.Collections;
using System.Drawing;



namespace MeshPoints
{
    public class DeconstructMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructMesh class.
        /// </summary>
        public DeconstructMesh()
          : base("DeconstructMesh", "decM",
              "Deconstructing the mesh",
              "MyPlugIn", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Base mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Vertices", "v", "Verticies of a face", GH_ParamAccess.list);
            pManager.AddGenericParameter("Face", "f", "Faces of a mesh", GH_ParamAccess.list);
            pManager.AddGenericParameter("Colours", "c", "Mesh vertex colour", GH_ParamAccess.list);
            pManager.AddGenericParameter("Normals", "n", "Mesh normals", GH_ParamAccess.item);
            pManager.AddGenericParameter("Quality_AR", "q", "Quality of mesh", GH_ParamAccess.list);
            pManager.AddGenericParameter("Quality_SK", "q", "Quality of mesh", GH_ParamAccess.list);
            pManager.AddGenericParameter("MeshAR", "q", "Qualtity of mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Input
            Mesh m = new Mesh();
            DA.GetData(0, ref m);

            //Variables
            MeshVertexList verticies = m.Vertices; //Vertices of the mesh
            MeshFaceList faces = m.Faces; //Faces of the mesh
            //MeshVertexColorList colors = m.VertexColors;
            MeshVertexNormalList normals = m.Normals; //FaceNormals of the mesh

            //_var for mesh quality
            MeshFace face = new MeshFace(); // might delete later.
            MeshQuality quality = new MeshQuality();
            List<Point3f> pts = new List<Point3f>(); //list of vertices of a mesh face

            //_var for AR
            double AR = 0; 
            List<double> dist = new List<double>(); //list distacens between vertices in a mesh face, following mesh edges CCW
            List<double> qualityAR = new List<double>(); //Metrics of mesh quality_aspectRatio

            //_var for SK
            double SK = 0;
            List<double> angle = new List<double>();
            double angleIdeal = 90;
            double angleRad = 0;
            int neigbourPt = 0;
            List<double> qualitySK = new List<double>(); //Metrics of mesh quality_aspectRatio

            m.VertexColors.CreateMonotoneMesh(Color.White);
            int test = m.VertexColors.Count - 1;
            Random rnd = new Random(5);
            //int colorIndex = 0; slett?

            Point3f[] myArray = verticies.ToPoint3fArray();
            List<Point3f> vertList = new List<Point3f>(myArray);

            Mesh meshColor = new Mesh();
            Mesh singleMesh = new Mesh();

            #region Code
            for (int i = 0; i < faces.Count; i++) //find the distances between vertices in a face and then calculates AR
            {
                #region Get Vertices 
                face = faces.GetFace(i);
                if (face.IsQuad)
                {
                    faces.GetFaceVertices(i, out Point3f p1, out Point3f p4, out Point3f p3, out Point3f p2); // Vertices CCW of meshface //wanna insert in list pts right away
                    pts.Add(p1);
                    pts.Add(p2);
                    pts.Add(p3);
                    pts.Add(p4);

                    pts.Add(p1); //dublicate list, wanna do this more efficient
                    pts.Add(p2);
                    pts.Add(p3);
                    pts.Add(p4);
                    neigbourPt = 3;

                }
                else if (face.IsTriangle)
                {
                    //faces.GetFaceVertices(i, out Point3f p1, out Point3f p3, out Point3f p2); //Need to change
                    neigbourPt = 2;
                }
                #endregion

                #region AspectRatio
                for (int n = 0; n < pts.Count/2; n++)
                {
                    dist.Add(pts[n].DistanceTo(pts[n + 1])); //Add the distance between the points, following mesh edges CCW
                }
                dist.Sort();
                AR = (dist[0] / dist[3]); //calculates AR
                qualityAR.Add(AR);  //wanna add AR to the property quality.Aspectratio
                #endregion

                #region Skewness
                for (int n = 0; n < pts.Count/2; n++) //gjelder quads
                {
                    Vector3f a = new Vector3f(pts[n].X - pts[n+1].X, pts[n].Y - pts[n+1].Y, pts[n].Z - pts[n+1].Z); //creat a vector from a vertice to a neighbour vertice
                    Vector3f b = new Vector3f(pts[n].X - pts[n + neigbourPt].X, pts[n].Y - pts[n + neigbourPt].Y, pts[n].Z - pts[n + neigbourPt].Z); //creat a vector from a vertice to the other neighbour vertice
                    angleRad = Math.Abs( Math.Acos(Vector3f.Multiply(a, b) / (a.Length * b.Length))); //calc angles in radians between vectors
                    angle.Add( angleRad * 180/Math.PI); //convert from rad to deg
                }
                angle.Sort();
                SK = 1 - Math.Max((angle[3] - angleIdeal) / (180 - angleIdeal), (angleIdeal - angle[0]) / (angleIdeal));
                qualitySK.Add(SK);
                #endregion

                #region Color
                //Create single mesh
                List<Point3d> pts3d = new List<Point3d>();
                for (int n = 0; n < 4; n++) //change 4 to a genertic parameter (triangle/quad)
                {
                    //singleMesh.Vertices.Ad
                    Point3d pt3d = new Point3d(pts[n]);
                    pts3d.Add(pt3d);
                    singleMesh.Vertices.Add(pt3d); //add vertices to a single mesh
                   
                }
                MeshFace mf = new MeshFace();
                mf.Set(0, 1, 2, 3);
                singleMesh.Faces.AddFace(mf);
                //singleMesh.Faces.AddFace(face);
                singleMesh.FaceNormals.ComputeFaceNormals();
                singleMesh.Compact();

                //Color AR
                /*
                if (AR > 0.9)
                {
                    //singleMesh.VertexColors.CreateMonotoneMesh(Color.Green)
                    singleMesh.VertexColors.Add(Color.Green);
                    singleMesh.VertexColors.Add(Color.Green);
                    singleMesh.VertexColors.Add(Color.Green);
                    singleMesh.VertexColors.Add(Color.Green);
                }
                
                else if (AR > 0.7)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                }
                else if (AR > 0.6)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                }
                else if (AR > 0)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Red);
                }
                */

                meshColor.Append(singleMesh);

                #region slett?
                /*
                for (int n = 0; n < 3; n++)
                {
                    // get the verticies of face i, update the color... will overwrite all verticies exept the boundaries
                    
                    var indexColor = vertList.FindIndex(Predicate<pts[0]> match);
                    if (AR > 0.75)
                    {
                        m.VertexColors[colorIndex] = Color.FromArgb(0, 204, 0);
                    }
                    else if (AR > 0.5)
                    {
                        m.VertexColors[colorIndex] = Color.FromArgb(255, 255, 204);
                    }
                    else if (AR > 0.25)
                    {
                        m.VertexColors[colorIndex] = Color.FromArgb(255, 128, 0);
                    }
                    else if (AR > 0)
                    {
                        m.VertexColors[colorIndex] = Color.FromArgb(204, 0, 0);
                    }
                    colorIndex++;

                }

                */

                /*m.VertexColors.CreateMonotoneMesh(Color.White);
                Random rnd = new Random(5);
                for (int n = 0; n < m.VertexColors.Count - 1; n++)
                {
                    m.VertexColors[n] = Color.FromArgb(rnd.Next(256), rnd.Next(256), rnd.Next(256));
                }*/
                //m.VertexColors.SetColor()
                //colors.Add( m.VertexColors.SetColor(faces[0], Color.Green));

                /*
                if (AR > 0.75)
                {
                    m.VertexColors.SetColor(faces[i], Color.Green );
                }
                else if (AR > 0.5)
                {
                    m.VertexColors.SetColor(faces[i], Color.Yellow);
                }
                else if (AR > 0.25)
                {
                    m.VertexColors.SetColor(faces[i], Color.Orange);
                }
                else if (AR > 0)
                {
                    m.VertexColors.SetColor(faces[i], Color.Red);
                }
                */
                #endregion
                #endregion

                dist.Clear();
                pts.Clear();
                singleMesh = new Mesh();
            }
            meshColor.Weld(0.1);

            #endregion


            #region slett?
            /*
            for (int n = 0; n < m.VertexColors.Count - 1; n++)
            {
                if (n < 20)
                {
                    m.VertexColors[n] = Color.FromArgb(0, 0, 0);
                }
                else
                {
                    m.VertexColors[n] = Color.FromArgb(rnd.Next(256), rnd.Next(256), rnd.Next(256));
                }
            }

            if (AR > 0.75)
            {
                m.VertexColors[i] = Color.FromArgb(0, 204, 0);
            }
            else if (AR > 0.5)
            {
                m.VertexColors[i] = Color.FromArgb(255, 255, 204);
            }
            else if (AR > 0.25)
            {
                m.VertexColors[i] = Color.FromArgb(255, 128, 0);
            }
            else if (AR > 0)
            {
                m.VertexColors[i] = Color.FromArgb(204, 0, 0);
            }*/
            #endregion

            MeshVertexColorList colors = meshColor.VertexColors;




            //Output
            DA.SetDataList(0, verticies); 
            DA.SetDataList(1, faces);
            DA.SetDataList(2, colors);
            DA.SetDataList(3, normals);
            DA.SetDataList(4, qualityAR);
            DA.SetDataList(5, qualitySK);
            DA.SetData(6, singleMesh);
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
            get { return new Guid("d5276674-6b3e-499d-a3ba-ad05e2205ccb"); }
        }
    }
}
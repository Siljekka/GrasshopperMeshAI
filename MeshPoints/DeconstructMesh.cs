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
            pManager.AddGenericParameter("Quality check", "qc", "AspectRatio(1) or Skewness(2)", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Vertices", "v", "Verticies of a face", GH_ParamAccess.list); //0
            pManager.AddGenericParameter("Face", "f", "Faces of a mesh", GH_ParamAccess.list); //1
            pManager.AddGenericParameter("Normals", "n", "Mesh normals", GH_ParamAccess.list); //2
            pManager.AddGenericParameter("Quality Values", "q", "Quality of mesh", GH_ParamAccess.list); //3
            pManager.AddGenericParameter("Qualiy Colours", "c", "Mesh vertex colour", GH_ParamAccess.item); //4
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Input
            Mesh m = new Mesh();
            double check = 0;
            DA.GetData(0, ref m);
            DA.GetData(1, ref check);

            //Variables
            MeshVertexList verticies = m.Vertices; //Vertices of the mesh
            MeshFaceList faces = m.Faces; //Faces of the mesh
            MeshVertexNormalList normals = m.Normals; //FaceNormals of the mesh

            //_var for mesh quality
            MeshFace face = new MeshFace(); // might delete later. Used to desiced if quad/triangle
            MeshQuality quality = new MeshQuality();
            List<Point3f> pts = new List<Point3f>(); //list of vertices of a mesh face

            //_var for Quality Check
            List<double> dist = new List<double>(); //list distacens between vertices in a mesh face, following mesh edges CCW
            List<double> qualityValueList = new List<double>();
            double qualityValue = 0;
            double angleIdeal = 0; //ideal angle in degrees
            double angleRad = 0; //angle in radians
            List<double> angle = new List<double>(); //list of angles in a element
            int neigbourPt = 0; //variable used in skweness calcualtion
            
            //_var for coloring
            m.VertexColors.CreateMonotoneMesh(Color.White); //slett?
            Mesh meshColor = new Mesh();
            Mesh singleMesh = new Mesh();
            MeshFace mf = new MeshFace();

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
                    angleIdeal = 90;
                }
                else if (face.IsTriangle)
                {
                    //Get the vertices
                    //faces.GetFaceVertices(i, out Point3f p1, out Point3f p3, out Point3f p2); //Need to change
                    neigbourPt = 2; 
                    angleIdeal = 60;
                }
                #endregion

                #region Single Mesh
                //Create single mesh
                for (int n = 0; n < pts.Count / 2; n++) 
                {
                    singleMesh.Vertices.Add(pts[n]); //add vertices to a single mesh
                }
                mf.Set(0, 1, 2, 3);
                singleMesh.Faces.AddFace(mf);
                #endregion

                #region Quality Check
                if (check == 1) //Aspect ratio
                {
                    for (int n = 0; n < pts.Count / 2; n++)
                    {
                        dist.Add(pts[n].DistanceTo(pts[n + 1])); //Add the distance between the points, following mesh edges CCW
                    }
                    dist.Sort();
                    qualityValue = (dist[0] / dist[3]); //calculates AR
                    qualityValueList.Add(qualityValue);  //wanna add AR to the property quality.Aspectratio
                }
                else if (check == 2)
                {
                    //Only for quads
                    for (int n = 0; n < pts.Count / 2; n++) //gjelder quads
                    {
                        Vector3f a = new Vector3f(pts[n].X - pts[n + 1].X, pts[n].Y - pts[n + 1].Y, pts[n].Z - pts[n + 1].Z); //creat a vector from a vertice to a neighbour vertice
                        Vector3f b = new Vector3f(pts[n].X - pts[n + neigbourPt].X, pts[n].Y - pts[n + neigbourPt].Y, pts[n].Z - pts[n + neigbourPt].Z); //creat a vector from a vertice to the other neighbour vertice
                        angleRad = Math.Abs(Math.Acos(Vector3f.Multiply(a, b) / (a.Length * b.Length))); //calc angles in radians between vectors
                        angle.Add(angleRad * 180 / Math.PI); //convert from rad to deg
                    }
                    angle.Sort();
                    qualityValue = 1 - Math.Max((angle[3] - angleIdeal) / (180 - angleIdeal), (angleIdeal - angle[0]) / (angleIdeal)); //for quads
                    qualityValueList.Add(qualityValue);
                }
                else 
                {
                    qualityValueList.Add(0);
                }
                #endregion

                #region Color

                //Color
                //can change to switch-case
                if (qualityValue > 0.9)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Green);
                }
                else if (qualityValue > 0.7)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                }
                else if (qualityValue > 0.6)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                }
                else if (qualityValue > 0)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Red);
                }
                
                meshColor.Append(singleMesh);

                dist.Clear();
                pts.Clear();
                mf = new MeshFace();
                singleMesh = new Mesh();
            }

            #endregion


            MeshVertexColorList colors = meshColor.VertexColors;

            #endregion

            //Output
            DA.SetDataList(0, verticies); 
            DA.SetDataList(1, faces);
            DA.SetDataList(2, normals);
            DA.SetDataList(3, qualityValueList);
            DA.SetData(4, meshColor);
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
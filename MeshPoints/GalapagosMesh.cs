using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using System.Drawing;

namespace MeshPoints
{
    public class GalapagosMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GalapagosMesh class.
        /// </summary>
        public GalapagosMesh()
          : base("GalapagosMesh", "gM",
              "Optimize mesh with gene pool",
              "MyPlugIn", "Evolutionary Solving")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Base mesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Gene Pool X-dir", "qp", "Gene pool list", GH_ParamAccess.list);
            pManager.AddGenericParameter("Gene Pool Y-dir", "qp", "Gene pool list", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh quality", "mq", "Quality of total mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Input
            Mesh m = new Mesh();
            List<double> genesX = new List<double>();
            List<double> genesY = new List<double>();

            DA.GetData(0, ref m);
            DA.GetDataList(1, genesX);
            DA.GetDataList(2, genesY);

            //Variables
            double check = 1; //checking AspectRatio --------------------------------------CHANGE!
            MeshVertexList vertices = m.Vertices; //Vertices of the mesh
            MeshFaceList faces = m.Faces; //Faces of the mesh
            MeshVertexNormalList normals = m.Normals; //FaceNormals of the mesh

            //_var for mesh quality
            MeshQuality mq = new MeshQuality();
            List<Point3d> pts = new List<Point3d>(); //list of vertices of a mesh face

            //_var for Quality Check
            List<double> dist = new List<double>(); //list distacens between vertices in a mesh face, following mesh edges CCW
            List<double> qualityValueList = new List<double>();
            double angleIdeal = 0; //ideal angle in degrees
            double angleRad = 0; //angle in radians
            List<double> angle = new List<double>(); //list of angles in a element
            int neigbourPt = 0; //variable used in skweness calcualtion

            //_var for coloring
            m.VertexColors.CreateMonotoneMesh(Color.White); //slett?
            Mesh meshColor = new Mesh();
            Mesh singleMesh = new Mesh();
            MeshFace mf = new MeshFace();

            //Galapagos
            double devX = 0;
            double devY = 0;
            bool BCright = false; //temporary
            bool BCleft = false;  //temporary
            bool BCtop = false;  //temporary
            bool BCbottom = false;  //temporary
            int nx = 3; //temporary
            List<Point3d> verNew = new List<Point3d>(); //list with updraged vertices


            #region Code
            //upgrade the vertices of mesh
            //_ Check if DistanceTo > 0
            //_change from BC... to node.BC.... when node class is ready
            // add warning message if length of geneX or geneY does not match vertices.Count
            //_change to methods

            for(int i = 0; i < vertices.Count; i++)
            {
                //Deviation X
                if (genesX[i] > 0 & !BCright)
                {
                    devX = Math.Abs(vertices[i].X - vertices[i + 1].X) / 2 * genesX[i];
                }
                else if (genesX[i] < 0 & !BCleft)
                {
                    devX = Math.Abs(vertices[i].X - vertices[i - 1].X) / 2 * genesX[i];
                }
                else { devX = 0; }

                
                //Deviation Y
                if (genesY[i] > 0 & !BCtop)
                {
                    devY = Math.Abs(vertices[i].Y - vertices[i + nx].Y) / 2 * genesY[i]; //nx...
                }
                else if (genesY[i] < 0 & !BCbottom)
                {
                    devY = Math.Abs(vertices[i].Y - vertices[i - nx].Y) / 2 * genesY[i]; //nx...                
                }
                else { devY = 0; }

                //Update vertices
                verNew.Add( new Point3d(vertices[i].X + devX, vertices[i].Y + devY, vertices[i].Z + 0));              
            }


            for (int i = 0; i < faces.Count; i++) //find the distances between vertices in a face and then calculates AR
            {
                #region Get Vertices 
                mq.MeshFace = faces.GetFace(i);
                if (mq.MeshFace.IsQuad)
                {
                    pts.Add(verNew[mq.MeshFace.A]); //same index in verNew and vertices list
                    pts.Add(verNew[mq.MeshFace.B]);
                    pts.Add(verNew[mq.MeshFace.C]);
                    pts.Add(verNew[mq.MeshFace.D]);

                    pts.Add(verNew[mq.MeshFace.A]);
                    pts.Add(verNew[mq.MeshFace.B]);
                    pts.Add(verNew[mq.MeshFace.C]);
                    pts.Add(verNew[mq.MeshFace.D]);

                    neigbourPt = 3;
                    angleIdeal = 90;
                }
                else if (mq.MeshFace.IsTriangle)
                {
                    pts.Add(verNew[mq.MeshFace.A]);
                    pts.Add(verNew[mq.MeshFace.B]);
                    pts.Add(verNew[mq.MeshFace.C]);

                    pts.Add(verNew[mq.MeshFace.A]);
                    pts.Add(verNew[mq.MeshFace.B]);
                    pts.Add(verNew[mq.MeshFace.C]);

                    neigbourPt = 2;
                    angleIdeal = 60;
                }
                #endregion

                #region Single Mesh  
                //Create single mesh //Change single mesh to class..
                for (int n = 0; n < pts.Count / 2; n++)
                {
                    singleMesh.Vertices.Add(pts[n]); //add vertices to a single mesh
                }
                if (mq.MeshFace.IsQuad)
                {
                    mf.Set(0, 1, 2, 3);
                }
                else if (mq.MeshFace.IsTriangle)
                {
                    mf.Set(0, 1, 2);
                }

                singleMesh.Faces.AddFace(mf);
                singleMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                singleMesh.Normals.ComputeNormals(); //Control if needed
                singleMesh.Compact(); //to ensure that it calculate
                #endregion

                #region Quality Check
                if (check == 1) //Aspect ratio
                {
                    for (int n = 0; n < pts.Count / 2; n++)
                    {
                        dist.Add(pts[n].DistanceTo(pts[n + 1])); //Add the distance between the points, following mesh edges CCW
                    }
                    dist.Sort();
                    mq.AspectRatio = (dist[0] / dist[dist.Count - 1]);
                    qualityValueList.Add(mq.AspectRatio);
                }
                else if (check == 2)
                {
                    for (int n = 0; n < pts.Count / 2; n++)
                    {
                        Vector3d a = new Vector3d(pts[n].X - pts[n + 1].X, pts[n].Y - pts[n + 1].Y, pts[n].Z - pts[n + 1].Z); //creat a vector from a vertice to a neighbour vertice
                        Vector3d b = new Vector3d(pts[n].X - pts[n + neigbourPt].X, pts[n].Y - pts[n + neigbourPt].Y, pts[n].Z - pts[n + neigbourPt].Z); //creat a vector from a vertice to the other neighbour vertice
                        angleRad = Math.Abs(Math.Acos(Vector3d.Multiply(a, b) / (a.Length * b.Length))); //calc angles in radians between vectors
                        angle.Add(angleRad * 180 / Math.PI); //convert from rad to deg
                    }
                    angle.Sort();
                    mq.Skewness = 1 - Math.Max((angle[angle.Count - 1] - angleIdeal) / (180 - angleIdeal), (angleIdeal - angle[0]) / (angleIdeal));
                    qualityValueList.Add(mq.Skewness);
                }
                else
                {
                    qualityValueList.Add(0);
                }
                #endregion

                #region Color
                //Color
                if (qualityValueList[i] > 0.9)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Green);
                }
                else if (qualityValueList[i] > 0.7)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                }
                else if (qualityValueList[i] > 0.6)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                }
                else if (qualityValueList[i] > 0)
                {
                    singleMesh.VertexColors.CreateMonotoneMesh(Color.Red);
                }

                meshColor.Append(singleMesh);

                dist.Clear();
                angle.Clear();
                pts.Clear();
                mf = new MeshFace();
                singleMesh = new Mesh();
            }

            #endregion

            //MeshVertexColorList colors = meshColor.VertexColors;

            #endregion
            //Output
            DA.SetDataList(0, vertices);
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
            get { return new Guid("219e8033-a05c-473a-8219-f7a6c96c7256"); }
        }
    }
}
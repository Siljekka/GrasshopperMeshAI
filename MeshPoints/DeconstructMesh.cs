using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;

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
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Input
            Mesh m = new Mesh();
            MeshQuality quality = new MeshQuality();

            DA.GetData(0, ref m);
            #endregion


            #region Code
            var verticies = m.Vertices; //Vertices of the mesh
            var faces = m.Faces; //Faces of the mesh
            var normals = m.FaceNormals; //FaceNormals of the mesh
            var colours = m.VertexColors;


            // Quality: Aspect Ratio
            List<double> distances = new List<double>(); //list with distacens between points in mesh face
            List<double> qualityAR = new List<double>();

            for (int i = 0; i < faces.Count; i++) //find the distances between vertices in a face and then calculates AR
            {
                distances.Clear();  //emties the distances list
                faces.GetFaceVertices(i, out Point3f p1, out Point3f p2, out Point3f p3, out Point3f p4);

                double dist1 = p1.DistanceTo(p2);
                double dist2 = p2.DistanceTo(p3);
                double dist3 = p3.DistanceTo(p4);
                double dist4 = p4.DistanceTo(p1);
                distances.Add(dist1);
                distances.Add(dist2);
                distances.Add(dist3);
                distances.Add(dist4);
                distances.Sort();  //sorterer listen med distances

                double AR = (distances[0] / distances[3]); // calculates AR
                qualityAR.Add(AR);  // Puts the ARs in one list

                /*
                if (AR > 0.75)
                {
                    m.VertexColors.SetColor(faces[i], Color.Green);
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


            }

            //quality.aspectRatio = qualityAR;

            // Warp angle:



            #endregion


            #region Output
            DA.SetDataList(0, verticies); //Vertices
            DA.SetDataList(1, faces);
            DA.SetDataList(2, colours);
            DA.SetDataList(3, normals);
            DA.SetDataList(4, qualityAR);
            #endregion


        }

        //Methods:


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
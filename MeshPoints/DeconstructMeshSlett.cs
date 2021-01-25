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
          : base("DeconstructMesh", "Nickname",
              "Description",
              "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Mesh", "m", "Base mesh", GH_ParamAccess.item); 
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Vertices", "v", "Verticies of a face", GH_ParamAccess.list);
            pManager.AddGenericParameter("Face", "f", "Faces of a mesh", GH_ParamAccess.list);
            pManager.AddGenericParameter("Colours", "c", "Mesh vertex colour", GH_ParamAccess.list);
            pManager.AddGenericParameter("Normals", "n", "Mesh normals", GH_ParamAccess.list);
            pManager.AddGenericParameter("Quality", "q", "Quality of mesh", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Variables
            Mesh m = new Mesh();
            var faces = m.Faces;
            List<Point3d> verticies = new List<Point3d>();
            var normals = m.Normals;
            //var colours;
            //MeshQuality quality;


            //Input
            DA.GetData(0, ref m);




            //Output
            DA.SetDataList(0, verticies); //Vertices
            DA.SetDataList(1, faces); //Vertices
            //DA.SetDataList(2, colours); //Vertices
            DA.SetDataList(3, normals); //Vertices
            //DA.SetData(4, quality); //Vertices


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
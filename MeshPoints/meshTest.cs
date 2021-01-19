using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MeshPoints
{
    public class meshTest : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the meshTest class.
        /// </summary>
        public meshTest()
          : base("MeshTest", "MeshPts",
              "Create mesh between four given points",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "p", "Insert 4 points", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Point3d> pts = new List<Point3d>();
            Mesh m = new Mesh();
            
            Point3d p1 = new Point3d();
            Point3d p2 = new Point3d();
            Point3d p3 = new Point3d();
            Point3d p4 = new Point3d();

            //Input
            DA.GetDataList(0, pts);

            p1 = pts[0];
            p2 = pts[1];
            p3 = pts[2];
            p4 = pts[3]; 
            
            m.Vertices.Add(p1);
            m.Vertices.Add(p2);
            m.Vertices.Add(p3);
            m.Vertices.Add(p4);

            //mf.Set(0, 1, 2, 3);
            
            // test 1


            MeshFace mf = new MeshFace(0,1,2,3);
            m.Faces.AddFace(mf);
            m.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            m.Compact(); //to ensure that it calculate

            // Output
            DA.SetData(0, m);

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
            get { return new Guid("03d95de3-30bd-4549-80bc-0b44d557910b"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;

namespace MeshPoints
{
    public class MeshPointsComponentVer1 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public MeshPointsComponentVer1()
          : base("MeshPoints", "MeshPts",
              "Create mesh between given points",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "pts", "Insert list of points", GH_ParamAccess.tree); //Point??
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Mesh between points", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Variables
            Mesh m = new Mesh();
            Mesh allMesh = new Mesh();
            MeshFace mf = new MeshFace();
            GH_Point p1 = new GH_Point();
            GH_Point p2 = new GH_Point();
            GH_Point p3 = new GH_Point();
            GH_Point p4 = new GH_Point();

            int nx = 0;
            int ny = 0;
            
            GH_Structure<GH_Point> pts = new GH_Structure<GH_Point>();

            //Input
            DA.GetDataTree(0, out pts);

            #region Loop
            nx = pts.Branches.Count;
            ny = pts.Branches[0].Count;

            for (int i = 0; i < ny - 1; i++)
            {

                for (int j = 0; j < nx - 1; j++)
                {
                    p1 = pts[i][j];
                    p2 = pts[i+1][j];
                    p3 = pts[i+1][j + 1];
                    p4 = pts[i][j + 1];

                    p1.CastTo<Point3d>(out Point3d p1proxy);
                    p2.CastTo<Point3d>(out Point3d p2proxy);
                    p3.CastTo<Point3d>(out Point3d p3proxy);
                    p4.CastTo<Point3d>(out Point3d p4proxy);

                    m.Vertices.Add(p1proxy);
                    m.Vertices.Add(p2proxy);
                    m.Vertices.Add(p3proxy);
                    m.Vertices.Add(p4proxy);

                    mf.Set(0, 1, 2, 3);
                    m.Faces.AddFace(mf);
                    m.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
                    m.Compact(); //to ensure that it calculate
                    allMesh.Append(m);
                    m = new Mesh();
                }


            }
            #endregion

            allMesh.Weld(0.1);

            // Output
            DA.SetData(0, allMesh);
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
            get { return new Guid("235955ae-889b-4281-a7cf-6df443f17d74"); }
        }
    }
}
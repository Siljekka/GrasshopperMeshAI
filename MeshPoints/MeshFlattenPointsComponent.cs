using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;


// In order to load the result of this wizard, you will also need to
// add the output bin/ folder of this project to the list of loaded
// folder in Grasshopper.
// You can use the _GrasshopperDeveloperSettings Rhino command for that.

namespace MeshPoints
{
    public class MeshFlattenPointsComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public MeshFlattenPointsComponent()
          : base("MeshFlattenPoints", "MeshPts",
              "Create mesh between given points",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            //pManager.AddGenericParameter("Points", "pts", "Insert list of points", GH_ParamAccess.tree); //AddPointParameter?
            pManager.AddGenericParameter("Points", "pts", "Insert list of points", GH_ParamAccess.list); //change to point
            pManager.AddIntegerParameter("nx", "nx", "", GH_ParamAccess.item); //meantime 
            pManager.AddIntegerParameter("ny", "ny", "", GH_ParamAccess.item); //meantime
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
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Variables
            //List<Mesh> meshList = new List<Mesh>();
            Mesh m = new Mesh();
            Mesh allMesh = new Mesh();
            MeshFace mf = new MeshFace();
            Point3d p1 = new Point3d();
            Point3d p2 = new Point3d();
            Point3d p3 = new Point3d();
            Point3d p4 = new Point3d();
            int nx = 0;
            int ny = 0;
            int counter = 0;
            //GH_Structure<GH_Point> pts = new GH_Structure<GH_Point>();
            List<Point3d> pts = new List<Point3d>();

            //Input
            DA.GetDataList(0, pts);
            DA.GetData(1, ref nx);
            DA.GetData(2, ref ny);
            
            #region Loop
            for (int j = 0; j < ny - 1; j++)
            { 
                for (int i = 0; i < nx -1; i++)
                {
                    p1 = pts[counter];
                    p2 = pts[counter +1];
                    p3 = pts[counter + nx + 1];
                    p4 = pts[counter + nx];

                    m.Vertices.Add(p1);
                    m.Vertices.Add(p2);
                    m.Vertices.Add(p3);
                    m.Vertices.Add(p4);

                    mf.Set(0, 1, 2, 3);
                    m.Faces.AddFace(mf);

                    m.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
                    m.Compact(); //to ensure that it calculate
                    allMesh.Append(m);
                    m = new Mesh();
                    counter++;
                }
                counter++;
            }
            #endregion

            // Output
            DA.SetData(0, allMesh);
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("b421404b-70ec-4cf5-b4a7-263c781d5b51"); }
        }
    }
}

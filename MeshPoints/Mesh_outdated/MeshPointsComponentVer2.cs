using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;

//Created: 18.01.21
//Name: Silje Knutsvik Kalleberg
//File: MeshPointsComponentVer2.cs
//Operating system: Win10
//Compiler & Version: Visual Studio 2019 C# for Rhino/Grasshopper v6
//
//Component meshes between input points and return a consistant mesh.
//The version meshes without welding.

namespace MeshPoints
{
    public class MeshPointsComponentVer2 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public MeshPointsComponentVer2()
          : base("MeshPointsVer2", "MeshPts",
              "Create mesh between given points w0 welding",
              "MyPlugIn", "Outdated")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "pts", "Insert list of points", GH_ParamAccess.tree);
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
            GH_Point pt = new GH_Point();
            GH_Structure<GH_Point> pts = new GH_Structure<GH_Point>();
            int nx = 0;
            int ny = 0;
            int counter = 0;

            //Input
            DA.GetDataTree(0, out pts);

            #region Loop Vertices
            nx = pts.Branches.Count; //number points in x-dir; Control
            ny = pts.Branches[0].Count; //number points in y-dir; Control

            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    pt = pts[i][j]; //Get a point from input
                    pt.CastTo<Point3d>(out Point3d ptProxy); //From GH_Point to Point3d
                    m.Vertices.Add(ptProxy); //Add point as mesh vertice
                }
            }
            #endregion

            #region Loop MeshFace
            for (int i = 0; i < nx - 1; i++)
            {
                for (int j = 0; j < ny - 1; j++)
                {
                    m.Faces.AddFace(counter, counter + 1, counter + ny + 1, counter + ny); //Add MeshFace; Constrol nx vs ny
                    counter++;
                }
                counter++; //Skip when new "row" to mesh
            }

            #endregion

            m.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            m.Normals.ComputeNormals(); //Control if needed
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
            get { return new Guid("5ea3f97e-8d41-4d40-b0ad-a4b3b3819a61"); }
        }
    }
}


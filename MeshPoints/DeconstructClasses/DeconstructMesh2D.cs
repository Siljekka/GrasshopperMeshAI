using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

//Deconstruct the class Mesh2D

namespace MeshPoints
{
    public class Deconstruct_Mesh2D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Deconstruct_Mesh2D class.
        /// </summary>
        public Deconstruct_Mesh2D()
          : base("Deconstruct SurfaceMesh", "decSurf",
              "Deconstruct SurfaceMesh class",
              "MyPlugIn", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh2D", "m2D", "Mesh2D class", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "e", "List of elements", GH_ParamAccess.list); //0
            pManager.AddGenericParameter("Nodes", "n", "List of nodes", GH_ParamAccess.list); //1
            pManager.AddGenericParameter("Mesh", "m", "Mesh", GH_ParamAccess.item); //2
            pManager.AddGenericParameter("Geometry", "g", "Geometry information", GH_ParamAccess.item); //3
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            Mesh2D m = new Mesh2D();
            DA.GetData(0, ref m);

            /*
            List<Point3d> points = new List<Point3d>();
            NurbsSurface surface = m.Geometry.Brep.Faces[0].ToNurbsSurface();
            foreach (Node node in m.Nodes)
            {
                surface.ClosestPoint(node.Coordinate, out double PointU, out double PointV);
                Point3d newPoint = new Point3d(PointU/surface.Domain(0).T1, PointV/surface.Domain(1).T1, 0);
                points.Add(newPoint);
            }*/
            
            


            //output
            DA.SetDataList(0, m.Elements);
            DA.SetDataList(1, m.Nodes);
            DA.SetData(2, m.mesh);
            DA.SetData(3, m.Geometry);
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
                return Properties.Resources.Icon_DeconstructSurfaceMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6cb2b0ff-e396-4fb5-a0b9-33ff368636ac"); }
        }
    }
}
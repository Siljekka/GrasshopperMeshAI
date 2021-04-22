using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.DeconstructClasses
{
    public class DeconstructSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructMesh3d class.
        /// </summary>
        public DeconstructSmartMesh()
          : base("Deconstruct SmartMesh", "decMesh",
              "Deconstructing SmartMesh class",
              "MyPlugIn", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "mesh", "SmartMesh class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "e", "List of elements", GH_ParamAccess.list); 
            pManager.AddGenericParameter("Nodes", "n", "List of nodes", GH_ParamAccess.list); 
            pManager.AddGenericParameter("Geometry", "geo", "Geometry information", GH_ParamAccess.item); 
            pManager.AddGenericParameter("Mesh", "m", "Mesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Normalized points", "", "", GH_ParamAccess.list);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            SmartMesh mesh = new SmartMesh();
            DA.GetData(0, ref mesh);

            List<Point3d> points = new List<Point3d>();
            NurbsSurface surface = mesh.Geometry.Brep.Faces[0].ToNurbsSurface();
            surface.SetDomain(0, new Interval(0, 1));
            surface.SetDomain(1, new Interval(0, 1));
            foreach (Node node in mesh.Nodes)
            {
                surface.ClosestPoint(node.Coordinate, out double PointU, out double PointV);
                Point3d newPoint = new Point3d(PointU, PointV, 0);
                points.Add(newPoint);
            }

            // Output
            DA.SetDataList(0, mesh.Elements);
            DA.SetDataList(1, mesh.Nodes);
            DA.SetData(2, mesh.Geometry);
            DA.SetData(3, mesh.Mesh);
            DA.SetDataList(4, points);
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
                return Properties.Resources.Icon_DeconstructSolidMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("97c30c27-48c9-41ac-b09d-d02f80e806f6"); }
        }
    }
}
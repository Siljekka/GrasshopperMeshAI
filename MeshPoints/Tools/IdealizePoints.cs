using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.Tools
{
    public class IdealizePoints : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MapPoints class.
        /// </summary>
        public IdealizePoints()
          : base("Idealize Points", "map",
              "Map points of a surface to a perfect quad/cube.",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "sm", "SmartMesh Class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "np", "SmartMesh Class with normalized node coordinates.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SmartMesh mesh = new SmartMesh();
            DA.GetData(0, ref mesh);


            NurbsSurface surface = mesh.Geometry.Brep.Faces[0].ToNurbsSurface();
            surface.SetDomain(0, new Interval(0, 1));
            surface.SetDomain(1, new Interval(0, 1));
            foreach (Node node in mesh.Nodes)
            {
                surface.ClosestPoint(node.Coordinate, out double PointU, out double PointV);
                Point3d newPoint = new Point3d(PointU, PointV, 0);
                int id = mesh.Nodes.IndexOf(node);
                mesh.Nodes[id].Coordinate = newPoint;
            }

            DA.SetData(0, mesh);
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
            get { return new Guid("14c464b4-b62c-4dec-9420-8bf9826578f0"); }
        }
    }
}
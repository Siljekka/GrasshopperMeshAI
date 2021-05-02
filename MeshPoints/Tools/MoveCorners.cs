using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;

namespace MeshPoints.Tools
{
    public class MoveCorners : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MoveCorners class.
        /// </summary>
        public MoveCorners()
          : base("Move Corners", "mc",
              "Move corners of a surface with foure edges",
              "MyPlugIn", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "srf", "Surface to modify", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Gene pool for translation in u direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("v genes", "qp", "Gene pool for translation in v direction", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "srf", "¨Modified Surface ", GH_ParamAccess.item);
            pManager.AddGenericParameter("New corner points", "cp", "Modified corner points", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            List<double> genesU = new List<double>();
            List<double> genesV = new List<double>();
            Brep refSrf = new Brep();
            DA.GetData(0, ref refSrf);
            DA.GetDataList(1, genesU);
            DA.GetDataList(2, genesV);

            // Code
            // 1. Check input OK.
            if (!DA.GetData(0, ref refSrf)) return;
            if (genesU.Count < 4) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Increase genes in u-dir."); return; }
            if (genesV.Count < 4) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Increase genes in v-dir."); return; }

            List<Point3d> cornerPts = GetCornerPoints(refSrf);
            List<Point3d> newLocation = ModifyCornerLocation(cornerPts, genesU, genesV);
            NurbsSurface srf = NurbsSurface.CreateFromCorners(newLocation[0], newLocation[1], newLocation[2], newLocation[3]);

            // Output
            DA.SetData(0, srf);
            DA.SetDataList(1, newLocation);
        }

        // Methods
        private List<Point3d> ModifyCornerLocation(List<Point3d> refPoints, List<double> genesU, List<double> genesV)
        {
            Vector3d uVector1 = refPoints[1] - refPoints[0];
            Vector3d uVector2 = refPoints[2] - refPoints[3];

            Vector3d vVector1 = refPoints[3] - refPoints[0];
            Vector3d vVector2 = refPoints[2] - refPoints[1];

            // first point
            Point3d newPoint1 = new Point3d
               (
               refPoints[0].X + (uVector1.X * genesU[0] + vVector1.X * genesV[0]) * 0.5 * 0.875, // 0.5 for half the length, 0.875 for a space free of nodes
               refPoints[0].Y + (uVector1.Y * genesU[0] + vVector1.Y * genesV[0]) * 0.5 * 0.75,
               0
               );

            Point3d newPoint2 = new Point3d
                (
                refPoints[1].X + (uVector1.X * genesU[1] + vVector2.X * genesV[1]) * 0.5 * 0.875, // 0.5 for half the length, 0.875 for a space free of nodes
                refPoints[1].Y + (uVector1.Y * genesU[1] + vVector2.Y * genesV[1]) * 0.5 * 0.75,
                0
                );

            Point3d newPoint3 = new Point3d
                (
                refPoints[2].X + (uVector2.X * genesU[2] + vVector2.X * genesV[2]) * 0.5 * 0.875, // 0.5 for half the length, 0.875 for a space free of nodes
                refPoints[2].Y + (uVector2.Y * genesU[2] + vVector2.Y * genesV[2]) * 0.5 * 0.75,
                0
                );
            Point3d newPoint4 = new Point3d
                (
                refPoints[3].X + (uVector2.X * genesU[3] + vVector1.X * genesV[3]) * 0.5 * 0.875, // 0.5 for half the length, 0.875 for a space free of nodes
                refPoints[3].Y + (uVector2.Y * genesU[3] + vVector1.Y * genesV[3]) * 0.5 * 0.75,
                0
                );
            List<Point3d> newLocation = new List<Point3d>() { newPoint1, newPoint2, newPoint3, newPoint4 };
            return newLocation;

        }
        private List<Point3d> GetCornerPoints(Brep refSrf)
        {
            List<Point3d> cornerPts = new List<Point3d>();
            Brep[] planarBrep = Brep.CreatePlanarBreps(refSrf.Edges, 0.0001); // make planar brep on floor i     
            NurbsSurface refNurbsSrf = NurbsSurface.CreateNetworkSurface(planarBrep[0].Edges, 0, 0.0001, 0.0001, 0.0001, out int error); // make planar brep to nurbssurface

            var uDomain = refNurbsSrf.Domain(0);
            var vDomain = refNurbsSrf.Domain(1);

            cornerPts.Add(refNurbsSrf.PointAt(uDomain[0], vDomain[0]));  // make point on surface
            cornerPts.Add(refNurbsSrf.PointAt(uDomain[1], vDomain[0]));  // make point on surface
            cornerPts.Add(refNurbsSrf.PointAt(uDomain[1], vDomain[1]));  // make point on surface
            cornerPts.Add(refNurbsSrf.PointAt(uDomain[0], vDomain[1]));  // make point on surface

            return cornerPts;
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
            get { return new Guid("b316cbf7-b649-473d-8818-17279bac573d"); }
        }
    }
}
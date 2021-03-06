﻿using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry.Intersect;

namespace MeshPoints.CreateMesh
{
    public class SolidSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateSolidMesh_Sweep class.
        /// </summary>
        public SolidSmartMesh()
          : base("Solid SmartMesh", "SolidSM",
              "Creates a SmartMesh of hexahedral-elements.",
              "SmartMesh", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "brep", "Brep.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Index", "i", "Index of the bottom face of the brep.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("u", "u", "Number element in u-direction.", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("v", "v", "Number element in v-direction.", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("w", "w", "Number element in w-direction.", GH_ParamAccess.item, 4);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "mesh", "Mesh (hexahedral-elements).", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Brep brep = new Brep();
            int bottomFace = 0;
            int u = 0;
            int v = 0;
            int w = 0;
            DA.GetData(0, ref brep);
            DA.GetData(1, ref bottomFace);
            DA.GetData(2, ref u);
            DA.GetData(3, ref v);
            DA.GetData(4, ref w);


            // 1. Check input OK
            if (!DA.GetData(0, ref brep)) return;
            if (!DA.GetData(1, ref bottomFace)) return;
            if (u == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "u can not be zero."); return; }
            if (v == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "v can not be zero."); return; }
            if (w == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "w can not be zero."); return; }

            // 2. Assign geometrical properties to mesh
            SmartMesh smartMesh = new SmartMesh();
            Geometry brepGeometry = new Geometry(brep, bottomFace);
            smartMesh.nu = u + 1;
            smartMesh.nv = v + 1;
            smartMesh.nw = w + 1;
            smartMesh.Type = "Solid";
            smartMesh.Geometry = brepGeometry;
            
            // 3. Find Rails
            List<Curve> rails = new List<Curve>() { brepGeometry.Edges[8], brepGeometry.Edges[9], brepGeometry.Edges[10], brepGeometry.Edges[11] };

            // 4. Divide each brep edge in w direction (rail) into nw points.
            DataTree<Point3d> railPoints = DivideRailIntoWPoints(rails, brepGeometry.Faces[0], w);

            // 5. Create NurbsSurface for each nw-floor
            DataTree<Curve> intersectionCurve = GetIntersectionCurveBrepAndRailPoints(railPoints, brep);
            if (intersectionCurve == null) return;
            List<NurbsSurface> surfaceAtW = CreateNurbSurfaceAtEachFloor(intersectionCurve);

            // 6. Make grid of points in u and v direction at leven nw
            List<Point3d> meshPoints = CreateGridOfPointsAtEachFloor(u, v, surfaceAtW, railPoints);
            
            // 7. Create nodes 
            smartMesh.CreateNodes(meshPoints, u, v, w); // assign Coordiantes, GlobalId and Boundary Conditions

            //8. create elements
            smartMesh.CreateHexElements();

            //9. Create global mesh
            smartMesh.CreateMesh();

            // Output
            DA.SetData(0, smartMesh);
            DA.SetData(1, smartMesh.Mesh);
        }

        #region Methods
        /// <summary>
        /// Divide each brep edge w-direction into nw points. The brep edges in w-direction are named rail.
        /// </summary>
        /// <returns> DataTree with points on each rail. Branch: floor level.</returns>
        private DataTree<Point3d> DivideRailIntoWPoints(List<Curve> rails, BrepFace brepBottomFace, int w)
        {
            DataTree<Point3d> railPoints = new DataTree<Point3d>(); 

            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(w, true, out Point3d[] pt);
                List<Point3d> point = pt.ToList();
                brepBottomFace.ClosestPoint(point[0], out double PointOnCurveU, out double PointOnCurveV);
                Point3d testPoint = brepBottomFace.PointAt(PointOnCurveU, PointOnCurveV);
                Vector3d distanceToFace = testPoint - point[0];

                if (distanceToFace.Length > 0.001) { point.Reverse(); }
                for (int j = 0; j < point.Count; j++)
                {
                    railPoints.Add(point[j], new GH_Path(j)); // tree with points on each rail. Branch: floor
                }
            }

            // Check if the rails must be re-oredered to generate elements with nodes counting ccw
            Curve testCurve = Curve.CreateControlPointCurve(railPoints.Branch(0), 1);
            Vector3d direction = railPoints.Branch(w)[0] - railPoints.Branch(0)[0];
            string curveOrientation = testCurve.ClosedCurveOrientation(direction).ToString();
            if (curveOrientation == "Clockwise")
            {
                for (int i = 0; i < railPoints.BranchCount; i++)
                {
                    railPoints.Branch(i).Reverse();
                }
            }
            return railPoints;
        }

        /// <summary>
        /// Get intersectionCurve between Brep and plane made from RailPoints. Curve is closed.
        /// </summary>
        /// <returns> DataTree with closed curves at each Floor. Branch: Floor Level </returns>
        private DataTree<Curve> GetIntersectionCurveBrepAndRailPoints(DataTree<Point3d> railPoints, Brep brep)
        {
            DataTree<Curve> intersectionCurve = new DataTree<Curve>();
            List<Plane> planes = new List<Plane>();

            for (int i = 0; i < railPoints.BranchCount; i++)
            {
                Vector3d vec1 = railPoints.Branch(i)[1] - railPoints.Branch(i)[0];
                Vector3d vec2 = railPoints.Branch(i)[3] - railPoints.Branch(i)[0];
                Plane plane = new Plane(railPoints.Branch(i)[0], vec1, vec2);
                Intersection.BrepPlane(brep, plane, 0.0001, out Curve[] iCrv, out Point3d[] iPt); // make intersection curve between brep and plane on floor i
                List<Curve> intCrv = iCrv.ToList();
                
                for (int j = 0; j < intCrv.Count; j++) { intCrv[j].MakeClosed(0.0001); intersectionCurve.Add(intCrv[j], new GH_Path(i)); }  // make curve closed and add to intersectionCurve
                if (intersectionCurve.Branch(i).Count != 1) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Brep input is not OK."); return null; }
                if (intersectionCurve.Branch(i)[0].ClosedCurveOrientation(plane).ToString() == "Clockwise") { intersectionCurve.Branch(i)[0].Reverse(); }

                planes.Add(plane);
                intCrv.Clear();
            }
            return intersectionCurve;
        }

        /// <summary>
        /// Create NurbSurface on each Floor with the contour of intersectionCurves
        /// </summary>
        /// <returns> List with NurbSurface</returns>
        private List<NurbsSurface> CreateNurbSurfaceAtEachFloor(DataTree<Curve> intersectionCurve)
        {
            List<NurbsSurface> surfaceAtNw = new List<NurbsSurface>();
            for (int i = 0; i < intersectionCurve.BranchCount; i++)
            {
                List<Brep> planarBrep = Brep.CreatePlanarBreps(intersectionCurve.Branch(i), 0.0001).ToList(); // make planar brep on floor i     

                for (int j = 0; j < planarBrep.Count; j++)
                {
                    NurbsSurface nurbsSurface = NurbsSurface.CreateNetworkSurface(planarBrep[j].Edges, 0, 0.0001, 0.0001, 0.0001, out int error); // make planar brep to nurbssurface
                    surfaceAtNw.Add(nurbsSurface);
                }
                planarBrep.Clear();
            }
            return surfaceAtNw;
        }

        /// <summary>
        /// Copy a grid of points in u and v direction onto every floor.
        /// </summary>
        /// <returns> DataTree with points on Brep. Branch: floor level.</returns>
        private List<Point3d> CreateGridOfPointsAtEachFloor(int u, int v, List<NurbsSurface> surfaceAtNw, DataTree<Point3d> railPoints)
        {
            List<Point3d> pointList = new List<Point3d>();
            Vector3d direction = Vector3d.Zero;
            for (int i = 0; i < surfaceAtNw.Count; i++) // loop floors
            {
                if (i == surfaceAtNw.Count - 1)
                {
                    direction = (railPoints.Branch(i-1)[0] - railPoints.Branch(i)[0]);
                    direction.Reverse();
                }
                else
                {
                    direction = (railPoints.Branch(i + 1)[0] - railPoints.Branch(i)[0]);
                }
                List<Point3d> pt = CreateGridOfPointsUV(u, v, surfaceAtNw[i], railPoints.Branch(i)[0], direction);
                pointList.AddRange(pt);
            }
            return pointList;
        }

        /// <summary>
        /// Makes grid of points in U and V direction
        /// </summary>
        /// <returns> List of points in U and V direction</returns>
        private List<Point3d> CreateGridOfPointsUV(int u, int v, NurbsSurface surface, Point3d railPoint, Vector3d direction)
        {
            var domU = surface.Domain(0);
            var domV = surface.Domain(1);
            surface.UVNDirectionsAt(0, 0, out Vector3d uDir, out Vector3d vDir, out Vector3d nDir);
            if (Point3d.Subtract(surface.PointAt(domU.T1, domV.T0), railPoint).Length < 0.1)
            {
                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                {
                    domU = new Interval(domU.T1, domU.T0);
                    List<Point3d> points = GeneratePoints(0, surface, domU, domV, u, v);
                    return points; 
                }
                else
                {
                    domU = new Interval(domU.T1, domU.T0);
                    List<Point3d> points = GeneratePoints(1, surface, domU, domV, u, v);
                    return points;
                }
            }
            else if (Point3d.Subtract(surface.PointAt(domU.T0, domV.T1), railPoint).Length < 0.1)
            {
                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                {
                    domV = new Interval(domV.T1, domV.T0);
                    List<Point3d> points = GeneratePoints(0, surface, domU, domV, u, v);
                    return points;
                }
                else
                {
                    domV = new Interval(domV.T1, domV.T0);
                    List<Point3d> points = GeneratePoints(1, surface, domU, domV, u, v);
                    return points;
                }

            }
            else if (Point3d.Subtract(surface.PointAt(domU.T1, domV.T1), railPoint).Length < 0.1)
            {
                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                {
                    domU = new Interval(domU.T1, domU.T0);
                    domV = new Interval(domV.T1, domV.T0);
                    List<Point3d> points = GeneratePoints(1, surface, domU, domV, u, v);
                    return points;
                }
                else 
                {
                    domU = new Interval(domU.T1, domU.T0);
                    domV = new Interval(domV.T1, domV.T0);
                    List<Point3d> points = GeneratePoints(0, surface, domU, domV, u, v);
                    return points;
                }
            }
            else 
            {
                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                {
                    List<Point3d> points = GeneratePoints(1, surface, domU, domV, u, v);
                    return points;
                }
                else
                {
                    List<Point3d> points = GeneratePoints(0, surface, domU, domV, u, v);
                    return points;
                }
            }
           
        }
        private List<Point3d> GeneratePoints(int Case, NurbsSurface surface, Interval domainU, Interval domainV, int u, int v)
        {
            List<Point3d> pt = new List<Point3d>();

            if (Case == 0)
            {
                double stepU = domainU.Length / (double)u;
                double stepV = domainV.Length / (double)v;

                double pointU = domainU.T0;
                double pointV = domainV.T0;
                for (double j = 0; j <= v; j++)
                {
                    for (double k = 0; k <= u; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV));  // make point on surface
                        pointU = pointU + stepU;
                    }
                    pointV = pointV + stepV;
                    pointU = domainU.T0;
                }
            }
            else if (Case == 1)
            {
                double stepU = domainU.Length / (double)v;
                double stepV = domainV.Length / (double)u;

                double pointU = domainU.T0;
                double pointV = domainV.T0;
                for (double j = 0; j <= v; j++)
                {
                    for (double k = 0; k <= u; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV));  // make point on surface
                        pointV = pointV + stepV;
                    }
                    pointU = pointU + stepU;
                    pointV = domainV.T0;
                }
            }
            else { pt = null; }
            return pt;
        }


        #endregion

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Icon_Solid;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("c7862868-a0a6-46f7-bb93-88c7ab55e3fb"); }
        }
    }
}
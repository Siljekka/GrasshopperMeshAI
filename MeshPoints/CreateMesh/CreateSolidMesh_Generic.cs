using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry.Intersect;
using Rhino.Geometry.Collections;

namespace MeshPoints.CreateMesh
{
    public class CreateSolidMesh_Generic : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateSolidMesh_Sweep class.
        /// </summary>
        public CreateSolidMesh_Generic()
          : base("CreateSolidMesh", "solid",
              "Generate solid mesh. Independent on how surface composing the brep is made.",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "bp", "Brep", GH_ParamAccess.item);
            pManager.AddIntegerParameter("BottomSurface", "bottom", "Index of bottom surface of Brep", GH_ParamAccess.item);
            pManager.AddIntegerParameter("u", "u", "division in u direction", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("v", "v", "division in v direction", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("w", "w", "division in w direction", GH_ParamAccess.item, 4);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SolidMesh", "solid", "Creates a SolidMesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("mesh", "m", "Mesh (solid elements).", GH_ParamAccess.item);
            //pManager.AddGenericParameter("test2", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            /* Todo:
             * Legg til/endre INP-valid.
            */
            // Input
            Brep brep = new Brep();
            int bottomFace = 0;
            int nu = 0;
            int nv = 0;
            int nw = 0;
            DA.GetData(0, ref brep);
            DA.GetData(1, ref bottomFace);
            DA.GetData(2, ref nu);
            DA.GetData(3, ref nv);
            DA.GetData(4, ref nw);

            #region Variables
            //Variables
            //Mesh3D solidMesh = new Mesh3D();
            Mesh globalMesh = new Mesh();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<NurbsSurface> surfaceAtNw = new List<NurbsSurface>();
            DataTree<Curve> intersectionCurve = new DataTree<Curve>();
            DataTree<Point3d> railPoints = new DataTree<Point3d>();
            DataTree<Point3d> meshPoints = new DataTree<Point3d>();
            List<Plane> planes = new List<Plane>();
            #endregion

            if (!DA.GetData(0, ref brep)) return;
            if (!DA.GetData(1, ref bottomFace)) return;
            if (nu == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nu can not be zero."); return; }
            if (nv == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nv can not be zero."); return; }
            if (nw == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nw can not be zero."); return; }



            /* Todo:
             * Fjern bottomFace og brep og ha heller brepGeometry som input
             * kan erstatte FindRails
             * vurder å lage en tuple for sortbrepProperties, evt legg metode inn i Geometry-klassen
            */


            // 0. Add properties to Geometry
            Geometry brepGeometry = new Geometry(brep, SortBrepFaces(brep, bottomFace), SortBrepEdges(brep, bottomFace), SortBrepVertex(brep, bottomFace));

            // 1. Find Rails
            List<Curve> rails = FindRails(brep, bottomFace);

            // 2. Divide each brep edge in w direction (rail) into nw points.
            railPoints = DivideRailIntoNwPoints(rails, brep, nw, bottomFace);

            // 3. Create NurbsSurface for each nw-floor
            intersectionCurve = GetIntersectionCurveBrepAndRailPoints(railPoints, brep);
            if (intersectionCurve == null) return;
            surfaceAtNw = CreateNurbSurfaceAtEachFloor(intersectionCurve);
         
            // 4. Make grid of points in u and v direction at leven nw
            meshPoints = CreateGridOfPointsAtEachFloor(nu, nv, surfaceAtNw, railPoints);
            
            //5. Create nodes 
            nodes = CreateNodes(meshPoints, nu, nv, nw); // assign Coordiantes, GlobalId and Boundary Conditions


            //6. create elements
            elements = CreateHexElements(meshPoints, nodes, nu, nv); // assign ElementId, ElementMesh and Nodes incl. Coordiantes, GlobalId, LocalId and Boundary Conditions), elementId, elementMesh.

            // 6. Check if brep can be interpret by Abaqus
            //IsBrepCompatibleWithAbaqus(elements[0], solidMesh);

            //7. Create global mesh
            globalMesh = CreateGlobalMesh(elements);

            // 8. Add properties to SolidMesh
            Mesh3D solidMesh = new Mesh3D(nu, nv, nw, nodes, elements, globalMesh);
            solidMesh.Geometry = brepGeometry;
            solidMesh.inp = true;

            // Output
            DA.SetData(0, solidMesh);
            DA.SetData(1, solidMesh.mesh);
        }

        #region Methods
        private List<BrepFace> SortBrepFaces(Brep brep, int bottomFace)
        {
            List<BrepFace> faceSorted = new List<BrepFace>();
            // Find top and bottom edge
            List<BrepFace> brepFace = brep.Faces.ToList();
            List<int> indexAdjecentFaces = (brepFace[bottomFace].AdjacentFaces()).ToList();
            indexAdjecentFaces.Add(bottomFace);
            for (int i = 0; i < brepFace.Count; i++)
            {
                if (!indexAdjecentFaces.Contains(brepFace.IndexOf(brepFace[i])))
                {
                    faceSorted.Add(brepFace[bottomFace]); // bottom face
                    faceSorted.Add(brepFace[i]); // top face
                    continue;
                }
            }
            indexAdjecentFaces.Remove(bottomFace);
            foreach (int index in indexAdjecentFaces) { faceSorted.Add(brepFace[index]); }
            return faceSorted;
        }
        private List<BrepEdge> SortBrepEdges(Brep brep, int bottomFace)
        {
            List<BrepFace> faceSorted = SortBrepFaces(brep, bottomFace);
            List<BrepEdge> edgeSorted = new List<BrepEdge>();
            List<int> indexAdjecentEdges = new List<int>();

            // Find edges connected to top and bottom face.
            indexAdjecentEdges.AddRange(faceSorted[0].AdjacentEdges().ToList());
            indexAdjecentEdges.AddRange(faceSorted[1].AdjacentEdges().ToList());

            // Add edges to list.
            foreach (int index in indexAdjecentEdges) { edgeSorted.Add(brep.Edges[index]); }

            // Find rest of edges
            List<BrepEdge> brepEdgesCopy = brep.Edges.ToList();
            List<Curve> rails = new List<Curve>(brepEdgesCopy);
            foreach (int index in indexAdjecentEdges) { rails.Remove(brepEdgesCopy[index]); }

            // Add rest of edges to list.
            foreach (BrepEdge edge in rails) { edgeSorted.Add(edge); }
            return edgeSorted;
        }
        private List<BrepVertex> SortBrepVertex(Brep brep, int bottomFace)
        {
            List<BrepVertex> vertexSorted = new List<BrepVertex>();
            List<BrepFace> brepFaces = SortBrepFaces(brep, bottomFace);
            List<BrepVertex> brepVertex = brep.Vertices.ToList();

            foreach (BrepVertex vertex in brepVertex)
            {
                bool isOnBottomFace = IsOnFace(vertex.Location, brepFaces[0]);
                if (isOnBottomFace) { vertexSorted.Add(vertex); }
            }

            foreach (BrepVertex vertex in brepVertex)
            {
                bool isOnTopFace = IsOnFace(vertex.Location, brepFaces[1]);
                if (isOnTopFace) { vertexSorted.Add(vertex); }
            }
            return vertexSorted;
        }
        private bool IsOnFace(Point3d point, BrepFace face)
        {
            bool isOnFace = false;

            face.ClosestPoint(point, out double PointOnCurveU, out double PointOnCurveV);
            Point3d testPoint = face.PointAt(PointOnCurveU, PointOnCurveV);  // make test point 
            double distanceToFace = (testPoint - point).Length; // calculate distance between testPoint and node
            if (distanceToFace <= 0.0001 & distanceToFace >= -0.0001) // if distance = 0: node is on edge
            {
                isOnFace = true;
            }
            return isOnFace;
        }
        private bool IsOnEdge(Point3d point, BrepEdge edge)
        {
            bool isOnEdge = false;

            edge.ClosestPoint(point, out double PointOnCurve);
            Point3d testPoint = edge.PointAt(PointOnCurve);  // make test point 
            double distanceToEdge = (testPoint - point).Length; // calculate distance between testPoint and node
            if (distanceToEdge <= 0.0001 & distanceToEdge >= -0.0001) // if distance = 0: node is on edge
            {
                isOnEdge = true;
            }
            return isOnEdge;
        }

        /// <summary>
        /// Find the edges of brep composing rails
        /// </summary>
        /// <returns> List with rails. </returns>
        private List<Curve> FindRails(Brep brep, int bottomFace)
        {
            // Sort brep edges.
            List<BrepEdge> brepEdges = SortBrepEdges(brep, bottomFace);

            // Find rails.
            List<Curve> rails = new List<Curve>() { brepEdges[8], brepEdges[9], brepEdges[10], brepEdges[11] };
            #region Old Code
            /*
            List<BrepFace> brepFace = brep.Faces.ToList();
            List<int> indexAdjecentFaces = (brepFace[bottomFace].AdjacentFaces()).ToList();
            List<int> indexAdjecentEdges = (brepFace[bottomFace].AdjacentEdges()).ToList();
            indexAdjecentFaces.Add(bottomFace);
            for (int i = 0; i < brepFace.Count; i++)
            {
                if (!indexAdjecentFaces.Contains(brepFace.IndexOf(brepFace[i])))
                {
                    BrepFace brepBottomFace = brepFace[bottomFace];
                    BrepFace brepTopFace = brepFace[i]; // top face
                    indexAdjecentEdges.AddRange(brepBottomFace.AdjacentEdges());
                    indexAdjecentEdges.AddRange(brepTopFace.AdjacentEdges());
                    continue;
                }
            }*/

            // Find rails
            /*
            List<int> indexAdjecentEdges = new List<int>();
            indexAdjecentEdges.AddRange(brepFace[0].AdjacentEdges().ToList());
            indexAdjecentEdges.AddRange(brepFace[1].AdjacentEdges().ToList());
            List<BrepEdge> brepEdges = brep.Edges.ToList();
            List<Curve> rails = new List<Curve>(brepEdges);
            foreach (int index in indexAdjecentEdges) { rails.Remove(brepEdges[index]); }*/
            #endregion
            #region Old Code
            /*
            foreach (BrepEdge edge in brepEdges) // check if node is on edge
            {
                List<Point3d> edgePoints = new List<Point3d> { edge.StartVertex.Location, edge.EndVertex.Location };
                bool pointOnFace4 = false;
                bool pointOnFace5 = false;
                foreach (Point3d point in edgePoints)
                {
                    brepFace[0].ClosestPoint(point, out double PointOnCurveUFace4, out double PointOnCurveVFace4);
                    brepFace[1].ClosestPoint(point, out double PointOnCurveUFace5, out double PointOnCurveVFace5);
                    Point3d testPointFace4 = brepFace[0].PointAt(PointOnCurveUFace4, PointOnCurveVFace4);  // make test point
                    Point3d testPointFace5 = brepFace[1].PointAt(PointOnCurveUFace5, PointOnCurveVFace5);  // make test point
                    double distanceToFace4 = testPointFace4.DistanceTo(point); // calculate distance between testPoint and node
                    double distanceToFace5 = testPointFace5.DistanceTo(point); // calculate distance between testPoint and node
                    if ((distanceToFace4 <= 0.0001 & distanceToFace4 >= -0.0001)) // if distance = 0: node is on edge
                    {
                        pointOnFace4 = true;
                    }
                    else if ((distanceToFace5 <= 0.0001 & distanceToFace5 >= -0.0001))
                    {
                        pointOnFace5 = true;
                    }
                }
                if (pointOnFace4 & pointOnFace5)
                {
                    rails.Add(edge);  //get edge1 of brep = rail 1
                }
            }*/
            #endregion
            return rails;
        }

        /// <summary>
        /// Check if mesh is compatible with Abaqus
        /// </summary>
        /// <returns> Nothing. Assign propertie to solidMesh. </returns>
        private void IsBrepCompatibleWithAbaqus(Element element, Mesh3D solidMesh)
        {
            List<Point3d> nodes = new List<Point3d> { element.Node1.Coordinate, element.Node2.Coordinate, element.Node3.Coordinate, element.Node4.Coordinate };
            NurbsCurve curve = PolyCurve.CreateControlPointCurve(nodes, 2).ToNurbsCurve();
            Vector3d direction = element.Node5.Coordinate - element.Node1.Coordinate;
            curve.MakeClosed(0.001);
            string curveOrientation = curve.ClosedCurveOrientation(direction).ToString();
            if (curveOrientation == "Clockwise")
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Abaqus can not interpret order of nodes. Error in brep input. ");
                solidMesh.inp = false;
            }
            else { solidMesh.inp = true; }
        }

        /// <summary>
        /// Divide each brep edge w-direction into nw points. The brep edges in w-direction are named rail.
        /// </summary>
        /// <returns> DataTree with points on each rail. Branch: floor level.</returns>
        private DataTree<Point3d> DivideRailIntoNwPoints(List<Curve> rails, Brep brep, int nw, int bottomFace)
        {
            DataTree<Point3d> railPoints = new DataTree<Point3d>(); // todo: erstatt brep og bottomFace med geometry
            BrepFace brepBottomFace = brep.Faces[bottomFace]; //todo: fix input.

            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(nw, true, out Point3d[] pt);
                List<Point3d> point = pt.ToList();
                brepBottomFace.ClosestPoint(point[0], out double PointOnCurveU, out double PointOnCurveV);
                Point3d testPoint = brepBottomFace.PointAt(PointOnCurveU, PointOnCurveV);
                Vector3d distanceToFace = testPoint - point[0];

                if (distanceToFace.Length > 0.001) { point.Reverse(); }
                //if (!IsOnFace(point[0], brepBottomFace)) { point.Reverse(); } //todo: test if this works instead of over 
                for (int j = 0; j < point.Count; j++)
                {
                    railPoints.Add(point[j], new GH_Path(j)); //tree with nw points on each rail. Branch: floor
                }
            }

            // Check if the rails must be re-oredered to generate elements with nodes counting ccw
            Curve testCurve = Curve.CreateControlPointCurve(railPoints.Branch(0), 1);
            Vector3d direction = railPoints.Branch(nw)[0] - railPoints.Branch(0)[0];
            string curveOrientation = testCurve.ClosedCurveOrientation(direction).ToString();
            if (curveOrientation == "Clockwise")
            {
                for (int i = 0; i < railPoints.BranchCount; i++)
                {
                    railPoints.Branch(i).Reverse();
                }
            }
            #region Old code
            /*
                Point3d[] nwPt;
                List<Point3d> nwPoints = new List<Point3d>();
                //DataTree<Point3d> railPoints = new DataTree<Point3d>();

                // Check if the rails must be re-oredered to generate elements with nodes counting ccw
                NurbsSurface endSurface = brep.Surfaces[5].ToNurbsSurface();
                var u = endSurface.Domain(0);
                var v = endSurface.Domain(1);

                Vector3d vector1 = (rails[0].PointAtEnd - rails[0].PointAtStart);
                Vector3d vector2 = endSurface.NormalAt(u.T1 * 0.5, v.T1 * 0.5); //(0, 0);
                Vector3d normal = Vector3d.CrossProduct(vector1, vector2);
                double angle = Vector3d.VectorAngle(vector1, vector2, normal);
                if (angle > Math.PI / 2)
                {
                    rails.Reverse();
                }

                //Divide each rail into nw points.
                for (int i = 0; i < rails.Count; i++)
                {
                    rails[i].DivideByCount(nw, true, out nwPt);  //divide each rail in nw number of points
                    nwPoints = nwPt.ToList();
                    for (int j = 0; j < nwPoints.Count; j++)
                    {
                        railPoints.Add(nwPoints[j], new GH_Path(j)); //tree with nw points on each rail. Branch: floor
                    }
                }*/
            #endregion
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
                //Vector3d normal = Vector3d.CrossProduct(vec1, vec2);
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
        private DataTree<Point3d> CreateGridOfPointsAtEachFloor(int nu, int nv, List<NurbsSurface> surfaceAtNw, DataTree<Point3d> railPoints)
        {
            List<Point3d> pt = new List<Point3d>();
            DataTree<Point3d> points = new DataTree<Point3d>();
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
                pt = CreateGridOfPointsUV(nu, nv, surfaceAtNw[i], railPoints.Branch(i)[0], direction);
                points.AddRange(pt, new GH_Path(i)); // add points to datatree. Branch: floor level
                pt.Clear();
            }
            return points;
        }


        /// <summary>
        /// Makes grid of points in U and V direction
        /// </summary>
        /// <returns> List of points in U and V direction</returns>
        private List<Point3d> CreateGridOfPointsUV(int nu, int nv, NurbsSurface surface, Point3d railPoint, Vector3d direction)
        {
            var u = surface.Domain(0);
            var v = surface.Domain(1);
            surface.UVNDirectionsAt(0, 0, out Vector3d uDir, out Vector3d vDir, out Vector3d nDir);
            if (Point3d.Subtract(surface.PointAt(u.T1, v.T0), railPoint).Length < 0.1)
            {
             /*   Vector3d vec1 = surface.PointAt(u.T1, v.T0) - surface.PointAt(u.T0, v.T0);
                Vector3d vec2 = surface.PointAt(u.T0, v.T1) - surface.PointAt(u.T0, v.T0);*/

                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                {
                    u = new Interval(u.T1, u.T0);
                    List<Point3d> points = GeneratePoints(0, surface, u, v, nu, nv);
                    return points; 
                }
                else
                {
                    u = new Interval(u.T1, u.T0);
                    List<Point3d> points = GeneratePoints(1, surface, u, v, nu, nv);
                    return points;
                }
            }
            else if (Point3d.Subtract(surface.PointAt(u.T0, v.T1), railPoint).Length < 0.1)
            {
                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                {
                    v = new Interval(v.T1, v.T0);
                    List<Point3d> points = GeneratePoints(0, surface, u, v, nu, nv);
                    return points;
                }
                else
                {
                    v = new Interval(v.T1, v.T0);
                    List<Point3d> points = GeneratePoints(1, surface, u, v, nu, nv);
                    return points;
                }

            }
            else if (Point3d.Subtract(surface.PointAt(u.T1, v.T1), railPoint).Length < 0.1)
            {
                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                { 
                    u = new Interval(u.T1, u.T0);
                    v = new Interval(v.T1, v.T0);
                    List<Point3d> points = GeneratePoints(1, surface, u, v, nu, nv);
                    return points;
                }
                else 
                {
                    u = new Interval(u.T1, u.T0);
                    v = new Interval(v.T1, v.T0);
                    List<Point3d> points = GeneratePoints(0, surface, u, v, nu, nv);
                    return points;
                }
            }
            else 
            {
                if (Vector3d.VectorAngle(Vector3d.CrossProduct(uDir, vDir), direction) > Math.PI / 2)
                {
                    List<Point3d> points = GeneratePoints(1, surface, u, v, nu, nv);
                    return points;
                }
                else
                {
                    List<Point3d> points = GeneratePoints(0, surface, u, v, nu, nv);
                    return points;
                }
            }
            #region Old Code
            /*
            var u = surface.Domain(0);
            var v = surface.Domain(1);

            double stepU = 1 / ((double)nu) * u.Length;
            double stepV = 1 / ((double)nv) * v.Length;

            string curveOrientation = (intersectionCurve[0].ToNurbsCurve()).ClosedCurveOrientation(plane).ToString();
            if (curveOrientation == "CounterClockwise")
            {

                double pointU = 0;
                double pointV = 0;
                for (double j = 0; j <= nv; j++)
                {
                    for (double k = 0; k <= nu; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV));  // make point on surface
                        pointU = pointU + stepU;
                    }
                    pointV = pointV + stepV;
                    pointU = 0;
                }
            }
            else
            {
                double pointU = 0; //u.Length;
                double pointV = v.Length;
                for (double j = 0; j <= nv; j++)
                {
                    for (double k = 0; k <= nu; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV)); // make point on surface
                        pointU = pointU + stepU;
                    }
                    pointV = pointV - stepV;
                    pointU = 0;//u.Length;
                }
            }*/
            #endregion
        }
        private List<Point3d> GeneratePoints(int Case, NurbsSurface surface, Interval u, Interval v, int nu, int nv)
        {
            List<Point3d> pt = new List<Point3d>();

            if (Case == 0)
            {
                double stepU = u.Length / (double)nu;
                double stepV = v.Length / (double)nv;

                double pointU = u.T0;
                double pointV = v.T0;
                for (double j = 0; j <= nv; j++)
                {
                    for (double k = 0; k <= nu; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV));  // make point on surface
                        pointU = pointU + stepU;
                    }
                    pointV = pointV + stepV;
                    pointU = u.T0;
                }
            }
            else if (Case == 1)
            {
                double stepU = u.Length / (double)nv;
                double stepV = v.Length / (double)nu;

                double pointU = u.T0;
                double pointV = v.T0;
                for (double j = 0; j <= nv; j++)
                {
                    for (double k = 0; k <= nu; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV));  // make point on surface
                        pointV = pointV + stepV;
                    }
                    pointU = pointU + stepU;
                    pointV = v.T0;
                }
            }
            else { pt = null; }
            return pt;
        }

        /// <summary>
        /// Create Nodes: assign Coordiantes, GlobalId and Boundary Conditions
        /// </summary>
        /// <returns> List with nodes incl properties</returns>
        private List<Node> CreateNodes(DataTree<Point3d> meshPoints, int nu, int nv, int nw)
        {
            List<Node> nodes = new List<Node>();
            int id = 0;
            nu = nu + 1; //input nu = nu - 1. Exs: nu = 3, total points in u-direction is 4;
            nv = nv + 1; //input nv = nv - 1. Exs: nv = 3, total points in v-direction is 4;
            nw = nw + 1; //input nw = nw - 1. Exs: nw = 3, total points in w-direction is 4;

            for (int i = 0; i < nw; i++)
            {
                int row = 0;
                int column = 0;
                for (int j = 0; j < meshPoints.Branch(i).Count; j++)
                {
                    Node node = new Node(id, meshPoints.Branch(i)[j]); // assign Global ID and cooridinates
                    id++;

                    if (column == 0 | column == nu - 1) { node.BC_U = true; } // assign BCU
                    if (row == 0 | row == nv - 1) { node.BC_V = true; } // assign BCV
                    if (i == 0 | i == nw - 1) { node.BC_W = true; } // assign BCW

                    column++;
                    if (column == nu)
                    {
                        row++;
                        column = 0;
                    }
                    nodes.Add(node);
                }
            }
            return nodes;
        }


        /*
        /// <summary>
        /// Create Elements: assign ElementId, ElementMesh and Nodes incl. Coordiantes, GlobalId, LocalId and Boundary Conditions), elementId, elementMesh.
        /// </summary>
        /// <returns>List with elements incl properties</returns>
        private List<Element> CreateHexElements(List<Node> nodes, int nu, int nv, int nw)
        {
            List<Element> elements = new List<Element>();

            nu = nu + 1; //input nu = nu - 1. Exs: nu = 3, total points in u-direction is 4;
            nv = nv + 1; //input nv = nv - 1. Exs: nv = 3, total points in v-direction is 4;
            nw = nw + 1; //input nv = nv - 1. Exs: nv = 3, total points in v-direction is 4; // Hilde: riktig ?
            int counter = 0;
            int jumpInPlane = 0;
            int numNodeInLevel = nu*nv;
            for (int i = 0; i < (nu - 1) * (nv - 1) * (nw - 1); i++)  // loop levels
            {
                if (jumpInPlane < (nu - 1))
                {
                    List<int> connectivity = new List<int>();

                    connectivity.Add(counter);
                    connectivity.Add(counter + 1);
                    connectivity.Add(counter + 1 + nu);
                    connectivity.Add(counter + nu);
                    connectivity.Add(counter + numNodeInLevel);
                    connectivity.Add(counter + 1 + numNodeInLevel);
                    connectivity.Add(counter + 1 + nu + numNodeInLevel);
                    connectivity.Add(counter + nu + numNodeInLevel);

                    counter++;
                    jumpInPlane++;


                    List<Node> elementNodes = new List<Node>();
                    foreach (int id in connectivity)
                    {
                        elementNodes.Add(nodes[id]);
                    }

                    Element element = new Element(i, elementNodes, connectivity);

                    // create local mesh
                    Mesh localMesh = new Mesh();
                    foreach (Node node in elementNodes)
                    {
                        localMesh.Vertices.Add(node.Coordinate); //0
                    }

                    localMesh.Faces.AddFace(0, 1, 5, 4);
                    localMesh.Faces.AddFace(1, 2, 6, 5);
                    localMesh.Faces.AddFace(2, 3, 7, 6);
                    localMesh.Faces.AddFace(3, 0, 4, 7);
                    localMesh.Faces.AddFace(0, 1, 2, 3);
                    localMesh.Faces.AddFace(4, 5, 6, 7);
                    /*
                    mesh.Faces.AddFace(connectivity[0], connectivity[1], connectivity[2], connectivity[3]);
                    mesh.Faces.AddFace(connectivity[4], connectivity[5], connectivity[6], connectivity[7]);
                    mesh.Faces.AddFace(connectivity[0], connectivity[1], connectivity[5], connectivity[4]);
                    mesh.Faces.AddFace(connectivity[1], connectivity[2], connectivity[6], connectivity[5]);
                    mesh.Faces.AddFace(connectivity[2], connectivity[3], connectivity[7], connectivity[6]);
                    mesh.Faces.AddFace(connectivity[3], connectivity[0], connectivity[4], connectivity[7]);
                    
                    localMesh.Normals.ComputeNormals();  //Control if needed
                    localMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                    localMesh.Compact(); //to ensure that it calculate

                    //add element and mesh to element list
                    element.mesh = localMesh;
                    elements.Add(element);
                }
                else
                {
                    jumpInPlane = 0;
                    counter++;
                }
            }                
            
            return elements;
        }
        */

        private List<Element> CreateHexElements(DataTree<Point3d> meshPoints, List<Node> nodes, int nu, int nv)
        {
            Element e = new Element();
            Mesh mesh = new Mesh();
            List<Element> elements = new List<Element>();
            List<Point3d> ptsBot = new List<Point3d>();
            List<Point3d> ptsTop = new List<Point3d>();
            int elemId = 0;

            nu = nu + 1; //input nu = nu - 1. Exs: nu = 3, total points in u-direction is 4;
            nv = nv + 1; //input nv = nv - 1. Exs: nv = 3, total points in v-direction is 4;

            for (int i = 0; i < meshPoints.BranchCount - 1; i++)  // loop levels
            {
                ptsBot = meshPoints.Branch(i);
                ptsTop = meshPoints.Branch(i + 1);
                int count2 = 0;
                int counter = meshPoints.Branch(0).Count * i;

                for (int j = 0; j < meshPoints.Branch(0).Count - nu - 1; j++) // loop elements in a level
                {
                    e.Id = elemId;
                    e.IsCube = true;
                    List<Node> elementNodes = new List<Node>();
                    List<int> connectivity = new List<int>();

                    if (count2 < nu - 1)
                    {

                        connectivity.Add(counter);
                        connectivity.Add(counter + 1);
                        connectivity.Add(counter + nu + 1);
                        connectivity.Add(counter + nu);

                        connectivity.Add(counter + nu * nv);
                        connectivity.Add(counter + 1 + nu * nv);
                        connectivity.Add(counter + nu + 1 + nu * nv);
                        connectivity.Add(counter + nu + nu * nv);

                        foreach (int id in connectivity)
                        {
                            elementNodes.Add(nodes[id]);
                        }

                        Element element = new Element(i, elementNodes, connectivity);

                        // create local mesh
                        Mesh localMesh = new Mesh();
                        foreach (Node node in elementNodes)
                        {
                            localMesh.Vertices.Add(node.Coordinate); //0
                        }

                        localMesh.Faces.AddFace(0, 1, 5, 4);
                        localMesh.Faces.AddFace(1, 2, 6, 5);
                        localMesh.Faces.AddFace(2, 3, 7, 6);
                        localMesh.Faces.AddFace(3, 0, 4, 7);
                        localMesh.Faces.AddFace(0, 1, 2, 3);
                        localMesh.Faces.AddFace(4, 5, 6, 7);

                        localMesh.Normals.ComputeNormals();  //Control if needed
                        localMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                        localMesh.Compact(); //to ensure that it calculate

                        //add element and mesh to element list
                        element.mesh = localMesh;
                        elements.Add(element);

                        count2++;
                        elemId++;
                        counter++;
                    }
                    else { count2 = 0; counter++; }
                }

            }
            return elements;
        }

        /// <summary>
        /// Create Global mesh
        /// </summary>
        /// <returns>Global mesh</returns>
        /*
        private Mesh CreateGlobalMesh(List<Node> nodes, int nu, int nv, int nw)
        {
            nu++;
            nv++;
            nw++; // to do: sjekk disse


            Mesh mesh = new Mesh();

            foreach (Node node in nodes)
            {
                mesh.Vertices.Add(node.Coordinate);
            }

            // mesh planes in w -dir
            int counter = 0;
            int jump = 0; // new v sequence
            for (int i = 0; i < nw; i++)
            {
                for (int j = 0; j < (nu-1)*(nv-1); j++)
                {
                    mesh.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);

                    counter++;
                    jump++;
                    if (jump == (nu - 1)) // check if done with a v sequence
                    {
                        counter++;
                        jump = 0; // new v sequence
                    }
                }
                counter++;
            }


            // mesh planes in u - dir
            counter = 0; 
            jump = 0; // new w sequence
            for (int i = 0; i < nu; i++)
            {
                for (int j = 0; j < (nv - 1) * (nw - 1); j++)
                {
                    mesh.Faces.AddFace(counter, counter + 1, counter + nv + 1, counter + nv);

                    counter++;
                    jump++;
                    if (jump == (nv - 1)) // check if done with a w sequence
                    {
                        counter++;
                        jump = 0; // new w sequence
                    }
                }
                counter++;
            }


            // mesh planes in v - dir
            counter++;
            jump = 0; // new u sequence
            for (int i = 0; i < nv; i++)
            {
                for (int j = 0; j < (nw - 1) * (nu - 1); j++)
                {
                    mesh.Faces.AddFace(counter, counter + 1, counter + nw + 1, counter + nw);

                    counter++;
                    jump++;
                    if (jump == (nw - 1)) // check if done with a u sequence
                    {
                        counter++;
                        jump = 0; // new u sequence
                    }
                }
                counter++;
            }

            mesh.Normals.ComputeNormals();  //Control if needed
            mesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            mesh.Compact(); //to ensure that it calculate

            return mesh;
        }*/

        private Mesh CreateGlobalMesh(List<Element> elements)
        {
            Mesh allMesh = new Mesh();
            foreach (Element el in elements)
            {
                allMesh.Append(el.mesh);
            }
            allMesh.Weld(0.01);

            return allMesh;
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
                return Properties.Resources.Icon_SolidMesh;
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
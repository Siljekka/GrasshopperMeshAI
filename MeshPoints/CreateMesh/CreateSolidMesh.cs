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
    public class CreateSolidMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateSolidMesh_Sweep class.
        /// </summary>
        public CreateSolidMesh()
          : base("CreateSolidMesh", "solid",
              "Generate solid mesh. Independent on how surface composing the brep is made.",
              "SmartMesh", "Mesh")
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
            pManager.AddGenericParameter("SmartMesh", "solid", "Creates a SolidMesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "m", "Mesh (solid elements).", GH_ParamAccess.item);
            //pManager.AddGenericParameter("test2", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // to do: Hilde, slett?
            /* Todo:
             * Legg til/endre INP-valid.
            */
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



            /* Todo:
            * Fjern bottomFace og brep og ha heller brepGeometry som input
            * kan erstatte FindRails
'            */

            // 1. Check input OK
            if (!DA.GetData(0, ref brep)) return;
            if (!DA.GetData(1, ref bottomFace)) return;
            if (u == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "u can not be zero."); return; }
            if (v == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "v can not be zero."); return; }
            if (w == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "w can not be zero."); return; }

            SmartMesh smartMesh = new SmartMesh();

            // 2. Assign geometrical properties to mesh
            Geometry brepGeometry = new Geometry(brep, bottomFace);
            smartMesh.nu = u + 1;
            smartMesh.nv = v + 1;
            smartMesh.nw = w + 1;
            smartMesh.Type = "Solid";
            smartMesh.Geometry = brepGeometry;
            
            // 3. Find Rails
            List<Curve> rails = FindRails(brepGeometry.Edges);

            // 4. Divide each brep edge in w direction (rail) into nw points.
            DataTree<Point3d> railPoints = DivideRailIntoWPoints(rails, brepGeometry.Faces[0], w);

            // 5. Create NurbsSurface for each nw-floor
            DataTree<Curve> intersectionCurve = GetIntersectionCurveBrepAndRailPoints(railPoints, brep);
            if (intersectionCurve == null) return;
            List<NurbsSurface> surfaceAtW = CreateNurbSurfaceAtEachFloor(intersectionCurve);

            // 6. Make grid of points in u and v direction at leven nw
            DataTree<Point3d> meshPoints = CreateGridOfPointsAtEachFloor(u, v, surfaceAtW, railPoints);
            
            // 7. Create nodes 
            smartMesh.Nodes = CreateNodes(meshPoints, u, v, w); // assign Coordiantes, GlobalId and Boundary Conditions

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
        /// Find the edges of brep composing rails
        /// </summary>
        /// <returns> List with rails. </returns>
        private List<Curve> FindRails(List<BrepEdge> brepEdges)
        {
            // Find rails.
            List<Curve> rails = new List<Curve>() { brepEdges[8], brepEdges[9], brepEdges[10], brepEdges[11] };
            // to do: Hilde, slett old code.
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
            #endregion // to do:
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
                //if (!IsOnFace(point[0], brepBottomFace)) { point.Reverse(); } //todo: Hilde: test if this works instead of over 
                for (int j = 0; j < point.Count; j++)
                {
                    railPoints.Add(point[j], new GH_Path(j)); //tree with nw points on each rail. Branch: floor
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
            // to do: Hilde, slett
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
        private DataTree<Point3d> CreateGridOfPointsAtEachFloor(int u, int v, List<NurbsSurface> surfaceAtNw, DataTree<Point3d> railPoints)
        {
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
                List<Point3d> pt = CreateGridOfPointsUV(u, v, surfaceAtNw[i], railPoints.Branch(i)[0], direction);
                points.AddRange(pt, new GH_Path(i)); // add points to datatree. Branch: floor level
                pt.Clear();
            }
            return points;
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
             /*   Vector3d vec1 = surface.PointAt(u.T1, v.T0) - surface.PointAt(u.T0, v.T0);
                Vector3d vec2 = surface.PointAt(u.T0, v.T1) - surface.PointAt(u.T0, v.T0);*/

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

        /// <summary>
        /// Create Nodes: assign Coordiantes, GlobalId and Boundary Conditions
        /// </summary>
        /// <returns> List with nodes incl properties</returns>
        private List<Node> CreateNodes(DataTree<Point3d> meshPoints, int u, int v, int w)
        {
            List<Node> nodes = new List<Node>();
            int id = 0;
            int nu = u + 1; // number nodes in u dir
            int nv = v + 1; // number nodes in v dir 
            int nw = w + 1; // number nodes in w dir 

            for (int i = 0; i < nw; i++)
            {
                int row = 0;
                int column = 0;
                for (int j = 0; j < meshPoints.Branch(i).Count; j++)
                {
                    bool BC_U = false;
                    bool BC_V = false;
                    bool BC_W = false;

                    if (column == 0 | column == nu - 1) { BC_U = true; } // assign BCU
                    if (row == 0 | row == nv - 1) { BC_V = true; } // assign BCV
                    if (i == 0 | i == nw - 1) { BC_W = true; } // assign BCW

                    Node node = new Node(id, meshPoints.Branch(i)[j], BC_U, BC_V, BC_W);
                    id++;

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

                        Element element = new Element(elemId, elementNodes, connectivity);

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
                        element.Mesh = localMesh;
                        elements.Add(element);

                        count2++;
                        elemId++;
                        counter++;
                    }
                    else { count2 = 0; counter++; }
                }

            }
            return elements;
        }// to do: slett

        /// <summary>
        /// Create Global mesh
        /// </summary>
        /// <returns>Global mesh</returns>

        private Mesh CreateGlobalMesh(List<Element> elements)
        {
            Mesh allMesh = new Mesh();
            foreach (Element el in elements)
            {
                allMesh.Append(el.Mesh);
            }
            allMesh.Weld(0.01);

            return allMesh;
        } // to do: slett

        private Mesh CreateGlobalMeshNew(List<Node> nodes, int nu, int nv, int nw)
        {
            // to do: fix senere. mesh w- dir er riktig
            nu++;
            nv++;
            nw++; 

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
               for (int j = 0; j < (nu - 1) * (nv - 1 ); j++)
               {
                   if (jump < nu - 1)
                   {

                       mesh.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);

                       counter++;
                       jump++;
                   }
                   else { counter++; jump = 0; j--; }
               }
               counter = (i + 1) * nu * nv ;
               jump = 0;
           }

           // mesh planes in u - dir
           counter = 0;
           jump = 0; // new w sequence
           for (int i = 0; i < nu; i++)
           {
               for (int j = 0; j < (nv - 1) * (nw - 1); j++)
               {
                   if (jump < nv - 1)
                   {

                       mesh.Faces.AddFace(counter, counter + (nu*nv), counter + (nu * nv) + nu, counter + nu );

                       counter += nu;
                       jump++;
                   }
                   else { counter += nu; jump = 0; j--; }
               }
               counter = (i + 1);
               jump = 0;
           }
       

            // mesh planes in v - dir
            counter=0;
            int counterU = 0;
            jump = 0; // new u sequence
            for (int i = 0; i < nv; i++)
            {
                for (int j = 0; j < (nu - 1) * (nw - 1); j++)
                {
                    if (jump < nw - 1)
                    {

                        mesh.Faces.AddFace(counter, counter + 1, counter + (nu * nv) + 1, counter + nu*nv);

                        counter += nu*nv;
                        jump++;
                    }
                    else { counterU++; counter = counterU; jump = 0; j--; }
                }
                counter = (i + 1)*nu;
                jump = 0;
                counterU = counter;
            }

            
            
            mesh.Normals.ComputeNormals();  //Control if needed
            mesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            mesh.Compact(); //to ensure that it calculate

            return mesh;
        }// to do: slett

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
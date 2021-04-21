using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry.Intersect;

namespace MeshPoints.CreateMesh
{
    public class CreateSolidMesh_OldSweep : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public CreateSolidMesh_OldSweep()
          : base("Create SolidMesh (OldSweep)", "solid",
              "Creates a solid mesh (brep can not be made from box",
              "MyPlugIn", "Outdated")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "bp", "Brep", GH_ParamAccess.item);
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
            pManager.AddGenericParameter("Vertices", "vertices", "Mesh (solid elements).", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            /*
            // Input
            Brep brep = new Brep();
            int nu = 0;
            int nv = 0;
            int nw = 0;
            DA.GetData(0, ref brep);
            DA.GetData(1, ref nu);
            DA.GetData(2, ref nv);
            DA.GetData(3, ref nw);

            #region Variables
            //Variables
            SmartMesh solidMesh = new SmartMesh();
            Mesh allMesh = new Mesh();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<NurbsSurface> surfaceAtNw = new List<NurbsSurface>();
            DataTree<Curve> intersectionCurve = new DataTree<Curve>();
            DataTree<Point3d> railPoints = new DataTree<Point3d>();
            DataTree<Point3d> meshPoints = new DataTree<Point3d>();
            List<Plane> planes = new List<Plane>();
            #endregion

            if (!brep.IsValid) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No valid brep input found."); return; } //todo: is this one needed?
            if (nu == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nu can not be zero."); return; }
            if (nv == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nv can not be zero."); return; }
            if (nw == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nw can not be zero."); return; }


            // 1. Assign properties to SolidMesh
            solidMesh.nu = nu;
            solidMesh.nv = nv;
            solidMesh.nw = nw;

            //2. Divide each brep edge in w direction (rail) into nw points.
            railPoints = DivideRailIntoNwPoints(brep, solidMesh.nw);

            //3. Create NurbsSurface for each nw-floor
            intersectionCurve = GetIntersectionCurveBrepAndRailPoints(railPoints, brep).Item1;
            planes = GetIntersectionCurveBrepAndRailPoints(railPoints, brep).Item2;
            surfaceAtNw = CreateNurbSurfaceAtEachFloor(intersectionCurve);

            //4. Make grid of points in u and v direction at leven nw
            meshPoints = CreateGridOfPointsAtEachFloor(solidMesh.nu, solidMesh.nv, surfaceAtNw, intersectionCurve, planes);

            //5. Create nodes and elements
            nodes = CreateNodes(meshPoints, solidMesh.nu, solidMesh.nv, solidMesh.nw); // assign Coordiantes, GlobalId and Boundary Conditions
            elements = CreateHexElements(meshPoints, nodes, solidMesh.nu, solidMesh.nv); // assign ElementId, ElementMesh and Nodes incl. Coordiantes, GlobalId, LocalId and Boundary Conditions), elementId, elementMesh.

            // Check if brep can be interpret by Abaqus
            IsBrepCompatibleWithAbaqus(elements[0], solidMesh);

            //6. Create global mesh
            allMesh = CreateGlobalMesh(elements);

            //7. Add properties to SolidMesh
            solidMesh.Nodes = nodes;
            solidMesh.Elements = elements;
            solidMesh.mesh = allMesh;
           
            // Output
            DA.SetData(0, solidMesh);
            DA.SetData(1, solidMesh.mesh);
            */
        }
        /*
        #region Methods

        /// <summary>
        /// Check if mesh is compatible with Abaqus
        /// </summary>
        /// <returns> Nothing. Assign propertie to solidMesh. </returns>
        private void IsBrepCompatibleWithAbaqus(Element element, SmartMesh solidMesh)
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
        private DataTree<Point3d> DivideRailIntoNwPoints(Brep brep, int nw)
        {
            Point3d[] nwPt;
            List<Point3d> nwPoints = new List<Point3d>();
            DataTree<Point3d> railPoints = new DataTree<Point3d>();

            // Find edges composing the rails and add into list
            Curve rail1 = brep.Edges[0];  //get edge1 of brep = rail 1
            Curve rail2 = brep.Edges[11];  //get edge2 of brep = rail 2
            Curve rail3 = brep.Edges[10]; //get edge3 of brep = rail 3
            Curve rail4 = brep.Edges[9]; //get edge4 of brep = rail 4

            List<Curve> rails = new List<Curve>() { rail1, rail2, rail3, rail4 };
            rails.Reverse();

            //Divide each rail into nw points.
            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(nw, true, out nwPt);  //divide each rail in nw number of points
                nwPoints = nwPt.ToList(); //.OrderBy(f => f.Z).ToList();     //order the points according to z-value
                for (int j = 0; j < nwPoints.Count; j++)
                {
                    railPoints.Add(nwPoints[j], new GH_Path(j)); //tree with nw points on each rail. Branch: floor
                }
            }
            return railPoints;
        }


        /// <summary>
        /// Get intersectionCurve between Brep and plane made from RailPoints. Curve is closed.
        /// </summary>
        /// <returns> DataTree with closed curves at each Floor. Branch: Floor Level </returns>
        private Tuple<DataTree<Curve>, List<Plane>> GetIntersectionCurveBrepAndRailPoints(DataTree<Point3d> railPoints, Brep brep)
        {
            DataTree<Curve> intersectionCurve = new DataTree<Curve>();
            List<Plane> planes = new List<Plane>();

            for (int i = 0; i < railPoints.BranchCount; i++)
            {
                Vector3d vec1 = railPoints.Branch(i)[1] - railPoints.Branch(i)[0];
                Vector3d vec2 = railPoints.Branch(i)[3] - railPoints.Branch(i)[0];
                Vector3d normal = Vector3d.CrossProduct(vec1, vec2);
                Plane plane = new Plane(railPoints.Branch(i)[0], normal);
                //Plane.FitPlaneToPoints(railPoints.Branch(i), out Plane plane); // make plane on floor i
                Intersection.BrepPlane(brep, plane, 0.00001, out Curve[] iCrv, out Point3d[] iPt); // make intersection curve between brep and plane on floor i
                List<Curve> intCrv = iCrv.ToList();
                planes.Add(plane);

                for (int j = 0; j < intCrv.Count; j++) { intCrv[j].MakeClosed(0.0001); intersectionCurve.Add(intCrv[j], new GH_Path(i)); }  // make curve closed and add to intersectionCurve
                if (intersectionCurve.Branch(i).Count != 1) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Brep input is not OK."); }

                intCrv.Clear();
            }
            return new Tuple<DataTree<Curve>, List<Plane>>(intersectionCurve, planes);
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
        private DataTree<Point3d> CreateGridOfPointsAtEachFloor(int nu, int nv, List<NurbsSurface> surfaceAtNw, DataTree<Curve> intersectionCurve, List<Plane> planes)
        {
            List<Point3d> pt = new List<Point3d>();
            DataTree<Point3d> points = new DataTree<Point3d>();

            for (int i = 0; i < surfaceAtNw.Count; i++) // loop floors
            {
                pt = CreateGridOfPointsUV(nu, nv, surfaceAtNw[i], intersectionCurve.Branch(i), planes[i]);
                points.AddRange(pt, new GH_Path(i)); // add points to datatree. Branch: floor level
                pt.Clear();
            }
            return points;
        }

        /// <summary>
        /// Makes grid of points in U and V direction
        /// </summary>
        /// <returns> List of points in U and V direction</returns>
        private List<Point3d> CreateGridOfPointsUV(int nu, int nv, NurbsSurface surfaceAtNw, List<Curve> intersectionCurve, Plane plane)
        {
            List<Point3d> pt = new List<Point3d>();
            NurbsSurface surface = surfaceAtNw;

            var u = surface.Domain(0);
            var v = surface.Domain(1);

            double stepU = 1 / ((double)nu) * u.Length;
            double stepV = 1 / ((double)nv) * v.Length;


            /*
            surface.SetDomain(0, new Interval(0, 1)); // set domain for surface 0-direction
            surface.SetDomain(1, new Interval(0, 1)); // set domain for surface 1-direction
            double stepU = 1 / ((double)nu - 1);
            double stepV = 1 / ((double)nv - 1);*/
        /*
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
                double pointU = 0;//u.Length;
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
            }
            return pt;
        }


        /// <summary>
        /// Create Nodes: assign Coordiantes, GlobalId and Boundary Conditions
        /// </summary>
        /// <returns> List with nodes incl properties</returns>
        private List<Node> CreateNodes(DataTree<Point3d> meshPoints, int nu, int nv, int nw)
        {
            List<Node> nodes = new List<Node>();
            int count1 = 0;
            nu = nu + 1; //input nu = nu - 1. Exs: nu = 3, total points in u-direction is 4;
            nv = nv + 1; //input nv = nv - 1. Exs: nv = 3, total points in v-direction is 4;
            nw = nw + 1; //input nw = nw - 1. Exs: nw = 3, total points in w-direction is 4;

            for (int i = 0; i < nw ; i++)
            {
                int row = 0;
                int column = 0;
                for (int j = 0; j < meshPoints.Branch(i).Count; j++)
                {
                    Node node = new Node(count1, meshPoints.Branch(i)[j]); // assign Global ID and cooridinates
                    count1++;

                    if (column == 0 | column == nu - 1) { node.BC_U = true; } // assign BCU
                    if (row == 0 | row == nv - 1) { node.BC_V = true; } // assign BCV
                    if (i == 0 | i == nw - 1 ) { node.BC_W = true; } // assign BCW

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

        /// <summary>
        /// Create Elements: assign ElementId, ElementMesh and Nodes incl. Coordiantes, GlobalId, LocalId and Boundary Conditions), elementId, elementMesh.
        /// </summary>
        /// <returns>List with elements incl properties</returns>
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
                    //e.IsCube = true;
                    if (count2 < nu - 1)
                    {
                        Node n1 = new Node(1, nodes[counter].GlobalId, ptsBot[j], nodes[counter].BC_U, nodes[counter].BC_V, nodes[counter].BC_W);
                        e.Node1 = n1;

                        Node n2 = new Node(2, nodes[counter + 1].GlobalId, ptsBot[j + 1], nodes[counter + 1].BC_U, nodes[counter + 1].BC_V, nodes[counter + 1].BC_W);
                        e.Node2 = n2;

                        Node n3 = new Node(3, nodes[counter + nu + 1].GlobalId, ptsBot[j + nu + 1], nodes[counter + nu + 1].BC_U, nodes[counter + nu + 1].BC_V, nodes[counter + nu + 1].BC_W);
                        e.Node3 = n3;

                        Node n4 = new Node(4, nodes[counter + nu].GlobalId, ptsBot[j + nu], nodes[counter + nu].BC_U, nodes[counter + nu].BC_V, nodes[counter + nu].BC_W);
                        e.Node4 = n4;

                        Node n5 = new Node(5, nodes[counter + nu * nv].GlobalId, ptsTop[j], nodes[counter + nu * nv].BC_U, nodes[counter + nu * nv].BC_V, nodes[counter + nu * nv].BC_W);
                        e.Node5 = n5;

                        Node n6 = new Node(6, nodes[counter + 1 + nu * nv].GlobalId, ptsTop[j + 1], nodes[counter + 1 + nu * nv].BC_U, nodes[counter + 1 + nu * nv].BC_V, nodes[counter + 1 + nu * nv].BC_W);
                        e.Node6 = n6;

                        Node n7 = new Node(7, nodes[counter + nu + 1 + nu * nv].GlobalId, ptsTop[j + nu + 1], nodes[counter + nu + 1 + nu * nv].BC_U, nodes[counter + nu + 1 + nu * nv].BC_V, nodes[counter + nu + 1 + nu * nv].BC_W);
                        e.Node7 = n7;

                        Node n8 = new Node(8, nodes[counter + nu + nu * nv].GlobalId, ptsTop[j + nu], nodes[counter + nu + nu * nv].BC_U, nodes[counter + nu + nu * nv].BC_V, nodes[counter + nu + nu * nv].BC_W);
                        e.Node8 = n8;

                        
                        mesh.Vertices.Add(e.Node1.Coordinate); //0
                        mesh.Vertices.Add(e.Node2.Coordinate); //1
                        mesh.Vertices.Add(e.Node3.Coordinate); //2
                        mesh.Vertices.Add(e.Node4.Coordinate); //3
                        mesh.Vertices.Add(e.Node5.Coordinate); //4
                        mesh.Vertices.Add(e.Node6.Coordinate); //5
                        mesh.Vertices.Add(e.Node7.Coordinate); //6
                        mesh.Vertices.Add(e.Node8.Coordinate); //7

                        mesh.Faces.AddFace(0, 1, 5, 4);
                        mesh.Faces.AddFace(1, 2, 6, 5);
                        mesh.Faces.AddFace(2, 3, 7, 6);
                        mesh.Faces.AddFace(3, 0, 4, 7);
                        mesh.Faces.AddFace(0, 1, 2, 3);
                        mesh.Faces.AddFace(4, 5, 6, 7);

                        mesh.Normals.ComputeNormals();  //Control if needed
                        mesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                        mesh.Compact(); //to ensure that it calculate
                        e.mesh = mesh;

                        //add element and mesh to element list
                        elements.Add(e);

                        //clear
                        e = new Element();
                        mesh = new Mesh();

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
        */

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
            get { return new Guid("bf8907fb-fb39-41c7-aa44-c0af8111dccb"); }
        }
    }
}
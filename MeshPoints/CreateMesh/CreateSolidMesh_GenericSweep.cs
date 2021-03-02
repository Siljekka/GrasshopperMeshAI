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
    public class CreateSolidMesh_GenericSweep : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public CreateSolidMesh_GenericSweep()
          : base("Create Mesh3D (GenericSweep)", "mesh3DG",
              "Creates a solid mesh (more generic)",
              "MyPlugIn", "Mesh")
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
            pManager.AddGenericParameter("Mesh3D", "m3D", "Creates a Mesh3D", GH_ParamAccess.item);
            pManager.AddGenericParameter("test", "", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("test2", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Variables
            Brep bp = new Brep();

            Mesh3D m3D = new Mesh3D();
            Element e = new Element();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();

            Mesh mesh = new Mesh();
            Mesh allMesh = new Mesh();

            Point3d[] nwPt;

            List<Point3d> nwPoints = new List<Point3d>();
            List<Point3d> ptsBot = new List<Point3d>();
            List<Point3d> ptsTop = new List<Point3d>();
            DataTree<Point3d> railPoints = new DataTree<Point3d>();
            DataTree<Point3d> points = new DataTree<Point3d>();

            List<Curve> rails = new List<Curve>();

            int nu = 0;
            int nv = 0;
            int nw = 0;
            int row = 0;
            int column = 0;
            int count1 = 0;
            int count2 = 0;
            int counter = 0;
            int elemId = 0;

            string curveOrientation = null;
            Curve[] iCrv1;
            Point3d[] iPt1;
            List<Curve> intCrv = new List<Curve>();
            List<Plane> surfacePlanes = new List<Plane>();

            List<Brep> planarBrep = new List<Brep>();
            List<Curve> interCrv = new List<Curve>();
            Brep brep = new Brep();
            NurbsSurface nurbsSurface = null;
            List<NurbsSurface> nwSurface = new List<NurbsSurface>();
            DataTree<Curve> intersectionCurve = new DataTree<Curve>();
            NurbsSurface surface = null;
            List<Point3d> pt = new List<Point3d>();


            //Input
            DA.GetData(0, ref bp);
            DA.GetData(1, ref nu);
            DA.GetData(2, ref nv);
            DA.GetData(3, ref nw);

            if ((nu == 0) | (nv == 0) | (nw == 0)) { return; } // todo: add warning message


            // Code
            m3D.nu = nu;
            m3D.nv = nv;
            m3D.nw = nw;

            Curve rail1 = bp.Edges[0];  //get edge1 of brep = rail 1
            Curve rail2 = bp.Edges[9];  //get edge2 of brep = rail 2
            Curve rail3 = bp.Edges[10]; //get edge3 of brep = rail 3
            Curve rail4 = bp.Edges[11]; //get edge4 of brep = rail 4

            rails.Add(rail1);
            rails.Add(rail2);
            rails.Add(rail3);
            rails.Add(rail4);
            rails.Reverse();



            //Divide each rail into nw points.
            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(m3D.nw, true, out nwPt);  //divide each rail in nw number of points
                nwPoints = nwPt.OrderBy(f => f.Z).ToList();     //order the points according to z-value
                for (int j = 0; j < nwPoints.Count; j++)
                {
                    railPoints.Add(nwPoints[j], new GH_Path(j)); //tree with nw points on each rail. Branch: floor
                }
            }

            // Make edgecurves for each floor
            for (int i = 0; i < railPoints.BranchCount; i++)
            {
                Plane.FitPlaneToPoints(railPoints.Branch(i), out Plane plane);
                surfacePlanes.Add(plane);

                Intersection.BrepPlane(bp, plane, 0.00001, out iCrv1, out iPt1);
                intCrv = iCrv1.ToList();

                for (int j = 0; j < intCrv.Count; j++) { intCrv[j].MakeClosed(0.0001); intersectionCurve.Add(intCrv[j], new GH_Path(i)); } //
                //if (intersectionCurve.Branch(i).Count != 1) { return; } //todo: add error message "BrerailPointsis not OK."

                planarBrep = Brep.CreatePlanarBreps(intersectionCurve.Branch(i), 0.0001).ToList();

                for (int j = 0; j < planarBrep.Count; j++)
                {
                    nurbsSurface = NurbsSurface.CreateNetworkSurface(planarBrep[j].Edges, 0, 0.000001, 0.000001, 0.000001, out int error);
                    nwSurface.Add(surface);
                }
                planarBrep.Clear();
                intCrv.Clear();
                interCrv.Clear();
            }
            

            #region Makes grid of points in v and u direction at level nw

            double stepU = 1 / ((double)m3D.nu - 1);
            double stepV = 1 / ((double)m3D.nv - 1);
            
            for (int i = 0; i < nwSurface.Count; i++)
            { 
                curveOrientation = intersectionCurve.Branch(i)[0].ClosedCurveOrientation().ToString();

                surface = nwSurface[i];
                surface.SetDomain(0, new Interval(0, 1));
                surface.SetDomain(1, new Interval(0, 1));

                if (curveOrientation == "CounterClockwise")
                {
                    for (double j = 0; j <= 1; j += stepV)
                    {
                        for (double k = 0; k <= 1; k += stepU)
                        {
                            pt.Add(surface.PointAt(j, k));
                        }
                    }
                }
                else
                {
                    for (double j = 0; j <= 1; j += stepV)
                    {
                        for (double k = 1; k >= 0; k -= stepU)
                        {
                            pt.Add(surface.PointAt(j, k));
                        }
                    }
                }
                points.AddRange(pt, new GH_Path(i));
                pt.Clear();
            }
            #endregion


            #region Create Nodes

            // Create Nodes  
            count1 = 0;
            for (int i = 0; i < m3D.nw + 1; i++)
            {
                row = 0;
                column = 0;
                for (int j = 0; j < points.Branch(i).Count; j++)
                {
                    Node node = new Node(count1, points.Branch(i)[j]); //Assign Global ID and cooridinates
                    count1++;

                    if (column == 0 | column == m3D.nu - 1) { node.BC_U = true; }
                    if (row == 0 | row == m3D.nv - 1) { node.BC_V = true; }
                    if (i == 0 | i == nw) { node.BC_W = true; }

                    column++;

                    if (column == m3D.nu)
                    {
                        row++;
                        column = 0;
                    }
                    nodes.Add(node);
                }
            }
            #endregion

            
            #region Create Elements and Mesh
            // Create Elements and Mesh
            
            for (int i = 0; i < points.BranchCount - 1; i++)  // loop levels
            {
                
                ptsBot = points.Branch(i);
                ptsTop = points.Branch(i + 1);
                count2 = 0;
                counter = points.Branch(0).Count*i;

                for (int j = 0; j < points.Branch(0).Count - m3D.nu - 1; j++) // loop elements in a level
                {
                    e.Id = elemId;
                    if (count2 < m3D.nu - 1)
                    {
                        Node n1 = new Node(1, nodes[counter].GlobalId, ptsBot[j], nodes[counter].BC_U, nodes[counter].BC_V, nodes[counter].BC_W);
                        e.Node1 = n1;

                        Node n2 = new Node(2, nodes[counter + 1].GlobalId, ptsBot[j + 1], nodes[counter + 1].BC_U, nodes[counter + 1].BC_V, nodes[counter + 1].BC_W);
                        e.Node2 = n2;

                        Node n3 = new Node(3, nodes[counter + m3D.nu + 1].GlobalId, ptsBot[j + m3D.nu + 1], nodes[counter + m3D.nu + 1].BC_U, nodes[counter + m3D.nu + 1].BC_V, nodes[counter + m3D.nu + 1].BC_W);
                        e.Node3 = n3;

                        Node n4 = new Node(4, nodes[counter + m3D.nu].GlobalId, ptsBot[j + m3D.nu], nodes[counter + m3D.nu].BC_U, nodes[counter + m3D.nu].BC_V, nodes[counter + m3D.nu].BC_W);
                        e.Node4 = n4;

                        Node n5 = new Node(5, nodes[counter + m3D.nu * m3D.nv].GlobalId, ptsTop[j], nodes[counter + m3D.nu * m3D.nv].BC_U, nodes[counter + m3D.nu * m3D.nv].BC_V, nodes[counter + m3D.nu * m3D.nv].BC_W);
                        e.Node5 = n5;

                        Node n6 = new Node(6, nodes[counter + 1 + m3D.nu * m3D.nv].GlobalId, ptsTop[j + 1], nodes[counter + 1 + m3D.nu * m3D.nv].BC_U, nodes[counter + 1 + m3D.nu * m3D.nv].BC_V, nodes[counter + 1 + m3D.nu * m3D.nv].BC_W);
                        e.Node6 = n6;

                        Node n7 = new Node(7, nodes[counter + m3D.nu + 1 + m3D.nu * m3D.nv].GlobalId, ptsTop[j + m3D.nu + 1], nodes[counter + m3D.nu + 1 + m3D.nu * m3D.nv].BC_U, nodes[counter + m3D.nu + 1 + m3D.nu * m3D.nv].BC_V, nodes[counter + m3D.nu + 1 + m3D.nu * m3D.nv].BC_W);
                        e.Node7 = n7;

                        Node n8 = new Node(8, nodes[counter + m3D.nu + m3D.nu * m3D.nv].GlobalId, ptsTop[j + m3D.nu], nodes[counter + m3D.nu + m3D.nu * m3D.nv].BC_U, nodes[counter + m3D.nu + m3D.nu * m3D.nv].BC_V, nodes[counter + m3D.nu + m3D.nu * m3D.nv].BC_W);
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

                        //create global mesh
                        allMesh.Append(mesh);
                        allMesh.Weld(0.01);

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
            #endregion
            

            //Add properties to Mesh3D
            m3D.Nodes = nodes;
            m3D.Elements = elements;
            m3D.mesh = allMesh;

               
            // Output
            DA.SetData(0, m3D);
            DA.SetDataList(1, ptsTop);
            DA.SetDataTree(2, points);

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
            get { return new Guid("bf8907fb-fb39-41c7-aa44-c0af8111dccb"); }
        }
    }
}
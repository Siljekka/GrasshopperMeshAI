﻿using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel.Data;

namespace MeshPoints.CreateMesh
{
    public class CreateMesh3D_Sweep : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public CreateMesh3D_Sweep()
          : base("Create Mesh3D (sweep)", "mesh3D",
              "Creates a solid mesh",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "bp", "Brep", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh2D", "m2D", "Mesh2D for bottom surface of Brep", GH_ParamAccess.item);
            pManager.AddIntegerParameter("w", "w", "division in w direction", GH_ParamAccess.item, 4);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh3D", "m3D", "Creates a Mesh3D", GH_ParamAccess.item);
           // pManager.AddGenericParameter("test", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Variables
            Brep bp = new Brep();
            
            Mesh2D m2D = new Mesh2D();
            Mesh3D m3D = new Mesh3D();
            Element e = new Element();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();

            Mesh mesh = new Mesh();
            Mesh allMesh = new Mesh();

            Point3d[] nwPt;
            Point3d p1 = new Point3d();
            Point3d p2 = new Point3d();
            Point3d p3 = new Point3d();
            Point3d p4 = new Point3d();

            List<Point3d> nwPts = new List<Point3d>();
            List<Point3d> ptsBot = new List<Point3d>();
            List<Point3d> ptsTop = new List<Point3d>();
            DataTree<Point3d> p = new DataTree<Point3d>();
            DataTree<Point3d> pts = new DataTree<Point3d>();

            List<Curve> rails = new List<Curve>();

            int nw = 0;
            int row = 0;
            int column = 0;
            int count1 = 0;
            int count2 = 0;
            int elemId = 0;


            //Input
            DA.GetData(0, ref bp);
            DA.GetData(1, ref m2D);
            DA.GetData(2, ref nw);


            // Code
            m3D.nu = m2D.nu;
            m3D.nv = m2D.nv;
            m3D.nw = nw;

            Curve rail1 = bp.Edges[0];  //get edge1 of brep = rail 1
            Curve rail2 = bp.Edges[9];  //get edge2 of brep = rail 2
            Curve rail3 = bp.Edges[10]; //get edge3 of brep = rail 3
            Curve rail4 = bp.Edges[11]; //get edge4 of brep = rail 4

            rails.Add(rail1);
            rails.Add(rail2);
            rails.Add(rail3);
            rails.Add(rail4);


            //Divide each rail into nw points.
            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(nw, true, out nwPt);  //divide each rail in nw number of points
                nwPts = nwPt.OrderBy(f => f.Z).ToList();     //order the points according to z-value
                for (int j = 0; j < nwPts.Count; j++)
                {
                    p.Add(nwPts[j], new GH_Path(i)); //tree with nw points on each rail. Branch: rail
                }
            }


            #region Makes grid of points in v and u direction at level nw
            //Makes grid of points in v and u direction at level with nw.
            for (int i = 0; i < p.Branch(0).Count; i++) //Decide which w-level the grid are made
            {
                p1 = p.Branch(0)[i]; //points on rail 1
                p2 = p.Branch(1)[i]; //points on rail 2
                p3 = p.Branch(2)[i]; //points on rail 4 
                p4 = p.Branch(3)[i]; //points on rail 3 

                double spanV1 = (p4 - p1).Length / (m2D.nv - 1);  //distance between points on v1 direction
                double spanV2 = (p3 - p2).Length / (m2D.nv - 1);  //distance between points on v2 direction
                Vector3d vecV1 = new Vector3d((p4 - p1) / (p4 - p1).Length);  //vector in v1 direction
                Vector3d vecV2 = new Vector3d((p3 - p2) / (p3 - p2).Length);  //vector in v2 direction

                for (int j = 0; j < m2D.nv; j++) //Loop for v-direction
                {
                    p1 = new Point3d(p1.X + spanV1 * j * vecV1.X, p1.Y + spanV1 * j * vecV1.Y, p1.Z + spanV1 * j * vecV1.Z);  //makes point in v1 direction
                    p2 = new Point3d(p2.X + spanV2 * j * vecV2.X, p2.Y + spanV2 * j * vecV2.Y, p2.Z + spanV2 * j * vecV2.Z);  //makes point in v2 direction

                    double spanU = (p2 - p1).Length / (m2D.nu - 1); //distance between points on u direction
                    Vector3d vecU = new Vector3d((p2 - p1) / (p2 - p1).Length);  //vector in u direction

                    for (int k = 0; k < m2D.nu; k++) //Loop for u-direction
                    {
                        Point3d pt = new Point3d(p1.X + spanU * k * vecU.X, p1.Y + spanU * k * vecU.Y, p1.Z + spanU * k * vecU.Z);  //make points in u direction (between v1 and v2)
                        pts.Add(pt, new GH_Path(i));  //Tree with points in v and u direction. Branch: level nw
                    }
                    p1 = p.Branch(0)[i]; //reset point
                    p2 = p.Branch(1)[i]; //reset point
                }

                #region old colde
                //Brep berep = Brep.CreateFromMesh(m2D.mesh, true);
                //BrepEdgeList bpEdge = berep.Edges;
                //NurbsSurface srf;
                //List<NurbsSurface> srfs = new List<NurbsSurface>();
                //DataTree<Mesh> mes = new DataTree<Mesh>();
                //int y = 1;
                //srf = NurbsSurface.CreateFromCorners(pts.Branch(0)[i], pts.Branch(1)[i], pts.Branch(2)[i], pts.Branch(3)[i]);
                //srfs.Add(srf);
                //me = Mesh.CreateFromSurface(srfs[i]);
                //me.CopyFrom(m2D.mesh);
                //mes.Add(me, new GH_Path(y));
                //y++;
                #endregion
            }
            #endregion


            #region Create Nodes
            // Create Nodes   
            for (int i = 0; i < nw + 1; i++)
            {
                for (int j = 0; j < pts.Branch(i).Count; j++)
                {
                    Node node = new Node(count1, pts.Branch(i)[j]); //Assign Global ID and cooridinates
                    count1++;

                    if (row == 0 | row == m3D.nv - 1) { node.BC_V = true; }
                    if (column == 0 | column == m3D.nu - 1) { node.BC_U = true; }

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

            for (int i = 0; i < pts.BranchCount - 1; i++)
            {
                count2 = 0;

                ptsBot = pts.Branch(i);
                ptsTop = pts.Branch(i + 1);

                for (int j = 0; j < pts.Branch(i).Count - m3D.nu - 1; j++)
                {
                    e.Id = elemId;
                    if (count2 < m3D.nu - 1)
                    {
                        Node n1 = new Node(1, nodes[j].GlobalId, ptsBot[j], nodes[j].BC_U, nodes[j].BC_V);
                        e.Node1 = n1;

                        Node n2 = new Node(2, nodes[j + 1].GlobalId, ptsBot[j + 1], nodes[j + 1].BC_U, nodes[j + 1].BC_V);
                        e.Node2 = n2;

                        Node n3 = new Node(3, nodes[j + m3D.nu + 1].GlobalId, ptsBot[j + m3D.nu + 1], nodes[j + m3D.nu + 1].BC_U, nodes[j + m3D.nu + 1].BC_V);
                        e.Node3 = n3;

                        Node n4 = new Node(4, nodes[j + m3D.nu].GlobalId, ptsBot[j + m3D.nu], nodes[j + m3D.nu].BC_U, nodes[j + m3D.nu].BC_V);
                        e.Node4 = n4;

                        Node n5 = new Node(5, nodes[j + m3D.nu * m3D.nv].GlobalId, ptsTop[j], nodes[j + m3D.nu * m3D.nv].BC_U, nodes[j + m3D.nu * m3D.nv].BC_V);
                        e.Node5 = n5;

                        Node n6 = new Node(6, nodes[j + 1 + m3D.nu * m3D.nv].GlobalId, ptsTop[j + 1], nodes[j + 1 + m3D.nu * m3D.nv].BC_U, nodes[j + 1 + m3D.nu * m3D.nv].BC_V);
                        e.Node6 = n6;

                        Node n7 = new Node(7, nodes[j + m3D.nu + 1 + m3D.nu * m3D.nv].GlobalId, ptsTop[j + m3D.nu + 1], nodes[j + m3D.nu + 1 + m3D.nu * m3D.nv].BC_U, nodes[j + m3D.nu + 1 + m3D.nu * m3D.nv].BC_V);
                        e.Node7 = n7;

                        Node n8 = new Node(8, nodes[j + m3D.nu + m3D.nu * m3D.nv].GlobalId, ptsTop[j + m3D.nu], nodes[j + m3D.nu + m3D.nu * m3D.nv].BC_U, nodes[j + m3D.nu + m3D.nu * m3D.nv].BC_V);
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
                    }
                    else { count2 = 0; }

                }
            }
            #endregion

            //Add properties to Mesh3D
            m3D.Nodes = nodes;
            m3D.Elements = elements;
            m3D.mesh = allMesh;


            // Output
            DA.SetData(0, m3D);
            //DA.SetDataList(1, nodes);


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
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;
using Rhino.Geometry.Intersect;
using Grasshopper;
using Grasshopper.Kernel.Data;

namespace MeshPoints.CreateMesh
{
    public class CreateMesh3D_Box : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateMesh3D_Box class.
        /// </summary>
        public CreateMesh3D_Box()
          : base("Create Mesh3D (box)", "mesh3DBox",
              "Creates a solid mesh for a box-brep",
              "SmartMesh", "Outdated")
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
            //pManager.AddGenericParameter("test", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            /*
            // Variables

            //Point3d[] nwPt;
            Brep bp = new Brep();

            SmartMesh m2D = new SmartMesh();
            SmartMesh m3D = new SmartMesh();
            Element e = new Element();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();

            Mesh m = new Mesh();
            Mesh allMesh = new Mesh();
            Mesh mesh = new Mesh();
            List<Mesh> mBot = new List<Mesh>();
            List<Mesh> mTop = new List<Mesh>();
            List<Mesh> ms = new List<Mesh>();
            DataTree<Mesh> mSrf = new DataTree<Mesh>();

            Point3d[] nwPt;
            List<Point3d> nwPts = new List<Point3d>();
            List<Point3d> ptsBot = new List<Point3d>();
            List<Point3d> ptsTop = new List<Point3d>();
            List<Point3d> pts = new List<Point3d>();


            int nw = 0;
            int row = 0;
            int column = 0;
            int count1 = 0;
            int count2 = 0;
            int elemId = 0;


            DA.GetData(0, ref bp);
            DA.GetData(1, ref m2D);
            DA.GetData(2, ref nw);



            // Code
            m3D.nu = m2D.nu;
            m3D.nv = m2D.nv;
            m3D.nw = nw;


            Curve rail = bp.Edges[1]; // Get rail
            rail.DivideByCount(nw, true, out nwPt); //Divide rail into nw points
            nwPts = nwPt.OrderBy(f => f.Z).ToList(); //Sort points after Z-coordinate


            // Copy the mesh in w-direction and store the mesh in data tree
            mSrf.Add(m2D.mesh, new GH_Path(0)); //Mesh at 0 "floor"

            for (int i = 0; i < nwPts.Count-1; i++)
            {
                Vector3d vec1 = nwPts[i + 1] - nwPts[i];
                double dist = vec1.Length;
                m = m2D.mesh.Offset(dist*(i+1)); //Offset the floor
                mSrf.Add(m, new GH_Path(i+1));   //Mesh at i "floor" (do not catch 0 "floor")
            }

       
            #region Create Vertices and Nodes   
            for (int i = 0; i < nw+1; i++)
            {
                ms= mSrf.Branch(i);  //Get mesh on i "floor"
                pts = (ms[0].Vertices.ToPoint3dArray()).ToList();  //Get Vertices on i "floor"

                for (int j = 0; j < pts.Count; j++) 
                {
                    Node node = new Node(count1, pts[j]); //Assign Global ID and cooridinates
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
            for (int i = 0; i < mSrf.BranchCount-1; i++)
            {
                
                count2 = 0;

                mBot = mSrf.Branch(i);
                mTop = mSrf.Branch(i+1);
                ptsBot = (mBot[0].Vertices.ToPoint3dArray()).ToList();
                ptsTop = (mTop[0].Vertices.ToPoint3dArray()).ToList();


                for (int j = 0; j < pts.Count - m3D.nu-1; j++)
                {
                    e.Id = elemId;
                    if (count2 < m3D.nu-1)
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
                        allMesh.Weld(0.1);

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

            //Add properties to Mesh3D
            m3D.Nodes = nodes;
            m3D.Elements = elements;
            m3D.mesh = allMesh;

            #endregion


            #region old code
            /*
            for (int i = 0; i < nwPts.Count; i++)
            {
                Plane plane = new Plane(nwPts[i], vec);
                planes.Add(plane);
            }

            Interval a = new Interval(0, 20);
            foreach (Plane p in planes)
            {
                Curve[] crv;
                Point3d[] pt;
                Intersection.BrepPlane(bp, p, 0.01, out crv, out pt);
                Brep[] srf = Brep.CreatePlanarBreps(crv, 0.001);
                srfs.Add(srf[0]);
            }
            */
            /*
            List<Point3d> nwPts = new List<Point3d>();

            rail.DivideByCount(nw, true, out nwPt);
            nwPts = nwPt.OrderBy(f => f.Z).ToList();


            Mesh m = meshBot.mesh;
            Mesh m2D = meshBot.mesh;
            Mesh allmesh = new Mesh();

            for (int i = 0; i < nwPts.Count - 1; i++)
            {
                Vector3d vec = nwPts[i + 1] - nwPts[i];
                double dist = vec.Length;
                m = m2D.Offset(dist, true, vec);
                m2D = m2D.Offset(dist, false, vec);

                allmesh.Append(m);
            }
            List<Point3f> vert = new List<Point3f>();
            List<Point3f> pts = new List<Point3f>();

            vert = allmesh.Vertices.OrderByDescending(f => f.X).ToList();
            vert = vert.OrderByDescending(f => f.Y).ToList();
            vert = vert.OrderByDescending(f => f.Z).ToList();


            for (int i = 0; i < vert.Count-1; i++)
            {
                if (vert[i] != vert[i + 1])
                {
                    pts.Add(vert[i]);
                }
            }
            pts.Add(vert[vert.Count-1]);
            */
            /*
            #endregion

            // Output
            DA.SetData(0, m3D);
            //DA.SetDataList(1, nodes);
            */

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
            get { return new Guid("ab9128dd-4d53-449a-a03f-26070f75d9e7"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;
using MeshPoints.Classes;

namespace MeshPoints
{
    public class MeshPoints3D_outdated : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Mesh3D class.
        /// </summary>
        public MeshPoints3D_outdated()
          : base("Mesh3D", "M3D",
              "Creates mesh for solid elements",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "pt", "Points on a Brep", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Mesh on Brep", GH_ParamAccess.item);
            pManager.AddGenericParameter("test", "m", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Variables
            List<Point3d> pts = new List<Point3d>();

            Mesh allMesh = new Mesh();
            Element e = new Element();
            List<Element> elements = new List<Element>();
            List<Point3d> nodePts = new List<Point3d>();

            MeshPoints3D_outdated m = new MeshPoints3D_outdated();


            int u = 2;
            int v = 2;
            int w = 2;

            DA.GetDataList(0, pts);


            //Code

            #region Finds u, v and w
            bool completeU = false;
            bool completeV = false;
            bool completeW = false;

            #region finds u
            for (int i = 0; i < pts.Count - 2; i++)
            {
                Vector3d vec1 = (pts[i + 1] - pts[i]); //distance from start point to point[i+1]
                Vector3d vec2 = (pts[i + 2] - pts[i + 1]); //distance from start point to point[i+2]
                if (vec1.IsParallelTo(vec2) == 1)
                {
                    if (!completeU)
                    {
                        u++; //count inner points in first row
                    }
                }
                else
                {
                    completeU = true;
                }
            }
            #endregion

            #region finds v
            for (int i = 0; i < pts.Count - 2 * u; i = i + u)
            {
                Vector3d vec1 = (pts[i + u] - pts[i]); //distance from start point to point[i+1]
                Vector3d vec2 = (pts[i + 2 * u] - pts[i + u]); //distance from start point to point[i+2]
                if (vec1.IsParallelTo(vec2) == 1)
                {
                    if (!completeV)
                    {
                        v++; //count inner points in first row
                    }
                }
                else
                {
                    completeV = true;
                }
            }
            #endregion

            #region finds w
            for (int i = 0; i < pts.Count - 2 * u * v; i = i + u * v)
            {
                Vector3d vec1 = (pts[i + u * v] - pts[i]); //distance from start point to point[i+1]
                Vector3d vec2 = (pts[i + 2 * u * v] - pts[i + u * v]); //distance from start point to point[i+2]
                if (vec1.IsParallelTo(vec2) == 1)
                {
                    if (!completeW)
                    {
                        w++; //count inner points in first row
                    }
                }
                else
                {
                    completeW = true;
                }
            }
            #endregion

            //List<int> o = new List<int>();
            //o.Add(u);
            //o.Add(v);
            //o.Add(w);
            #endregion

            /*
            #region Nodes and Vertices
            int count1 = 0;
            int count2 = 0;
            int count3 = 0;
            int count = 1;
            int global = 0;

            // Arranging the nodes
            #region Arranging the nodes
            for (int i = 0; i <= (pts.Count - (u+u*v)); i++)
            {
                if (count1 <= (u - 2) )
                {
                    if (count2 <= v - 2)
                    {
                        e.Id = i;
                        e.Node1 = new Node(1, global, pts[i]);
                        e.Node2 = new Node(2, global, pts[(i + 1)]);
                        e.Node3 = new Node(3, global, pts[i + u]);
                        e.Node4 = new Node(4, global, pts[(i + 1) + u]);
                        e.Node5 = new Node(5, global, pts[i + u * v]);
                        e.Node6 = new Node(6, global, pts[(i + 1) + u * v]);
                        e.Node7 = new Node(7, global, pts[i + u + u * v]);
                        e.Node8 = new Node(8, global, pts[(i + 1) + u + u * v]);
                        elements.Add(e);

                        List<Point3d> n = new List<Point3d>()
                        {
                        e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                            e.Node5.Coordinate, e.Node6.Coordinate, e.Node7.Coordinate, e.Node8.Coordinate,
                        };
                        nodePts.AddRange(n);
                        count1++;
                        global++;
                    }
                    count2 = 0;
                }
                else
                {
                    if (i == (v * u)*count)
                    {
                        count++;
                        count3++;
                        count2 = 0;
                    }
                    count1 = 0;
                    count2++;
                }
            }
            #endregion

            // Make Vertices
            #region Make Vertices
            for (int i = 0; i < nodes.Count; i++)
            {
                allMesh.Vertices.Add(nodes[i]);
            }
            #endregion
            
            #endregion
            */




            //Output
            //DA.SetDataList(0, nodes);
            //DA.SetData(1, m);



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
                return Properties.Resources.test;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6ac80e83-5e46-4a2b-b43a-b25931f7d307"); }
        }
    }
}
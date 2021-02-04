using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;

namespace MeshPoints
{
    public class MeshPoints2D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MeshPoints2D class.
        /// </summary>
        public MeshPoints2D()
          : base("MeshPoints2D", "mp2D",
              "Mesh list with flatten points in 2D",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "pts", "Insert flatten list of points", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Mesh between points", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Variables
            Mesh2D m = new Mesh2D();
            Element e = new Element();
            Mesh mesh = new Mesh();
            Mesh allMesh = new Mesh();

            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<Point3d> pts = new List<Point3d>();

            m.Nx = 2;//Number points in x-dir, start by adding first/last point in a row
            m.Ny = 1; //Number points in y-dir, start by adding first point in a colomn
            int counter = 0;

            double dist1 = 0;
            double dist2 = 0;
            Boolean completeRow = false;

            //Input
            DA.GetDataList(0, pts);

            #region Find nx and ny
            for (int i = 0; i < pts.Count - 2; i++)
            {
                dist1 = Math.Abs(pts[0].X - pts[i + 1].X); //distance from start point to point[i+1]
                dist2 = Math.Abs(pts[0].X - pts[i + 2].X);  //distance from start point to point[i+2]
                if (dist1 < dist2)
                {
                    if (!completeRow)
                    {
                        m.Nx++; //count inner points in first row
                    }
                }
                else
                {
                    m.Ny++; //count each end of row
                    completeRow = true;
                }
            }
            #endregion

            #region Create Vertices and Nodes   
            int row = 0;
            int column = 0;
            //OBS: use nodes instead of vertices??
            for (int i = 0; i < pts.Count; i++)
            {
                allMesh.Vertices.Add(pts[i]);
                Node node = new Node(i, pts[i]); //Assign Global ID and cooridinates
                if (row == 0 | row == m.Ny - 1) {node.BC_Y = true;}
                if (column == 0 | column == m.Nx - 1) {node.BC_X = true;}

                column++;

                if (column == m.Nx)
                {
                    row++;
                    column = 0;
                }
                nodes.Add(node);
            }
            #endregion

            #region Create Elements and Mesh
            int newRow = 0;
            for (int i = 0; i < (m.Nx - 1) * (m.Ny - 1); i++)
            {
                //add properties
                e.Id = i;

                e.Node1 = nodes[counter]; //OBS: bug with LocalID... changes nodes when changing e.Node1.LocalId.. wanna make a copy ????????????
                e.Node1.LocalId = 1;

                e.Node2 = nodes[counter + 1];
                e.Node2.LocalId = 2;

                e.Node3 = nodes[counter + m.Nx + 1];
                e.Node3.LocalId = 3;

                e.Node4 = nodes[counter + m.Nx];
                e.Node4.LocalId = 4;

                //create local mesh for element
                mesh.Vertices.Add(e.Node1.Coordinate);
                mesh.Vertices.Add(e.Node2.Coordinate);
                mesh.Vertices.Add(e.Node3.Coordinate);
                mesh.Vertices.Add(e.Node4.Coordinate);
                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.Normals.ComputeNormals();  //Control if needed
                mesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                mesh.Compact(); //to ensure that it calculate
                e.mesh = mesh;

                //add element and mesh to element list
                elements.Add(e);

                //create global mesh
                allMesh.Faces.AddFace(counter, counter + 1, counter + m.Nx + 1, counter + m.Nx);

                //clear
                e = new Element();
                mesh = new Mesh();

                //element counter
                counter++;
                newRow++; ;
                if (newRow == (m.Nx - 1)) //new row
                {
                    counter++;
                    newRow = 0;
                }
            }
            #endregion


            //OBS: burde finne en annen løsning for meshingen...
            allMesh.Normals.ComputeNormals();  //Control if needed
            allMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            allMesh.Compact(); //to ensure that it calculate

            //Add properties to Mesh2D
            m.Nodes = nodes;
            m.Elements = elements;
            m.mesh = allMesh;

            // Output
            DA.SetData(0, m); 
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
            get { return new Guid("c43287d0-d99d-40d8-81d4-b967ec6f8263"); }
        }
    }
}
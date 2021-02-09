using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;

//Create a Mesh2D from list of flatten points

namespace MeshPoints.CreateMesh
{
    public class MeshPoints2D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MeshPoints2D class.
        /// </summary>
        public MeshPoints2D()
          : base("Create Mesh2D", "mp2D",
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
            pManager.AddGenericParameter("Mesh2D", "m", "Mesh2D from given points", GH_ParamAccess.item);
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

            int row = 0;
            int column = 0;
            int newRow = 0;
            int counter = 0;
            double dist1 = 0;
            double dist2 = 0;
            Boolean completeRow = false;

            //Input
            DA.GetDataList(0, pts);

            
            #region Find nu and nv
            m.nu = 2;//Number points in u-dir, start by adding first/last point in a row
            m.nv = 2; //Number points in v-dir, start by adding first point in a colomn
            for (int i = 0; i < pts.Count - 2; i++)
            {
                Vector3d vec1 = (pts[i + 1] - pts[i]) / (pts[i + 1] - pts[i]).Length;
                Vector3d vec2 = (pts[i + 2] - pts[i+1]) / (pts[i + 2] - pts[i+1]).Length;
                if (vec1.IsParallelTo(vec2) == 1)
                {
                    if (!completeRow)
                    {
                        m.nu++; //count inner points in first row
                    }
                }
                else
                {
                    m.nv++; //count each end of row
                    completeRow = true;
                }
            }
            m.nv = m.nv/2;
            #endregion
            
            #region Create Vertices and Nodes   
            //OBS: use nodes instead of vertices??
            for (int i = 0; i < pts.Count; i++)
            {
                allMesh.Vertices.Add(pts[i]);
                Node node = new Node(i, pts[i]); //Assign Global ID and cooridinates
                if (row == 0 | row == m.nv - 1) {node.BC_V = true;}
                if (column == 0 | column == m.nu - 1) {node.BC_U = true;}

                column++;

                if (column == m.nu)
                {
                    row++;
                    column = 0;
                }
                nodes.Add(node);
            }
            #endregion

            #region Create Elements and Mesh
            for (int i = 0; i < (m.nu - 1) * (m.nv - 1); i++)
            {
                //create nodes
                e.Id = i;
                Node n1 = new Node(1, nodes[counter].GlobalId, nodes[counter].Coordinate, nodes[counter].BC_U, nodes[counter].BC_V);
                e.Node1 =  n1;

                Node n2 = new Node(2, nodes[counter+1].GlobalId, nodes[counter+1].Coordinate, nodes[counter+1].BC_U, nodes[counter+1].BC_V);
                e.Node2 = n2;

                Node n3 = new Node(3, nodes[counter + m.nu + 1].GlobalId, nodes[counter + m.nu + 1].Coordinate, nodes[counter + m.nu + 1].BC_U, nodes[counter + m.nu + 1].BC_V);
                e.Node3 = n3;

                Node n4 = new Node(4, nodes[counter + m.nu].GlobalId, nodes[counter + m.nu].Coordinate, nodes[counter + m.nu].BC_U, nodes[counter + m.nu].BC_V);
                e.Node4 = n4;

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

                //create global mesh
                allMesh.Faces.AddFace(counter, counter + 1, counter + m.nu + 1, counter + m.nu);
                
                //add element and mesh to element list
                elements.Add(e);

                //clear
                e = new Element();
                mesh = new Mesh();

                //element counter
                counter++;
                newRow++; ;
                if (newRow == (m.nu - 1)) //new row
                {
                    counter++;
                    newRow = 0;
                }
            }
            
            //OBS: burde finne en annen løsning for meshingen...
            allMesh.Normals.ComputeNormals();  //Control if needed
            allMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            allMesh.Compact(); //to ensure that it calculate

            //Add properties to Mesh2D
            m.Nodes = nodes;
            m.Elements = elements;
            m.mesh = allMesh;
            #endregion
            
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
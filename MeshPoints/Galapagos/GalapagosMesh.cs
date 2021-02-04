using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using System.Drawing;
using MeshPoints.Classes;

// Modify Mesh2D with evolutionary solver (Galapagos) to optimize the mesh quality

namespace MeshPoints
{
    public class GalapagosMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GalapagosMesh class.
        /// </summary>
        public GalapagosMesh()
          : base("GalapagosMesh", "gM",
              "Optimize mesh with gene pool",
              "MyPlugIn", "Evolutionary Solving")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh2D", "m2d", "Input Mesh2D", GH_ParamAccess.item);
            pManager.AddGenericParameter("Gene Pool U-dir", "qp", "Gene pool list", GH_ParamAccess.list);
            pManager.AddGenericParameter("Gene Pool V-dir", "qp", "Gene pool list", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh2D", "m2d", "Updated mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Variables
            Mesh2D m = new Mesh2D();
            List<double> genesX = new List<double>();
            List<double> genesY = new List<double>();

            Element e = new Element();
            Mesh mesh = new Mesh();
            Mesh allMesh = new Mesh();
            Mesh2D mUpdated = new Mesh2D();

            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<Point3d> newPts = new List<Point3d>(); //list with upgraded vertices
            double devX = 0;
            double devY = 0;
            int newRow = 0;
            int counter = 0;

            //Input
            DA.GetData(0, ref m);
            DA.GetDataList(1, genesX);
            DA.GetDataList(2, genesY);

            #region Update nodes
            for (int i = 0; i < m.Nodes.Count; i++)
            {
                //Deviation U
                if (genesX[i] > 0 & !m.Nodes[i].BC_U)
                {
                    devX = Math.Abs(m.Nodes[i].Coordinate.X - m.Nodes[i + 1].Coordinate.X) / 2 * genesX[i];
                }
                else if (genesX[i] < 0 & !m.Nodes[i].BC_U)
                {
                    devX = Math.Abs(m.Nodes[i].Coordinate.X - m.Nodes[i - 1].Coordinate.X) / 2 * genesX[i];
                }
                else { devX = 0; }


                //Deviation V
                if (genesY[i] > 0 & !m.Nodes[i].BC_V)
                {
                    devY = Math.Abs(m.Nodes[i].Coordinate.Y - m.Nodes[i + 1].Coordinate.Y) / 2 * genesY[i]; //nx...
                }
                else if (genesY[i] < 0 & !m.Nodes[i].BC_V)
                {
                    devY = Math.Abs(m.Nodes[i].Coordinate.Y - m.Nodes[i - 1].Coordinate.Y) / 2 * genesY[i]; //nx...                
                }
                else { devY = 0; }

                //Update vertices
                Point3d pt = new Point3d(m.Nodes[i].Coordinate.X + devX, m.Nodes[i].Coordinate.Y + devY, m.Nodes[i].Coordinate.Z + 0);
                Node n = new Node(i, pt, m.Nodes[i].BC_U, m.Nodes[i].BC_V);
                
                nodes.Add(n);
                allMesh.Vertices.Add(pt);
            }
            #endregion

            #region Element and mesh
            for (int i = 0; i < (m.nu - 1) * (m.nv - 1); i++)
            {
                e.Id = i;

                e.Node1 = nodes[counter]; //OBS: bug with LocalID... changes nodes when changing e.Node1.LocalId.. wanna make a copy ????????????
                e.Node1.LocalId = 1;

                e.Node2 = nodes[counter + 1];
                e.Node2.LocalId = 2;

                e.Node3 = nodes[counter + m.nu + 1];
                e.Node3.LocalId = 3;

                e.Node4 = nodes[counter + m.nu];
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
            #endregion

            //OBS: burde finne en annen løsning for meshingen...
            allMesh.Normals.ComputeNormals();  //Control if needed
            allMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            allMesh.Compact(); //to ensure that it calculate

            //Add properties to Mesh2D
            mUpdated.Nodes = nodes;
            mUpdated.Elements = elements;
            mUpdated.mesh = allMesh;

            //Output
            DA.SetData(0, mUpdated);
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
            get { return new Guid("219e8033-a05c-473a-8219-f7a6c96c7256"); }
        }
    }
}
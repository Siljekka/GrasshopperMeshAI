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
            // variables
            Mesh2D m = new Mesh2D();
            Element e = new Element();
            Mesh mesh = new Mesh();
            Mesh allMesh = new Mesh();

            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<Point3d> meshPts = new List<Point3d>();
            List<Point3d> meshPtsSorted = new List<Point3d>();

            int row = 0;
            int column = 0;
            int newRow = 0;
            int counter = 0;
            int numPtsDirection1 = 2; //Number points in 1. direction, start by adding first and last point in the 1. direction
            int numPtsDirection2 = 0;
            Boolean completeCountDirection1 = false;

            // input
            DA.GetDataList(0, meshPts);

            if (meshPts.Count < 4) { return; }// add warning message

            #region Find nu and nv
            // count the number of points in 1. direction and 2. direction
            for (int i = 0; i < meshPts.Count - 2; i++)
            {
                Vector3d vec1 = (meshPts[i + 1] - meshPts[i]);
                Vector3d vec2 = (meshPts[i + 2] - meshPts[i]);
                double dot = Vector3d.Multiply(vec1, vec2);

                if (dot > 0)
                {
                    if (!completeCountDirection1)
                    {
                        numPtsDirection1++; // count points in direction 1
                    }
                }
                else
                {
                    completeCountDirection1 = true;
                }
            }
            numPtsDirection2 = meshPts.Count / numPtsDirection1;

            // assign directions to u and v axis
            bool switchAxis = false;
            Vector3d vecDirection1 = (meshPts[1] - meshPts[0]) / (meshPts[1] - meshPts[0]).Length;
            Vector3d vecDirection2 = (meshPts[numPtsDirection1]-meshPts[0]) / (meshPts[numPtsDirection1] - meshPts[0]).Length;
            if (vecDirection1.X < vecDirection2.X)
            {
                // meshPts has "Column format"
                m.nu = numPtsDirection2; // u is direction 2 
                m.nv = numPtsDirection1; // v is direction 1
                switchAxis = true;
            }
            else
            {
                // meshPts has "Row format"
                m.nu = numPtsDirection1; // u is direction 1
                m.nv = numPtsDirection2; // v is direction 2
            }
            #endregion

            // change meshPts from "Column format" to "Row format":
            if (switchAxis)
            {
                for (int k = 0; k < m.nv; k++)
                {
                    for (int n = 0; n < m.nu; n++)
                    {
                        meshPtsSorted.Add(meshPts[n * m.nv + k]); 
                    }
                }
                meshPts = meshPtsSorted;
            }
            
            #region Create Vertices and Nodes   
            //OBS: should use nodes instead of vertices
            for (int i = 0; i < meshPts.Count; i++)
            {
                allMesh.Vertices.Add(meshPts[i]);
                Node node = new Node(i, meshPts[i]); // assign Global ID and cooridinates
                if (row == 0 | row == m.nv - 1) { node.BC_V = true; } // assign BC v-dir
                if (column == 0 | column == m.nu - 1) { node.BC_U = true; } // assign BC u-dir
                
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
            
            //OBS: should find a better solution for meshing
            allMesh.Normals.ComputeNormals();  //Control if needed
            allMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            allMesh.Compact(); //to ensure that it calculate

            // add properties to Mesh2D
            m.Nodes = nodes;
            m.Elements = elements;
            m.mesh = allMesh;
            m.Nodes = nodes;
            #endregion

            // output
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
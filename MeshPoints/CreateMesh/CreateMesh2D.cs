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
            Mesh globalMesh = new Mesh();

            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<Point3d> meshPts = new List<Point3d>();

            int nu = 0;
            int nv = 0;
            int uSequence = 0;
            int vSequence = 0;
            int counter = 0;
            int numPtsDirectionU = 2; // number points in 1. direction, start by adding first and last point in the 1. direction
            int numPtsDirectionV = 0;
            double dotProduct = 0;
            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;
            Boolean completeCountOfPtsDirectionU = false;

            // input
            DA.GetDataList(0, meshPts);

            if (meshPts.Count < 4) { return; } // add warning message

            #region Find nu and nv
            // count the number of points in 1. direction and 2. direction
            for (int i = 0; i < meshPts.Count - 2; i++)
            {
                vec1 = (meshPts[i + 1] - meshPts[i]);
                vec2 = (meshPts[i + 2] - meshPts[i]);
                dotProduct = Vector3d.Multiply(vec1, vec2);

                if (dotProduct > 0)
                {
                    if (!completeCountOfPtsDirectionU)
                    {
                        numPtsDirectionU++; // count points in direction 1
                    }
                }
                else
                {
                    completeCountOfPtsDirectionU = true;
                }
            }
            numPtsDirectionV = meshPts.Count / numPtsDirectionU;
            nu = numPtsDirectionU;
            nv = numPtsDirectionV; 
            #endregion

            #region Create Vertices and Nodes   
            for (int i = 0; i < meshPts.Count; i++)
            {
                globalMesh.Vertices.Add(meshPts[i]);
                Node node = new Node(i, meshPts[i]); // assign global id and cooridinates
                if (uSequence == 0 | uSequence == nu - 1) { node.BC_U = true; } // assign BC u-dir
                if (vSequence == 0 | vSequence == nv - 1) { node.BC_V = true; } // assign BC v-dir
                
                uSequence++;
                if (uSequence == nu)
                {
                    vSequence++;
                    uSequence = 0;
                }
                nodes.Add(node);
            }
            #endregion

            #region Create Elements and Mesh
            uSequence = 0;
            for (int i = 0; i < (nu - 1) * (nv - 1); i++) // loop elements
            {
                e = CreateElement(i, nodes, counter, nu, nv);
                elements.Add(e); // add element to list of elements
                globalMesh = CreateGlobalMesh(globalMesh, counter, nu, nv);
                
                counter++;
                uSequence++;
                if (uSequence == (nu - 1)) // check if done with a v sequence
                {
                    counter++;
                    uSequence = 0; // new v sequence
                }
            }

            globalMesh = MakeConsistent(globalMesh);
            m = new Mesh2D(nu, nv, nodes, elements, globalMesh);
            #endregion

            // output
            DA.SetData(0, m);
        }

        Mesh MakeConsistent(Mesh m)
        {
            m.Normals.ComputeNormals();  // todo: control if needed
            m.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            m.Compact(); // to ensure that it calculate
            return m;
        }

        Mesh CreateGlobalMesh(Mesh m, int counter, int nu, int nv)
        {
            m.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);
            return m;
        }

        Element CreateElement(int id, List<Node> nodes, int counter, int nu, int nv)
        {
            Element e = new Element();
            e.Id = id;

            Node n1 = new Node(1, nodes[counter].GlobalId, nodes[counter].Coordinate, nodes[counter].BC_U, nodes[counter].BC_V);
            e.Node1 = n1;

            Node n2 = new Node(2, nodes[counter + 1].GlobalId, nodes[counter + 1].Coordinate, nodes[counter + 1].BC_U, nodes[counter + 1].BC_V);
            e.Node2 = n2;

            Node n3 = new Node(3, nodes[counter + nu + 1].GlobalId, nodes[counter + nu+ 1].Coordinate, nodes[counter + nu + 1].BC_U, nodes[counter + nu + 1].BC_V);
            e.Node3 = n3;

            Node n4 = new Node(4, nodes[counter + nu].GlobalId, nodes[counter + nu].Coordinate, nodes[counter + nu].BC_U, nodes[counter + nu].BC_V);
            e.Node4 = n4;

            Mesh m = new Mesh();
            m.Vertices.Add(e.Node1.Coordinate);
            m.Vertices.Add(e.Node2.Coordinate);
            m.Vertices.Add(e.Node3.Coordinate);
            m.Vertices.Add(e.Node4.Coordinate);
            m.Faces.AddFace(0, 1, 2, 3);

            e.mesh = m;
            return e;
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
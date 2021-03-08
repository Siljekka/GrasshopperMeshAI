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
            #region Variables
            Mesh2D m = new Mesh2D();
            Element e = new Element();
            Mesh globalMesh = new Mesh();

            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<Point3d> meshPts = new List<Point3d>();
            int nu = 0; // number of grids in u direction
            int nv = 0; // number of grids in v direction
            #endregion

            // input
            DA.GetDataList(0, meshPts);

            #region Controll if input ok
            if (meshPts.Count < 3)
            {
                throw new ArgumentOutOfRangeException(
                    message: "To few points to mesh",
                    paramName: "meshPt");
            }

            /* Todo: Add warning if the list is not flatten
            if (meshPts == null)
            {
                throw new ArgumentNullException(
                    message: "List of points is invalid",
                    paramName: "meshPt");
            }
            */
            #endregion

            var numberOfGrids = GetNumberOfGrids(meshPts);
            nu = numberOfGrids.Item1;
            nv = numberOfGrids.Item2;

            nodes = CreateNodes(meshPts, nu, nv);
            elements = CreateQuadElements(nodes, nu, nv);
            globalMesh = CreateGlobalMesh(meshPts, nu, nv);

            m = new Mesh2D(nu, nv, nodes, elements, globalMesh);
            // output
            DA.SetData(0, m);
        }

        /// <summary>
        /// count the number of grids in u direction and v direction
        /// </summary>
        /// <param name="meshPts"></param>
        /// <returns></returns>
        Tuple<int, int> GetNumberOfGrids(List<Point3d> meshPts)
        {
            int numGridsInU = 2; // number points in U direction, start by adding first and last point in the u direction
            int numGridsInV = 0;
            double dotProduct = 0;
            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;
            Boolean completeCountOfGridsInU = false;
            for (int i = 0; i < meshPts.Count - 2; i++)
            {
                vec1 = (meshPts[i + 1] - meshPts[i]);
                vec2 = (meshPts[i + 2] - meshPts[i]);
                dotProduct = Vector3d.Multiply(vec1, vec2);

                if (dotProduct > 0)
                {
                    if (!completeCountOfGridsInU) {numGridsInU++; } // count grids in u direction
                }
                else { completeCountOfGridsInU = true;}
            }
            numGridsInV = meshPts.Count / numGridsInU;
            return Tuple.Create(numGridsInU, numGridsInV);
        }

        /// <summary>
        /// Create global nodes by assigning global id, coordinate, boundary condiditon in u and v direction
        /// </summary>
        /// <param name="meshPts"></param>
        /// <returns></returns>
        List<Node> CreateNodes(List<Point3d> meshPts, int nu, int nv)
        {
            List<Node> nodes = new List<Node>();
            int uSequence = 0;
            int vSequence = 0;
            for (int i = 0; i < meshPts.Count; i++)
            {
                Node node = new Node(i, meshPts[i]); // assign global id and cooridinates

                // assign boundary condition
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
            return nodes;
        }

        /// <summary>
        /// Creates quad elements by assigning id, local nodes and mesh.
        /// </summary>
        /// <param name="nodes"></param>
        /// <param name="nu"></param>
        /// <param name="nv"></param>
        /// <returns></returns>
        List<Element> CreateQuadElements(List<Node> nodes, int nu, int nv)
        {
            Element e = new Element();
            List<Element> elements = new List<Element>();
            int uSequence = 0;
            int counter = 0;

            for (int i = 0; i < (nu - 1) * (nv - 1); i++) // loop elements
            {
                Node n1 = new Node(1, nodes[counter].GlobalId, nodes[counter].Coordinate, nodes[counter].BC_U, nodes[counter].BC_V);
                Node n2 = new Node(2, nodes[counter + 1].GlobalId, nodes[counter + 1].Coordinate, nodes[counter + 1].BC_U, nodes[counter + 1].BC_V);
                Node n3 = new Node(3, nodes[counter + nu + 1].GlobalId, nodes[counter + nu + 1].Coordinate, nodes[counter + nu + 1].BC_U, nodes[counter + nu + 1].BC_V);
                Node n4 = new Node(4, nodes[counter + nu].GlobalId, nodes[counter + nu].Coordinate, nodes[counter + nu].BC_U, nodes[counter + nu].BC_V);
                Mesh m = new Mesh();
                m.Vertices.Add(n1.Coordinate);
                m.Vertices.Add(n2.Coordinate);
                m.Vertices.Add(n3.Coordinate);
                m.Vertices.Add(n4.Coordinate);
                m.Faces.AddFace(0, 1, 2, 3);

                e = new Element(i, n1, n2, n3, n4, m);
                elements.Add(e); // add element to list of elements

                counter++;
                uSequence++;
                if (uSequence == (nu - 1)) // check if done with a v sequence
                {
                    counter++;
                    uSequence = 0; // new v sequence
                }
            }
            return elements;
        }

        /// <summary>
        /// Creates a global mesh for the geometry.
        /// </summary>
        /// <param name="meshPt"></param>
        /// <param name="nu"></param>
        /// <param name="nv"></param>
        /// <returns></returns>
        Mesh CreateGlobalMesh(List<Point3d> meshPt, int nu, int nv)
        {
            Mesh globalMesh = new Mesh();
            int counter = 0;
            int uSequence = 0;

            foreach (Point3d pt in meshPt)
            {
                globalMesh.Vertices.Add(pt);
            }
            for (int i = 0; i < (nu - 1) * (nv - 1); i++) // loop elements
            {
                globalMesh.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);

                counter++;
                uSequence++;
                if (uSequence == (nu - 1)) // check if done with a v sequence
                {
                    counter++;
                    uSequence = 0; // new v sequence
                }
            }
            globalMesh.Normals.ComputeNormals();  // todo: control if needed
            globalMesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            globalMesh.Compact(); // to ensure that it calculate

            return globalMesh;
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
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using System.Drawing;
using MeshPoints.Classes;
using Rhino.Geometry.Intersect;

// Move mesh vertices of a Mesh2D with gene pools. Use with evolutionary solver to optimize mesh quality.

namespace MeshPoints
{
    public class GalapagosMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GalapagosMesh class.
        /// </summary>
        public GalapagosMesh()
          : base("Move Mesh Vertices", "mmv",
              "Move mesh vertices with gene pools",
              "MyPlugIn", "Modify Mesh2D")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "srf", "Input source surface", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh2D", "m2d", "Input Mesh2D", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Gene pool for translation in u direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("v genes", "qp", "Gene pool for translation in v direction", GH_ParamAccess.list);
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
            // variables
            Mesh2D m = new Mesh2D();
            Brep srf = new Brep();
            List<double> genesU = new List<double>();
            List<double> genesV = new List<double>();

            double overlapTolerance = 0.99; // ensure no collision of vertices, reduce number to avoid "the look of triangles".
            int newRow = 0;
            int counter = 0;
            Point3d meshPoint = new Point3d();
            Node n = new Node();
            Element e = new Element();
            Mesh mesh = new Mesh();
            Mesh allMesh = new Mesh();
            Mesh2D meshUpdated = new Mesh2D();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
    
            // input
            DA.GetData(0, ref srf);
            DA.GetData(1, ref m);
            DA.GetDataList(2, genesU);
            DA.GetDataList(3, genesV);

            if ( (genesU.Count < m.Nodes.Count) | (genesV.Count < m.Nodes.Count)) { return; }// add warning message

            #region Update nodes
            Vector3d translationVectorUDirection = (m.Nodes[1].Coordinate - m.Nodes[0].Coordinate) * 0;  //dummy-vector: only to be able to assign value later
            Vector3d translationVectorVDirection = (m.Nodes[1].Coordinate - m.Nodes[0].Coordinate) * 0;  //dummy-vector: only to be able to assign value later
            for (int i = 0; i < m.Nodes.Count; i++)
            {
                // translation in u direction
                if (genesU[i] > 0 & !m.Nodes[i].BC_U)
                {
                    translationVectorUDirection = 0.5 * (m.Nodes[i+1].Coordinate - m.Nodes[i].Coordinate) * genesU[i];
                }
                else if (genesU[i] < 0 & !m.Nodes[i].BC_U)
                {
                    translationVectorUDirection = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - 1].Coordinate) * genesU[i];
                }
                else { translationVectorUDirection = translationVectorUDirection * 0;  }

                // translation in v direction
                if (genesV[i] > 0 & !m.Nodes[i].BC_V)
                {
                    translationVectorVDirection = 0.5 * (m.Nodes[i + m.nu].Coordinate - m.Nodes[i].Coordinate) * genesV[i];

                }
                else if (genesV[i] < 0 & !m.Nodes[i].BC_V)
                {
                    translationVectorVDirection = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i-m.nu].Coordinate) * genesV[i];
                }
                else { translationVectorVDirection = translationVectorVDirection * 0; }

                meshPoint = new Point3d(m.Nodes[i].Coordinate.X + (translationVectorUDirection.X + translationVectorVDirection.X)*overlapTolerance ,
                    m.Nodes[i].Coordinate.Y + (translationVectorUDirection.Y + translationVectorVDirection.Y)* overlapTolerance,
                    m.Nodes[i].Coordinate.Z + 0);

                // project meshPoint to surface
                var meshPointProjected = Intersection.ProjectPointsToBreps(
                    new List<Brep> { srf }, // brep on which to project
                    new List<Point3d> { meshPoint }, // some random points to project
                    new Vector3d(0, 0, 1), // project on Z axis
                    0.01);

                n = new Node(i, meshPointProjected[0], m.Nodes[i].BC_U, m.Nodes[i].BC_V);
                nodes.Add(n);
                allMesh.Vertices.Add(meshPointProjected[0]);
            }
            #endregion

            #region Element and mesh
            for (int i = 0; i < (m.nu - 1) * (m.nv - 1); i++)
            {
                e.Id = i;

                e.Node1 = nodes[counter];
                e.Node1.LocalId = 1;

                e.Node2 = nodes[counter + 1];
                e.Node2.LocalId = 2;

                e.Node3 = nodes[counter + m.nu + 1];
                e.Node3.LocalId = 3;

                e.Node4 = nodes[counter + m.nu];
                e.Node4.LocalId = 4;

                // create local mesh for element
                mesh.Vertices.Add(e.Node1.Coordinate);
                mesh.Vertices.Add(e.Node2.Coordinate);
                mesh.Vertices.Add(e.Node3.Coordinate);
                mesh.Vertices.Add(e.Node4.Coordinate);
                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.Normals.ComputeNormals();  // control if needed
                mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
                mesh.Compact(); // to ensure that it calculate
                e.mesh = mesh;

                // create global mesh
                allMesh.Faces.AddFace(counter, counter + 1, counter + m.nu + 1, counter + m.nu);

                // add element to list of elements
                elements.Add(e);

                // clear
                e = new Element();
                mesh = new Mesh();

                // element counter
                counter++;
                newRow++; ;
                if (newRow == (m.nu - 1)) //new row
                {
                    counter++;
                    newRow = 0;
                }
            }
            #endregion

            // OBS: should find a better way to mesh
            allMesh.Normals.ComputeNormals();  // control if needed
            allMesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            allMesh.Compact(); // to ensure that it calculate

            // add properties to Mesh2D
            meshUpdated.Nodes = nodes;
            meshUpdated.Elements = elements;
            meshUpdated.mesh = allMesh;

            // output
            DA.SetData(0, meshUpdated);
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
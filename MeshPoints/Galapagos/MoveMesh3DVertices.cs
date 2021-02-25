using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using System.Drawing;
using MeshPoints.Classes;
using Rhino.Geometry.Intersect;


namespace MeshPoints.Galapagos
{
    public class GalapagosMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MoveMesh3DVertices class.
        /// </summary>
        public GalapagosMesh()
          : base("Move Mesh3D Vertices", "m3dv",
              "Move mesh vertices with gene pools",
              "MyPlugIn", "Modify Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "bp", "Input source brep", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh3D", "m2d", "Input Mesh2D", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Gene pool for translation in u direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("v genes", "qp", "Gene pool for translation in v direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("w genes", "qp", "Gene pool for translation in w direction", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh3D", "m3d", "Updated mesh", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // variables
            Mesh3D m = new Mesh3D();
            Brep brep = new Brep();
            List<double> genesU = new List<double>();
            List<double> genesV = new List<double>();
            List<double> genesW = new List<double>();

            double overlapTolerance = 0.99; // ensure no collision of vertices, reduce number to avoid "the look of triangles".
            int newRow = 0;
            int counter = 0;
            Point3d meshPoint = new Point3d();
            Node n = new Node();
            Element e = new Element();
            Mesh mesh = new Mesh();
            Mesh globalMesh = new Mesh();
            Mesh3D meshUpdated = new Mesh3D();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();

            // input
            DA.GetData(0, ref brep);
            DA.GetData(1, ref m);
            DA.GetDataList(2, genesU);
            DA.GetDataList(3, genesV);
            DA.GetDataList(4, genesW);

            if ((genesU.Count < m.Nodes.Count) | (genesV.Count < m.Nodes.Count) | (genesW.Count < m.Nodes.Count)) { return; }// add warning message


            Point3d testPoint = new Point3d();
            Point3d meshPointProjected = new Point3d();
            
            double distanceToFace = 1;
            double distanceToCurve = 1;
            
            BrepFace face = null;
            BrepEdge edge = null;
            

            #region Update nodes
            Vector3d translationVectorUDirection = (m.Nodes[1].Coordinate - m.Nodes[0].Coordinate) * 0;  //dummy-vector: only to be able to assign value later
            Vector3d translationVectorVDirection = (m.Nodes[1].Coordinate - m.Nodes[0].Coordinate) * 0;  //dummy-vector: only to be able to assign value later
            Vector3d translationVectorWDirection = (m.Nodes[1].Coordinate - m.Nodes[0].Coordinate) * 0;  //dummy-vector: only to be able to assign value later
            
            BrepFaceList brepFace = brep.Faces; //add edges of brep to brepEdge
            BrepEdgeList brepEdge = brep.Edges;
            Curve edgeCurve;

            for (int i = 0; i < m.Nodes.Count; i++)
            {
                bool IsOnFace = false;
                bool IsOnEdge = false;
                foreach (BrepFace bFace in brepFace) // check if node is on edge
                {
                    bFace.ClosestPoint(m.Nodes[i].Coordinate, out double PointOnCurveU, out double PointOnCurveV);
                    testPoint = bFace.PointAt(PointOnCurveU, PointOnCurveV);  // make test point 
                    distanceToFace = testPoint.DistanceTo(m.Nodes[i].Coordinate); // calculate distance between testPoint and node
                    if (distanceToFace <= 0.000001 & distanceToFace >= -0.000001) // if distance = 0: node is on edge
                    {
                        if (m.Nodes[i].BC_U & m.Nodes[i].BC_V & m.Nodes[i].BC_W) // cornerpoints
                        {
                            IsOnFace = false;
                            IsOnEdge = false;
                        }
                        else if ((!m.Nodes[i].BC_U & !m.Nodes[i].BC_V) | (!m.Nodes[i].BC_U & !m.Nodes[i].BC_W) | (!m.Nodes[i].BC_V & !m.Nodes[i].BC_W))
                        {
                            IsOnFace = true;
                            IsOnEdge = false;
                            face = bFace;
                        }
                        else if ((m.Nodes[i].BC_U & m.Nodes[i].BC_V) | (m.Nodes[i].BC_U & m.Nodes[i].BC_W) | (m.Nodes[i].BC_V & m.Nodes[i].BC_W))
                        {
                            IsOnFace = false;
                            IsOnEdge = true;
                            foreach (BrepEdge bEdge in brepEdge) // check if node is on edge
                            {
                                bEdge.ClosestPoint(m.Nodes[i].Coordinate, out double PointOnCurve);
                                testPoint = bEdge.PointAt(PointOnCurve);  // make test point 
                                distanceToCurve = testPoint.DistanceTo(m.Nodes[i].Coordinate); // calculate distance between testPoint and node
                                if (distanceToCurve <= 0.000001 & distanceToCurve >= -0.000001) { edge = bEdge; } // if distance = 0: node is on edge
                            }
                        }
                    }
                }



                // translation in u direction
                if (genesU[i] >= 0 & !m.Nodes[i].BC_U) // not restrained in U
                {
                    translationVectorUDirection = 0.5 * (m.Nodes[i + 1].Coordinate - m.Nodes[i].Coordinate) * genesU[i];
                    if (IsOnFace) //if nodes is on edge, set new meshPoint
                    {
                        translationVectorUDirection = 0.5 * (m.Nodes[i + 1].Coordinate - m.Nodes[i].Coordinate) * genesU[i];
                        meshPoint = new Point3d(m.Nodes[i].Coordinate.X + (translationVectorUDirection.X + translationVectorVDirection.X + translationVectorWDirection.X) * overlapTolerance,
                            m.Nodes[i].Coordinate.Y + (translationVectorUDirection.Y + translationVectorVDirection.Y + translationVectorWDirection.Y) * overlapTolerance,
                            m.Nodes[i].Coordinate.Z + (translationVectorUDirection.Z + translationVectorVDirection.Z + translationVectorWDirection.Z) * overlapTolerance);
                    }
                    else if (IsOnEdge)
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(m.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(m.Nodes[i + 1].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(0.49 * genesU[i]); // move node along edgeCurve
                    }


                }
                else if (genesU[i] <= 0 & !m.Nodes[i].BC_U)
                {
                    translationVectorUDirection = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - 1].Coordinate) * genesU[i];
                    if (IsOnFace) //if nodes is on edge, set new meshPoint
                    {

                    }
                    else if (IsOnEdge)
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(m.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(m.Nodes[i - 1].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(-0.49 * genesU[i]); // move node along edgeCurve
                    }
                }
                else { translationVectorUDirection = translationVectorUDirection * 0; } // restrained in U

                // translation in v direction
                if (genesV[i] >= 0 & !m.Nodes[i].BC_V) // not restrained in V
                {
                    translationVectorVDirection = 0.5 * (m.Nodes[i + m.nu].Coordinate - m.Nodes[i].Coordinate) * genesV[i];
                    if (IsOnFace) //if nodes is on edge, set new meshPoint
                    {

                    }
                    else if (IsOnEdge)
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(m.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(m.Nodes[i + m.nu].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(0.49 * genesV[i]); // move node along edgeCurve
                    }
                }
                else if (genesV[i] <= 0 & !m.Nodes[i].BC_V)
                {
                    translationVectorVDirection = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - m.nu].Coordinate) * genesV[i];
                    if (IsOnFace) //if nodes is on edge, set new meshPoint
                    {

                    }
                    else if (IsOnEdge)
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(m.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(m.Nodes[i - m.nu].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(-0.49 * genesV[i]); // move node along edgeCurve
                    }
                }
                else { translationVectorVDirection = translationVectorVDirection * 0; } // restrained in V

                // translation in w direction
                if (genesW[i] >= 0 & !m.Nodes[i].BC_W) // not restrained in W
                {
                    translationVectorWDirection = 0.5 * (m.Nodes[i + m.nu * m.nv].Coordinate - m.Nodes[i].Coordinate) * genesW[i];
                    if (IsOnFace) //if nodes is on edge, set new meshPoint
                    {

                    }
                    else if (IsOnEdge)
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(m.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(m.Nodes[i + m.nu * m.nv].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(0.49 * genesW[i]); // move node along edgeCurve
                    }
                }
                else if (genesW[i] <= 0 & !m.Nodes[i].BC_W)
                {
                    translationVectorWDirection = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - m.nu * m.nv].Coordinate) * genesW[i];
                    if (IsOnFace) //if nodes is on edge, set new meshPoint
                    {

                    }
                    else if (IsOnEdge)
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(m.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(m.Nodes[i - m.nu * m.nv].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(-0.49 * genesW[i]); // move node along edgeCurve
                    }
                }
                else { translationVectorWDirection = translationVectorWDirection * 0; } // restrained in W

                if (IsOnFace)
                {
                    meshPoint = new Point3d(m.Nodes[i].Coordinate.X + (translationVectorUDirection.X + translationVectorVDirection.X + translationVectorWDirection.X) * overlapTolerance,
                        m.Nodes[i].Coordinate.Y + (translationVectorUDirection.Y + translationVectorVDirection.Y + translationVectorWDirection.Y) * overlapTolerance,
                        m.Nodes[i].Coordinate.Z + (translationVectorUDirection.Z + translationVectorVDirection.Z + translationVectorWDirection.Z) * overlapTolerance);

                    Brep srf = face.DuplicateFace(false);
                    meshPoint = srf.ClosestPoint(meshPoint); // "Project" meshPoint to surface.
                }

                if (!IsOnEdge & !IsOnFace)
                {
                    meshPoint = new Point3d(m.Nodes[i].Coordinate.X + (translationVectorUDirection.X + translationVectorVDirection.X + translationVectorWDirection.X) * overlapTolerance,
                        m.Nodes[i].Coordinate.Y + (translationVectorUDirection.Y + translationVectorVDirection.Y + translationVectorWDirection.Y) * overlapTolerance,
                        m.Nodes[i].Coordinate.Z + (translationVectorUDirection.Z + translationVectorVDirection.Z + translationVectorWDirection.Z) * overlapTolerance);
                }

                meshPointProjected = brep.ClosestPoint(meshPoint); // "Project" meshPoint to surface.
                
                // todo: fix projecting onto brep
                // project meshPoint to brep
                // var meshPointProjected = Intersection.ProjectPointsToBreps(
                //   new List<Brep> { brep }, // brep on which to project
                // new List<Point3d> { meshPoint }, // some random points to project
                //new Vector3d(0, 0, 1), // project on Z axis
                // 0.01);
                //n = new Node(i, meshPointProjected[0], m.Nodes[i].BC_U, m.Nodes[i].BC_V, m.Nodes[i].BC_W);

                n = new Node(i, meshPointProjected, m.Nodes[i].BC_U, m.Nodes[i].BC_V, m.Nodes[i].BC_W); // todo: fix local id;
                nodes.Add(n);
                globalMesh.Vertices.Add(meshPointProjected);
            }
            #endregion

            #region Element and mesh3D
            counter = 0;
            for (int j = 0; j < m.nw ; j++) // loop trough levels
            {
                counter = m.nu * m.nv*j;
                for (int i = 0; i < (m.nu - 1) * (m.nv - 1); i++) // mesh a level
                {
                    int id = i * (j + 1); // element id
                    e = CreateElement(id, nodes, counter, m.nu, m.nv);
                    elements.Add(e); // add element and mesh to element list
                    globalMesh = CreateGlobalMesh(globalMesh, counter, m.nu, m.nv);

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
            }
            #endregion
            // todo: should find a better way to mesh
            globalMesh = MakeConsistent(globalMesh);
            meshUpdated = new Mesh3D(m.nu, m.nv, m.nw, nodes, elements, globalMesh);

            // Output
            DA.SetData(0, meshUpdated);
        }

        Mesh MakeConsistent(Mesh mesh)
        {
            mesh.Normals.ComputeNormals();  // todo: control if needed
            mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            mesh.Compact(); // to ensure that it calculate
            return mesh;
        }

        Mesh CreateGlobalMesh(Mesh m, int counter, int nu, int nv)
        {
            int meshPtsAtLevel = nu * nv;
            m.Faces.AddFace(counter, counter + 1, counter + meshPtsAtLevel + 1, counter + meshPtsAtLevel);
            m.Faces.AddFace(counter + 1, counter + nu + 1, counter + meshPtsAtLevel + nu + 1, counter + meshPtsAtLevel + 1);
            m.Faces.AddFace(counter + nu + 1, counter + nu, counter + meshPtsAtLevel + nu, counter + meshPtsAtLevel + nu + 1);
            m.Faces.AddFace(counter + nu, counter, counter + meshPtsAtLevel, counter + meshPtsAtLevel + nu);
            m.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);
            m.Faces.AddFace(counter + meshPtsAtLevel, counter + meshPtsAtLevel + 1, counter + meshPtsAtLevel + nu + 1, counter + meshPtsAtLevel + nu);
            return m;
        }

        Element CreateElement(int id, List<Node> nodes, int counter, int nu, int nv)
        {
            Element e = new Element();
            int meshPtsAtLevel = nu * nv;
            
            e.Id = id;
            e.Node1 = nodes[counter];
            e.Node1.LocalId = 1;
            e.Node2 = nodes[counter + 1];
            e.Node2.LocalId = 2;
            e.Node3 = nodes[counter + nu + 1];
            e.Node3.LocalId = 3;
            e.Node4 = nodes[counter + nu];
            e.Node4.LocalId = 4;
            e.Node5 = nodes[counter + meshPtsAtLevel];
            e.Node5.LocalId = 5;
            e.Node6 = nodes[counter + meshPtsAtLevel + 1];
            e.Node6.LocalId = 6;
            e.Node7 = nodes[counter + meshPtsAtLevel + nu + 1];
            e.Node7.LocalId = 7;
            e.Node8 = nodes[counter + meshPtsAtLevel + nu];
            e.Node8.LocalId = 8;

            Mesh mesh = new Mesh();
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
            mesh = MakeConsistent(mesh);
            e.mesh = mesh;
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
            get { return new Guid("6de864c8-742d-4bba-96f7-54f3577082c0"); }
        }
    }
}
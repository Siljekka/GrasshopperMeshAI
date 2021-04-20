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
        /// Initializes a new instance of the MoveMesh2DVertices class.
        /// </summary>
        public GalapagosMesh()
          : base("Move Mesh Vertices", "mmv",
              "Move mesh vertices with gene pools",
              "MyPlugIn", "Modify Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Geometry", "geo", "Input source geometry", GH_ParamAccess.item);
            pManager.AddGenericParameter("SmartMesh", "sm", "Input a SmartMesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Gene pool for translation in u direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("v genes", "qp", "Gene pool for translation in v direction", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "sm", "Updated mesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "m", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Variables
            Brep srf = new Brep();
            Mesh3D inputMesh = new Mesh3D();
            List<double> genesU = new List<double>();
            List<double> genesV = new List<double>();

            // Input
            DA.GetData(0, ref srf);
            DA.GetData(1, ref inputMesh);
            DA.GetDataList(2, genesU);
            DA.GetDataList(3, genesV);

            if (!DA.GetData(0, ref srf)) return;
            if (!DA.GetData(1, ref inputMesh)) return;
            if (!DA.GetDataList(2, genesU)) return;
            if (!DA.GetDataList(3, genesV)) return;

            Mesh3D meshUpdated = new Mesh3D();
            Mesh allMesh = new Mesh();
            Node n = new Node();
            Element e = new Element();
            BrepEdge edge = null;
            Curve edgeCurve;
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            Point3d testPoint = new Point3d();
            Point3d meshPoint = new Point3d();
            Point3d meshPointProjected = new Point3d();
            double overlapTolerance = 0.99; // ensure no collision of vertices, reduce number to avoid "the look of triangles".
            double distanceToCurve = 1;
            int newRow = 0;
            int counter = 0;


            if ( (genesU.Count < inputMesh.Nodes.Count) | (genesV.Count < inputMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Too few genes"); return; }// add warning message
            
            /* Todo:
             *  Fikse inn i metoder, finn ut hvilke metoder som kan vær i klasse
             *  Fiks sjekk av input
             */

            #region Update nodes
            BrepEdgeList brepEdge = srf.Edges;  //add edges of surface to brepEdge
            Vector3d translationVectorUDirection = Vector3d.Zero; //dummy-vector: only to be able to assign value later
            Vector3d translationVectorVDirection = Vector3d.Zero; //dummy-vector: only to be able to assign value later
           

            for (int i = 0; i < inputMesh.Nodes.Count; i++)
            {
                bool IsOnCurve = false;
                foreach (BrepEdge bEdge in brepEdge) // check if node is on edge
                {
                    if (inputMesh.Nodes[i].BC_U & inputMesh.Nodes[i].BC_V) { IsOnCurve = false; }
                    else 
                    {
                        IsOnCurve = inputMesh.Nodes[i].IsOnEdge(bEdge);
                        edge = bEdge;
                    }
                    /*
                    bEdge.ClosestPoint(inputMesh.Nodes[i].Coordinate, out double PointOnCurve);
                    testPoint = bEdge.PointAt(PointOnCurve);  // make test point 
                    distanceToCurve = testPoint.DistanceTo(inputMesh.Nodes[i].Coordinate); // calculate distance between testPoint and node
                    if (distanceToCurve <= 0.000001 & distanceToCurve >= -0.000001) // if distance = 0: node is on edge
                    {
                        if (inputMesh.Nodes[i].BC_U & inputMesh.Nodes[i].BC_V) { IsOnCurve = false; } // cornerpoints: IsOnCurve must be false
                        else { IsOnCurve = true; edge = bEdge; }
                    }*/
                }
                
                // translation in u direction
                if (genesU[i] >= 0 & !inputMesh.Nodes[i].BC_U) // not restrained in U
                {
                translationVectorUDirection = 0.5 * (inputMesh.Nodes[i + 1].Coordinate - inputMesh.Nodes[i].Coordinate) * genesU[i]; // for all nodes not on edge

                    if (IsOnCurve) //if nodes is on edge, set new meshPoint
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(inputMesh.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(inputMesh.Nodes[i + 1].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(0.49 * genesU[i]); // move node along edgeCurve
                    }
                }
                else if (genesU[i] <= 0 & !inputMesh.Nodes[i].BC_U) // not restrained in U
                {
                    translationVectorUDirection = 0.5 * (inputMesh.Nodes[i].Coordinate - inputMesh.Nodes[i - 1].Coordinate) * genesU[i]; // for all nodes not on edge

                    if (IsOnCurve) //if node is on edge, set new meshPoint
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(inputMesh.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(inputMesh.Nodes[i - 1].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(-0.49 * genesU[i]); // move node along edgeCurve
                    }
                }
                else { translationVectorUDirection = translationVectorUDirection * 0; } // restrained in U

                // translation in v direction
                if (genesV[i] >= 0 & !inputMesh.Nodes[i].BC_V) // not restrained in V
                { 
                    translationVectorVDirection = 0.5 * (inputMesh.Nodes[i + inputMesh.nu].Coordinate - inputMesh.Nodes[i].Coordinate) * genesV[i]; // for all nodes not on edge

                    if (IsOnCurve) //if node is on edge, set new meshPoint
                    {
                        edgeCurve = edge.DuplicateCurve(); 
                        edgeCurve.SetStartPoint(inputMesh.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(inputMesh.Nodes[i + inputMesh.nu].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(0.49 * genesV[i]); // move node along edgeCurve
                    }
                }
                else if (genesV[i] <= 0 & !inputMesh.Nodes[i].BC_V) // not restrained in V
                {
                    translationVectorVDirection = 0.5 * (inputMesh.Nodes[i].Coordinate - inputMesh.Nodes[i - inputMesh.nu].Coordinate) * genesV[i]; // for all nodes not on edge

                    if (IsOnCurve) //if point is on edge, set new meshPoint
                    {
                        edgeCurve = edge.DuplicateCurve();
                        edgeCurve.SetStartPoint(inputMesh.Nodes[i].Coordinate); //forces start point of edgeCurve
                        edgeCurve.SetEndPoint(inputMesh.Nodes[i - inputMesh.nu].Coordinate); //forces end point of edgeCurve
                        meshPoint = edgeCurve.PointAtNormalizedLength(-0.49 * genesV[i]); // move node along edgeCurve
                    }
                }
                else { translationVectorVDirection = translationVectorVDirection * 0; } // restrained in V
                
                if (!IsOnCurve) // if point is NOT on edge, set new meshPoint
                {
                   meshPoint = new Point3d(inputMesh.Nodes[i].Coordinate.X + (translationVectorUDirection.X + translationVectorVDirection.X) * overlapTolerance,
                       inputMesh.Nodes[i].Coordinate.Y + (translationVectorUDirection.Y + translationVectorVDirection.Y) * overlapTolerance,
                       inputMesh.Nodes[i].Coordinate.Z + 0);
                }
                    
                meshPointProjected = srf.ClosestPoint(meshPoint); // "Project" meshPoint to surface.

                // Old code - replaced with line above: project meshPoint to surface
                //List<Point3d> testP = new List<Point3d>();
                /*
                meshPointProjected = Intersection.ProjectPointsToBreps(
                    new List<Brep> { srf }, // brep on which to project
                    new List<Point3d> { meshPoint }, // some random points to project
                    new Vector3d(0, 0, 1), // project on Z axis
                    0.01);
                */

                n = new Node(i, meshPointProjected, inputMesh.Nodes[i].BC_U, inputMesh.Nodes[i].BC_V);
                nodes.Add(n);
                allMesh.Vertices.Add(meshPointProjected);
                
            }
            #endregion
            
            #region Element and mesh
            for (int i = 0; i < (inputMesh.nu - 1) * (inputMesh.nv - 1); i++)
            {
                Mesh mesh = new Mesh();
                List<Node> elementNodes = new List<Node>();
                List<int> connectivity = new List<int>();
                connectivity.Add(counter);
                connectivity.Add(counter + 1);
                connectivity.Add(counter + inputMesh.nu + 1);
                connectivity.Add(counter + inputMesh.nu);

                foreach (int id in connectivity)
                {
                    elementNodes.Add(nodes[id]);
                    mesh.Vertices.Add(nodes[id].Coordinate);
                };

                Element element = new Element(i, elementNodes, connectivity);

                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.Normals.ComputeNormals();  // control if needed
                mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
                mesh.Compact();
                element.mesh = mesh;

                elements.Add(element); // add element to list of elements

                // create global mesh
                allMesh.Faces.AddFace(counter, counter + 1, counter + inputMesh.nu + 1, counter + inputMesh.nu);

                // element counter
                counter++;
                newRow++; ;
                if (newRow == (inputMesh.nu - 1)) //new row
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
            meshUpdated.Geometry = inputMesh.Geometry;
            
            // output
            DA.SetData(0, meshUpdated);
            DA.SetData(1, meshUpdated.mesh);
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
                return Properties.Resources.Icon_MoveSurfaceMeshVertices;
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
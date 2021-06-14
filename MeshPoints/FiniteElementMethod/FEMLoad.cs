﻿using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MeshPoints.Classes;


namespace MeshPoints.FiniteElementMethod
{
    public class FEMLoad : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMLoad class.
        /// </summary>
        public FEMLoad()
          : base("FEM Load", "Load",
              "Create load for the FEM Solver.",
              "SmartMesh", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {        
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh.", GH_ParamAccess.item); 
            pManager.AddIntegerParameter("Type", "type", "Load type: Point load = 1, Surface load = 2.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Position", "pos", "If type = 1: List of coordinates for point loads.", GH_ParamAccess.list); 
            pManager.AddIntegerParameter("Surface", "srf", "If type = 2: Face index of geometry to apply surface load.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Vector", "vec", "List of vectors of the loads. If surface load, only one vector.", GH_ParamAccess.list);

            pManager[0].Optional = true;
            pManager[1].Optional = true;
            pManager[2].Optional = true;
            pManager[3].Optional = true;

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
             pManager.AddGenericParameter("Load", "load", "List of residual forces (R).", GH_ParamAccess.list);
             pManager.AddGenericParameter("Points", "pts", "List of points subjected to load.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            // Input

            SmartMesh smartMesh = new SmartMesh(); 
            int loadType = 0;
            List<Vector3d> loadVectors = new List<Vector3d>();
            List<Point3d> loadPosition = new List<Point3d>();
            List<int> surfaceIndex = new List<int>();

            DA.GetData(0, ref smartMesh);
            DA.GetData(1, ref loadType);
            DA.GetDataList(2, loadPosition);
            DA.GetDataList(3, surfaceIndex);
            DA.GetDataList(4, loadVectors);

            // Code

            if (!DA.GetData(0, ref smartMesh)) return;
            List<Point3d> pointsWithLoad = new List<Point3d>();

            double[] residualForces = new double[smartMesh.Nodes.Count * 3];

            // Assign external load
            if (loadType == 1) // point load
            {           
                for (int i = 0; i < loadVectors.Count; i++)
                {
                    int nodeIndex = GetClosestNodeIndex(smartMesh.Nodes, loadPosition[i]);

                    // Deconstruct load vector
                    double xLoad = loadVectors[i].X;
                    double yLoad = loadVectors[i].Y;
                    double zLoad = loadVectors[i].Z;

                    // Construct residual force list
                    residualForces[nodeIndex * 3] = residualForces[nodeIndex * 3] + xLoad; 
                    residualForces[nodeIndex * 3 + 1] = residualForces[nodeIndex * 3 + 1] + yLoad; 
                    residualForces[nodeIndex * 3 + 2] = residualForces[nodeIndex * 3 + 2] + zLoad;
                    pointsWithLoad.Add(smartMesh.Nodes[nodeIndex].Coordinate);
                }
            }
            else if (loadType == 2) // surface load
            {
                for (int p = 0; p < loadVectors.Count; p++)
                {
                    BrepFace surface = smartMesh.Geometry.Faces[surfaceIndex[p]];
                    List<int> nodeIndexOnSurface = GetNodeIndexOnSurface(smartMesh.Nodes, surface);

                    // Prepare load lumping
                    foreach (Element element in smartMesh.Elements)
                    {
                        for (int i = 0; i < element.Connectivity.Count; i++)
                        {
                            for (int j = 0; j < nodeIndexOnSurface.Count; j++)
                            {
                                if (element.Connectivity[i] == nodeIndexOnSurface[j])
                                {
                                    List<List<Node>> faceList = element.GetFaces();

                                    for (int k = 0; k < faceList.Count; k++)
                                    {
                                        List<Node> face = faceList[k];
                                        int counter = 0;
                                        for (int n = 0; n < face.Count; n++)
                                        {
                                            if (nodeIndexOnSurface.Contains(face[n].GlobalId))
                                            {
                                                counter++;
                                            }
                                        }
                                        if (counter == 4)
                                        {
                                            Point3d n0 = face[0].Coordinate;
                                            Point3d n1 = face[1].Coordinate;
                                            Point3d n2 = face[2].Coordinate;
                                            Point3d n3 = face[3].Coordinate;

                                            double area1 = Math.Abs(0.5 * Vector3d.CrossProduct(n0 - n3, n1 - n3).Length);
                                            double area2 = Math.Abs(0.5 * Vector3d.CrossProduct(n1 - n3, n2 - n3).Length);

                                            // check area
                                            Vector3d normal = element.Mesh.FaceNormals[k];
                                            if (Vector3d.VectorAngle(n1 - n0, n3 - n0, normal) >= Math.PI)
                                            {
                                                area1 = Math.Abs(0.5 * Vector3d.CrossProduct(n1 - n0, n2 - n0).Length);
                                                area2 = Math.Abs(0.5 * Vector3d.CrossProduct(n3 - n0, n2 - n0).Length);
                                            }
                                            if (Vector3d.VectorAngle(n2 - n1, n3 - n1, normal) >= Math.PI)
                                            {
                                                area1 = Math.Abs(0.5 * Vector3d.CrossProduct(n2 - n1, n3 - n1).Length);
                                                area2 = Math.Abs(0.5 * Vector3d.CrossProduct(n0 - n1, n3 - n1).Length);
                                            }

                                            if (Vector3d.VectorAngle(n3 - n2, n1 - n2, normal) >= Math.PI)
                                            {
                                                area1 = Math.Abs(0.5 * Vector3d.CrossProduct(n3 - n2, n0 - n2).Length);
                                                area2 = Math.Abs(0.5 * Vector3d.CrossProduct(n1 - n2, n0 - n2).Length);
                                            }

                                            if (Vector3d.VectorAngle(n0 - n3, n2 - n3, normal) >= Math.PI)
                                            {
                                                area1 = Math.Abs(0.5 * Vector3d.CrossProduct(n0 - n3, n1 - n3).Length);
                                                area2 = Math.Abs(0.5 * Vector3d.CrossProduct(n2 - n3, n1 - n3).Length);
                                            }

                                            double faceArea = area1 + area2;

                                            foreach (Node node in face)
                                            {
                                                residualForces[3 * node.GlobalId + 0] = residualForces[node.GlobalId * 3 + 0] + loadVectors[p].X * faceArea / (double)4;
                                                residualForces[3 * node.GlobalId + 1] = residualForces[node.GlobalId * 3 + 1] + loadVectors[p].Y * faceArea / (double)4;
                                                residualForces[3 * node.GlobalId + 2] = residualForces[node.GlobalId * 3 + 2] + loadVectors[p].Z * faceArea / (double)4;

                                                if (!pointsWithLoad.Contains(node.Coordinate))
                                                {
                                                    pointsWithLoad.Add(node.Coordinate);
                                                }
                                            }

                                            i = element.Connectivity.Count; // break
                                            j = nodeIndexOnSurface.Count; // break
                                            k = faceList.Count; // break
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }

            // Output

            DA.SetDataList(0, residualForces);
            DA.SetDataList(1, pointsWithLoad);
        }

        #region Methods

        /// <summary>
        /// Get the index of the node closest to the load position.
        /// </summary>
        /// <returns> Index of node closest to load position.</returns>
        private int GetClosestNodeIndex(List<Node> nodes, Point3d loadPosition)
        {
            double minDistance = 10000000;
            int nodeIndex = 0;

            for (int i = 0; i < nodes.Count; i++)
            {
                Vector3d vec = nodes[i].Coordinate - loadPosition;
                double distance = vec.Length;
                if (distance < minDistance)
                {
                    minDistance = distance;
                    nodeIndex = i;
                }
                if (distance < 0.0001) { break; }
            }
            return nodeIndex;
        }

        /// <summary>
        /// Get the index of the nodes on the surface.
        /// </summary>
        /// <returns> List of index of nodes on the surface.</returns>
        private List<int> GetNodeIndexOnSurface(List<Node> nodes, BrepFace surface)
        {
            List<int> nodeIndexOnSurface = new List<int>();
            for (int i = 0; i < nodes.Count; i++)
            {
                if (nodes[i].IsOnFace(surface))
                {
                    nodeIndexOnSurface.Add(i);
                }
            }
            return nodeIndexOnSurface;
        }


        #endregion

            /// <summary>
            /// Provides an Icon for the component.
            /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Icon_Load;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("8972a393-9603-486a-bf29-a436a72d2c8d"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.FiniteElementMethod
{
    public class FEMLoad : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMLoad class.
        /// </summary>
        public FEMLoad()
          : base("FEM Load", "FEM load",
              "Create load for FEM solver.",
              "SmartMesh", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {        
            pManager.AddGenericParameter("SmartMesh", "SM", "Input a SmartMesh.", GH_ParamAccess.item); 
            pManager.AddIntegerParameter("Type", "type", "Load type: Point load = 1, Surface load = 2", GH_ParamAccess.item);
            pManager.AddGenericParameter("Position", "pos", "List of coordinates for point loads.", GH_ParamAccess.list); 
            pManager.AddIntegerParameter("Surface", "srf", "Input surface index of geometry to apply load to.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Vector", "vec", "List of vectors for the loads. If surface load, only one vector.", GH_ParamAccess.list);

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
             pManager.AddGenericParameter("Load", "load", "List of residual forces(R).", GH_ParamAccess.list);
             pManager.AddGenericParameter("Points", "pts", "List of node coordinates applied load to", GH_ParamAccess.list);
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
            int surfaceIndex = 0;

            DA.GetData(0, ref smartMesh);
            DA.GetData(1, ref loadType);
            DA.GetDataList(2, loadPosition);
            DA.GetData(3, ref surfaceIndex);
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
                // Assumption: only one load vector for surface load

                BrepFace surface = smartMesh.Geometry.Faces[surfaceIndex];
                List<int> nodeIndexOnSurface = GetNodeIndexOnSurface(smartMesh.Nodes, surface);

                // Prepare load lumping
                //Brep surfaceAsBerep = surface.ToBrep();
                //double area = Rhino.Geometry.AreaMassProperties.Compute(surfaceAsBerep).Area; // to do: slett?
                //int loadCounter = 0; // to do: slett?

                foreach (Element element in smartMesh.Elements)
                {
                    foreach (int connectivity in element.Connectivity)
                    {
                        foreach (int id in nodeIndexOnSurface)
                        {
                            if (connectivity == id)
                            {
                                List<List<Node>> faceList = element.GetFaces();
                                for (int i = 0; i < faceList.Count; i++)
                                {
                                    List<Node> face = faceList[i];
                                    for (int j = 0; j < face.Count; j++)
                                    {
                                        Node node = face[j];
                                        if (node.GlobalId == id)
                                        {
                                            Point3d A1 = face[0].Coordinate;
                                            Point3d B1 = face[1].Coordinate;
                                            Point3d C1 = face[3].Coordinate;

                                            Point3d A2 = face[1].Coordinate;
                                            Point3d B2 = face[2].Coordinate;
                                            Point3d C2 = face[3].Coordinate;

                                            // check area
                                            double area1 = Math.Abs(0.5 * (A1.X * (B1.Y - C1.Y) + B1.X * (C1.Y - A1.Y) + C1.X * (A1.Y - B1.Y)));
                                            double area2 = Math.Abs(0.5 * (A2.X * (B2.Y - C2.Y) + B2.X * (C2.Y - A2.Y) + C2.X * (A2.Y - B2.Y)));
                                            double faceArea = area1 + area2;

                                            residualForces[3 * id + 0] = residualForces[id * 3 + 0] + loadVectors[0].X * faceArea / (double)4;
                                            residualForces[3 * id + 1] = residualForces[id * 3 + 1] + loadVectors[0].Y * faceArea / (double)4;
                                            residualForces[3 * id + 2] = residualForces[id * 3 + 2] + loadVectors[0].Z * faceArea / (double)4;

                                            if (!pointsWithLoad.Contains(smartMesh.Nodes[id].Coordinate))
                                            {
                                                pointsWithLoad.Add(smartMesh.Nodes[id].Coordinate);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                /*
            List<List<Node>> faceList = element.GetFaces();
            for (int i = 0; i < faceList.Count; i++)
            {
                List<Node> face = faceList[i];
                for(int j = 0; j < face.Count; j++)
                {
                    Node node = face[j];
                    foreach (int id in nodeIndexOnSurface)
                    {
                        if (node.GlobalId == id)
                        {
                                Point3d A1 = face[0].Coordinate;
                                Point3d B1 = face[1].Coordinate;
                                Point3d C1 = face[3].Coordinate;

                                Point3d A2 = face[1].Coordinate;
                                Point3d B2 = face[2].Coordinate;
                                Point3d C2 = face[3].Coordinate;

                                // check area
                                double area1 = Math.Abs(0.5 * (A1.X * (B1.Y - C1.Y) + B1.X * (C1.Y - A1.Y) + C1.X * (A1.Y - B1.Y)));
                                double area2 = Math.Abs(0.5 * (A2.X * (B2.Y - C2.Y) + B2.X * (C2.Y - A2.Y) + C2.X * (A2.Y - B2.Y)));
                                double faceArea = area1 + area2;

                                residualForces[3 * id + 0] = residualForces[id * 3 + 0] + loadVectors[0].X * faceArea / (double)4;
                                residualForces[3 * id + 1] = residualForces[id * 3 + 1] + loadVectors[0].Y * faceArea / (double)4;
                                residualForces[3 * id + 2] = residualForces[id * 3 + 2] + loadVectors[0].Z * faceArea / (double)4;

                                if (!pointsWithLoad.Contains(smartMesh.Nodes[id].Coordinate))
                                {
                                    pointsWithLoad.Add(smartMesh.Nodes[id].Coordinate);
                                }

                            }
                    }
                }
            }
        }*/

                // to do: slett?
                /*
                foreach (int nodeIndex in nodeIndexOnSurface)
                {
                    if (smartMesh.Nodes[nodeIndex].Type == "Corner") { loadCounter++; } // corner node
                    else if (smartMesh.Nodes[nodeIndex].Type == "Edge") { loadCounter += 2; } // edge node
                    else { loadCounter += 4; } // face node
                }

                List<double> nodalLoad = new List<double>();
                nodalLoad.Add(loadVectors[0].X * area / (double) loadCounter);
                nodalLoad.Add(loadVectors[0].Y * area / (double) loadCounter); 
                nodalLoad.Add(loadVectors[0].Z * area / (double) loadCounter); 

                // Load lumping
                double addLoad = 0;
                foreach (int nodeIndex in nodeIndexOnSurface)
                {
                    for (int loadDir = 0; loadDir < 3; loadDir++)
                    {
                        if (smartMesh.Nodes[nodeIndex].Type == "Corner") { addLoad = nodalLoad[loadDir]; } // corner node
                        else if (smartMesh.Nodes[nodeIndex].Type == "Edge") { addLoad = nodalLoad[loadDir] * 2; } // edge node
                        else { addLoad = nodalLoad[loadDir] * 4; } // face node

                        residualForces[ 3 * nodeIndex + loadDir] = residualForces[nodeIndex * 3 + loadDir] + addLoad;
                    }
                    pointsWithLoad.Add(smartMesh.Nodes[nodeIndex].Coordinate);
                }
                */
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
                return null;
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
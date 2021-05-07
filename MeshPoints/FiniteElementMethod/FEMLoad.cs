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
            pManager.AddGenericParameter("SmartMesh", "SmartMesh", "Input a SmartMesh", GH_ParamAccess.item); 
            pManager.AddIntegerParameter("Load type", "load type", "Point load = 1, Surface load = 2", GH_ParamAccess.item);
            pManager.AddGenericParameter("Position", "pos", "Coordinate for point load", GH_ParamAccess.list); 
            pManager.AddIntegerParameter("Surface index", "Input surface index of geometry to apply load to", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("Load vectors", "vec", "List of vectors for the loads. If surface load, only one vector", GH_ParamAccess.list);

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
             pManager.AddGenericParameter("Load", "load", "List of residual (R)", GH_ParamAccess.list);
             pManager.AddGenericParameter("Nodes", "nodes", "List of node coordinates applied load to", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // assume only perpendicular negativ load

            #region Input
            SmartMesh mesh = new SmartMesh(); // to do: change to MeshGeometry elns
            int loadType = 0;
            List<Vector3d> loadVectors = new List<Vector3d>();
            List<Point3d> loadPosition = new List<Point3d>();
            int surfaceIndex = 0;

            DA.GetData(0, ref mesh);
            DA.GetData(1, ref loadType);
            DA.GetDataList(2, loadPosition);
            DA.GetData(3, ref surfaceIndex);
            DA.GetDataList(4, loadVectors);
            #endregion

            #region Code
            List<Node> nodes = mesh.Nodes;
            List<Element> element = mesh.Elements;
            int nodeDOFS = 3;
            List<Point3d> pointsWithLoad = new List<Point3d>();

            List<double> globalCoordinateLoadList = new List<double>();
            for (int i = 0; i < mesh.Nodes.Count * nodeDOFS; i++)
            {
                globalCoordinateLoadList.Add(0);
            }

            // Assign external load
            if (loadType == 1)
            {
                // Point load
                for (int i = 0; i < loadVectors.Count; i++)
                {
                    int nodeIndex = GetClosestNodeIndex(mesh.Nodes, loadPosition[i]);

                    // Deconstruct load to x,y,z
                    double xLoad = loadVectors[i].X;
                    double yLoad = loadVectors[i].Y;
                    double zLoad = loadVectors[i].Z;

                    globalCoordinateLoadList[nodeIndex * nodeDOFS] = globalCoordinateLoadList[nodeIndex * nodeDOFS] + xLoad; 
                    globalCoordinateLoadList[nodeIndex * nodeDOFS + 1] = globalCoordinateLoadList[nodeIndex * nodeDOFS + 1] + yLoad; 
                    globalCoordinateLoadList[nodeIndex * nodeDOFS + 2] = globalCoordinateLoadList[nodeIndex * nodeDOFS + 2] + zLoad;
                    pointsWithLoad.Add(nodes[nodeIndex].Coordinate);
                }

            }
            else if (loadType == 2)
            {
                // Surface load
                
                BrepFace surface = mesh.Geometry.Faces[surfaceIndex];
                List<int> nodeIndexOnSurface = GetNodeIndexOnSurface(mesh.Nodes, surface);

                Brep surfaceAsBerep = surface.ToBrep();
                double area = Rhino.Geometry.AreaMassProperties.Compute(surfaceAsBerep).Area;
                int loadCounter = 0;
                foreach (int nodeIndex in nodeIndexOnSurface)
                {
                    int numBC = 0; // to do: gjør om til methode i klasse
                    if (nodes[nodeIndex].BC_U){ numBC++; }
                    if (nodes[nodeIndex].BC_V) { numBC++; }
                    if (nodes[nodeIndex].BC_W) { numBC++; }
                    switch (numBC)
                    {
                        case 3: // corner node
                            loadCounter++; break;
                        case 2: // edge node
                            loadCounter += 2; break;
                        case 1: // midle node
                            loadCounter += 4; break;
                    }
                }

                List<double> nodalLoad = new List<double>();
                nodalLoad.Add(loadVectors[0].X * area / (double) loadCounter); // assume one load vector
                nodalLoad.Add(loadVectors[0].Y * area / (double) loadCounter); // assume one load vector
                nodalLoad.Add(loadVectors[0].Z * area / (double) loadCounter); // assume one load vector

                double addLoad = 0;
                foreach (int nodeIndex in nodeIndexOnSurface)
                {
                    for (int loadDir = 0; loadDir < 3; loadDir++)
                    {
                        int numBC = 0;
                        if (nodes[nodeIndex].BC_U) { numBC++; }
                        if (nodes[nodeIndex].BC_V) { numBC++; }
                        if (nodes[nodeIndex].BC_W) { numBC++; }
                        switch (numBC)
                        {
                            case 3: // corner node
                                addLoad = nodalLoad[loadDir]; break;
                            case 2: // edge node
                                addLoad = nodalLoad[loadDir] * 2; break;
                            case 1: // midle node
                                addLoad = nodalLoad[loadDir] * 4; ; break;
                        }
                        globalCoordinateLoadList[ nodeDOFS * nodeIndex + loadDir] = globalCoordinateLoadList[nodeIndex * nodeDOFS + loadDir] + addLoad;
                    }
                    pointsWithLoad.Add(nodes[nodeIndex].Coordinate);

                }

                /*
                List<int> elementIndexOnSurface = GetElementIndexConnectedToNodes(mesh.Elements, nodeIndexOnSurface);

                List<Node> elementNodesToLump = new List<Node>();

                foreach (int i in elementIndexOnSurface)
                {
                    List<Node> elementNodes = mesh.Elements[i].Nodes;
                    if (mesh.Type == "shell")
                    {
                        elementNodesToLump = elementNodes;
                    }
                    else
                    { 
                        // get nodes on side of element
                        switch(surfaceIndex)
                        {
                            case 0:
                                elementNodesToLump = new List<Node>() { elementNodes[0], elementNodes[1], elementNodes[2], elementNodes[3] };
                                break;
                            case 1:
                                elementNodesToLump = new List<Node>() { elementNodes[4], elementNodes[5], elementNodes[6], elementNodes[7] };
                                break;
                            case 2: 
                                elementNodesToLump = new List<Node>() { elementNodes[0], elementNodes[1], elementNodes[5], elementNodes[4] };
                                break;
                            case 3:
                                elementNodesToLump = new List<Node>() { elementNodes[1], elementNodes[2], elementNodes[6], elementNodes[5] };
                                break;
                            case 4:
                                elementNodesToLump = new List<Node>() { elementNodes[2], elementNodes[3], elementNodes[7], elementNodes[6] };
                                break;
                            case 5:
                                elementNodesToLump = new List<Node>() { elementNodes[3], elementNodes[0], elementNodes[4], elementNodes[7] };
                                break;
                        }

                    }
                    LoadLumping(elementNodesToLump, loadSize, R);
                }*/



            }
            #endregion

            DA.SetDataList(0, globalCoordinateLoadList) ;
            DA.SetDataList(1, pointsWithLoad);

        }

        #region Methods

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

        private double CalculateArea(List<Node> nodes)
        {
            Point3d A = nodes[0].Coordinate;
            Point3d B = nodes[1].Coordinate;
            Point3d C = nodes[2].Coordinate;
            Point3d D = nodes[3].Coordinate;
            double areaABC = 0.5 * (A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y));
            double areaCDA = 0.5 * (C.X * (D.Y - A.Y) + D.X * (A.Y - C.Y) + A.X * (C.Y - D.Y));
            return (areaABC + areaCDA);
        }
        private void LoadLumping(List<Node> elementNodes, double loadSize, List<double> R)
        {
           //double area = CalculateArea()
        }

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

        private List<int> GetElementIndexConnectedToNodes(List<Element> elements, List<int> nodeIndexToCheck)
        {
            List<int> elementIndexOnSurface = new List<int>();

            for (int i = 0; i < elements.Count; i++)
            {
                for (int j = 0; j < elements[0].Connectivity.Count; j++)
                {
                    if (nodeIndexToCheck.Contains(elements[i].Connectivity[j]))
                    {
                        elementIndexOnSurface.Add(i);
                        break;
                    }
                }
            }
            return elementIndexOnSurface;
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
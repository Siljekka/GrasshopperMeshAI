using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;

namespace MeshPoints
{
    public class FEMLoad : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMLoad class.
        /// </summary>
        public FEMLoad()
          : base("FEM Load", "FEM load",
              "Create load for FEM solver.",
              "MyPlugIn", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {        
            pManager.AddGenericParameter("MeshGeometry", "MeshFeometry", "Input a MeshGeometry", GH_ParamAccess.item); // to do: change name
            pManager.AddGenericParameter("Load type", "load type", "Point load = 1, Surface load = 2", GH_ParamAccess.item);
            pManager.AddGenericParameter("Position", "pos", "Coordinate for point load", GH_ParamAccess.item); 
            pManager.AddGenericParameter("Surface index", "Input surface index of geometry to apply load to", "", GH_ParamAccess.item);
            //pManager.AddGenericParameter("Load size", "size", "Load size (kN for point load, kN/mm^2 for surface)", GH_ParamAccess.item); // to do: se om dette skal være en liste ift input i FEMsolver
            pManager.AddGenericParameter("Load vectors", "vec", "List of vectors for the loads", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
             pManager.AddGenericParameter("Load", "load", "List of residual (R)", GH_ParamAccess.list);
        }

    /// <summary>
    /// This is the method that actually does the work.
    /// </summary>
    /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
    protected override void SolveInstance(IGH_DataAccess DA)
        {
            // assume only perpendicular negativ load

            #region Input
            Mesh3D mesh = new Mesh3D(); // to do: change to MeshGeometry elns
            int loadType = 0;
            List<Vector3d> loadVectors = new List<Vector3d>(); 
            Point3d loadPosition = new Point3d();
            int surfaceIndex = 0;

            DA.GetData(0, ref mesh);
            DA.GetData(1, ref loadType);
            DA.GetData(2, ref loadPosition);
            DA.GetData(3, ref surfaceIndex);
            DA.GetDataList(4, loadVectors);
            #endregion

            // to do: må gjøre noe med vektor til load... 
            #region Code
            List<Node> nodes = mesh.Nodes;
            List<Element> element = mesh.Elements;
            int nodeDOFS = 3;
            if (mesh.Type == "shell") { nodeDOFS = 2; }

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
                    int nodeIndex = GetClosestNodeIndex(mesh.Nodes, loadPosition);
                    // Deconstruct load to x,y,z
                    double xLoad = loadVectors[i].X;
                    double yLoad = loadVectors[i].Y;
                    double zLoad = loadVectors[i].Z;

                    globalCoordinateLoadList[nodeIndex * nodeDOFS] = xLoad; // to do: change if for shell
                    globalCoordinateLoadList[nodeIndex * nodeDOFS + 1] = yLoad; // to do: change if for shell
                    globalCoordinateLoadList[nodeIndex * nodeDOFS + 2] = zLoad; // to do: change if for shell
                }

            }
            else if (loadType == 2)
            {
                // Surface load
                /*
                BrepFace surface = mesh.GeometryInformation.Faces[surfaceIndex];
                List<int> nodeIndexOnSurface = GetNodeIndexOnSurface(mesh.Nodes, surface);
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
                }
  
                */

            }
            #endregion

            DA.SetDataList(0, globalCoordinateLoadList) ;
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

        private bool IsPointOnSurface(Point3d point, BrepFace face)
        {
            bool nodeIsOnGeometry = false;
            face.ClosestPoint(point, out double PointOnCurveU, out double PointOnCurveV);
            Point3d testPoint = face.PointAt(PointOnCurveU, PointOnCurveV);  // make test point 
            double distanceToFace = testPoint.DistanceTo(point); // calculate distance between testPoint and node
            if (distanceToFace <= 0.0001 & distanceToFace >= -0.0001) // if distance = 0: node is on edge
            {
                nodeIsOnGeometry = true;
            }
            return nodeIsOnGeometry;
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
                if (IsPointOnSurface(nodes[i].Coordinate, surface))
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
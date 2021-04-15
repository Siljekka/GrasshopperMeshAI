using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
namespace MeshPoints
{
    public class FEMBoundaryCondition : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMBoundaryCondition class.
        /// </summary>
        public FEMBoundaryCondition()
          : base("FEM Boundary condtion", "BC",
              "Create boundary consition for FEM solver.",
              "MyPlugIn", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("MeshGeometry", "MeshFeometry", "Input a MeshGeometry", GH_ParamAccess.item); // to do: change name
            pManager.AddGenericParameter("Edge/Face indices with BC", "", "", GH_ParamAccess.list); // to do: change name
            pManager.AddGenericParameter("Tx", "", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("Ty", "", "", GH_ParamAccess.item); 
            pManager.AddGenericParameter("Tw", "", "", GH_ParamAccess.item);
           // pManager.AddGenericParameter("Rx", "", "", GH_ParamAccess.item); 
            //pManager.AddGenericParameter("Rv", "", "", GH_ParamAccess.item); 
            //pManager.AddGenericParameter("Rz", "", "", GH_ParamAccess.item);

            pManager[2].Optional = true;
            pManager[3].Optional = true;
            pManager[4].Optional = true; 
           // pManager[5].Optional = true;
            //pManager[6].Optional = true;
            //pManager[7].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Boundary conditions", "BC", "List of DOFS that are fixed", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Input
            Mesh3D mesh = new Mesh3D(); // to do: change to MeshGeometry elns
            List<int> indicesOfGeometryWithBC = new List<int>();
            bool Tx = false;
            bool Ty = false;
            bool Tz = false;
            //bool Rx = false;
            //bool Ry = false;
            //bool Rz = false;

            DA.GetData(0, ref mesh);
            DA.GetDataList(1, indicesOfGeometryWithBC);
            DA.GetData(2, ref Tx);
            DA.GetData(3, ref Ty);
            DA.GetData(4, ref Tz);
            //DA.GetData(5, ref Rx);
            //DA.GetData(6, ref Ry);
            //DA.GetData(7, ref Rz);
            #endregion

            #region Code
            List<bool> applyBC = new List<bool>() { Tx, Ty, Tz }; // alternatively add rotation

            // Settings
            List<List<int>> applyBCToDOF = new List<List<int>>();
            bool nodeIsOnGeometry = false;
            int nodeDOFS = 2;
            if (String.Equals(mesh.Type, "solid")) { nodeDOFS = 3; }

            // Loop each dof for each node
            for (int i = 0; i < mesh.Nodes.Count; i++)
            {
                Node node = mesh.Nodes[i];

                // Check if node is on geometry to apply BC to
                if (String.Equals(mesh.Type, "shell"))
                {
                    for (int n = 0; n < indicesOfGeometryWithBC.Count; n++)
                    {
                        BrepEdge edge = mesh.Geometry.Edges[indicesOfGeometryWithBC[n]];
                        if (IsPointOnEdge(node.Coordinate, edge)) 
                        {
                            nodeIsOnGeometry = true;
                            break;
                        }
                    }
                }
                else
                {
                    for (int n = 0; n < indicesOfGeometryWithBC.Count; n++)
                    {
                        BrepFace face = mesh.Geometry.Faces[indicesOfGeometryWithBC[n]];
                        if (IsPointOnSurface(node.Coordinate, face))
                        {
                            nodeIsOnGeometry = true;
                            break;
                        }
                    }
                }

                // Loop each DOF of node. Add 0 is no BC and 1 if BC.
                List<int> BCNode = new List<int>();
                for (int j = 0; j < nodeDOFS; j++)
                {
                    if (nodeIsOnGeometry & applyBC[j] == true)
                    {
                        BCNode.Add(1);
                        continue;
                    }
                    BCNode.Add(0); // add information of a dof of the node
                }
                applyBCToDOF.Add(BCNode); // add information of the node
                nodeIsOnGeometry = false; // reset
            }

            #endregion


            #region Output
            DA.SetDataList(0, applyBCToDOF);
            #endregion
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

        private bool IsPointOnEdge(Point3d point, BrepEdge edge)
        {
            bool nodeIsOnGeometry = false;
            edge.ClosestPoint(point, out double PointOnCurve);
            Point3d testPoint = edge.PointAt(PointOnCurve);  // make test point 
            double distanceToEdge = testPoint.DistanceTo(point); // calculate distance between testPoint and node
            if (distanceToEdge <= 0.0001 & distanceToEdge >= -0.0001) // if distance = 0: node is on edge
            {
                nodeIsOnGeometry = true;
            }
            return nodeIsOnGeometry;
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
            get { return new Guid("7268f076-dbf3-4dbd-a363-41cd9966ad96"); }
        }
    }
}
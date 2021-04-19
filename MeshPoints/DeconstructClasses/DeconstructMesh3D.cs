using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.DeconstructClasses
{
    public class DeconstructMesh3D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructMesh3d class.
        /// </summary>
        public DeconstructMesh3D()
          : base("Deconstruct SolidMesh", "decSolid",
              "Deconstructing SolidMesh class",
              "MyPlugIn", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh3D", "m3D", "Mesh3D class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "e", "List of elements", GH_ParamAccess.list); //0
            pManager.AddGenericParameter("Nodes", "n", "List of nodes", GH_ParamAccess.list); //1
            pManager.AddGenericParameter("Mesh", "m", "Mesh", GH_ParamAccess.item); //2
            pManager.AddGenericParameter("Geometry", "geo", "Geometry information", GH_ParamAccess.item); //3
            pManager.AddGenericParameter("NormalizedCoordinates", "coord", "Coordinates of nodes, normalized", GH_ParamAccess.item); //4
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            Mesh3D m = new Mesh3D();
            DA.GetData(0, ref m);

            List<Node> nodes = m.Nodes;
            List<Point3d> points = new List<Point3d>();
            List<Point3d> newPoints = new List<Point3d>();
            double maxX = 0;
            double maxY = 0;
            double maxZ = 0;
            double minX = 0;
            double minY = 0;
            double minZ = 0;


            foreach (Node node in nodes) { points.Add(node.Coordinate); }

            for (int i = 0; i < points.Count; i++)
            {
                if (points[i].X > maxX) { maxX = points[i].X; }
                if (points[i].X < minX) { minX = points[i].X; }
                if (points[i].Y > maxY) { maxY = points[i].Y; }
                if (points[i].Y < minY) { minY = points[i].Y; }
                if (points[i].X > maxZ) { maxZ = points[i].Z; }
                if (points[i].X < minZ) { minZ = points[i].Z; }

            }
            
            for (int i = 0; i < nodes.Count; i++)
            {
                double lengthofdiagonal = Math.Sqrt((maxX - minX) * (maxX - minX) + (maxY - minY) * (maxY - minY) + (maxZ - minZ) * (maxZ - minZ));
                double normalizedX = (nodes[i].Coordinate.X - minX) / lengthofdiagonal;
                double normalizedY = (nodes[i].Coordinate.Y - minY) / lengthofdiagonal;
                double normalizedZ = (nodes[i].Coordinate.Z - minZ) / lengthofdiagonal;
                Point3d newPoint = new Point3d(normalizedX, normalizedY, normalizedZ);
                newPoints.Add(newPoint);
            }
           
            //output
            DA.SetDataList(0, m.Elements);
            DA.SetDataList(1, m.Nodes);
            DA.SetData(2, m.mesh);
            DA.SetData(3, m.Geometry);
            DA.SetData(4, newPoints);
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
                return Properties.Resources.Icon_DeconstructSolidMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("97c30c27-48c9-41ac-b09d-d02f80e806f6"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MeshPoints.CreateMesh
{
    public class CreateTriangleMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateDelauneyTriangleMesh class.
        /// </summary>
        public CreateTriangleMesh()
          : base("Triangle mesh", "TriMesh",
              "Creates a triangle mesh on a (2D) brep using built-in Delaunay method",
              "MyPlugin", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Brep surface", "b", "Insert brep of surface.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Max distance", "l", "Insert maximum distance between edge points.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Inner count", "c", "Insert amount of inner points.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Triangle mesh", "m", "Triangle mesh (Delauney)", GH_ParamAccess.item);
            pManager.AddGenericParameter("Bounding rectangle", "r", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("edgepoints", "e", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("edgepoints", "f", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Brep meshSurface = new Brep();
            double edgeNodeDistance = 0.0;
            double totalEdgeNodeCount = 0;
            DA.GetData(0, ref meshSurface);
            DA.GetData(1, ref edgeNodeDistance);
            DA.GetData(2, ref totalEdgeNodeCount);


            // 1. Generate a (BRep) rectangle containing the meshSurface
            // 2. Populate with points (based on sth)
            // 3. Remove points outside meshSurface
            Brep boundingRectangle = CreateBoundingRectangle(meshSurface);

            // 1. Populate meshSurface edge with points (based on input)
            List<Point3d> edgePointsSurface = CreateEdgePointsByCount(meshSurface, totalEdgeNodeCount);
            List<Point3d> edgePointsBoundingRectangle = CreateEdgePointsByCount(boundingRectangle, totalEdgeNodeCount);
            //List<Point3d> innerPoints = new List<Point3d>();
            

            

            // Output
            Mesh triangleMesh = new Mesh();
            DA.SetData(0, triangleMesh);
            DA.SetDataList(1, boundingRectangle.Edges);
            DA.SetDataList(2, edgePointsSurface);
            DA.SetDataList(3, edgePointsBoundingRectangle);

        }
        
        /// <summary>
        /// Creates a <see cref="Brep"/> rectangle around an arbitrary plane <see cref="Brep"/> geometry.
        /// </summary>
        private Brep CreateBoundingRectangle(Brep meshSurface)
        {
            var xList = new List<double>();
            var yList = new List<double>();
            // var zList = new List<double>();
            foreach (BrepVertex vertex in meshSurface.Vertices)
            {
                xList.Add(vertex.Location.X);
                yList.Add(vertex.Location.Y);
                // zList.Add(vertex.Location.Z);
            }
            var p1 = new Point3d(xList.Min(), yList.Min(), 0);
            var p2 = new Point3d(xList.Max(), yList.Min(), 0);
            var p3 = new Point3d(xList.Max(), yList.Max(), 0);
            var p4 = new Point3d(xList.Min(), yList.Max(), 0);

            Brep boundingRectangle = Brep.CreateFromCornerPoints(p1, p2, p3, p4, RhinoMath.DefaultDistanceToleranceMillimeters);

            return boundingRectangle;  
        }

        /// <summary>
        /// Creates points on the <see cref="BrepEdge"/>s of a <see cref="Brep"/> based on a given total edge node count.
        /// </summary>
        /// <returns>Returnds a list of <see cref="Point3d"/> along the edge of a <see cref="Brep"/>.</returns>
        private List<Point3d> CreateEdgePointsByCount(Brep meshSurface, double totalEdgeNodeCount)
        {
            List<Point3d> edgePoints = new List<Point3d>();
            double totalEdgeLength = 0;
            foreach (Curve edge in meshSurface.Edges)
            {
                totalEdgeLength += edge.GetLength();
            }
            foreach (Curve edge in meshSurface.Edges)
            {
                double[] tValues;
                double edgeLength = edge.GetLength();
                var edgeNodeCount = Convert.ToInt32(totalEdgeNodeCount * (edgeLength / totalEdgeLength) + 1);
                
                tValues = edge.DivideByCount(edgeNodeCount, true);
                foreach (double t in tValues)
                {
                    edgePoints.Add(edge.PointAt(t));
                }
            }
            return edgePoints;
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
            get { return new Guid("a07a01a6-a771-4d75-9f96-e87ece274885"); }
        }
    }
}
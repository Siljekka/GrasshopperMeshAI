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
        /// Initializes a new instance of the CreateTriangleMesh class.
        /// </summary>
        public CreateTriangleMesh()
          : base("Triangle Mesh", "TriMesh",
              "Creates a triangle mesh on a (2D) Brep surface using built-in Delaunay method",
              "MyPlugin", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Brep surface", "b", "Insert brep of surface.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Edge Node Count", "c", "Insert wanted amount of edge nodes.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Inner Node Count", "n", "Insert wanted amount of inner nodes.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Triangle Mesh", "m", "Triangle mesh (Delauney)", GH_ParamAccess.item);
            //pManager.AddGenericParameter("Edge Nodes", "e", "Edge nodes of triangle mesh", GH_ParamAccess.list);
            //pManager.AddGenericParameter("Inner Nodes", "i", "Inner nodes of triangle mesh", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Inputs
            Brep meshSurface = new Brep();
            double totalEdgeNodeCount = 0;
            double totalInnerNodeCount = 0;
            DA.GetData(0, ref meshSurface);
            DA.GetData(1, ref totalEdgeNodeCount);
            DA.GetData(2, ref totalInnerNodeCount);

            #region Main
            // 1. Create nodes along the edges of the surface.
            List<List<Point3d>> edgeNodesSurface = CreateEdgePointsByCount(meshSurface, totalEdgeNodeCount);

            // 2. Create points inside the surface by creating a bounding box, populating it with points, and culling all points not inside the surface.
            Brep boundingBoxSurface = CreateBoundingBoxBrep(meshSurface);
            List<Point3d> nodeGridBoundingBox = CreatePointGridRectangle(boundingBoxSurface, meshSurface, totalInnerNodeCount);
            List<Point3d> nodesInsideSurface = CullPointsOutsideSurface(nodeGridBoundingBox, meshSurface);

            // 3. Flatten list of points to use for triangle meshing and cast to compatible data structure (Node2List) for Delaunay method.
            List<Point3d> nodeCollection = new List<Point3d>();
            List<Point3d> flattenedEdgeNodes = new List<Point3d>(); // collect separately for output
            foreach(var list in edgeNodesSurface)
            {
                flattenedEdgeNodes.AddRange(list);
            }
            nodeCollection.AddRange(flattenedEdgeNodes);
            nodeCollection.AddRange(nodesInsideSurface);
            var meshNodes = new Grasshopper.Kernel.Geometry.Node2List(nodeCollection);
            
            // 4. Throw all our points into the Delaunay mesher. Adjust jitter_amount as needed.
            var meshFaces = new List<Grasshopper.Kernel.Geometry.Delaunay.Face>();
            var triangleMesh = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Mesh(meshNodes, 0.01, ref meshFaces); // todo: what is "double jitter_amount"?

            // 5. Sometimes the mesh acts up; in these cases it is necessary to cull mesh faces that are outside the surface.
            Mesh culledTriangleMesh = CullMeshFacesOutsideSurface(triangleMesh, meshSurface);
            #endregion

            // Outputs
            DA.SetData(0, culledTriangleMesh);
            //DA.SetDataList(1, flattenedEdgeNodes);
            //DA.SetDataList(2, nodesInsideSurface);

        }
        #region Methods
        /// <summary>
        /// Cull unwanted mesh faces by checking if their center points are outside the actual surface of the mesh.
        /// </summary>
        /// <returns>A <see cref="Mesh"/> with (hopefully) no outside mesh faces.</returns>
        private Mesh CullMeshFacesOutsideSurface(Mesh meshSurface, Brep brep)
        {
            Mesh insideFaces = meshSurface.DuplicateMesh();
            for (int i = meshSurface.Faces.Count-1; i>0; i--) // reverse iteration to maintain indices
            {
                if (!IsPointOnBrepSurface(meshSurface.Faces.GetFaceCenter(i), brep))
                {
                    insideFaces.Faces.RemoveAt(i);
                }
            }
            return insideFaces;
        }

        /// <summary>
        /// Takes an input point and a Brep surface. If the distance between input point 
        /// and the closest point on the Brep ~ 0, the point is deemed on the surface.
        /// </summary>
        /// <returns>True if point is on Brep.</returns>
        private bool IsPointOnBrepSurface(Point3d point, Brep brep)
        {
            var testPointSurfaceDistance = point.DistanceTo(brep.ClosestPoint(point));

            if (testPointSurfaceDistance < RhinoMath.SqrtEpsilon) 
            { 
                return true; 
            }
            else 
            { 
                return false; 
            }
        }

        /// <summary>
        /// Takes a list of <see cref="Point3d"/> and checks if points are inside an input <see cref="Brep"/> surface.
        /// </summary>
        /// <returns>A list of <see cref="Point3d"/> containing points inside the input <see cref="Brep"/> surface.</returns>
        private List<Point3d> CullPointsOutsideSurface(List<Point3d> pointGrid, Brep meshSurface)
        {
            var insidePoints = new List<Point3d>();

            foreach (Point3d point in pointGrid)
            {
                if (IsPointOnBrepSurface(point, meshSurface))
                {
                    insidePoints.Add(point);
                }
            }
            return insidePoints;
        }

        /// <summary>
        /// Creates a <see cref="Brep"/> rectangle around an arbitrary plane <see cref="Brep"/> geometry and fills it with points.
        /// </summary>
        private Brep CreateBoundingBoxBrep(Brep meshSurface)
        {
            BoundingBox boundingBox = meshSurface.GetBoundingBox(false);
            Brep boundingBoxBrep = Brep.CreateFromCornerPoints(
                new Point3d(boundingBox.Min.X, boundingBox.Min.Y, 0),
                new Point3d(boundingBox.Max.X, boundingBox.Min.Y, 0),
                new Point3d(boundingBox.Max.X, boundingBox.Max.Y, 0),
                new Point3d(boundingBox.Min.X, boundingBox.Max.Y, 0),
                RhinoMath.ZeroTolerance
                );

            return boundingBoxBrep;
        }

        /// <summary>
        /// Creates points on the <see cref="BrepEdge"/>s of a <see cref="Brep"/> based on a given total edge node count.
        /// </summary>
        /// <returns>Returnds a list of <see cref="Point3d"/> along the edge of a <see cref="Brep"/>.</returns>
        private List<List<Point3d>> CreateEdgePointsByCount(Brep meshSurface, double totalEdgeNodeCount)
        {
            var edgePoints = new List<List<Point3d>>();
            double totalEdgeLength = 0;
            foreach (Curve edge in meshSurface.Edges)
            {
                totalEdgeLength += edge.GetLength();
            }
            foreach (Curve edge in meshSurface.Edges)
            {
                double[] tValues;
                var innerEdgePoints = new List<Point3d>();
                double edgeLength = edge.GetLength();
                var edgeNodeCount = Convert.ToInt32(totalEdgeNodeCount * (edgeLength / totalEdgeLength) + 1);
                
                tValues = edge.DivideByCount(edgeNodeCount, true);
                foreach (double t in tValues)
                {
                    innerEdgePoints.Add(edge.PointAt(t));
                }
                edgePoints.Add(innerEdgePoints);
            }
            return edgePoints;
        }

        /// <summary>
        /// Creates a grid of points in a bounding box <see cref="Brep"/> based on a given number of wanted internal nodes in a surface.
        /// </summary>
        /// <returns>A List of <see cref="Point3d"/> describing a grid of points in a rectangle.</returns>
        private List<Point3d> CreatePointGridRectangle(Brep boundingBox, Brep meshSurface, double nodeCount)
        {
            var gridPoints = new List<Point3d>(); // output

            // boundingBoxEdgeNodeCount is crudely implemented and is only accurate at a high number of nodes && a quadratic bounding box.
            // It tries to calculate how many evenly spaced nodes we need on the edge of the bounding box to 
            // achieve the wanted amount of inner nodes. Send an e-mail to magnus@kunnas.no for explanation.
            var boundingBoxEdgeNodeCount = Math.Sqrt(nodeCount * boundingBox.GetArea() / meshSurface.GetArea()) * 4 - 4;

            List<List<Point3d>> edgeGrid = CreateEdgePointsByCount(boundingBox, boundingBoxEdgeNodeCount); 

            var edge1 = edgeGrid[0];
            var edge2 = edgeGrid[1];
            var edge3 = edgeGrid[2];
            // var edge4 = edgeGrid[3]; // not used
            edge3.Reverse();

            var pointsInUdirection = edge1.Count();
            var pointsInVdirection = edge2.Count();

            for(int i = 0; i<pointsInUdirection; i++)
            {
                double[] tValues;
                var line = new LineCurve(edge1[i], edge3[i]); // draw lines between points on opposite edges
                tValues = line.DivideByCount(pointsInVdirection - 1, true); // # of divisions is one less than # of nodes along edge
                foreach (double t in tValues)
                {
                    gridPoints.Add(line.PointAt(t));
                }
            }
            return gridPoints;
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
            get { return new Guid("a07a01a6-a771-4d75-9f96-e87ece274885"); }
        }
    }
}
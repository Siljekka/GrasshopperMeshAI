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
            pManager.AddMeshParameter("prunedTriangle mesh", "m", "Triangle mesh (Delauney) pruned", GH_ParamAccess.item);
            pManager.AddGenericParameter("Bounding rectangle", "r", "", GH_ParamAccess.list);
            pManager.AddGenericParameter("meshpoints", "e", "", GH_ParamAccess.list);
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
            // todo: refactor 
            Brep boundingRectangle = CreateBoundingRectangle(meshSurface);
            BoundingBox boundingBox = meshSurface.GetBoundingBox(false);
            Brep boundingBoxBrep = Brep.CreateFromCornerPoints(
                new Point3d(boundingBox.Min.X, boundingBox.Min.Y, 0),
                new Point3d(boundingBox.Max.X, boundingBox.Min.Y, 0),
                new Point3d(boundingBox.Max.X, boundingBox.Max.Y, 0),
                new Point3d(boundingBox.Min.X, boundingBox.Max.Y, 0),
                RhinoMath.ZeroTolerance
                );


            var edgePointsOfBoundingRectangle = CreateEdgePointsByCount(boundingBoxBrep, totalEdgeNodeCount);
            var pointGridBoundingRectangle = CreatePointGridRectangle(edgePointsOfBoundingRectangle);
            var gridPointsInsideMeshSurface = PrunePointsOutsideSurface(pointGridBoundingRectangle, meshSurface);

            // 1. Populate meshSurface edge with points (based on input)
            List<List<Point3d>> edgePointsSurface = CreateEdgePointsByCount(meshSurface, totalEdgeNodeCount);

            //List<Point3d> innerPoints = new List<Point3d>();
            var pointCollection = new List<Point3d>();
            foreach(var list in edgePointsSurface)
            {
                pointCollection.AddRange(list);
            }
            pointCollection.AddRange(gridPointsInsideMeshSurface);

            var meshPoints = new Grasshopper.Kernel.Geometry.Node2List(pointCollection);
            
            var meshFaces = new List<Grasshopper.Kernel.Geometry.Delaunay.Face>();
            var triangleMesh = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Mesh(meshPoints, 0.01, ref meshFaces);

            Mesh prunedTriangleMesh = PruneMeshFacesOutsideSurface(triangleMesh, meshSurface);

            // Output
            DA.SetData(0, triangleMesh);
            DA.SetData(1, prunedTriangleMesh);
            DA.SetDataList(2, gridPointsInsideMeshSurface);
            DA.SetDataList(3, pointCollection);

        }
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

        private Mesh PruneMeshFacesOutsideSurface(Mesh mesh, Brep brep)
        {
            Mesh insideFaces = mesh;
            for (int i = mesh.Faces.Count-1; i>0; i--) // reverse iteration to maintain indices
            {
                if (!IsPointOnBrepSurface(mesh.Faces.GetFaceCenter(i), brep))
                {
                    insideFaces.Faces.RemoveAt(i);
                }
            }
            return insideFaces;
        }

        private List<Point3d> PrunePointsOutsideSurface(List<Point3d> pointGrid, Brep meshSurface)
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

            

            //// Based on Solution 4 (3D) from:
            //// http://www.eecs.umich.edu/courses/eecs380/HANDOUTS/PROJ2/InsidePoly.html
            //bool IsPointInside(Point3d testPoint, List<Point3d> boundingPoints)
            //{
            //    double totalAngle = 0;
            //    for (int i = 0; i<boundingPoints.Count-1; i++)
            //    {
            //        var nextIndex = i + 1;
            //        var vec1 = new Vector3d(
            //            boundingPoints[i].X - testPoint.X,
            //            boundingPoints[i].Y - testPoint.Y,
            //            boundingPoints[i].Z - testPoint.Z
            //            );
            //        var vec2 = new Vector3d(
            //            boundingPoints[nextIndex].X - testPoint.X,
            //            boundingPoints[nextIndex].Y - testPoint.Y,
            //            boundingPoints[nextIndex].Z - testPoint.Z
            //            );
            //        totalAngle += Vector3d.VectorAngle(vec1, vec2);
            //    }

            //    var twoPi = RhinoMath.TwoPI;
            //    var e = RhinoMath.SqrtEpsilon;
            //    if (twoPi < totalAngle+e && totalAngle-e < twoPi)
            //    {
            //        return true;
            //    } 
            //    else
            //    {
            //        return false;
            //    }
            //}
        }

        /// <summary>
        /// Creates a <see cref="Brep"/> rectangle around an arbitrary plane <see cref="Brep"/> geometry and fills it with points.
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
            // List<Point3d> cornerPoints = new List<Point3d>() { p1, p2, p3, p4 };

            return boundingRectangle;  
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

        private List<Point3d> CreatePointGridRectangle(List<List<Point3d>> edgeGrid)
        {
            var edge1 = edgeGrid[0];
            //var edge2 = edgeGrid[1];
            var edge3 = edgeGrid[2];
            //var edge4 = edgeGrid[3];
            edge3.Reverse();
            //edge4.Reverse();

            var gridPoints = new List<Point3d>();

            var pointsInUDirection = edgeGrid[0].Count();
            var pointsInVDirection = edgeGrid[1].Count();
            for(int i = 0; i<pointsInUDirection; i++)
            {
                double[] tValues;
                var line = new LineCurve(edge1[i], edge3[i]);
                tValues = line.DivideByCount(pointsInVDirection - 1, true);
                foreach (double t in tValues)
                {
                    gridPoints.Add(line.PointAt(t));
                }
            }
            return gridPoints;
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
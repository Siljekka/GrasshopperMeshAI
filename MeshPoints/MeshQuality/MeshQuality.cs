using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino;
using System;
using System.Linq;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Drawing;

namespace MeshPoints
{
    public class MeshQuality : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Mesh_Quality class.
        /// </summary>
        public MeshQuality()
          : base("Mesh Quality", "mq",
              "Mesh Quality",
              "MyPlugIn", "Quality")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh2D", "m", "Insert Mesh2D class", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Quality metric", "q", "Aspect Ratio = 1, Skewness = 2, Jacobian = 3", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality", "mq", "Mesh Quality for elements", GH_ParamAccess.list);
            pManager.AddGenericParameter("Avg. Aspect Ratio", "ar", "Average aspect ratio of all elements.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Avg. Skewness", "sk", "Average skewness of all elements.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Avg. Jacobian", "jb", "Average Jacobian ratio of all elements", GH_ParamAccess.item);
            pManager.AddGenericParameter("Color Mesh", "cm", "Color map of quality check", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Variables
            //variables
            Mesh2D mesh2D = new Mesh2D();
            Mesh colorMesh = new Mesh();

            List<Quality> qualityList = new List<Quality>(); // list of Quality for each element in the mesh
            List<double> vertexDistance = new List<double>(); //list distances between vertices in a mesh face, following mesh edges CCW (counter-clockwise)
            List<double> elementAngles = new List<double>(); //list of angles in a element

            // determines which quality check type to color mesh with
            // 1 = aspect ratio, 2 = skewness, 3 = jacobian
            int qualityCheckType = 0; 

            double idealAngle = 90; //ideal angle in degrees
            int neighborPoint = 3; //variable used in skewness calculation

            double sumAspectRatio = 0;
            double sumSkewness = 0;
            double sumJacobianRatio = 0;
            #endregion

            //input
            DA.GetData(0, ref mesh2D);
            DA.GetData(1, ref qualityCheckType);

            #region Calculate mesh quality
            foreach (Element e in mesh2D.Elements)
            {
                Quality elementQuality = new Quality();

                List < Point3d > pts = new List<Point3d>()
                { 
                        e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                        e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                };

                for (int n = 0; n < pts.Count / 2; n++)
                {   
                    //Aspect Ratio
                    vertexDistance.Add(pts[n].DistanceTo(pts[n + 1])); //Add the distance between the points, following mesh edges CCW

                    //Skewness
                    //create a vector from a vertex to a neighbouring vertex
                    Vector3d a = new Vector3d(pts[n].X - pts[n + 1].X, pts[n].Y - pts[n + 1].Y, pts[n].Z - pts[n + 1].Z); 
                    Vector3d b = new Vector3d(pts[n].X - pts[n + neighborPoint].X, pts[n].Y - pts[n + neighborPoint].Y, pts[n].Z - pts[n + neighborPoint].Z);

                    //calculate angles in radians between vectors
                    double angleRad = Math.Abs(Math.Acos(Vector3d.Multiply(a, b) / (a.Length * b.Length))); 
                    elementAngles.Add(angleRad * 180 / Math.PI); //convert from rad to deg
                }
                
                vertexDistance.Sort();
                elementAngles.Sort();

                elementQuality.AspectRatio = (vertexDistance[0] / vertexDistance[vertexDistance.Count - 1]);
                elementQuality.Skewness = 1 - Math.Max((elementAngles[elementAngles.Count - 1] - idealAngle) / (180 - idealAngle), (idealAngle - elementAngles[0]) / (idealAngle));
                elementQuality.JacobianRatio = CalculateJacobianRatioOfPlaneQuadElement(e);
                
                elementQuality.element = e;
                e.MeshQuality = elementQuality;

                sumAspectRatio += elementQuality.AspectRatio;
                sumSkewness += elementQuality.Skewness;
                sumJacobianRatio += elementQuality.JacobianRatio;

                qualityList.Add(elementQuality);
        
                vertexDistance.Clear();
                elementAngles.Clear();
            }

            double avgAspectRatio = sumAspectRatio / mesh2D.Elements.Count;
            double avgSkewness = sumSkewness / mesh2D.Elements.Count;
            double avgJacobianRatio = sumJacobianRatio / mesh2D.Elements.Count;
            #endregion

            #region Color the mesh based on quality type
            // 1 = aspect ratio
            if (qualityCheckType == 1) 
            {
                foreach (Quality q in qualityList)
                {
                    if (q.AspectRatio > 0.9)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                    }
                    else if (q.AspectRatio > 0.7)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                    }
                    else if (q.AspectRatio > 0.6)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                    }
                    else if (q.AspectRatio > 0)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                    }
                    colorMesh.Append(q.element.mesh);
                }
            }
            // 2 = skewness
            else if (qualityCheckType == 2)
            {
                foreach (Quality q in qualityList)
                {
                    if (q.Skewness > 0.9)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                    }
                    else if (q.Skewness > 0.7)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                    }
                    else if (q.Skewness > 0.6)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                    }
                    else if (q.Skewness > 0)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                    }
                    colorMesh.Append(q.element.mesh);
                }
            }
            // Todo implement mesh coloring for jacobian ratio
            #endregion

            #region Outputs
            DA.SetDataList(0, qualityList);
            DA.SetData(1, avgAspectRatio);
            DA.SetData(2, avgSkewness);
            DA.SetData(3, avgJacobianRatio);
            DA.SetData(4, colorMesh);
            #endregion
        }

        #region Component methods

        List<Point3d> CalculateLocalCoordinatesOfPlaneElement(Element meshFace)
        {
            /*
            1. Calculate surface normal.
            2. Calculate transformation matrix between surface normal and unit z-vector [0,0,1].
            3. Extract global corner points of meshFace.
            4. Transform points to local coordinate system (x', y', z'=0)
            output: list of Point3d
             */

            var elementPoints = new List<Point3d>(){
                meshFace.Node1.Coordinate, meshFace.Node2.Coordinate, meshFace.Node3.Coordinate, meshFace.Node4.Coordinate
            };

            var elementNormal = meshFace.mesh.NormalAt(0, 1, 0, 0, 0); // // for planar element normal is the same in every point
            var elementPlane = new Plane(meshFace.Node1.Coordinate, elementNormal);
            var xyPlane = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));

            var elementTransform = Transform.PlaneToPlane(elementPlane, xyPlane); // generates a transformation matrix for mapping points to xy-plane

            var transformedPoints = new List<Point3d>();
            foreach (Point3d point in elementPoints)
            {
                transformedPoints.Add(elementTransform * point);
            }

            return transformedPoints;

            // todo remove old code
            //var elementPoints = new List<Point3d>(){
            //    meshFace.Node1.Coordinate, meshFace.Node2.Coordinate, meshFace.Node3.Coordinate, meshFace.Node4.Coordinate
            //};
            //// Check if points are co-planar to determine if we want to continue
            //if (!Point3d.ArePointsCoplanar(elementPoints, RhinoMath.ZeroTolerance))
            //{
            //    throw new ArgumentException("Points of element are not co-planar.");
            //}

            //var meshBrep = Brep.CreateFromCornerPoints(
            //    elementPoints[0],
            //    elementPoints[1],
            //    elementPoints[2],
            //    elementPoints[3],
            //    0
            //    );

            //if (!meshBrep.IsSurface)
            //{
            //    throw new ArgumentNullException("Brep does not contain a surface.");
            //}
            //var meshBrepSurface = meshBrep.Surfaces[0]; // this should only have one element
            //var surfaceNormal = meshBrepSurface.NormalAt(0, 0); // for planar element normal is the same in every point

        }
        double CalculateJacobianRatioOfPlaneQuadElement(Element meshFace)
        {
            /*Alternate method:
            1. Define surface from corner points. It already is
            2. Find normal.
            3. Transform to plane element.
            4. Use local coordinates to calculate Jacobian.
             */
            var localPoints = CalculateLocalCoordinatesOfPlaneElement(meshFace);

            var gX = new List<Double>()
            {
                localPoints[0].X, localPoints[1].X, localPoints[2].X, localPoints[3].X,
            };
            var gY = new List<Double>()
            {
                localPoints[0].Y, localPoints[1].Y, localPoints[2].Y, localPoints[3].Y,
            };
            /* 
            1. Collect corner points.
            2. Define corner of natural quad element.
            3. Calculate the Jacobian determinant of each corner.
            4. The Jacobian ratio for mesh quality is the ratio of the minimum and maximum jacobian determinants of the element.
               A value of 1 is a perfect rectangular mesh.
             */

            // NodeN.Coordinate is a reference to the inherited Point3d object of the Node-class
            // Todo refactor Node-class, why is Point3d named Coordinate???
            //var gX = new List<Double>() // global x-coordinates of corner points
            //{
            //    meshFace.Node1.Coordinate.X, meshFace.Node2.Coordinate.X, meshFace.Node3.Coordinate.X, meshFace.Node4.Coordinate.X,
            //};
            //var gY = new List<Double>() // global y-coordinates of corner points
            //{
            //    meshFace.Node1.Coordinate.Y, meshFace.Node2.Coordinate.Y, meshFace.Node3.Coordinate.Y, meshFace.Node4.Coordinate.Y,
            //};

            var naturalCornerPoints = new List<List<Double>>
            {
                new List<double>{ -1, -1}, new List<double> { 1, -1 }, new List<double> { 1, 1 }, new List<double> { -1, 1 }
            };

            var jacobiansOfElement = new List<Double>();

            // Calculate the Jacobian determinant of each corner point
            for (int n=0; n<naturalCornerPoints.Count; n++)
            {
                double nX = naturalCornerPoints[n][0];
                double nY = naturalCornerPoints[n][1];

                // See documentation for derivation of formula
                double jacobianOfCorner = 0.0625 *
                    (
                    ((1 - nY) * (gX[1] - gX[0]) + (1  +nY) * (gX[2] - gX[3]))*
                    ((1 - nX) * (gY[3] - gY[0]) + (1 + nX) * (gY[2] - gY[1]))
                    -
                    ((1 - nY) * (gY[1] - gY[0]) + (1 + nY) * (gY[2] - gY[3])) *
                    ((1 - nX) * (gX[3] - gX[0]) + (1 + nX) * (gX[2] - gX[1]))
                    );

                jacobiansOfElement.Add(jacobianOfCorner);
            };

            // Minimum element divided by maximum element. A value of 1 is a perfect square.
            double jacobianRatio = jacobiansOfElement.Min() / jacobiansOfElement.Max();

            return jacobianRatio;
        }
        #endregion

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.MeshQuality;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("31a76520-576a-4b37-a64f-8d51178e4be7"); }
        }
    }
}
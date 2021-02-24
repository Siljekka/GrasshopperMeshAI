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
            
            pManager[1].Optional = true; // coloring the mesh is optional
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

            // Determines which quality check type to color mesh with
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

            // Color the mesh based on quality type
            ColorMesh(colorMesh, qualityList, qualityCheckType);

            #region Outputs
            DA.SetDataList(0, qualityList);
            DA.SetData(1, avgAspectRatio);
            DA.SetData(2, avgSkewness);
            DA.SetData(3, avgJacobianRatio);
            DA.SetData(4, colorMesh);
            #endregion
        }

        #region Component methods
        
        /// <summary>
        /// Transforms the corner points of an arbitrary 3D plane quad surface to a 2D plane.
        /// </summary>
        /// <param name="meshFace">An <see cref="Element"/> object describing a mesh face; see <see cref="Element"/> class for attributes.</param>
        /// <returns>A list of <see cref="Point3d"/> where the Z (third) coordinate is 0.</returns>
        List<Point3d> TransformPlaneQuadSurfaceTo2DPoints(Element meshFace)
        {
            /*
             * 1. Calculate surface normal.
             * 2. Define surface plane using surface normal and first point of element.
             * 3. Define xy-plane by origin: (0, 0, 0) and unit Z-vector: [0, 0, 1]
             * 4. Define transformation from surface plane to xy-plane using PlaneToPlane().
             * 5. Transform points to local coordinate system (X', Y', Z'=0)
             * output: list of Point3d where Z (third) coordinate = 0
             */

            var transformedPoints = new List<Point3d>(); // output
            var elementPoints = new List<Point3d>(){
                meshFace.Node1.Coordinate, meshFace.Node2.Coordinate, meshFace.Node3.Coordinate, meshFace.Node4.Coordinate
            };

            // Check that surface is planar
            if (!Point3d.ArePointsCoplanar(elementPoints, RhinoMath.ZeroTolerance))
            {
                throw new ArgumentException(
                    message: "Corner points of input surface are not co-planar.",
                    paramName: "elementPoints");
            }

            var elementVectors = new List<Vector3d>
            {
                new Vector3d(elementPoints[1].X - elementPoints[0].X, elementPoints[1].Y - elementPoints[0].Y, elementPoints[1].Z - elementPoints[0].Z),
                new Vector3d(elementPoints[3].X - elementPoints[0].X, elementPoints[3].Y - elementPoints[0].Y, elementPoints[3].Z - elementPoints[0].Z)
            };
            
            // Cross product of two linearily independent vectors is the normal to the plane containing them
            var surfaceNormal = Vector3d.CrossProduct(elementVectors[0], elementVectors[1]);

            var surfacePlane = new Plane(meshFace.Node1.Coordinate, surfaceNormal);
            var xyPlane = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));

            var elementTransformation = Transform.PlaneToPlane(surfacePlane, xyPlane); 

            foreach (Point3d point in elementPoints)
            {
                transformedPoints.Add(elementTransformation * point);
            }

            return transformedPoints;
        }

        /// <summary>
        /// Calculates the Jacobian ratio of a plane and simple quadrilateral mesh element.
        /// </summary>
        /// <param name="meshFace">An <see cref="Element"/> object describing a mesh face; see <see cref="Element"/> class for attributes.</param>
        /// <returns>A <see cref="double"/> between 0.0 and 1.0. A negative Jacobian might indicate a self-intersecting element.</returns>
        double CalculateJacobianRatioOfPlaneQuadElement(Element meshFace)
        {
            /*
             * This method utilizes the idea of shape functions and natural coordinate system to calculate the Jacobian
             * of given points on a plane, simple quadrilateral element.
             * 
             * 1. Transform the input 3D plane element (and specifically the corner points) to a 2D space (X', Y', Z'=0).
             * 2. Define natural coordinates we want to calculate the Jacobian in. This could be the corner points (or 
             *    alternatively the Gauss points) of the quad element. 
             * 3. Calculate the Jacobian of each point. 
             * 4. The ratio is the ratio of the minimum and maximum Jacobian calculated, given as a double from 0.0 to 1.0.
             *    A negative Jacobian indicates a self-intersecting element and should not happen.
                */

            List<Point3d> localPoints = TransformPlaneQuadSurfaceTo2DPoints(meshFace);

            var gX = new List<Double>()
            {
                localPoints[0].X, localPoints[1].X, localPoints[2].X, localPoints[3].X,
            };
            var gY = new List<Double>()
            {
                localPoints[0].Y, localPoints[1].Y, localPoints[2].Y, localPoints[3].Y,
            };

            var naturalPoints = new List<List<Double>> // natural coordinates of corner points
            {
                new List<double>{ -1, -1}, new List<double> { 1, -1 }, new List<double> { 1, 1 }, new List<double> { -1, 1 }
            };

            #region Todo: Implement which points we want to evaluate the jacobians for (corner vs 4 gauss integration points)
            //double s = 0.57735; // this represents the Gauss point of an isoparametric quadrilateral element: sqrt(1/3)
            //var naturalGaussPoints = new List<List<Double>> // natural coordinates of Gauss points
            //{
            //    new List<double>{ -s, -s}, new List<double> { s, -s }, new List<double> { s, s }, new List<double> { -s, s }
            //};
            #endregion

            var jacobiansOfElement = new List<Double>();

            // Calculate the Jacobian determinant of each corner point
            for (int n=0; n<naturalPoints.Count; n++)
            {
                double nX = naturalPoints[n][0];
                double nY = naturalPoints[n][1];

                // See documentation for derivation of formula
                double jacobianDeterminantOfPoint = 0.0625 *
                    (
                    ((1 - nY) * (gX[1] - gX[0]) + (1 + nY) * (gX[2] - gX[3]))*
                    ((1 - nX) * (gY[3] - gY[0]) + (1 + nX) * (gY[2] - gY[1]))
                    -
                    ((1 - nY) * (gY[1] - gY[0]) + (1 + nY) * (gY[2] - gY[3])) *
                    ((1 - nX) * (gX[3] - gX[0]) + (1 + nX) * (gX[2] - gX[1]))
                    );

                jacobiansOfElement.Add(jacobianDeterminantOfPoint);
            };

            // Minimum element divided by maximum element. A value of 1 denotes a rectangular element.
            double jacobianRatio = jacobiansOfElement.Min() / jacobiansOfElement.Max();

            // Throw error if Jacobian ratio is outside of range: [0.0, 1.0]
            if (jacobianRatio < 0 || 1.0 < jacobianRatio)
            {
                throw new ArgumentOutOfRangeException(
                    paramName: "jacobianRatio", jacobianRatio,
                    message: "The Jacobian ratio is outside of the valid range [0.0, 1.0]. A negative value might indicate a self-intersecting element.");
            }

            return jacobianRatio;
        }
        
        /// <summary>
        /// Color each mesh face based on a given quality type.
        /// </summary>
        /// <param name="colorMesh">The output colored <see cref="Mesh"/> object.</param>
        /// <param name="qualityList">List of <see cref="Quality"/> objects for each element in the mesh.</param>
        /// <param name="qualityCheckType">Which quality type to color the mesh with.</param>
        void ColorMesh(Mesh colorMesh, List<Quality> qualityList, int qualityCheckType)
        {
            switch (qualityCheckType)
            {
                // 1 = Aspect ratio
                case 1:
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
                    break;
                // 2 = Skewness
                case 2:
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
                    break;
                // 3 = Jacobian
                case 3:
                    foreach (Quality q in qualityList)
                    {
                        if (q.JacobianRatio > 0.8)
                        {
                            q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                        }
                        else if (q.JacobianRatio > 0.5)
                        {
                            q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                        }
                        else if (q.JacobianRatio > 0.03)
                        {
                            q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                        }
                        else if (q.JacobianRatio >= 0)
                        {
                            q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                        }
                        else
                        {
                            q.element.mesh.VertexColors.CreateMonotoneMesh(Color.HotPink);
                        }
                        colorMesh.Append(q.element.mesh);
                    }
                    break;
            }
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
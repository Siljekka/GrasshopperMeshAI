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
            pManager.AddGenericParameter("Mesh3D", "m", "Insert Mesh3D class", GH_ParamAccess.item);

            pManager.AddIntegerParameter("Quality metric", "q", "Aspect Ratio = 1, Skewness = 2, Jacobian = 3", GH_ParamAccess.item);

            pManager[0].Optional = true; // mesh2D is optional
            pManager[1].Optional = true; // mesh3D is optional
            pManager[2].Optional = true; // coloring the mesh is optional
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
            Mesh2D mesh2D = new Mesh2D();
            Mesh3D mesh3D = new Mesh3D();
            Mesh colorMesh = new Mesh();

            // Determines which quality check type to color mesh with
            // 1 = aspect ratio, 2 = skewness, 3 = jacobian
            int qualityCheckType = 0;

            List<Quality> qualityList = new List<Quality>(); // list of Quality for each element in the mesh
            Quality elementQuality = new Quality(); 
            List<Element> elements = new List<Element>();
            double sumAspectRatio = 0;
            double sumSkewness = 0;
            double sumJacobianRatio = 0;
            #endregion

            #region Inputs
            DA.GetData(0, ref mesh2D);
            DA.GetData(1, ref mesh3D);
            DA.GetData(2, ref qualityCheckType);
            #endregion

            #region Calculate mesh quality
            elements = GetElements(mesh2D, mesh3D);

            foreach (Element e in elements)
            {
                elementQuality.AspectRatio = CalculateAspectRatio(e);
                elementQuality.Skewness = CalculateSkewness(e);
                elementQuality.JacobianRatio = CalculateJacobianRatioOfPlaneQuadElement(e);

                elementQuality.element = e;
                e.MeshQuality = elementQuality;

                sumAspectRatio += elementQuality.AspectRatio;
                sumSkewness += elementQuality.Skewness;
                sumJacobianRatio += elementQuality.JacobianRatio;

                qualityList.Add(elementQuality); // todo: check if this is needed
                elementQuality = new Quality();
            }

            double avgAspectRatio = Math.Round(sumAspectRatio / elements.Count, 3);
            double avgSkewness = Math.Round(sumSkewness / elements.Count, 3);
            double avgJacobianRatio = Math.Round(sumJacobianRatio / elements.Count, 3);
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



        List<Element> GetElements(Mesh2D mesh2D, Mesh3D mesh3D)
        {
            List<Element> elements = new List<Element>();

            // check mesh type
            if (mesh2D.Elements != null) { elements = mesh2D.Elements; } // mesh2D input
            else { elements = mesh3D.Elements; } // mesh3D input

            // if no elements are found
            if (elements == null)
            {
                throw new ArgumentException(
                    message: "No elements found.",
                    paramName: "elements");
            }
            return elements;
        }


        #region Caculate Asprect Ratio
        double CalculateAspectRatio(Element e)
        {
            double AR = 0;
            double maxDistance = 0;
            double minDistance = 0;
            double idealAR = 0.5 / Math.Sqrt(Math.Pow(0.5, 2)* 3);

            List<double> cornerToCornerDistance = new List<double>();
            List<double> cornerToCentroidDistance = new List<double>();
            List<double> faceToCentroidDistance = new List<double>();
            
            List<Point3d> faceCenterPts = FindFaceCenter(e);
            Point3d centroidPt = FindCentroidOfElement(e);
            List<Point3d> elementCornerPtsDublicated = FindCoordinatesOfNodesDublicated(e);

            // find distances from corners to centroid
            for (int n = 0; n < elementCornerPtsDublicated.Count / 2; n++)
            {
                cornerToCentroidDistance.Add(elementCornerPtsDublicated[n].DistanceTo(centroidPt));
                cornerToCornerDistance.Add(elementCornerPtsDublicated[n].DistanceTo(elementCornerPtsDublicated[n + 1])); // add the distance between the points, following mesh edges CCW
            }

            // find distances from face center to centroid
            for (int n = 0; n < faceCenterPts.Count; n++)
            {
                faceToCentroidDistance.Add(faceCenterPts[n].DistanceTo(centroidPt));
            }

            cornerToCentroidDistance.Sort();
            cornerToCornerDistance.Sort();
            faceToCentroidDistance.Sort();

            if (!e.IsCube) // mesh2d
            {
                minDistance = cornerToCornerDistance[0];
                maxDistance = cornerToCornerDistance[cornerToCornerDistance.Count - 1];
                AR = minDistance / maxDistance;
            }
            else // mesh3d
            {
                minDistance = Math.Min(cornerToCentroidDistance[0], faceToCentroidDistance[0]);
                maxDistance = Math.Max(cornerToCentroidDistance[cornerToCentroidDistance.Count - 1], faceToCentroidDistance[faceToCentroidDistance.Count - 1]);
                AR = (minDistance / maxDistance)/idealAR; // normalized AR
            }
            return AR;
        }
        #endregion

        #region Caculate Skewness
        double CalculateSkewness(Element e)
        {
            double SK = 0;
            double maxAngle = 0;
            double minAngle = 0;
            double idealAngle = 90; //ideal angle in degrees
            int neighborPoint = 3; //variable used in skewness calculation
            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;
            List<double> elementAngles = new List<double>();

            int numFaces = 1; // mesh2D
            if (e.IsCube) { numFaces = 6; } // mesh3D

            for (int i = 0; i < numFaces; i++)
            {
                e.mesh.Faces.GetFaceVertices(i, out Point3f n1, out Point3f n2, out Point3f n3, out Point3f n4); // todo: sjekk rekkefølge // todo: rett opp i rekkefølge og point3f til point3d
                List<Point3d> faceCornersDublicated = new List<Point3d>()
                {
                        n1, n2, n3, n4,
                        n1, n2, n3, n4
                };

                for (int n = 0; n < faceCornersDublicated.Count / 2; n++)
                {
                    //create a vector from a vertex to a neighbouring vertex
                    vec1 = faceCornersDublicated[n] - faceCornersDublicated[n + 1];
                    vec2 = faceCornersDublicated[n] - faceCornersDublicated[n + neighborPoint];

                    //calculate angles between vectors
                    double angleRad = Math.Abs(Math.Acos(Vector3d.Multiply(vec1, vec2) / (vec1.Length * vec2.Length)));
                    double angleDegree = angleRad * 180 / Math.PI;//convert from rad to deg
                    elementAngles.Add(angleDegree); 
                }

            }
            elementAngles.Sort();
            minAngle = Math.Abs(elementAngles[0]); // todo: controll if abs ok
            maxAngle = Math.Abs(elementAngles[elementAngles.Count - 1]);
            SK = 1 - Math.Max((maxAngle - idealAngle) / (180 - idealAngle), (idealAngle - minAngle) / (idealAngle));
            return SK;
        }
        #endregion

        #region Find Points on a element
        List<Point3d> FindCoordinatesOfNodesDublicated(Element e)
        {
            List<Point3d> pts = new List<Point3d>();
            if (!e.IsCube) // mesh2D
            {
                pts = new List<Point3d>()
                {
                    e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                    e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                };
            }
            else // mesh3D
            {
                pts = new List<Point3d>() 
                {
                    e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                    e.Node5.Coordinate, e.Node6.Coordinate, e.Node7.Coordinate, e.Node8.Coordinate,
                    e.Node1.Coordinate, e.Node2.Coordinate, e.Node3.Coordinate, e.Node4.Coordinate,
                    e.Node5.Coordinate, e.Node6.Coordinate, e.Node7.Coordinate, e.Node8.Coordinate,
                };

            }
            return pts;
        }

        List<Point3d> FindFaceCenter(Element e)
        {
            List<Point3d> faceCenterPts = new List<Point3d>();
            int numFaces = 1; // mesh2D
            if (e.IsCube) { numFaces = 6; } // mesh3D

            for (int i = 0; i < numFaces; i++)
            {
                faceCenterPts.Add(e.mesh.Faces.GetFaceCenter(i));
            }
            return faceCenterPts;
        }

        Point3d FindCentroidOfElement(Element e)
        {
            double sx = 0;
            double sy = 0;
            double sz = 0;
            List<Point3d> pts = FindCoordinatesOfNodesDublicated(e);
            foreach (Point3d pt in pts)
            {
                sx = sx + pt.X;
                sy = sy + pt.Y;
                sz = sz + pt.Z;
            }
            int n = pts.Count;
            Point3d centroidPt = new Point3d(sx/n, sy/n, sz/n);

            return centroidPt;
        }
       
        #endregion

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
                        else if (q.Skewness > 0.75)
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
using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino;
using System;
using System.Linq;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Drawing;
using MathNet.Numerics.LinearAlgebra;

namespace MeshPoints.MeshQuality
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
            pManager.AddGenericParameter("SmartMesh", "sm", "Insert a SmartMesh class", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Quality to color", "q", "Color mesh with quality metric: Aspect Ratio = 1, Skewness = 2, Jacobian = 3", GH_ParamAccess.item);
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
            // Input
            SmartMesh mesh = new SmartMesh(); 
            int qualityCheckType = 0;
            DA.GetData(0, ref mesh);
            DA.GetData(1, ref qualityCheckType);

            if (!DA.GetData(0, ref mesh)) return;

            // Code
            List<Quality> qualityList = new List<Quality>(); // list of Quality for each element in the mesh
            Quality elementQuality = new Quality();
            List<Element> elements = mesh.Elements;
            Mesh colorMesh = new Mesh();

            double sumAspectRatio = 0;
            double sumSkewness = 0;
            double sumJacobianRatio = 0;

            // to do: sjekk map til 2D er nødvendig når vi har surface... 
            foreach (Element e in elements)
            {
                elementQuality.AspectRatio = CalculateAspectRatio(e);
                //elementQuality.AspectRatio = CalculateAspectRatioAnsys(e); // to do: slett, old
                elementQuality.Skewness = CalculateSkewness(e);
                elementQuality.JacobianRatio = CalculateJacobianRatio(e);                
                //elementQuality.JacobianRatio = CalculateJacobianOf8NodeElementOLD(e); // to do: slett, old          
                //elementQuality.JacobianRatio = CalculateJacobianOfQuadElementOLD(e);    // to do: slett, old             

                elementQuality.element = e;
                e.MeshQuality = elementQuality;

                sumAspectRatio += elementQuality.AspectRatio;
                sumSkewness += elementQuality.Skewness;
                sumJacobianRatio += elementQuality.JacobianRatio;

                qualityList.Add(elementQuality);
                elementQuality = new Quality();
            }

            double avgAspectRatio = Math.Round(sumAspectRatio / (double)elements.Count, 3);
            double avgSkewness = Math.Round(sumSkewness / (double)elements.Count, 3);
            double avgJacobianRatio = Math.Round(sumJacobianRatio / (double )elements.Count, 3);

            // Color the mesh based on quality type
            ColorMesh(colorMesh, qualityList, qualityCheckType);

            // Output
            DA.SetDataList(0, qualityList);
            DA.SetData(1, avgAspectRatio);
            DA.SetData(2, avgSkewness);
            DA.SetData(3, avgJacobianRatio);
            DA.SetData(4, colorMesh);
        }

        #region Component methods
        double CalculateAspectRatio(Element element)
        {
            // Calculate Aspact Ratio like Abaqus

            // List of nodes
            List<Point3d> nodeCoordinates = new List<Point3d>(); ;
            foreach (Node node in element.Nodes)
            {
                nodeCoordinates.Add(node.Coordinate);
            }

            // New AR:
            // Find distances from corners to centroid (Abaqus)
            List<double> nodeToNodeDistance = new List<double>();
            if (element.Type != "Hex")
            {
                for (int n = 0; n < nodeCoordinates.Count - 1; n++)
                {
                    if (n == nodeCoordinates.Count - 1)
                    {
                        nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n - 3]));
                        continue;
                    }
                    nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n + 1])); // add the distance between the points, following mesh edges CCW
                }
            }
            else
            {
                for (int n = 0; n < nodeCoordinates.Count - 1; n++)
                {
                    if (n < 0.5 * nodeCoordinates.Count)
                    {
                        nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n + nodeCoordinates.Count / 2]));
                    }
                    if (n == (0.5 * nodeCoordinates.Count - 1) | n == (nodeCoordinates.Count - 1))
                    {
                        nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n - 3]));
                        continue;
                    }
                    nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n + 1]));
                }
            }
            nodeToNodeDistance.Sort();

            double minDistance = nodeToNodeDistance[0];
            double maxDistance = nodeToNodeDistance[nodeToNodeDistance.Count - 1];
            double AR = minDistance / maxDistance;
            return AR;
        }
        /// <summary>
        /// Calculates the Aspect Ratio of an element.
        /// </summary>
        double CalculateAspectRatioAnsys(Element element)
        {
            // Calculate Aspect Ratio like Ansys
            double AR = 0;
            double maxDistance = 0;
            double minDistance = 0;
            double idealAR = 0.5 / Math.Sqrt(Math.Pow(0.5, 2) * 3);

            List<double> nodeToNodeDistance = new List<double>();
            List<double> nodeToCentroidDistance = new List<double>();
            List<double> faceToCentroidDistance = new List<double>();

            List<Point3d> faceCenterPts = GetFaceCenter(element);
            Point3d centroidPt = GetCentroidOfElement(element);
            List<Node> nodes = new List<Node>(element.Nodes);

            // List of nodes 
            List<Point3d> nodeCoordinates = new List<Point3d>(); ;
            foreach (Node node in nodes)
            {
                nodeCoordinates.Add(node.Coordinate);
            }
            nodeCoordinates.Add(nodeCoordinates[0]);


            // Find distances from corners to centroid
            for (int n = 0; n < nodeCoordinates.Count - 1 ; n++)
            {
                nodeToCentroidDistance.Add(nodeCoordinates[n].DistanceTo(centroidPt));
                nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n + 1])); // add the distance between the points, following mesh edges CCW
            }

            // Find distances from face center to centroid
            for (int n = 0; n < faceCenterPts.Count; n++)
            {
                faceToCentroidDistance.Add(faceCenterPts[n].DistanceTo(centroidPt));
            }

            nodeToCentroidDistance.Sort();
            nodeToNodeDistance.Sort();
            faceToCentroidDistance.Sort();

            if (element.Type == "Quad") // surface
            {
                minDistance = nodeToNodeDistance[0];
                maxDistance = nodeToNodeDistance[nodeToNodeDistance.Count - 1];
                AR = minDistance / maxDistance;
            }
            else // solid
            {
                minDistance = Math.Min(nodeToCentroidDistance[0], faceToCentroidDistance[0]);
                maxDistance = Math.Max(nodeToCentroidDistance[nodeToCentroidDistance.Count - 1], faceToCentroidDistance[faceToCentroidDistance.Count - 1]);
                AR = (minDistance / maxDistance) / idealAR; // normalized AR
            }
            return AR;
        }

        /// <summary>
        /// Calculates the Skewness of an element.
        /// </summary>
        double CalculateSkewness(Element element)
        {
            double idealAngle = 90; // ideal angle in degrees
            int neighborPoint = 3; // variable used in skewness calculation, assume quad
            List<double> elementAngles = new List<double>();
            List<List<Node>> faces = element.GetFaces();

            for (int i = 0; i < faces.Count; i++)
            {
                // Create dublicated list of node
                List<Node> nodesOfFace = new List<Node>(faces[i]);
                int numNodesOfFace = nodesOfFace.Count;
                for (int n = 0; n < numNodesOfFace; n++)
                {
                    nodesOfFace.Add(nodesOfFace[n]);
                }

                for (int n = 0; n < numNodesOfFace ; n++)
                {
                    // Create a vector from a vertex to a neighbouring vertex
                    Vector3d vec1 = nodesOfFace[n].Coordinate - nodesOfFace[n + 1].Coordinate;
                    Vector3d vec2 = nodesOfFace[n].Coordinate - nodesOfFace[n + neighborPoint].Coordinate;
                    Vector3d normal = Vector3d.CrossProduct(vec1, vec2);

                    // Calculate angle
                    double angleRad = Vector3d.VectorAngle(vec1, vec2, normal); 
                    double angleDegree = angleRad * 180 / Math.PI; //convert from rad to deg
                    elementAngles.Add(angleDegree);
                }

            }

            elementAngles.Sort();
            double minAngle = elementAngles[0];
            double maxAngle = elementAngles[elementAngles.Count - 1];
            double SK;

            if (minAngle < 0) 
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Skewness: one or more angles are negative");
            }
            if (maxAngle > 180) // if chevron element
            { 
                SK = -1; return SK;
            }
            else
            {
                SK = 1 - Math.Max((maxAngle - idealAngle) / (180 - idealAngle), (idealAngle - minAngle) / (idealAngle)); return SK;
            }
        }

        /// <summary>
        /// Determines the Jacobian ratio based on if the element is from a <see cref="Mesh2D"/> or a <see cref="SmartMesh"/>.
        /// </summary>
        /// <param name="element">A single <see cref="Element"/> object from a mesh.</param>
        /// <returns>A <see cref="double"/> between 0.0 and 1.0 describing the ratio between the min and max values of the determinants of the Jacobian matrix of the element, evaluated in the corner nodes.</returns>
        private double CalculateJacobianRatioOLD(Element element) // to do: slett
        {
            // to do: slett
            double jacobian = 0;
            if (element.Type == "Quad") // surface
            {
                // jacobian = CalculateJacobianOfQuadElementOLD(element); magnus
            }
            else // solid
            {
                // jacobian = CalculateJacobianOf8NodeElementOLD(element); magnus
            }
            return jacobian;
        } 


        /// <summary>
        /// Returns a list of the face centroid to an element
        /// </summary>
        private List<Point3d> GetFaceCenter(Element element)
        {
            List<Point3d> faceCenterPts = new List<Point3d>();
            int numFaces = 1; // surface
            if (element.Type == "Hex") { numFaces = 6; } // solid

            for (int i = 0; i < numFaces; i++)
            {
                faceCenterPts.Add(element.Mesh.Faces.GetFaceCenter(i));
            }
            return faceCenterPts;
        }

        /// <summary>
        /// Return the centroid of an element
        /// </summary>
        Point3d GetCentroidOfElement(Element element)
        {
            double sx = 0;
            double sy = 0;
            double sz = 0;
            List<Node> nodes = element.Nodes;
            List<Point3d> pts = new List<Point3d>(); ;

            foreach (Node node in nodes)
            {
                Point3d pt = node.Coordinate;
                sx = sx + pt.X;
                sy = sy + pt.Y;
                sz = sz + pt.Z; 
            }
            int n = nodes.Count;
            Point3d centroidPt = new Point3d(sx / n, sy / n, sz / n);

            return centroidPt;
        }

        /// <summary>
        /// Calculate the Jacobian Ratio
        /// </summary>
        private double CalculateJacobianRatio(Element element) 
        {
            FEM _FEM = new FEM();
            int nodeDOFS = 2;
            if (element.Type == "Hex") { nodeDOFS = 3; }

            //List<Point3d> localPoints = TransformQuadSurfaceTo2DPoints(cornerPoints); to do: implement non-planar elements

            Matrix<double> globalCoordinates = Matrix<double>.Build.Dense(element.Nodes.Count, nodeDOFS);
            for (int i = 0; i < element.Nodes.Count; i++)
            {
                globalCoordinates[i, 0] = element.Nodes[i].Coordinate.X; // column of x coordinates
                globalCoordinates[i, 1] = element.Nodes[i].Coordinate.Y; // column of y coordinates
                if (nodeDOFS == 3)
                {
                    globalCoordinates[i, 2] = element.Nodes[i].Coordinate.Z; // colum of z coordinates
                }
            }

            // Calculate the Jacobian determinant of each node
            List<double> jacobiansOfElement = new List<double>();
            Matrix<double> gaussNodes = _FEM.GetGaussPoints(1, nodeDOFS);

            for (int n = 0; n < gaussNodes.RowCount; n++)  // loop gauss nodes
            {
                // Substitute the natural coordinates into the symbolic expression
                var r = gaussNodes.Row(n)[0];
                var s = gaussNodes.Row(n)[1];
                double t = 0;
                if (nodeDOFS == 3) { t = gaussNodes.Row(n)[2]; }

                // Partial derivatives of the shape functions
                Matrix<double> shapeFunctionsDerivatedNatural = _FEM.DerivateWithNatrualCoordinates(r, s, t, nodeDOFS);

                // Calculate Jacobian determinant
                Matrix<double> jacobianMatrix = shapeFunctionsDerivatedNatural.Multiply(globalCoordinates);
                double jacobianDeterminant = jacobianMatrix.Determinant();
                jacobiansOfElement.Add(jacobianDeterminant);
            }

            double jacobianRatio = 0;
            // If any of the determinants are negative, we have to divide the maximum with the minimum
            if (jacobiansOfElement.Any(x => x < 0))
            {
                jacobianRatio = jacobiansOfElement.Max() / jacobiansOfElement.Min();
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, $"One or more Jacobian determinants of element {element.Id} is negative.");
                if (jacobianRatio < 0)
                {
                    jacobianRatio = -1;
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, $"The Jacobian Ratio of element {element.Id} is negative.");
                }
            }
            else
            {
                jacobianRatio = jacobiansOfElement.Min() / jacobiansOfElement.Max();
            }

            return jacobianRatio;
        }

        /// <summary>
        /// Calculates the Jacobian Ratio of a hexahedral 8-node 3D element.
        /// </summary>
        /// <param name="element">An 8 node <see cref="Element"/> that is part of a <see cref="SmartMesh"/></param>
        /// <returns>A <see cref="double"/> between 0.0 and 1.0.</returns>
        private double CalculateJacobianOf8NodeElementOLD(Element element)
        {
                var naturalNodes = new List<List<Double>>
            {
                new List<double> { -1, -1, -1 }, new List<double> { 1, -1, -1}, new List<double> { 1, 1, -1 }, new List<double> { -1, 1, -1 },
                new List<double> { -1, -1, 1 }, new List<double> { 1, -1, 1 }, new List<double> { 1, 1, 1 }, new List<double> { -1, 1, 1 }
            };

            // Global X, Y, and Z-coordinates of the corner nodes of the actual element
            List<double> gX = new List<double>();
            List<double> gY = new List<double>();
            List<double> gZ = new List<double>();
            foreach (Node node in element.Nodes)
            {
                gX.Add(node.Coordinate.X);
                gY.Add(node.Coordinate.Y);
                gZ.Add(node.Coordinate.Z);


            }

            // Evaluate in each of the corner nodes.
            List<double> jacobiansOfElement = new List<double>();
            foreach (List<Double> node in naturalNodes)
            {
                // Substitute the natural coordinates into the symbolic expression
                var r = node[0];
                var s = node[1];
                var t = node[2];

                // Partial derivatives of the shape functions
                var N1Dr = -0.125 * (s - 1) * (t - 1);
                var N1Ds = -0.125 * (r - 1) * (t - 1);
                var N1Dt = -0.125 * (r - 1) * (s - 1);
                var N2Dr = 0.125 * (s - 1) * (t - 1);
                var N2Ds = 0.125 * (r + 1) * (t - 1);
                var N2Dt = 0.125 * (r + 1) * (s - 1);
                var N3Dr = -0.125 * (s + 1) * (t - 1);
                var N3Ds = -0.125 * (r + 1) * (t - 1);
                var N3Dt = -0.125 * (r + 1) * (s + 1);
                var N4Dr = 0.125 * (s + 1) * (t - 1);
                var N4Ds = 0.125 * (r - 1) * (t - 1);
                var N4Dt = 0.125 * (r - 1) * (s + 1);
                var N5Dr = 0.125 * (s - 1) * (t + 1);
                var N5Ds = 0.125 * (r - 1) * (t + 1);
                var N5Dt = 0.125 * (r - 1) * (s - 1);
                var N6Dr = -0.125 * (s - 1) * (t + 1);
                var N6Ds = -0.125 * (r + 1) * (t + 1);
                var N6Dt = -0.125 * (r + 1) * (s - 1);
                var N7Dr = 0.125 * (s + 1) * (t + 1);
                var N7Ds = 0.125 * (r + 1) * (t + 1);
                var N7Dt = 0.125 * (r + 1) * (s + 1);
                var N8Dr = -0.125 * (s + 1) * (t + 1);
                var N8Ds = -0.125 * (r - 1) * (t + 1);
                var N8Dt = -0.125 * (r - 1) * (s + 1);

                var sfDr = new List<double> // shape functions differentiated on r
                    {
                        N1Dr, N2Dr, N3Dr, N4Dr, N5Dr, N6Dr, N7Dr, N8Dr
                    };
                var sfDs = new List<double>
                    {
                        N1Ds, N2Ds, N3Ds, N4Ds, N5Ds, N6Ds, N7Ds, N8Ds
                    };
                var sfDt = new List<double>
                    {
                        N1Dt, N2Dt, N3Dt, N4Dt, N5Dt, N6Dt, N7Dt, N8Dt
                    };

                // Evaluates each partial derivative in the isoparametric node
                var calcDerivs = new List<Double>
                    {
                        MultiplyLists(gX, sfDr),
                        MultiplyLists(gX, sfDs),
                        MultiplyLists(gX, sfDt),

                        MultiplyLists(gY, sfDr),
                        MultiplyLists(gY, sfDs),
                        MultiplyLists(gY, sfDt),

                        MultiplyLists(gZ, sfDr),
                        MultiplyLists(gZ, sfDs),
                        MultiplyLists(gZ, sfDt)
                    };

                // Helper function to piecewise multiply elements of two lists of length 8 and add them together
                double MultiplyLists(List<double> a, List<double> b)
                {
                    double sum = 0.0;
                    for (int i = 0; i < 8; i++)
                    {
                        sum += (a[i] * b[i]);
                    }
                    return sum;
                }

                // Structure data in the form of a Jacobian matrix
                Matrix<double> jacobianMatrix = DenseMatrixModule.ofArray2(new double[,]
                {
                        {calcDerivs[0], calcDerivs[3], calcDerivs[6] },
                        {calcDerivs[1], calcDerivs[4], calcDerivs[7] },
                        {calcDerivs[2], calcDerivs[5], calcDerivs[8] },
                });

                var jacobianDeterminant = jacobianMatrix.Determinant();
                jacobiansOfElement.Add(jacobianDeterminant);

            }

            double jacobianRatio = 0;
            // If any of the determinants are negative, we have to divide the maximum with the minimum
            if (jacobiansOfElement.Any(x => x < 0))
            {
                jacobianRatio = jacobiansOfElement.Max() / jacobiansOfElement.Min();

                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, $"One or more Jacobian determinants of element {element.Id} is negative.");
                if (jacobianRatio < 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"The Jacobian Ratio of element {element.Id} is negative.");
                }
            }
            else
            {
                jacobianRatio = jacobiansOfElement.Min() / jacobiansOfElement.Max();
            }

            return jacobianRatio;
        } // to do: slett

        /// <summary>
        /// Calculates the Jacobian ratio of a simple quadrilateral mesh element.
        /// </summary>
        /// <param name="element">An <see cref="Element"/> object describing a mesh face; see <see cref="Element"/> class for attributes.</param>
        /// <returns>A <see cref="double"/> between 0.0 and 1.0.</returns>
        
        
        double CalculateJacobianOfQuadElementOLD(Element e)
        {
            /*
             * This method uses the idea of shape functions and natural coordinate system to 
             * calculate the Jacobian at given points on a simple quadrilateral element.
             * 
             * 1. Transform the input 3D quad element (and specifically the corner points) to a 2D space (X', Y', Z'=0).
             * 2. Define natural coordinates we want to calculate the Jacobian in. This could be the corner points (or 
             *    alternatively the Gauss points) of the quad element. 
             * 3. Calculate the Jacobian determinant of each point. 
             * 4. The ratio is the ratio of the minimum and maximum Jacobian calculated, given as a double from 0.0 to 1.0.
             *    A negative Jacobian indicates a self-intersecting or convex element and should not happen.
                */
            List<Point3d> cornerPoints = new List<Point3d>();
            foreach (Node node in e.Nodes)
            {
                cornerPoints.Add(node.Coordinate);
            }

            List<Point3d> localPoints = TransformQuadSurfaceTo2DPoints(cornerPoints);

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

            // If any of the determinants are negative, we have to divide the maximum with the minimum
            if (jacobiansOfElement.Any(x => x < 0))
            {
                jacobianRatio = jacobiansOfElement.Max() / jacobiansOfElement.Min();

                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, $"One or more Jacobian determinants of element {e.Id} is negative.");
                if (jacobianRatio < 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"The Jacobian Ratio of element {e.Id} is negative.");
                }
            }
            else
            {
                jacobianRatio = jacobiansOfElement.Min() / jacobiansOfElement.Max();
            }

            return jacobianRatio;
        } // to do: slett


        /// <summary>
        /// Transforms the corner points of an arbitrary 3D quad surface to a 2D plane.
        /// </summary>
        /// <param name="cornerPoints">A list of <see cref="Point3d"/> containing the corner points of a mesh element.</param>
        /// <returns>A list of <see cref="Point3d"/> where the Z (third) coordinate is 0.</returns>
        List<Point3d> TransformQuadSurfaceTo2DPoints(List<Point3d> cornerPoints)
        {
            List<Point3d> planeCornerPoints = new List<Point3d>();
            if (Point3d.ArePointsCoplanar(cornerPoints, RhinoMath.ZeroTolerance))
            {
                planeCornerPoints = cornerPoints;
            } 
            else
            {
                // Calculate average surface normal based on corner points
                List<Vector3d> pointNormals = new List<Vector3d>();
                for (int i = 0; i<cornerPoints.Count(); i++)
                {
                    var n = (i + 1) % 4; // next point in quad
                    var p = (i + 3) % 4; // previous point in quad

                    var pointVectors = new List<Vector3d>()
                    {
                        new Vector3d(cornerPoints[n].X - cornerPoints[0].X, cornerPoints[n].Y - cornerPoints[0].Y, cornerPoints[n].Z - cornerPoints[0].Z),
                        new Vector3d(cornerPoints[p].X - cornerPoints[0].X, cornerPoints[p].Y - cornerPoints[0].Y, cornerPoints[p].Z - cornerPoints[0].Z)
                    };
                    pointNormals.Add(Vector3d.CrossProduct(pointVectors[0], pointVectors[1]));
                }
                Vector3d averageNormal = pointNormals[0] + pointNormals[1] + pointNormals[2] + pointNormals[3];

                // Calculate average point of all corner locations.
                double avgX, avgY, avgZ;
                avgX = avgY = avgZ = 0;

                foreach (Point3d point in cornerPoints)
                {
                    avgX += 0.25 * point.X;
                    avgY += 0.25 * point.Y;
                    avgZ += 0.25 * point.Z;
                }
                Point3d averagePoint = new Point3d(avgX, avgY, avgZ);

                // Create a plane: through the "average point" & along the average normal.
                Plane projectPlane = new Plane(averagePoint, averageNormal);

                // Project all points onto the plane along the average normal
                Transform planeProjectionTransformation = Transform.ProjectAlong(projectPlane, averageNormal);
                foreach (Point3d point in cornerPoints)
                {
                    planeCornerPoints.Add(planeProjectionTransformation * point);
                }
            }
            
            // Finally, transform the input points to the xy-plane. 
            Transform elementTransformation = GetTransformationPlanePointsToXYPlane(planeCornerPoints);
            List<Point3d> transformedPoints = new List<Point3d>();
            foreach (Point3d point in planeCornerPoints)
            {
                transformedPoints.Add(elementTransformation * point);
            }
            
            return transformedPoints;


            // Inner method for getting a Transform object/matrix for mapping plane points in 3D space to the xy-plane (z=0).
            Transform GetTransformationPlanePointsToXYPlane(List<Point3d> points)
            {
                var elementVectors = new List<Vector3d>
                {
                    new Vector3d(points[1].X - points[0].X, points[1].Y - points[0].Y, points[1].Z - points[0].Z),
                    new Vector3d(points.Last().X - points[0].X, points.Last().Y - points[0].Y, points.Last().Z - points[0].Z)
                };

                // Cross product of two linearily independent vectors is the normal to the plane containing them
                var surfaceNormal = Vector3d.CrossProduct(elementVectors[0], elementVectors[1]);

                var surfacePlane = new Plane(points[0], surfaceNormal);
                var xyPlane = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));

                return Transform.PlaneToPlane(surfacePlane, xyPlane);
            }
        }
        
        /// <summary>
        /// Color each mesh face based on a given quality type.
        /// </summary>
        /// <param name="colorMesh">The output colored <see cref="Mesh"/> object.</param>
        /// <param name="qualityList">List of <see cref="Quality"/> objects for each element in the mesh.</param>
        /// <param name="qualityCheckType">Which quality type to color the mesh with; 1 = AR, 2 = SK, 3 = J.</param>
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
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                        }
                        else if (q.AspectRatio > 0.5)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                        }
                        else if (q.AspectRatio > 0.1)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                        }
                        else if (q.AspectRatio > 0)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                        }
                        colorMesh.Append(q.element.Mesh);
                    }
                    break;
                // 2 = Skewness
                case 2:
                    foreach (Quality q in qualityList)
                    {
                        if (q.Skewness > 0.75)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                        }
                        else if (q.Skewness > 0.5)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                        }
                        else if (q.Skewness > 0.1)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                        }
                        else if (q.Skewness > 0)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                        }
                        else 
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.HotPink); // invalid mesh
                        }
                        colorMesh.Append(q.element.Mesh);
                    }
                    break;
                // 3 = Jacobian
                case 3:
                    foreach (Quality q in qualityList)
                    {
                        if (q.JacobianRatio > 0.8)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                        }
                        else if (q.JacobianRatio > 0.5)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                        }
                        else if (q.JacobianRatio > 0.03)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                        }
                        else if (q.JacobianRatio >= 0)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                        }
                        else
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.HotPink); // invalid mesh
                        }
                        colorMesh.Append(q.element.Mesh);
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
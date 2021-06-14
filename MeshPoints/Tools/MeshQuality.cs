using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino;
using System;
using System.Linq;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Drawing;
using MathNet.Numerics.LinearAlgebra;

namespace MeshPoints.Tools
{
    public class MeshQuality : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Mesh_Quality class.
        /// </summary>
        public MeshQuality()
          : base("Mesh Quality", "quality",
              "Calculates quality of a SmartMesh.",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh Class.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Color", "col", "Color the mesh with quality metric: Aspect Ratio = 1, Skewness = 2, Jacobian = 3.", GH_ParamAccess.item);
            pManager[1].Optional = true; // coloring the mesh is optional
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality", "mq", "Quality Class.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Avg. Aspect Ratio", "AR", "Average aspect ratio of the mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Avg. Skewness", "SK", "Average skewness of the mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Avg. Jacobian Ratio", "JR", "Average jacobian ratio of the mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Color Mesh", "mesh", "Colored mesh with decided quality metric.", GH_ParamAccess.item);
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

        #region Methods
        /// <summary>
        /// Calculates the Aspect Ratio of an element.
        /// </summary>
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
            if (element.Nodes.Count == 3)
            {
                for (int n = 0; n < nodeCoordinates.Count; n++)
                {
                    if (n == nodeCoordinates.Count - 1)
                    {
                        nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n - 2]));
                        continue;
                    }
                    nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n + 1])); // add the distance between the points, following mesh edges CCW
                }
            }
            else if (element.Type != "Hex")
            {
                for (int n = 0; n < nodeCoordinates.Count; n++)
                {
                    if (n == nodeCoordinates.Count - 1)
                    {
                        nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n - 3]));
                        continue;
                    }
                    nodeToNodeDistance.Add(nodeCoordinates[n].DistanceTo(nodeCoordinates[n + 1])); // add the distance between the points, following mesh edges CCW
                }
            }
            else if (element.Type == "Hex")
            {
                for (int n = 0; n < nodeCoordinates.Count; n++)
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
            else { nodeToNodeDistance = null; }
            nodeToNodeDistance.Sort();

            double minDistance = nodeToNodeDistance[0];
            double maxDistance = nodeToNodeDistance[nodeToNodeDistance.Count - 1];
            double AR = minDistance / maxDistance;
            return AR;
        }

        /// <summary>
        /// Calculates the Skewness of an element.
        /// </summary>
        double CalculateSkewness(Element element)
        {
            double idealAngle = 90; // ideal angle for quad and hex in degrees
            int neighborPoint = 3; // variable used in skewness calculation, assume quad
            if (element.Nodes.Count == 3)
            {
                idealAngle = 60; // ideal angle for triangles in degrees
                neighborPoint = 2;
            }
            
            List<double> elementAngles = new List<double>();
            List<List<Node>> faces = element.GetFaces();

            for (int i = 0; i < faces.Count; i++)
            {
                // Create dublicated list of node
                List<Node> nodesOfFace = new List<Node>(faces[i]);
                nodesOfFace.AddRange(faces[i]);
                Vector3d faceNormal = element.Mesh.FaceNormals[i];
                //if (element.Type == "Hex") { faceNormal.Reverse(); }

                for (int n = 0; n < nodesOfFace.Count / 2; n++)
                {
                   // Create a vector from a vertex to a neighbouring vertex
                    Vector3d vec1 = nodesOfFace[n].Coordinate - nodesOfFace[n + 1].Coordinate;
                    Vector3d vec2 = nodesOfFace[n].Coordinate - nodesOfFace[n + neighborPoint].Coordinate;
                    Vector3d normal = Vector3d.CrossProduct(vec1, vec2);

                    // Calculate angle
                    double angleRad = Vector3d.VectorAngle(vec1, vec2, normal);
                    double testAngle = Vector3d.VectorAngle(Vector3d.CrossProduct(vec1, vec2), faceNormal);
                    if (testAngle >= Math.PI / (double)2) { normal.Reverse(); angleRad = Vector3d.VectorAngle(vec1, vec2, normal); }
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
        /// Calculate the Jacobian Ratio
        /// </summary>
        private double CalculateJacobianRatio(Element element) 
        {
            if (element.Nodes.Count == 3 | element.Nodes.Count == 6)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Component are not able to calculate Jacobian Ratio for tri and teth elements.");
                return 0;
            }
            FEM _FEM = new FEM();
            int nodeDOFS = 2;
            if (element.Type == "Hex") { nodeDOFS = 3; }

            // Get node coordinates
            List<Point3d> nodeCoordinates = new List<Point3d>();
            foreach (Node node in element.Nodes)
            {
                nodeCoordinates.Add(node.Coordinate);
            }

            // Project nodes to planar surface if quads
            if (element.Type == "Quad")
            {
                nodeCoordinates = TransformQuadSurfaceTo2DPoints(nodeCoordinates);
            }

            Matrix<double> globalCoordinates = Matrix<double>.Build.Dense(element.Nodes.Count, nodeDOFS);
            for (int i = 0; i < nodeCoordinates.Count; i++)
            {
                globalCoordinates[i, 0] = nodeCoordinates[i].X; // column of x coordinates
                globalCoordinates[i, 1] = nodeCoordinates[i].Y; // column of y coordinates
                if (nodeDOFS == 3)
                {
                    globalCoordinates[i, 2] = nodeCoordinates[i].Z; // colum of z coordinates
                }
            }

            // Calculate the Jacobian determinant of each node
            List<double> jacobiansOfElement = new List<double>();
            Matrix<double> gaussNodes = _FEM.GetNaturalCoordinate(1, nodeDOFS);

            for (int n = 0; n < gaussNodes.RowCount; n++)  // loop gauss nodes
            {
                // Substitute the natural coordinates into the symbolic expression
                var r = gaussNodes.Row(n)[0];
                var s = gaussNodes.Row(n)[1];
                double t = 0;
                if (nodeDOFS == 3) { t = gaussNodes.Row(n)[2];  }

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
                        if (q.AspectRatio > 0.75) //0.9
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
                        else if (q.Skewness >= 0)
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
                        if (q.JacobianRatio > 0.75)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                        }
                        else if (q.JacobianRatio > 0.5)
                        {
                            q.element.Mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                        }
                        else if (q.JacobianRatio > 0.1)
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
                return Properties.Resources.Icon_Quality;
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
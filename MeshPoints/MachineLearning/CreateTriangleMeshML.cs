using Grasshopper.Kernel;
using Keras.Models;
using MathNet.Numerics.LinearAlgebra;
using MeshPoints.Classes;
using Numpy;
using Rhino;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MeshPoints.MachineLearning
{
    public class CreateTriangleMeshML : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateTriangleMeshML class.
        /// </summary>
        public CreateTriangleMeshML()
          : base("Triangle Mesh (ML)", "NN2",
              "Creates a triangle mesh based on predictions from neural network.",
              "SmartMesh", "Machine Learning")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "Srf", "Input surface to be meshed.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Internal node count", "Int", "The number of internal nodes that should be predicted.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Model path", "Model", "The file path of the ML model to be used for prediction", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Triangle Mesh", "Mesh", "Triangle mesh (Delauney)", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep inputSurface = null;
            int internalNodeCount = 0;
            string modelPath = "";
            DA.GetData(0, ref inputSurface);
            DA.GetData(1, ref internalNodeCount);
            DA.GetData(2, ref modelPath);
            if (!DA.GetData(0, ref inputSurface)) { return; }
            if (!DA.GetData(1, ref internalNodeCount)) { return; }
            if (!DA.GetData(2, ref modelPath)) { return; }

            // Get transformation (and inverse transformation) to/from normalized surface
            Transform procrustesTransform = ProcrustesSuperimposition(inputSurface);
            if (!procrustesTransform.TryGetInverse(out Transform inverseProcrustes))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Procrustes transformation failed.");
                return;
            }

            // Get surface vertices in List<Point3d> (and not BrepVertexList) form.
            List<Point3d> surfaceVertices = new List<Point3d>();
            foreach (var point in inputSurface.Vertices)
            {
                surfaceVertices.Add(point.Location);
            }

            // Transform surface vertices with procrustes. Cast to list as output from transformList is datatype Array.
            List<Point3d> transformedSurfaceVertices = procrustesTransform.TransformList(surfaceVertices).ToList();

            // ARTIFICIAL INTELLIGENCE
            var predictedGrid = GridPrediction(transformedSurfaceVertices, modelPath);
            List<Point3d> predictedInternalNodes = GridPoint.InterpolateNodesFromGridScore(predictedGrid, internalNodeCount);

            // Inverse transform predicted nodes to fit with the input surface
            List<Point3d> inverseTransformedPredictedInternalNodes = inverseProcrustes.TransformList(predictedInternalNodes).ToList();

            // Mesh prediction together with input surface
            Mesh triangleMesh = DelaunayTriangulation(surfaceVertices, inverseTransformedPredictedInternalNodes, inputSurface);

            DA.SetData(0, triangleMesh);
        }

        public Mesh DelaunayTriangulation(List<Point3d> edgeNodes, List<Point3d> internalNodes, Brep meshSurface)
        {
            List<Point3d> nodeCollection = new List<Point3d>();
            nodeCollection.AddRange(edgeNodes);
            nodeCollection.AddRange(internalNodes);

            var meshNodes = new Grasshopper.Kernel.Geometry.Node2List(nodeCollection);

            var meshFaces = new List<Grasshopper.Kernel.Geometry.Delaunay.Face>();
            var triangleMesh = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Mesh(meshNodes, 0.01, ref meshFaces);

            Mesh culledTriangleMesh = CullMeshFacesOutsideSurface(triangleMesh, meshSurface);

            return culledTriangleMesh;
        }

        /// <summary>
        /// Cull unwanted mesh faces by checking if their center points are outside the actual surface of the mesh.
        /// </summary>
        /// <returns>A <see cref="Mesh"/> with (hopefully) no outside mesh faces.</returns>
        public Mesh CullMeshFacesOutsideSurface(Mesh meshSurface, Brep brep)
        {
            Mesh insideFaces = meshSurface.DuplicateMesh();
            for (int i = meshSurface.Faces.Count - 1; i > 0; i--) // reverse iteration to maintain indices
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
        public bool IsPointOnBrepSurface(Point3d point, Brep brep)
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

        public static Transform ProcrustesSuperimposition(Brep inputSurface)
        {
            Brep brepSurface = (Brep)inputSurface.Duplicate();
            int edgeCount = brepSurface.Vertices.Count;
            var referenceContour = CreateRegularNgon(edgeCount);

            // Parametrize surface for translation
            NurbsSurface parametrizedSurface = brepSurface.Faces[0].ToNurbsSurface();
            parametrizedSurface.SetDomain(0, new Interval(0, 1));
            parametrizedSurface.SetDomain(1, new Interval(0, 1));

            // Translation. Also perform the translation on brepSurface for Svd-stuff later
            Transform translationTransformation = Transform.Translation(Point3d.Origin - parametrizedSurface.PointAt(0.5, 0.5));
            brepSurface.Translate(Point3d.Origin - parametrizedSurface.PointAt(0.5, 0.5));

            // Create list of translated surfacePoints and build matrices for Svd-operations
            var surfacePoints = new List<List<double>>();
            foreach (var point in brepSurface.Vertices)
            {
                surfacePoints.Add(new List<double> { point.Location.X, point.Location.Y });
            }

            Matrix<double> referenceMatrix = Matrix<double>.Build.DenseOfRows(referenceContour);
            Matrix<double> surfaceMatrix = Matrix<double>.Build.DenseOfRows(surfacePoints);

            Matrix<double> H = surfaceMatrix.Transpose() * referenceMatrix;
            var svd = H.Svd();

            // Rotation; 2x2 rotation matrix: [[cosø, -sinø], [sinø, cosø]].
            var rotationMatrix = svd.U * svd.VT;
            var rotationInRadians = Math.Acos(rotationMatrix[0, 0]);
            Transform rotationTransformation = Transform.Rotation(-rotationInRadians, Point3d.Origin);

            // Scaling (based on Kabsch algorithm).
            var scalingFactor = surfaceMatrix.FrobeniusNorm() / Math.Sqrt(edgeCount);
            Transform scalingTransformation = Transform.Scale(Point3d.Origin, 1 / scalingFactor);

            // Order is important.
            Transform procrustesSuperimposition = rotationTransformation * scalingTransformation * translationTransformation;

            return procrustesSuperimposition;
        }

        public static List<List<double>> CreateRegularNgon(int edgeCount)
        {
            List<List<Double>> nGon = new List<List<Double>>();

            for (int i = 0; i < edgeCount; i++)
            {
                var coordinates = new List<Double>
                {
                    Math.Cos(2 * Math.PI * i / edgeCount),
                    Math.Sin(2 * Math.PI * i / edgeCount)
                };
                nGon.Add(coordinates);
            }

            return nGon;
        }

        private List<List<GridPoint>> GridPrediction(List<Point3d> contourPoints, string modelPath)
        {
            // Load model
            Keras.Keras.DisablePySysConsoleLog = true;

            var model = Sequential.LoadModel(modelPath);

            // Create empty grid
            var pointGrid = GridPoint.GeneratePointGrid();
            var patches = GridPoint.GeneratePatches(pointGrid);

            // Transform data to proper format for prediction
            double[] contourPointsArray = new double[contourPoints.Count * 2];
            int i = 0;
            foreach (var point in contourPoints)
            {
                contourPointsArray[i * 2] = point.X;
                contourPointsArray[i * 2 + 1] = point.Y;
                i++;
            }

            foreach (var patch in patches)
            {
                // Prepare patch data
                double[] patchCoordinates = new double[8];
                int j = 0;
                foreach (var point in patch)
                {
                    patchCoordinates[j * 2] = point.X;
                    patchCoordinates[j * 2 + 1] = point.Y;
                    j++;
                }

                // Write to numpy for Tensorflow
                var raw_features = np.append(contourPointsArray, patchCoordinates);
                var features = np.expand_dims(raw_features, axis: 0);

                // ARTIFICIAL INTELLIGENCE
                var patchPrediction = model.Predict(features)[0].GetData<float>();

                // Write scores to point grid
                int k = 0;
                foreach (var point in patch)
                {
                    point.Score = patchPrediction[k];
                    k++;
                }
            }
            return pointGrid;
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
                return Properties.Resources.Icon_MLTriangle;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("321310a7-a44e-4286-bc32-46cc3018d116"); }
        }
    }
}
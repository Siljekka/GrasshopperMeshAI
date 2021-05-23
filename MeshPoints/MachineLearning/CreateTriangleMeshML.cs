using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Keras.Models;
using Numpy;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace MeshPoints.MachineLearning
{
    public class CreateTriangleMeshML : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the PredictInternalNodesOfMesh class.
        /// </summary>
        public CreateTriangleMeshML()
          : base("Create Triangle Mesh (Machine Learning)", "TriMeshML",
              "Creates a triangle mesh based on predictions from neural network.",
              "SmartMesh", "Machine Learning")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "sf", "Brep", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Triangle Mesh", "tm", "Triangle mesh (Delauney)", GH_ParamAccess.item);
            pManager.AddGenericParameter("Internal nodes (ML)", "in", "Predicted internal nodes from NN-model.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Transformed contour", "tc", "debugging", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep inputSurface = null;
            DA.GetData(0, ref inputSurface);
            if (!DA.GetData(0, ref inputSurface)) { return; }
            if (inputSurface.Vertices.Count != 6) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input contour must be hexagon."); return; }

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

            // Create array suited for prediction [x1 y1 ... xn yn]
            var transformedSurfaceVerticesForPrediction = new double[inputSurface.Vertices.Count * 2];
            int i = 0;
            foreach (var point in transformedSurfaceVertices)
            {
                transformedSurfaceVerticesForPrediction[i * 2] = point.X;
                transformedSurfaceVerticesForPrediction[i * 2 + 1] = point.Y;
                i++;
            }

            // *** ARTIFICIAL INTELLIGENCE ***
            List<Point3d> predictedInternalNodes = NeuralNetworkPrediction(transformedSurfaceVerticesForPrediction);
            
            // Inverse transform predicted nodes to fit with the input surface
            List<Point3d> inverseTransformedPredictedInternalNodes = inverseProcrustes.TransformList(predictedInternalNodes).ToList();
            
            // Mesh prediction together with input surface
            Mesh triangleMesh = DelaunayTriangulation(surfaceVertices, inverseTransformedPredictedInternalNodes);

            DA.SetData(0, triangleMesh);
            DA.SetDataList(1, inverseTransformedPredictedInternalNodes);
            DA.SetDataList(2, transformedSurfaceVertices);

        }

        public Mesh DelaunayTriangulation(List<Point3d> edgeNodes, List<Point3d> internalNodes)
        {
            List<Point3d> nodeCollection = new List<Point3d>();
            nodeCollection.AddRange(edgeNodes);
            nodeCollection.AddRange(internalNodes);

            var meshNodes = new Grasshopper.Kernel.Geometry.Node2List(nodeCollection);

            var meshFaces = new List<Grasshopper.Kernel.Geometry.Delaunay.Face>();
            var triangleMesh = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Mesh(meshNodes, 0.01, ref meshFaces);

            return triangleMesh;
        }

        public Transform ProcrustesSuperimposition(Brep inputSurface)
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

        public List<List<Double>> CreateRegularNgon(int edgeCount)
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


        List<Point3d> NeuralNetworkPrediction(double[] brepCoordinates)
        {
            Keras.Keras.DisablePySysConsoleLog = true;
            // Load model
            var model = Sequential.LoadModel("C:\\Users\\magnus\\master\\GrasshopperMeshAI\\MeshPoints\\MachineLearning\\models\\direct-internal-nodes\\");

            var features = np.array(brepCoordinates);
            var predictionData = np.expand_dims(features, axis: 0);

            var predictionResult = model.Predict(predictionData);
            var prediction = predictionResult[0].GetData<float>();

            List<Point3d> predictedPoints = new List<Point3d> { 
                new Point3d(prediction[0], prediction[1], 0),
                new Point3d(prediction[2], prediction[3], 0)
            };

            
            return predictedPoints;
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
            get { return new Guid("321310a7-a44e-4286-bc32-46cc3018d116"); }
        }
    }
}
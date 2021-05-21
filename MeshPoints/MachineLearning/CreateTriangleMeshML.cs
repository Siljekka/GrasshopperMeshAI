using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Keras.Models;
using Numpy;

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
            pManager.AddGenericParameter("Internal Nodes", "in", "PointList containing the predicted internal nodes.", GH_ParamAccess.list);
            pManager.AddMeshParameter("Triangle Mesh", "tm", "Triangle mesh (Delauney)", GH_ParamAccess.item);
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

            // Turn surface data into wanted format for prediction (x1 y1 ... xn yn).
            var surfaceVerticesForPrediction = new double[inputSurface.Vertices.Count * 2];
            List<Point3d> surfaceVerticesForMeshing = new List<Point3d>();
            int i = 0;
            foreach (var point in inputSurface.Vertices)
            {
                surfaceVerticesForPrediction[i * 2] = point.Location.X;
                surfaceVerticesForPrediction[i * 2 + 1] = point.Location.Y;

                surfaceVerticesForMeshing.Add(new Point3d(point.Location.X, point.Location.Y, 0) );

                i++;
            }

            List<Point3d> predictedInternalNodes = NeuralNetworkPrediction(surfaceVerticesForPrediction);


            // Mesh edge nodes and internal nodes.
            List<Point3d> nodeCollection = new List<Point3d>();
            nodeCollection.AddRange(surfaceVerticesForMeshing);
            nodeCollection.AddRange(predictedInternalNodes);

            var meshNodes = new Grasshopper.Kernel.Geometry.Node2List(nodeCollection);

            // 4. Throw all our points into the Delaunay mesher. Adjust jitter_amount as needed.
            var meshFaces = new List<Grasshopper.Kernel.Geometry.Delaunay.Face>();
            var triangleMesh = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Mesh(meshNodes, 0.01, ref meshFaces);

            DA.SetDataList(0, predictedInternalNodes);
            DA.SetData(1, triangleMesh);

        }

        List<Point3d> NeuralNetworkPrediction(double[] brepCoordinates)
        {
            Keras.Keras.DisablePySysConsoleLog = true;
            // Load model
            var model = Sequential.LoadModel("C:\\Users\\mkunn\\skole\\master\\Mesh\\MeshPoints\\MachineLearning\\models\\direct-internal-nodes\\");

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
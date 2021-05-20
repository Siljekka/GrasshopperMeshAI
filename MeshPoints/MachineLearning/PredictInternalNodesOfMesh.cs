using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Keras.Models;
using Numpy;

namespace MeshPoints.MachineLearning
{
    public class PredictInternalNodesOfMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the PredictInternalNodesOfMesh class.
        /// </summary>
        public PredictInternalNodesOfMesh()
          : base("PredictInternalNodesOfMesh", "NNpred",
              "Predicts the number/position of the internal nodes of a polygon mesh.",
              "SmartMesh", "Machine Learning")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "sf", "Surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Internal Nodes", "in", "PointList containing the predicted internal nodes.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Surface inputSurface = null;
            DA.GetData(0, ref inputSurface);


            List<Point3d> predictedInternalNodes = NeuralNetworkPrediction(inputSurface);

            DA.SetData(0, predictedInternalNodes);

        }

        List<Point3d> NeuralNetworkPrediction(Surface surface)
        {
            // Load model
            var model = Sequential.LoadModel("models/direct-internal-nodes");

            // Turn surface data into wanted format (x1 y1 ... xn yn) in CCW order.
            NurbsSurface nurbsSurface = surface.ToNurbsSurface();
            List<double> surfacePointCoordinateList = new List<double>();
            foreach ( var point in nurbsSurface.Points)
            {
                surfacePointCoordinateList.Add(point.X);
                surfacePointCoordinateList.Add(point.Y);
            }

            var features = np.array(new[] { surfacePointCoordinateList });
            var predictionData = np.expand_dims(features, axis: 0);

            var predictionResult = model.Predict(predictionData);

            var prediction = predictionResult.GetData<double>();

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
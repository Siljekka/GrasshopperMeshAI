using Grasshopper.Kernel;
using Rhino.Geometry;
using MeshPoints.MachineLearning;
using System;
using System.Linq;
using System.Collections.Generic;
using Keras.Models;
using Numpy;

namespace MeshPoints.MachineLearning
{
    public class PredictInternalNodeCount : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the PredictInternalNodeCount class.
        /// </summary>
        public PredictInternalNodeCount()
          : base("Internal Node Count", "pinc",
              "Predict the number of internal nodes to be calculated by CreateTriangleMeshML",
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
            pManager.AddIntegerParameter("Internal node count", "inc", "The number of internal nodes that should be predicted.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep inputSurface = new Brep();
            DA.GetData(0, ref inputSurface);
            if (!DA.GetData(0, ref inputSurface)) { return; }
            if (inputSurface.Vertices.Count != 8) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Component currently only supports octagons."); return; }


            // Get transformation (and inverse transformation) to/from normalized surface 
            Transform procrustesTransform = CreateTriangleMeshML.ProcrustesSuperimposition(inputSurface);
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
            int incPrediction = NNPrediction(transformedSurfaceVertices);

            DA.SetData(0, incPrediction);
        }
        int NNPrediction(List<Point3d> contourPoints)
        {
            // Load model
            Keras.Keras.DisablePySysConsoleLog = true;
            var model = Sequential.LoadModel("C:\\Users\\mkunn\\skole\\master\\Mesh\\MeshPoints\\MachineLearning\\models\\nn1-8gon");

            // Transform data to proper format for prediction
            double[] contourPointsArray = new double[contourPoints.Count * 2];
            int i = 0;
            foreach (var point in contourPoints)
            {
                contourPointsArray[i * 2] = point.X;
                contourPointsArray[i * 2 + 1] = point.Y;
                i++;
            }

            // Write to numpy for Tensorflow
            var raw_features = np.array(contourPointsArray);
            var features = np.expand_dims(raw_features, axis: 0);

            // AI
            var nnPrediction = model.Predict(features)[0].GetData<float>()[0];
            
            // Round off and return
            return (int)nnPrediction;
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
            get { return new Guid("b5e9b9b1-38fd-4aa5-b44b-be0b5260eab4"); }
        }
    }
}
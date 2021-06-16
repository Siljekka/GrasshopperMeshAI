using Grasshopper.Kernel;
using MeshPoints.Classes;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshQualityMoveInternal;
using MeshPoints.Tools;

namespace MeshPoints.MachineLearning
{
    public class MeshQualityML : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MeshQualityML class.
        /// </summary>
        public MeshQualityML()
          : base("MeshQualityML", "mqml",
              "Mesh Quality calculator using ML.",
              "SmartMesh", "Machine Learning")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "m", "Insert SmartMesh class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Average Quality", "jb", "Average Jacobian ratio of all elements", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SmartMesh mesh = new SmartMesh();

            DA.GetData(0, ref mesh);

            // 0. input: 16 node, 4 internal nodes, square [-1,1]
            // 1. Extract nodes.
            List<Node> nodeList = mesh.Nodes;

            // 2. Transform nodes (procrustes) to get normalized input (can skip this for now).
            // We transform before input.
            // 3. Predict using model from ML.NET. Don't look at it! :|
            // NodeList indexes are weird because create surfaceMesh creates the mesh a bit weirdly.
            ModelInput nodeData = new ModelInput()
            {
                //X1 = Convert.ToSingle(nodeList[0].Coordinate.X),
                //Y1 = Convert.ToSingle(nodeList[0].Coordinate.Y),
                //X2 = Convert.ToSingle(nodeList[1].Coordinate.X),
                //Y2 = Convert.ToSingle(nodeList[1].Coordinate.Y),
                //X3 = Convert.ToSingle(nodeList[2].Coordinate.X),
                //Y3 = Convert.ToSingle(nodeList[2].Coordinate.Y),
                //X4 = Convert.ToSingle(nodeList[3].Coordinate.X),
                //Y4 = Convert.ToSingle(nodeList[3].Coordinate.Y),
                X5 = Convert.ToSingle(nodeList[5].Coordinate.X),
                Y5 = Convert.ToSingle(nodeList[5].Coordinate.Y),
                X6 = Convert.ToSingle(nodeList[9].Coordinate.X),
                Y6 = Convert.ToSingle(nodeList[9].Coordinate.Y),
                //X7 = Convert.ToSingle(nodeList[6].Coordinate.X),
                //Y7 = Convert.ToSingle(nodeList[6].Coordinate.Y),
                //X8 = Convert.ToSingle(nodeList[7].Coordinate.X),
                //Y8 = Convert.ToSingle(nodeList[7].Coordinate.Y),
                X9 = Convert.ToSingle(nodeList[6].Coordinate.X),
                Y9 = Convert.ToSingle(nodeList[6].Coordinate.Y),
                X10 = Convert.ToSingle(nodeList[10].Coordinate.X),
                Y10 = Convert.ToSingle(nodeList[10].Coordinate.Y),
                //X11 = Convert.ToSingle(nodeList[10].Coordinate.X),
                //Y11 = Convert.ToSingle(nodeList[10].Coordinate.Y),
                //X12 = Convert.ToSingle(nodeList[11].Coordinate.X),
                //Y12 = Convert.ToSingle(nodeList[11].Coordinate.Y),
                //X13 = Convert.ToSingle(nodeList[12].Coordinate.X),
                //Y13 = Convert.ToSingle(nodeList[12].Coordinate.Y),
                //X14 = Convert.ToSingle(nodeList[13].Coordinate.X),
                //Y14 = Convert.ToSingle(nodeList[13].Coordinate.Y),
                //X15 = Convert.ToSingle(nodeList[14].Coordinate.X),
                //Y15 = Convert.ToSingle(nodeList[14].Coordinate.Y),
                //X16 = Convert.ToSingle(nodeList[15].Coordinate.X),
                //Y16 = Convert.ToSingle(nodeList[15].Coordinate.Y)
            };
            //// 4. Output predicted quality.
            var predictedQuality = Convert.ToDouble(ConsumeModel.Predict(nodeData).Score);
            DA.SetData(0, predictedQuality);

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
                return Properties.Resources.Icon_MLQuality;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("f1dac36a-4de0-4233-b110-151a79795b72"); }
        }
    }
}
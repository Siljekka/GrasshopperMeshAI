using Grasshopper.Kernel;
using MathNet.Numerics.LinearAlgebra;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MeshPoints.Tools
{
    public class NormalizeSurface : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the NormalizeSurface class.
        /// </summary>
        public NormalizeSurface()
          : base("Normalize Surface", "normalize",
              "Translates a surface to origin and scales to [-1,1].",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "srf", "Surface.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Settings", "ts", "1 = simple normalization. 2 = procrustes superimposition.", GH_ParamAccess.item);
            pManager[1].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "brep", "Translated and scaled surface.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep inputSurface = null;
            double optionInput = 1;
            DA.GetData(0, ref inputSurface);
            DA.GetData(1, ref optionInput);
            if (!DA.GetData(0, ref inputSurface)) { return; }
            int option = (int)Math.Floor(optionInput);

            Brep surface = (Brep)inputSurface.Duplicate();
            Transform transformation = new Transform();
            if (option == 1)
            {
                transformation = SimpleNormalization(surface);
            }
            else if (option == 2)
            {
                transformation = ProcrustesSuperimposition(surface);
            }

            surface.Transform(transformation);

            DA.SetData(0, surface);
        }

        /// <summary>
        /// Creates a transformation that translates a surface to the origin and scales it so that it's bounding box fits between -1 and 1.
        /// </summary>
        /// <param name="inputSurface"></param>
        /// <returns>A <see cref="Transform" /> object consisting of scaling and translation.</returns>
        public Transform SimpleNormalization(Brep inputSurface)
        {
            Brep surface = (Brep)inputSurface.Duplicate();

            // Parametrize surface to find centroid.
            NurbsSurface parametrizedSurface = surface.Faces[0].ToNurbsSurface();
            parametrizedSurface.SetDomain(0, new Interval(0, 1));
            parametrizedSurface.SetDomain(1, new Interval(0, 1));

            // Translation.
            Transform translationTransformation = Transform.Translation(Point3d.Origin - parametrizedSurface.PointAt(0.5, 0.5));
            surface.Translate(Point3d.Origin - parametrizedSurface.PointAt(0.5, 0.5));

            // Scale surface to fit roughly between (-1, -1) and (1, 1)
            BoundingBox boundingBox = surface.GetBoundingBox(false);
            double[] boundingBoxPoints = { boundingBox.Max.MaximumCoordinate, boundingBox.Min.MaximumCoordinate };
            double scalingFactor = boundingBoxPoints.Select(x => Math.Abs(x)).Max();
            Transform scalingTransformation = Transform.Scale(Point3d.Origin, 1 / scalingFactor);

            return scalingTransformation * translationTransformation;
        }

        /// <summary>
        /// Creates a rotation and scale invariant transformation for an input <see cref="Brep"/> surface that best fits it
        /// to a regular polygon (with the same number of sides) inscribed in a unit circle. Simplified procrustes superimposition.
        /// </summary>
        /// <param name="inputSurface"></param>
        /// <returns>A <see cref="Transform" /> object consisting of rotation, scaling and translation.</returns>
        public Transform ProcrustesSuperimposition(Brep inputSurface)
        {
            Brep brepSurface = (Brep)inputSurface.Duplicate();
            int edgeCount = brepSurface.Vertices.Count;
            var referenceContour = CreateRegularUnitNgon(edgeCount);

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

            // Concatenation.
            Transform procrustesSuperimposition = rotationTransformation * scalingTransformation * translationTransformation;

            return procrustesSuperimposition;
        }

        public List<List<Double>> CreateRegularUnitNgon(int edgeCount)
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

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.Icon_NormalizeSurface;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("c738214a-5e14-46c9-a13a-126ff4ae912f"); }
        }
    }
}
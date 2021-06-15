using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace MeshPoints.Tools
{
    public class ProcrustesNormalization : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ProcrustesNormalization class.
        /// </summary>
        public ProcrustesNormalization()
          : base("ProcrustesNormalization", "Procrustes",
              "Transform a plane n-gon to best fit a unit n-gon",
              "SmartMesh", "Tools")
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
            pManager.AddGenericParameter("Surface", "sf", "Transformed Surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Surface inputSurface = null;
            DA.GetData(0, ref inputSurface);

            NurbsSurface surface = inputSurface.ToNurbsSurface();

            NurbsSurface transformedSurface = ProcrustesSuperimposition(surface);

            DA.SetData(0, transformedSurface);
        }

        public NurbsSurface ProcrustesSuperimposition(NurbsSurface surface)
        {
            int edgeCount = surface.Points.Count();
            var referenceContour = CreateRegularNgon(edgeCount);
            
            return surface;
        }

        public List<List<Double>> CreateRegularNgon(int edgeCount)
        {
            List<List<Double>> nGon = new List<List<Double>>();

            for (int i = 0; i<edgeCount; i++)
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
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("19464952-82db-4633-a673-7bdbdde86664"); }
        }
    }
}
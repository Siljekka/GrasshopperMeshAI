using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MeshPoints.CreateMesh
{
    public class CreateTriangleMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateDelauneyTriangleMesh class.
        /// </summary>
        public CreateTriangleMesh()
          : base("Triangle mesh", "TriMesh",
              "Creates a triangle mesh on a (2D) brep using built-in Delaunay method",
              "MyPlugin", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Brep surface", "b", "Insert brep of surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Triangle mesh", "m", "Triangle mesh (Delauney)", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Brep meshSurface = new Brep();
            DA.GetData(0, ref meshSurface);
            if (!meshSurface.IsSurface)
            {
                throw new ArgumentException("Input Brep must be a surface.", "meshSurface");
                //AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Input brep must be a surface.");
            }
            // Output
            Mesh triangleMesh = new Mesh();
            DA.SetData(0, triangleMesh);
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
            get { return new Guid("a07a01a6-a771-4d75-9f96-e87ece274885"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.CreateMesh
{
    public class CreateMesh3D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public CreateMesh3D()
          : base("Create Mesh3D (sweep)", "mesh3D",
              "Creates a solid mesh",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("bottomMesh2D", "bm2D", "Mesh2D for bottom surface", GH_ParamAccess.item);
            pManager.AddGenericParameter("topMesh2D", "bm2D", "Mesh2D for bottom surface", GH_ParamAccess.item);
            pManager.AddIntegerParameter("w", "w", "division in w direction", GH_ParamAccess.item, 4);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh3D", "m3D", "Creates a Mesh3D", GH_ParamAccess.item);
            pManager.AddGenericParameter("test", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Variables
            Mesh2D meshBottom = new Mesh2D(); 
            Mesh2D meshTop = new Mesh2D();
            int nw = 0;

            DA.GetData(0, ref meshBottom);
            DA.GetData(1, ref meshTop);
            DA.GetData(2, ref nw);


            // Code


            // Output
            DA.SetData(0, nw);
            //DA.SetDataList(1, );


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
            get { return new Guid("bf8907fb-fb39-41c7-aa44-c0af8111dccb"); }
        }
    }
}
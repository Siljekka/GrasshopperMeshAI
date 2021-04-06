using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.DeconstructClasses
{
    public class DeconstructMesh3D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructMesh3d class.
        /// </summary>
        public DeconstructMesh3D()
          : base("DeconstructMesh3D", "decM3D",
              "Deconstructing Mesh3D class",
              "MyPlugIn", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh3D", "m3D", "Mesh3D class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "e", "List of elements", GH_ParamAccess.list); //0
            pManager.AddGenericParameter("Nodes", "n", "List of nodes", GH_ParamAccess.list); //1
            pManager.AddGenericParameter("Mesh", "m", "Mesh", GH_ParamAccess.item); //2
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            Mesh3D m = new Mesh3D();
            DA.GetData(0, ref m);

            //output
            DA.SetDataList(0, m.Elements);
            DA.SetDataList(1, m.Nodes);
            DA.SetData(2, m.mesh);
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
                return Properties.Resources.Icon_DeconstructSolidMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("97c30c27-48c9-41ac-b09d-d02f80e806f6"); }
        }
    }
}
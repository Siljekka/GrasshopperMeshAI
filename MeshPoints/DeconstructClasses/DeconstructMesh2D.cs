using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

//Deconstruct the class Mesh2D

namespace MeshPoints
{
    public class Deconstruct_Mesh2D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Deconstruct_Mesh2D class.
        /// </summary>
        public Deconstruct_Mesh2D()
          : base("Deconstruct Mesh2D", "decM2D",
              "Deconstructing Mesh2D class",
              "MyPlugIn", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh2D", "m2D", "Mesh2D class", GH_ParamAccess.item);

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
            Mesh2D m = new Mesh2D();
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
                return Properties.Resources.Icon_DeconstructSurfaceMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6cb2b0ff-e396-4fb5-a0b9-33ff368636ac"); }
        }
    }
}
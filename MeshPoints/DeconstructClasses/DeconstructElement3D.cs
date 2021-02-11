using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.DeconstructClasses
{
    public class DeconstructElement3D : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        /// <summary>
        /// Initializes a new instance of the Deconstruct_Node class.
        /// </summary>
        public DeconstructElement3D()
        : base("Deconstruct Element 3D", "decE3D",
              "Deconstructing element class for 3D Mesh Class",
              "MyPlugIn", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Element", "e", "Element class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Node 1", "n", "Node 1", GH_ParamAccess.item); //0
            pManager.AddGenericParameter("Node 2", "n", "Node 2", GH_ParamAccess.item); //1
            pManager.AddGenericParameter("Node 3", "n", "Node 3", GH_ParamAccess.item); //2
            pManager.AddGenericParameter("Node 4", "n", "Node 4", GH_ParamAccess.item); //3
            pManager.AddGenericParameter("Node 5", "n", "Node 5", GH_ParamAccess.item); //4
            pManager.AddGenericParameter("Node 6", "n", "Node 6", GH_ParamAccess.item); //5
            pManager.AddGenericParameter("Node 7", "n", "Node 7", GH_ParamAccess.item); //6
            pManager.AddGenericParameter("Node 8", "n", "Node 8", GH_ParamAccess.item); //7

            pManager.AddGenericParameter("Id", "id", "Element Id", GH_ParamAccess.item); //8
            pManager.AddGenericParameter("Mesh", "m", "Element mesh", GH_ParamAccess.item); //9
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            Element e = new Element();
            DA.GetData(0, ref e);

            //output
            DA.SetData(0, e.Node1);
            DA.SetData(1, e.Node2);
            DA.SetData(2, e.Node3);
            DA.SetData(3, e.Node4);
            DA.SetData(4, e.Node5);
            DA.SetData(5, e.Node6);
            DA.SetData(6, e.Node7);
            DA.SetData(7, e.Node8);
            DA.SetData(8, e.Id);
            DA.SetData(9, e.mesh);
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
                return Properties.Resources.Image1;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("df3f3c5f-c32c-44ce-83f1-831a94edd1d8"); }
        }
    }
}
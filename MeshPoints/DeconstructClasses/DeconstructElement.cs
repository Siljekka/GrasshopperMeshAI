using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

// Deconstruct the class Element

namespace MeshPoints
{
    public class Deconstruct_Element : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Deconstruct_Element class.
        /// </summary>
        public Deconstruct_Element()
          : base("Deconstruct Element", "decE",
              "Deconstructing element class",
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

            pManager.AddGenericParameter("Id", "id", "Element Id", GH_ParamAccess.item); //4 
            pManager.AddGenericParameter("Mesh", "m", "Element mesh", GH_ParamAccess.item); //5
            //pManager.AddGenericParameter("Aspect Ratio", "ar", "Aspect ratio quality of element", GH_ParamAccess.item); //6
            //pManager.AddGenericParameter("Skewness", "ar", "Skewness quality of element", GH_ParamAccess.item); //7
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
            DA.SetData(4, e.Id);
            DA.SetData(5, e.mesh);
            //DA.SetData(6, e.quality.AspectRatio);
            //DA.SetData(7, e.quality.Skewness);
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
            get { return new Guid("149bdc58-2f8a-4b90-ae82-62389190e956"); }
        }
    }
}
﻿using Grasshopper.Kernel;
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
        : base("Deconstruct Element", "decE",
              "Deconstructing element class for SmartMesh Class",
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
            pManager.AddGenericParameter("Nodes", "nodes", "Nodes of element", GH_ParamAccess.list); //0
            pManager.AddGenericParameter("Connectivity", "con", "Connectivity of local to global nodes", GH_ParamAccess.list); //0
            pManager.AddGenericParameter("Type", "type", "Element type", GH_ParamAccess.item); //0
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
            DA.SetDataList(0, e.Nodes);
            DA.SetDataList(1, e.Connectivity);
            DA.SetData(2, e.Type);
            DA.SetData(3, e.Id);
            DA.SetData(4, e.mesh);
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
                return Properties.Resources.Icon_DeconstructSolidElement;
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
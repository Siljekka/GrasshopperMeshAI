using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MeshPoints
{
    public class GenerateINPfile : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GenerateINPfile class.
        /// </summary>
        public GenerateINPfile()
          : base("Generate inp-file", "inp",
              "Generate inp-file for .... Default material is STEEL.",
              "MyPlugIn", "inp-file")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "nodes", "Nodes from mesh", GH_ParamAccess.list);
            pManager.AddGenericParameter("ElementType (string)", "element", "String with element type from Abaqus", GH_ParamAccess.item);
            //pManager.AddGenericParameter("Material (string)", "material", "String with material from Abaqus", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("inp-file", "inp", "String containing inp-file", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
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
            get { return new Guid("1e9c1cd4-77ed-49a2-bd0c-34fb0871f992"); }
        }
    }
}
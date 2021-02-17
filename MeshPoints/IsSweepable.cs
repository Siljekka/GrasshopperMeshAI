using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;
using System.Linq;

namespace MeshPoints
{
    public class IsSweepable : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public IsSweepable()
          : base("IsSweepable", "Sweepable",
              "Checks if brep is Sweepable",
              "MyPlugIn", "Brep")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "bp", "Brep to check if sweepable", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("IsSweepable", "sweepable", "True if brep is sweepable", GH_ParamAccess.item);
            pManager.AddGenericParameter("testList", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Variables
            Brep bp = new Brep();
            List<Brep> bpFace1 = new List<Brep>();
            List<Curve> bpEdge1 = new List<Curve>();
            BrepFaceList bpFace;
            BrepEdgeList bpEdge;

            //Input
            DA.GetData(0, ref bp);

            //Code
            
            bpFace = bp.Faces;
            bpEdge = bp.Edges;

            foreach (BrepFace bf in bp.Faces)
            {
                bpFace1.Add(bf.Brep);

                foreach (BrepEdge be in bp.Edges)
                {
                    bpEdge1.Add(be.EdgeCurve);
                }
            }
            


            

            //Output
            //DA.SetData(0, );
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
            get { return new Guid("4b1bfdce-752e-4b21-8792-5704b036f4ea"); }
        }
    }
}
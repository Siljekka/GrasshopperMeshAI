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
              "Checks if brep is sweepable",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "brep", "Brep to check if sweepable.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("IsSweepable", "sweepable", "True if brep is sweepable.", GH_ParamAccess.item);
            pManager.AddGenericParameter("SweepableFaces", "faces", "List of sweepable faces of the brep.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Brep brep = new Brep();
            DA.GetData(0, ref brep);


            // Code
            bool sweepable = false;
            List<BrepFace> brepFaceDuplicate = new List<BrepFace>();
            List<BrepFace> sweepableFaces = new List<BrepFace>();
            BrepFaceList brepFace = brep.Faces;

            for (int i = 0; i < brepFace.Count; i++) // loop through every face of brep
            {
                List<int> indexAdjecentFaces = (brepFace[i].AdjacentFaces()).ToList();  // find index to faces adjacent to face i
                foreach (int j in indexAdjecentFaces)
                {
                    brepFaceDuplicate.Add(brepFace[j]); // make new list with faces adjacent to face i
                }
                brepFaceDuplicate.Add(brepFace[i]); // add face i to the list

                int countRemainingFaces = brepFace.Count - brepFaceDuplicate.Count; // count number of faces which are not adjacent to face i
                if (countRemainingFaces == 1) { sweepable = true; sweepableFaces.Add(brepFace[i]); } // check if brep is sweepable, and add sweepable face to list
                indexAdjecentFaces.Clear(); // clear list
                brepFaceDuplicate.Clear(); // clear list
            }

            if (!sweepable) { sweepableFaces.Add(null); } // if brep not sweepable; list of sweepable faces = null


            // Output
            DA.SetData(0, sweepable);
            DA.SetDataList(1, sweepableFaces);
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
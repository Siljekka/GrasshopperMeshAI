using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;

namespace MeshPoints
{
    public class FEMLoad : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMLoad class.
        /// </summary>
        public FEMLoad()
          : base("FEM Load", "FEM load",
              "Create load for FEM solver.",
              "MyPlugIn", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {        
            pManager.AddGenericParameter("MeshGeometry", "MeshFeometry", "Input a MeshGeometry", GH_ParamAccess.item); // to do: change name
            pManager.AddGenericParameter("Load type", "load type", "Point load = 1, Surface load = 2", GH_ParamAccess.item);
            pManager.AddGenericParameter("Position", "pos", "Coordinate for point load", GH_ParamAccess.item);  // to do: se om dette skal være en liste ift input i FEMsolver
            pManager.AddGenericParameter("Surface", "surface", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("Load size", "size", "Load size (kN for point load, kN/mm^2 for surface)", GH_ParamAccess.item); // to do: se om dette skal være en liste ift input i FEMsolver
        }
    

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
             pManager.AddGenericParameter("Load", "load", "String of load information", GH_ParamAccess.list);
        }

    /// <summary>
    /// This is the method that actually does the work.
    /// </summary>
    /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
    protected override void SolveInstance(IGH_DataAccess DA)
        {
            // assume only perpendicular negativ load

            #region Input
            Mesh3D mesh = new Mesh3D(); // to do: change to MeshGeometry elns
            int loadType = 0;
            Point3d loadPositions = new Point3d();
            double loadSize = 0;

            DA.GetData(0, ref mesh);
            DA.GetData(1, ref loadType);
            DA.GetData(2, ref loadPositions);
            DA.GetData(3, ref loadSize);
            #endregion

            #region Code
            List<Node> nodes = mesh.Nodes;
            List<Element> element = mesh.Elements;

            if (loadType == 1)
            {
                // get node closest to this point
            }
            else if (loadType == 2)
            { 
                // surface load
                /*
                 * loop nodes
                 * find elements connected
                 * calculate area
                */
            }
            #endregion

       // return a list of double : R

        }

        #region Methods
        private double CalculateArea(Element element)
        {
            double area = 0;
            

            return area;
        }

        #endregion

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
            get { return new Guid("8972a393-9603-486a-bf29-a436a72d2c8d"); }
        }
    }
}
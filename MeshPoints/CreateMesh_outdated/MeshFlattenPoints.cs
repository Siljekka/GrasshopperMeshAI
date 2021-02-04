using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;


namespace MeshPoints
{
    public class MeshFlattenPoints : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public MeshFlattenPoints()
          : base("MeshFlattenPoints", "MeshPts",
              "Create mesh between given a flatten list of points",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "pts", "Insert flatten list of points", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Mesh between points", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Variables
            Mesh m = new Mesh();
            List<Point3d> pts = new List<Point3d>();
           
            int nx = 2;//Number points in x-dir, start by adding first/last point in a row
            int ny = 1; //Number points in y-dir, start by adding first point in a colomn
            int counter = 0;

            double dist1 = 0;
            double dist2 = 0;
            Boolean completeRow = false;

            //Input
            DA.GetDataList(0, pts);

            #region Find nx and ny
            for (int i = 0; i < pts.Count-2; i++)
            {
                dist1 = Math.Abs(pts[0].X - pts[i + 1].X); //distance from start point to point[i+1]
                dist2 = Math.Abs(pts[0].X - pts[i + 2].X);  //distance from start point to point[i+2]
                if (dist1 < dist2)
                {
                    if (!completeRow) 
                    {   
                        nx++; //count inner points in first row
                    }   
                }
                else
                {
                    ny++; //count each end of row
                    completeRow = true; 
                }
            }
            #endregion
    
            #region Loop Vertices        
            for (int i = 0; i < pts.Count; i++)
            {
                m.Vertices.Add(pts[i]); //Add point as mesh vertice
            }
            #endregion

            int newRow = 0;
            #region Loop MeshFace
            for (int i = 0; i < (nx-1)*(ny-1); i++)
            {
                m.Faces.AddFace(counter, counter + 1, counter + nx + 1, counter + nx); //Add MeshFace; Constrol nx vs ny
                counter++;
                newRow++; ;
                if (newRow == (nx-1)) //new row
                {
                    counter++;
                    newRow = 0;
                }
            }
            #endregion

            m.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            m.Normals.ComputeNormals(); //Control if needed
            m.Compact(); //to ensure that it calculate

            // Output
            DA.SetData(0, m);
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
            get { return new Guid("32f15bcb-d828-4dcc-9575-7c5e258b6277"); }
        }
    }
}
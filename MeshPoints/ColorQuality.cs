using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Drawing;

namespace MeshPoints
{
    public class ColorQuality : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ColorQuality class.
        /// </summary>
        public ColorQuality()
          : base("Color Quality", "cq",
              "Add color map to mesh quality",
              "MyPlugIn", "Quality")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality list", "q", "List of mesh quality", GH_ParamAccess.list);
            pManager.AddGenericParameter("AR/SK", "q", "AR=1, SK=2", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Color Mesh", "cm", "Color map over quality check", GH_ParamAccess.item); //0
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //variables
            List<Quality> qList = new List<Quality>();
            double check = 0;
            Mesh colorMesh = new Mesh();

            //input
            DA.GetDataList(0, qList);
            DA.GetData(1, ref check);

            if (check == 1)
            {
                foreach (Quality q in qList)
                {
                    //Aspect Ratio
                    if (q.AspectRatio > 0.9)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                    }
                    else if (q.AspectRatio > 0.7)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                    }
                    else if (q.AspectRatio > 0.6)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                    }
                    else if (q.AspectRatio > 0)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                    }
                    colorMesh.Append(q.element.mesh);
                }
            }
            else if (check == 2)
            {
                foreach (Quality q in qList)
                {
                    if (q.Skewness > 0.9)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Green);
                    }
                    else if (q.Skewness > 0.7)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Yellow);
                    }
                    else if (q.Skewness > 0.6)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Orange);
                    }
                    else if (q.Skewness > 0)
                    {
                        q.element.mesh.VertexColors.CreateMonotoneMesh(Color.Red);
                    }
                    colorMesh.Append(q.element.mesh);
                }
            }
                
            //output
            DA.SetData(0, colorMesh);
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
            get { return new Guid("e8986d59-0c66-496e-9337-fc5044bf0c71"); }
        }
    }
}
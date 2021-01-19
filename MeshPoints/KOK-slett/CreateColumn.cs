using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace FirstPlugin
{
    public class CreateColumn : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the _4Plugin class.
        /// </summary>
        public CreateColumn()
          : base("CreateColumn", "cc",
              "create a column class object",
              "PC2021", "advanced")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
           pManager.AddPointParameter("basicPoint","bP","basic point of column", GH_ParamAccess.item);
pManager.AddNumberParameter("height", "h", "height of the column", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("column", "c", "our own class", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //variables
            Point3d pt = new Point3d(0, 0, 0);
            double h = 300;
            double rad = 10;

            //input
            DA.GetData(0, ref pt);
            DA.GetData(1, ref h);


            //code
            Point3d ept = new Point3d(pt.X,pt.Y,pt.Z+h);
            Line line = new Line(pt, ept);


            /*
             //Old version of constructor
            Column c = new Column();
            c.name = "straigth column";
            c.axis = line;
            c.brep = breps[0];
            */

            Curve rail = line.ToNurbsCurve(); //trick to read it 
            Circle ci = new Circle(pt, rad);
            Curve shape = ci.ToNurbsCurve();
            Brep[] breps = Brep.CreateFromSweep(rail, shape, true, 0.0001); //[] to make it a list


            //new constructur
            Column c = new Column("circular column", line, breps[0]);

            //output
            DA.SetData(0, c.name);
            DA.SetData(1, c.axis);
            DA.SetData(2, c.brep);
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
            get { return new Guid("ec646eab-26a8-4431-89d6-2598b9a7f912"); }
        }
    }
}
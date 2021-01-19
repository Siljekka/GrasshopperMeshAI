using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using BuildingCreator.Classes; //accessing the class
using System.Linq;
//using Rhino.Geometry.Intersect; //could add this, but since we use it only one time we do not use it


// In order to load the result of this wizard, you will also need to
// add the output bin/ folder of this project to the list of loaded
// folder in Grasshopper.
// You can use the _GrasshopperDeveloperSettings Rhino command for that.

namespace BuildingCreator
{
    public class BuildingCreatorComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public BuildingCreatorComponent()
          : base("BuildingCreator", "BuildingCreator",
              "Description",
              "PC2021", "Building")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Brep", "b", "Brep to create the builing inside", GH_ParamAccess.item); //0
            pManager.AddNumberParameter("Parameters", "param", "Secondary parameters to adjust the building", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Building class", "building", "", GH_ParamAccess.item); //only output
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Assign input data
            Brep boundary = new Brep(); //0
            List<double> settings = new List<double>(); //1

            if (!DA.GetData(0, ref boundary)) return;
            if (!DA.GetDataList(1, settings)) return;

            // Get settings later
            double floorHeight = settings[0];
            double colDiameter = settings[1];
            double slabThickness = settings[2];
            double dx = settings[3];
            double dy = settings[4];

            #endregion
            #region Create Builing

            Building building = new Building();
            building.Name = "MyBuilding";

            List<Brep> surfaces = CreateLevelSrufaces(boundary, floorHeight);

            if (building.Levels.Count < 2)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The number of levels of your building should be at least 2");
            }

            List<Level> levels = new List<Level>();

            for (int i = 0; i < surfaces.Count- 1; i++) 
            {
                //Create slab
                Slab slab = new Slab(surfaces[i], "concrete B30", slabThickness, "Slab: " + (i + 1).ToString());
                List<Point3d> columnPts = ColumnGrid(surfaces[0], dx, dy, colDiameter / 2);
                List<Column> columns = new List<Column>();

                for (int j = 0; j < columnPts.Count ; j++)
                {
                    double[] center = new double[3] { columnPts[j].X, columnPts[j].Y, columnPts[j].Z};
                    double cHeight = floorHeight - slabThickness;
                    string id = "Column: " +(i + 1).ToString();
                    string cMaterial= "concrete B30",
                    Column col = new Column(cHeight, colDiameter, center, cMaterial, id);
                    columns.Add(col);                
                }

                Level level = new Level(slab, column, "Slab: " + (i + 1).ToString());
                levels.Add(level);

             }
            Slab roof = new Slab(surfaces[surfaces.Count - 1], "Roof material", slabThickness, "Roof");
            levels.Add(new Level(roof, "Roof level"));

            building.Levels= levels;


            List<Point3d> columnPts = ColumnGrid(surfaces[0], dx, dy, colDiameter/2);
           
            #endregion
            //Set DAta
            DA.SetDataList(0, ref building);
        }

        private List<Brep> CreateLevelSrufaces(Brep b, double floorH)
        {
            List<Point3d> centers = new List<Point3d>();

            foreach (BrepFace f in b.Faces)
            {
                centers.Add(AreaMassProperties.Compute(f).Centroid);             
            }

            var sortedCenters = centers.OrderBy(pt => pt.Z).ToList(); // Sorting the center of points by their z-coord. 
            double height = sortedCenters[sortedCenters.Count - 1].Z - sortedCenters[0].Z;
            int numStories = (int) (height / floorH);


            List<Plane> planes = new List<Plane>();
            Vector3d unitZ = new Vector3d(0, 0, 1);

            double x0 = sortedCenters[0].X;
            double y0 = sortedCenters[0].Y;
            double z0 = sortedCenters[0].Z;

            for (int i = 0; i < numStories; i++)
            {
                Plane plane = new Plane(new Point3d(x0,y0,z0+(i*floorH)), unitZ);
            }

            List<Brep> srfs = new List<Brep>();

            foreach (Plane pln in planes)
            {
                Curve[] curves;
                Point3d[] pts;
                Rhino.Geometry.Intersect.Intersection.BrepPlane(b, pln, 0.001, out curves, out pts ); //use out to assign the lists directly (?)

                Brep[] srf = Brep.CreatePlanarBreps(curves, 0.001);
                srfs.Add(srf[0]);
            }

            return null;
        
        }

        private List<Point3d> ColumnGrid(Brep srfs, double dx, double dy, double colOffset)
        {
            Surface surface = srfs.Surfaces[0];
            surface.SetDomain(0, new Interval(0.0, 1.0));
            surface.SetDomain(1, new Interval(0.0, 1.0));

            Point3d p1 = surface.PointAt(0, 0);
            Point3d p2 = surface.PointAt(1, 0);
            Point3d p3 = surface.PointAt(1, 1);
            Point3d p4 = surface.PointAt(0, 1);

            List<Point3d> vertices = new List<Point3d>() { p1, p2, p3, p4 };

            double xMin = vertices.Min(pt => pt.X) + colOffset;
            double yMin = vertices.Min(pt => pt.Y) + colOffset;
            double xMax = vertices.Min(pt => pt.X);
            double yMax = vertices.Min(pt => pt.Y);

            double xDist = xMax - xMin;
            double yDist = yMax - yMin;

            int divX = (int)(xDist / dx) + 1;
            int divY = (int)(yDist / dy) + 1;

            List<Point3d> columnCenters = new List<Point3d>();

            BrepEdge[] edges = srfs.Edges.ToArray();
            List<Curve> edgeCurves = new List<Curve>();

            foreach (BrepEdge e in edges)
            {
                edgeCurves.Add( e.ToNurbsCurve());
            }

            Curve[] boundaryCurves = Curve.JoinCurves(edgeCurves);
            Plane plane = new Plane(p1, new Vector3d(0, 0, 1));

            for(int u = 0; u < divX; u++)
            {
                for (int v = 0; v < divY; v++)
                {
                    int inside = 0;
                    Point3d pt = new Point3d(xMin + (u * dx), yMin + (v * dy), p1.Z);
                    for (int i = 0; i < boundaryCurves.Length; i++)
                    {
                        if (boundaryCurves[i].Contains(pt, plane, 0.1) == PointContainment.Inside)
                        {
                            inside++;                       
                        }

                    }
                    if (inside == 1)
                    {

                        columnCenters.Add(pt);
                    }

                
                }
            }

            return null;





        }








        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("2ea70fdd-a570-4fdd-b5e9-713b3ace2297"); }
        }
    }
}

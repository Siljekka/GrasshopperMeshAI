using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;

namespace MeshPoints.CreateMesh
{
    public class SurfaceSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateSurfaceMesh class.
        /// </summary>
        public SurfaceSmartMesh()
          : base("Surface SmartMesh", "SurfaceSM",
              "Creates a SmartMesh of quadrilateral-elements.",
              "SmartMesh", "Mesh")
        {
        }


        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "srf", "Surface.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("u", "u", "Number element in u-direction.", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("v", "v", "Number element in v-direction.", GH_ParamAccess.item, 4);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "mesh", "Mesh (quadrilateral-elements).", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Brep brep = new Brep();
            int u = 0;
            int v = 0;
            DA.GetData(0, ref brep);
            DA.GetData(1, ref u);
            DA.GetData(2, ref v);

            // 1. Check input OK.
            if (!DA.GetData(0, ref brep)) return;
            if (u == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "u cannot be zero."); return; }
            if (v == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "v cannot be zero."); return; }

            // 2. Assign geometrical properties to mesh
            SmartMesh smartMesh = new SmartMesh();
            Geometry brepGeometry = new Geometry(brep, brep.Faces.ToList(), brep.Edges.ToList(), brep.Vertices.ToList());
            smartMesh.nu = u + 1;
            smartMesh.nv = v + 1;
            smartMesh.nw = 1;
            smartMesh.Type = "Surface";
            smartMesh.Geometry = brepGeometry;

            // 3. Generate grid of points on surface
            NurbsSurface nurbsSurface = NurbsSurface.CreateNetworkSurface(brep.Edges, 0, 0.0001, 0.0001, 0.0001, out int error); // make planar brep to nurbssurface
            List<Point3d> meshPoints = CreateGridOfPointsUV(nurbsSurface, u, v);

            // 4. Create nodes 
            smartMesh.CreateNodes(meshPoints, u, v, 0);

            // 5. Set elements
            smartMesh.CreateQuadElements();

            // 6. Set global mesh
            smartMesh.CreateMesh();

            // Output
            DA.SetData(0, smartMesh);
            DA.SetData(1, smartMesh.Mesh);
        }
        
        /// <summary>
        /// Makes grid of points in U and V direction
        /// </summary>
        /// <returns> List of points in U and V direction</returns>
        private List<Point3d> CreateGridOfPointsUV(NurbsSurface surface, int u, int v)
        {
            List<Point3d> pt = new List<Point3d>();

            var uDomain = surface.Domain(0);
            var vDomain = surface.Domain(1);
            double stepU = uDomain.Length / (double)u;
            double stepV = vDomain.Length / (double)v;

            double pointU = 0;
            double pointV = 0;
            for (double i = 0; i <= v; i++)
            {
                for (double j = 0; j <= u; j++)
                {
                    pt.Add(surface.PointAt(pointU, pointV));  // make point on surface
                    pointU = pointU + stepU;
                }
                pointV = pointV + stepV;
                pointU = 0;
            }
            return pt;
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
                return Properties.Resources.Icon_Surface;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("77485b0a-e12c-467e-8735-381d35f0f2ff"); }
        }
    }
}
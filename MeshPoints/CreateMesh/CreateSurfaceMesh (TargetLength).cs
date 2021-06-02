using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;

namespace MeshPoints.CreateMesh
{
    public class CreateSurfaceMesh__TargetLength_ : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateSurfaceMesh__TargetLength_ class.
        /// </summary>
        public CreateSurfaceMesh__TargetLength_()
          : base("CreateSurfaceMesh (TargetLength)", "Nickname",
              "Mesh a surface with a target length for the elements.",
              "SmartMesh", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "srf", "Surface", GH_ParamAccess.item);
            pManager.AddNumberParameter("TargetLength", "target", "Target Length for the elements", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "surface", "SmartMesh generated", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "m", "Mesh (surface elements).", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Brep brep = new Brep();
            double target = 0;
            DA.GetData(0, ref brep);
            DA.GetData(1, ref target);

            // 1. Check input OK.
            if (!DA.GetData(0, ref brep)) return;
            if (target == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Target length cannot be zero."); return; }

            // 2. Set number of divisions in u- and v-direction
            NurbsSurface surface = brep.Faces[0].ToNurbsSurface();
            int u = (int)Math.Ceiling(surface.Domain(0).Length / target);
            int v = (int)Math.Ceiling(surface.Domain(1).Length / target);

            // 2. Assign geometrical properties to mesh
            SmartMesh smartMesh = new SmartMesh();
            Geometry brepGeometry = new Geometry(brep, brep.Faces.ToList(), brep.Edges.ToList(), brep.Vertices.ToList());
            smartMesh.nu = u + 1;
            smartMesh.nv = v + 1;
            smartMesh.nw = 1;
            smartMesh.Type = "Surface";
            smartMesh.Geometry = brepGeometry;


            // 3. Generate grid of points on surface
            Brep[] planarBrep = Brep.CreatePlanarBreps(brep.Edges, 0.0001); // make planar brep on floor i     
            NurbsSurface nurbsSurface = NurbsSurface.CreateNetworkSurface(planarBrep[0].Edges, 0, 0.0001, 0.0001, 0.0001, out int error); // make planar brep to nurbssurface
            List<Point3d> meshPoints = CreateGridOfPointsUV(nurbsSurface, u, v); //brep.Faces[0].ToNurbsSurface()

            // 4. Create nodes 
            smartMesh.Nodes = CreateNodes(meshPoints, smartMesh.nu, smartMesh.nv);

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
            for (double j = 0; j <= v; j++)
            {
                for (double k = 0; k <= u; k++)
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
        /// Create global nodes by assigning global id, coordinate, boundary condiditon in u and v direction
        /// </summary>
        /// <returns></returns>
        List<Node> CreateNodes(List<Point3d> meshPoints, int nu, int nv)
        {
            List<Node> nodes = new List<Node>();
            int uSequence = 0;
            int vSequence = 0;
            for (int i = 0; i < meshPoints.Count; i++)
            {
                bool BC_U = false;
                bool BC_V = false;

                // assign boundary condition
                if (uSequence == 0 | uSequence == nu - 1) { BC_U = true; } // assign BC u-dir
                if (vSequence == 0 | vSequence == nv - 1) { BC_V = true; } // assign BC v-dir

                Node node = new Node(i, meshPoints[i], BC_U, BC_V); // assign global id and cooridinates

                uSequence++;
                if (uSequence == nu)
                {
                    vSequence++;
                    uSequence = 0;
                }
                nodes.Add(node);
            }
            return nodes;
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
                return Properties.Resources.Icon_TargetLength;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("fb93bcbe-0eaf-456b-bbb1-83c95361ad51"); }
        }
    }
}
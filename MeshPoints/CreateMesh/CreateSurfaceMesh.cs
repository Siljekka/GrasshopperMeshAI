﻿using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry.Intersect;
using Rhino.Geometry.Collections;

namespace MeshPoints.CreateMesh
{
    public class CreateSurfaceMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateSurfaceMesh class.
        /// </summary>
        public CreateSurfaceMesh()
          : base("CreateSurfaceMesh", "surface",
              "Mesh list with flatten points for planar breps",
              "MyPlugIn", "Mesh")
        {
        }


        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface", "srf", "Surface", GH_ParamAccess.item);
            pManager.AddIntegerParameter("u", "u", "division in u direction", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("v", "v", "division in v direction", GH_ParamAccess.item, 4);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SmartMesh", "SmartMesh generated", GH_ParamAccess.item);
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
            List<Point3d> meshPoints = CreateGridOfPointsUV(brep.Faces[0].ToNurbsSurface(), u, v);

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
                return Properties.Resources.Icon_SurfaceMesh;
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
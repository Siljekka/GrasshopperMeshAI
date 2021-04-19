using Grasshopper.Kernel;
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
            pManager.AddGenericParameter("SurfaceMesh", "surface", "SurfaceMesh from given points", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "m", "Mesh (surface elements).", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Surface surface = null;
            int nu = 0;
            int nv = 0;
            DA.GetData(0, ref surface);
            DA.GetData(1, ref nu);
            DA.GetData(2, ref nv);

            #region Variables
            Mesh2D surfaceMesh = new Mesh2D();
            Mesh globalMesh = new Mesh();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            List<Point3d> meshPoints = new List<Point3d>();
            #endregion


            // 0. Check input OK.
            if (!surface.IsValid) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No valid surface input found."); return; } //todo: is this one needed?
            if (nu == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nu can not be zero."); return; }
            if (nv == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nv can not be zero."); return; }

            // 1. Generate grid of points on surface
            meshPoints = CreateGridOfPointsUV(nu, nv, surface.ToNurbsSurface());

            // 2. Create nodes and elements
            nodes = CreateNodes(meshPoints, nu, nv);
            elements = CreateQuadElements(nodes, nu, nv);

            // 3. Create global mesh
            globalMesh = CreateGlobalMesh(meshPoints, nu, nv);

            // 4. Add properties to SolidMesh
            surfaceMesh = new Mesh2D(nu+1, nv+1, nodes, elements, globalMesh);

            // 5. Check if brep can be interpret by Abaqus
            IsBrepCompatibleWithAbaqus(surfaceMesh);

            // Output
            DA.SetData(0, surfaceMesh);
            DA.SetData(1, surfaceMesh.mesh);
        }
        /// <summary>
        /// Check if mesh is compatible with Abaqus
        /// </summary>
        /// <returns> Nothing. Assign propertie to surfaceMesh. </returns>
        private void IsBrepCompatibleWithAbaqus(Mesh2D surfaceMesh)
        {
            BoundingBox bb = surfaceMesh.Elements[0].mesh.GetBoundingBox(false);
            Vector3d vec1 = bb.GetCorners()[1] - bb.GetCorners()[0];
            Vector3d vec2 = bb.GetCorners()[3] - bb.GetCorners()[0];
            Vector3d vector1 = Vector3d.CrossProduct(vec1, vec2);
            Vector3d vector2 = surfaceMesh.Elements[0].mesh.FaceNormals[0];
            Vector3d normal = Vector3d.CrossProduct(vector1, vector2);
            double angle = Vector3d.VectorAngle(vector1, vector2, normal);
            if (angle < Math.PI / 2) { surfaceMesh.inp = true; }
            else { surfaceMesh.inp = false; }
        }
        /// <summary>
        /// Makes grid of points in U and V direction
        /// </summary>
        /// <returns> List of points in U and V direction</returns>
        private List<Point3d> CreateGridOfPointsUV(int nu, int nv, NurbsSurface surface)
        {
            List<Point3d> pt = new List<Point3d>();

            var u = surface.Domain(0);
            var v = surface.Domain(1);

            double stepU = 1 / ((double)nu) * u.Length;
            double stepV = 1 / ((double)nv) * v.Length;
            BoundingBox bb = surface.GetBoundingBox(true);
            Vector3d vec1 = bb.GetCorners()[1] - bb.GetCorners()[0];
            Vector3d vec2 = bb.GetCorners()[3] - bb.GetCorners()[0];
            Vector3d vector1 = Vector3d.CrossProduct(vec1, vec2);
            Vector3d vector2 = surface.NormalAt(u.T1 * 0.5, v.T1 * 0.5);
            
            Vector3d normal = Vector3d.CrossProduct(vector1, vector2);
            double angle = Vector3d.VectorAngle(vector1, vector2, normal);

            if (angle < Math.PI / 2) 
            { 
                double pointU = 0;
                double pointV = 0;
                for (double j = 0; j <= nv; j++)
                {
                    for (double k = 0; k <= nu; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV));  // make point on surface
                        pointU = pointU + stepU;
                    }
                    pointV = pointV + stepV;
                    pointU = 0;
                }
            }
            else
            {
                double pointU = 0;
                double pointV = v.Length;
                for (double j = 0; j <= nv; j++)
                {
                    for (double k = 0; k <= nu; k++)
                    {
                        pt.Add(surface.PointAt(pointU, pointV)); // make point on surface
                        pointU = pointU + stepU;
                    }
                    pointV = pointV - stepV;
                    pointU = 0; 
                }
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
            nu = nu + 1;
            nv = nv + 1;
            int uSequence = 0;
            int vSequence = 0;
            for (int i = 0; i < meshPoints.Count; i++)
            {
                Node node = new Node(i, meshPoints[i]); // assign global id and cooridinates

                // assign boundary condition
                if (uSequence == 0 | uSequence == nu - 1) { node.BC_U = true; } // assign BC u-dir
                if (vSequence == 0 | vSequence == nv - 1) { node.BC_V = true; } // assign BC v-dir

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
        /// Creates quad elements by assigning id, local nodes and mesh.
        /// </summary>
        /// <returns></returns>
        List<Element> CreateQuadElements(List<Node> nodes, int nu, int nv)
        {
            Element e = new Element();
            List<Element> elements = new List<Element>();
            int uSequence = 0;
            int counter = 0;
            nu = nu + 1;
            nv = nv + 1;

            for (int i = 0; i < (nu - 1) * (nv - 1); i++) // loop elements
            {
                Node n1 = new Node(1, nodes[counter].GlobalId, nodes[counter].Coordinate, nodes[counter].BC_U, nodes[counter].BC_V);
                Node n2 = new Node(2, nodes[counter + 1].GlobalId, nodes[counter + 1].Coordinate, nodes[counter + 1].BC_U, nodes[counter + 1].BC_V);
                Node n3 = new Node(3, nodes[counter + nu + 1].GlobalId, nodes[counter + nu + 1].Coordinate, nodes[counter + nu + 1].BC_U, nodes[counter + nu + 1].BC_V);
                Node n4 = new Node(4, nodes[counter + nu].GlobalId, nodes[counter + nu].Coordinate, nodes[counter + nu].BC_U, nodes[counter + nu].BC_V);
                Mesh mesh = new Mesh();
                mesh.Vertices.Add(n1.Coordinate);
                mesh.Vertices.Add(n2.Coordinate);
                mesh.Vertices.Add(n3.Coordinate);
                mesh.Vertices.Add(n4.Coordinate);
                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh

                e = new Element(i, n1, n2, n3, n4, mesh);
                elements.Add(e); // add element to list of elements
                
                counter++;
                uSequence++;
                if (uSequence == (nu - 1)) // check if done with a v sequence
                {
                    counter++;
                    uSequence = 0; // new v sequence
                }
            }
            return elements;
        }

        /// <summary>
        /// Creates a global mesh for the geometry.
        /// </summary>
        /// <returns></returns>
        Mesh CreateGlobalMesh(List<Point3d> meshPoints, int nu, int nv)
        {
            Mesh globalMesh = new Mesh();
            int counter = 0;
            int uSequence = 0;
            nu = nu + 1;
            nv = nv + 1;

            foreach (Point3d pt in meshPoints)
            {
                globalMesh.Vertices.Add(pt);
            }
            for (int i = 0; i < (nu - 1) * (nv - 1); i++) // loop elements
            {
                globalMesh.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);

                counter++;
                uSequence++;
                if (uSequence == (nu - 1)) // check if done with a v sequence
                {
                    counter++;
                    uSequence = 0; // new v sequence
                }
            }
            globalMesh.Normals.ComputeNormals();  // todo: control if needed
            globalMesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            globalMesh.Compact(); // to ensure that it calculate

            return globalMesh;
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
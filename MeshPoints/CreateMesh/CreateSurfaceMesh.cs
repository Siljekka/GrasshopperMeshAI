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
            int nu = 0;
            int nv = 0;
            DA.GetData(0, ref brep);
            DA.GetData(1, ref nu);
            DA.GetData(2, ref nv);


            /* Todo:
             *  Fikse hvordan punkt blir generert. Gjør som i solidMesh med PointAt(0,0) som referansepunkt.
             *  Fiks sånn at overflate alltid har u og v i samme rekkefølge uavhengig av hvordan den tegnes.
             */

            // 1. Check input OK.
            if (!DA.GetData(0, ref brep)) return;
            if (nu == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nu can not be zero."); return; }
            if (nv == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "nv can not be zero."); return; }


            // 1. Assign properties to Geometry Class
            Geometry brepGeometry = new Geometry(brep, brep.Faces.ToList(), brep.Edges.ToList(), brep.Vertices.ToList());

            // 2. Generate grid of points on surface
            List<Point3d> meshPoints = CreateGridOfPointsUV(brep.Faces[0].ToNurbsSurface(), nu, nv);

            // 2. Create nodes and elements
            List<Node> nodes = CreateNodes(meshPoints, nu, nv);
            List<Element> elements = CreateQuadElements(nodes, nu, nv);

            // 3. Create global mesh
            Mesh globalMesh = CreateGlobalMesh(meshPoints, nu, nv);

            //5. Add properties to SolidMesh
            Mesh3D smartMesh = new Mesh3D(nu+1, nv+1, nodes, elements, globalMesh);
            smartMesh.Geometry = brepGeometry;

            // Output
            DA.SetData(0, smartMesh);
            DA.SetData(1, smartMesh.mesh);
        }
        

        /// <summary>
        /// Makes grid of points in U and V direction
        /// </summary>
        /// <returns> List of points in U and V direction</returns>
        private List<Point3d> CreateGridOfPointsUV(NurbsSurface surface, int nu, int nv)
        {
            List<Point3d> pt = new List<Point3d>();

            var u = surface.Domain(0);
            var v = surface.Domain(1);
            double stepU = u.Length / (double)nu;
            double stepV = v.Length / (double)nv;

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
        List<Element> CreateQuadElements(List<Node> nodes, int nu, int nv) // to do: erstatt med Mesh3D metode
        {
            List<Element> elements = new List<Element>();
            int uSequence = 0;
            int counter = 0;
            nu = nu + 1;
            nv = nv + 1;

            for (int i = 0; i < (nu - 1) * (nv - 1); i++) // loop elements
            {
                Mesh mesh = new Mesh();
                List<Node> elementNodes = new List<Node>();
                List<int> connectivity = new List<int>();
                connectivity.Add(counter);
                connectivity.Add(counter + 1);
                connectivity.Add(counter + nu + 1);
                connectivity.Add(counter + nu);

                foreach (int id in connectivity)
                {
                    elementNodes.Add(nodes[id]);
                    mesh.Vertices.Add(nodes[id].Coordinate);
                };

                Element element = new Element(i, elementNodes, connectivity);

                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
                element.mesh = mesh;

                elements.Add(element); // add element to list of elements
                
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
        Mesh CreateGlobalMesh(List<Point3d> meshPoints, int nu, int nv) // to do: erstatt 
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
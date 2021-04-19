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

namespace MeshPoints.QuadRemesh
{
    public class CreateUnstructuredSolidMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateUnstructuredSolidMesh class.
        /// </summary>
        public CreateUnstructuredSolidMesh()
          : base("CopyMeshToPlanes", "cm",
              "Copy a mesh to planes",
              "MyPlugIn", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("BottomFace", "", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("nw", "", "", GH_ParamAccess.item, 4);
            pManager.AddGenericParameter("Mesh", "mb", "Mesh as base", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("output", "", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep brep = new Brep();
            int bottomFace = 0;
            int nw = 0;
            Mesh2D mesh = new Mesh2D();

            DA.GetData(0, ref brep);
            DA.GetData(1, ref bottomFace);
            DA.GetData(2, ref nw);
            DA.GetData(3, ref mesh);

            if (!DA.GetData(0, ref brep)) return;
            if (!DA.GetData(1, ref bottomFace)) return;
            if (nw == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "nw = 0"); return; }

            #region Code
            // to do: add error if not a planar mesh

            // 1. Assign properties to SolidMesh
            //solidMesh.inp = true;

            // 2. Find Rails
            List<Curve> rails = FindRails(brep, bottomFace);

            // 3. Divide each brep edge in w direction (rail) into nw points.
            DataTree<Point3d> railPoints = DivideRailIntoNwPoints(rails, brep.Faces[bottomFace], nw);
            
            // 4. Create Planes
            List<Plane> planes = GetPlanes(railPoints);

            // 5. Create elements

            // copy mesh to planes
            int globalNodeIndexAddition = mesh.Nodes.Count;
            List<Mesh2D> meshOnPlanes = new List<Mesh2D>();
            meshOnPlanes.Add(mesh);
            Plane basePlane = planes[0];

            for (int i = 1; i < planes.Count; i++)
            {
                Mesh2D meshToTransform = new Mesh2D();
                Transform tranformation = new Transform();
                meshToTransform = mesh;
                globalNodeIndexAddition = globalNodeIndexAddition * i;
                Point3d pointToTransform = new Point3d();
                tranformation = Transform.PlaneToPlane(basePlane, planes[i]);

                meshToTransform.mesh.Transform(tranformation); // transform mesh, to do: check if needed

                // transform global nodes
                foreach (Node nodeToTransform in meshToTransform.Nodes)
                {
                    nodeToTransform.GlobalId = nodeToTransform.GlobalId + globalNodeIndexAddition;
                    pointToTransform = nodeToTransform.Coordinate;
                    pointToTransform.Transform(tranformation);
                }

                // transform elements
                for (int j = 0; j < meshToTransform.Elements.Count; j++)
                {
                    Element elementToTransform = meshToTransform.Elements[j];
                    elementToTransform.mesh.Transform(tranformation); // fix mesh face of element, to do: check if needed

                    List<Node> nodesOfElementToTransform = new List<Node>() { elementToTransform.Node1, elementToTransform.Node2, elementToTransform.Node3 };
                    if (elementToTransform.IsQuad) { nodesOfElementToTransform.Add(elementToTransform.Node4); }

                    foreach (Node nodeToTransform in nodesOfElementToTransform)
                    {
                        nodeToTransform.GlobalId = nodeToTransform.GlobalId + globalNodeIndexAddition;
                        pointToTransform = nodeToTransform.Coordinate;
                        pointToTransform.Transform(tranformation);
                    }
                }
                meshOnPlanes.Add(meshToTransform); // add to list of mesh
            }
            /*
            // from SurfaceMesh to SolidMesh, Old method

            int elementIndex = 0;
            List<Element> solidElements = new List<Element>();
            for (int i = 0; i < meshOnPlanes.Count - 1; i++)
            {
                Mesh2D bottomMesh = meshOnPlanes[i];
                Mesh2D topMesh = meshOnPlanes[i+1];

                for (int j = 0; j < bottomMesh.Elements.Count; j++) // to do: make general for triangle mesh as well
                {
                    Mesh m = new Mesh();
                    Element e = new Element();
                    Element bottomElement = bottomMesh.Elements[j];
                    Element topElement = topMesh.Elements[j];
                    e.IsCube = true;
                    e.Id = elementIndex;

                    Node n1 = new Node();
                    n1 = bottomElement.Node1;
                    e.Node1 = n1;

                    Node n2 = new Node();
                    n2 = bottomElement.Node2;
                    e.Node2 = n2;

                    Node n3 = new Node();
                    n3 = bottomElement.Node3;
                    e.Node3 = n3;

                    Node n4 = new Node();
                    n4 = bottomElement.Node1;
                    e.Node4 = n4;

                    Node n5 = new Node();
                    n5 = topElement.Node1; n5.LocalId = 5;
                    e.Node5 = n5;

                    Node n6 = new Node();
                    n6 = topElement.Node2; n6.LocalId = 6;
                    e.Node6 = n6;

                    Node n7 = new Node();
                    n7 = topElement.Node3; n7.LocalId = 7;
                    e.Node7 = n7;

                    Node n8 = new Node();
                    n8 = topElement.Node4; n8.LocalId = 8;
                    e.Node8 = n8;

                    m.Vertices.Add(e.Node1.Coordinate); //0
                    m.Vertices.Add(e.Node2.Coordinate); //1
                    m.Vertices.Add(e.Node3.Coordinate); //2
                    m.Vertices.Add(e.Node4.Coordinate); //3
                    m.Vertices.Add(e.Node5.Coordinate); //4
                    m.Vertices.Add(e.Node6.Coordinate); //5
                    m.Vertices.Add(e.Node7.Coordinate); //6
                    m.Vertices.Add(e.Node8.Coordinate); //7

                    m.Faces.AddFace(0, 1, 5, 4);
                    m.Faces.AddFace(1, 2, 6, 5);
                    m.Faces.AddFace(2, 3, 7, 6);
                    m.Faces.AddFace(3, 0, 4, 7);
                    m.Faces.AddFace(0, 1, 2, 3);
                    m.Faces.AddFace(4, 5, 6, 7);

                    m.Normals.ComputeNormals();  //Control if needed
                    m.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                    m.Compact(); //to ensure that it calculate
                    e.mesh = m;

                    //add element and mesh to element list
                    solidElements.Add(e);

                    elementIndex++;
                }

                //  create global mesh
                Mesh allMesh = new Mesh();
                foreach (Element el in solidElements)
                {
                    allMesh.Append(el.mesh);
                }
                allMesh.Weld(0.01);
            }*/

            #endregion

            DA.SetDataList(0, meshOnPlanes);

        }


        #region Methods
        private List<Curve> FindRails(Brep brep, int bottomFaceIndex)
        {
            // Find top and bottom edge
            List<BrepFace> brepFace = brep.Faces.ToList();

            List<int> indexAdjecentFaces = (brepFace[bottomFaceIndex].AdjacentFaces()).ToList();
            List<int> indexAdjecentEdges = (brepFace[bottomFaceIndex].AdjacentEdges()).ToList();
            indexAdjecentFaces.Add(bottomFaceIndex);
            for (int i = 0; i < brepFace.Count; i++)
            {
                if (!indexAdjecentFaces.Contains(brepFace.IndexOf(brepFace[i])))
                {
                    BrepFace brepBottomFace = brepFace[bottomFaceIndex];
                    BrepFace brepTopFace = brepFace[i]; // top face
                    indexAdjecentEdges.AddRange(brepBottomFace.AdjacentEdges());
                    indexAdjecentEdges.AddRange(brepTopFace.AdjacentEdges());
                    continue;
                }
            }

            // Find rails
            List<BrepEdge> brepEdges = brep.Edges.ToList();
            List<Curve> rails = new List<Curve>(brepEdges);
            foreach (int index in indexAdjecentEdges) { rails.Remove(brepEdges[index]); }

            #region Old Code
            /*BrepEdgeList brepEdges = brep.Edges;
            List<Curve> rails = new List<Curve>();

            foreach (BrepEdge edge in brepEdges) // check if node is on edge
            {
                List<Point3d> edgePoints = new List<Point3d> { edge.StartVertex.Location, edge.EndVertex.Location };
                bool pointOnFace4 = false;
                bool pointOnFace5 = false;
                foreach (Point3d point in edgePoints)
                {
                    brep.Faces[bottomFace].ClosestPoint(point, out double PointOnCurveUFace4, out double PointOnCurveVFace4);
                    brep.Faces[topFace].ClosestPoint(point, out double PointOnCurveUFace5, out double PointOnCurveVFace5);
                    Point3d testPointFace4 = brep.Faces[bottomFace].PointAt(PointOnCurveUFace4, PointOnCurveVFace4);  // make test point
                    Point3d testPointFace5 = brep.Faces[topFace].PointAt(PointOnCurveUFace5, PointOnCurveVFace5);  // make test point
                    double distanceToFace4 = testPointFace4.DistanceTo(point); // calculate distance between testPoint and node
                    double distanceToFace5 = testPointFace5.DistanceTo(point); // calculate distance between testPoint and node
                    if ((distanceToFace4 <= 0.0001 & distanceToFace4 >= -0.0001)) // if distance = 0: node is on edge
                    {
                        pointOnFace4 = true;
                    }
                    else if ((distanceToFace5 <= 0.0001 & distanceToFace5 >= -0.0001))
                    {
                        pointOnFace5 = true;
                    }
                }
                if (pointOnFace4 & pointOnFace5)
                {
                    rails.Add(edge);  //get edge1 of brep = rail 1
                }
            }*/
            #endregion
            return rails;
        }
        private DataTree<Point3d> DivideRailIntoNwPoints(List<Curve> rails, BrepFace brepBottomFace, int nw)
        {
            DataTree<Point3d> railPoints = new DataTree<Point3d>();

            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(nw, true, out Point3d[] pt);
                List<Point3d> point = pt.ToList();
                
                brepBottomFace.ClosestPoint(point[0], out double PointOnCurveU, out double PointOnCurveV);
                Point3d testPoint = brepBottomFace.PointAt(PointOnCurveU, PointOnCurveV);
                Vector3d distanceToFace = testPoint - point[0];
                if (distanceToFace.Length > 0.001) { point.Reverse(); }

                for (int j = 0; j < point.Count; j++)
                {
                    railPoints.Add(point[j], new GH_Path(j)); //tree with nw points on each rail. Branch: floor
                }
            }

            // Check if the rails must be re-oredered to generate elements with nodes counting ccw
            Curve testCurve = Curve.CreateControlPointCurve(railPoints.Branch(0), 1);
            Vector3d direction = railPoints.Branch(nw)[0] - railPoints.Branch(0)[0];
            string curveOrientation = testCurve.ClosedCurveOrientation(direction).ToString();
            if (curveOrientation == "Clockwise")
            {
                for (int i = 0; i < railPoints.BranchCount; i++)
                {
                    railPoints.Branch(i).Reverse();
                }
            }
            #region Old Code:
            /*
            DataTree<Point3d> railPoints = new DataTree<Point3d>();
            if (rails.Count == 0) { return null; }

            //Divide each rail into nw points.
            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(nw, true, out Point3d[] nwPt);  //divide each rail in nw number of points
                List<Point3d> nwPoints = nwPt.ToList();
                for (int j = 0; j < nwPoints.Count; j++)
                {
                    railPoints.Add(nwPoints[j], new GH_Path(j)); //tree with nw points on each rail. Branch: floor
                }
            }

            // Check if the rails must be re-oredered to generate elements with nodes counting ccw
            Curve testCurve = Curve.CreateControlPointCurve(railPoints.Branch(0), 1);
            Vector3d direction = railPoints.Branch(nw)[0] - railPoints.Branch(0)[0];
            string curveOrientation = testCurve.ClosedCurveOrientation(direction).ToString();
            if (curveOrientation == "Clockwise")
            {
                rails.Reverse();
                for (int i = 0; i < railPoints.BranchCount; i++)
                {
                    railPoints.Branch(i).Reverse();
                }
            }
            */
            #endregion
            return railPoints;
        }
        private List<Plane> GetPlanes(DataTree<Point3d> railPoints)
        {
            if (railPoints == null) { return null; }
            List<Plane> planes = new List<Plane>();

            for (int i = 0; i < railPoints.BranchCount; i++)
            {
                Vector3d vec1 = railPoints.Branch(i)[1] - railPoints.Branch(i)[0];
                Vector3d vec2 = railPoints.Branch(i)[3] - railPoints.Branch(i)[0];
                Vector3d normal = Vector3d.CrossProduct(vec1, vec2);
                Plane plane = new Plane(railPoints.Branch(i)[0], vec1, vec2);
                planes.Add(plane);
            }
            return planes;
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
            get { return new Guid("edc7ac7c-86c7-4858-83ec-30a34dd92fe5"); }
        }
    }
}
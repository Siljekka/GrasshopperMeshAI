﻿using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel.Data;

namespace MeshPoints.CreateMesh
{
    public class SweepSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the SweepSmartMesh class.
        /// </summary>
        public SweepSmartMesh()
          : base("Sweep SmartMesh", "SweepSM",
              "Sweeps a surface SmartMesh.",
              "SmartMesh", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "brep", "Brep to be meshed with sweeping.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Index", "i", "Index of the bottom face of the brep.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("w", "w","Number element in w-direction.", GH_ParamAccess.item, 4);
            pManager.AddGenericParameter("SmartMesh", "SM", "Reference SmartMesh to sweep.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "mesh", "Mesh (hexahedral-elements).", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep brep = new Brep();
            int bottomFace = 0;
            int w = 0;
            SmartMesh refMesh = new SmartMesh();

            DA.GetData(0, ref brep);
            DA.GetData(1, ref bottomFace);
            DA.GetData(2, ref w);
            DA.GetData(3, ref refMesh);

            if (!DA.GetData(0, ref brep)) return;
            if (!DA.GetData(1, ref bottomFace)) return;
            if (w == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "w = 0"); return; }

            // to do: add error if not a planar mesh

            // 1. Find Rails
            List<Curve> rails = FindRails(brep, bottomFace);

            // 2. Divide each brep edge in w direction (rail) into w points.
            DataTree<Point3d> railPoints = DivideRailIntoWPoints(rails, brep.Faces[bottomFace], w);
            
            // 3. Create Planes
            List<Plane> planes = GetPlanes(railPoints);

            // 4. Create nodes
            List<Node> nodes = CreateNodes(refMesh, planes);

            // 5. Create elements
            List<Element> elements = CreateElements(refMesh, w, nodes);

            // 6. Create geometry info
            Geometry brepGeometry = new Geometry(brep, bottomFace);

            // 7. Create SmartMesh
            SmartMesh solidMesh = new SmartMesh(nodes, elements, "Solid");
            solidMesh.Geometry = brepGeometry;

            DA.SetData(0, solidMesh);
            DA.SetData(1, solidMesh.Mesh);
        }

        private List<Curve> FindRails(Brep brep, int bottomFaceIndex)
        {
            // Find top and bottom edge
            List<BrepFace> brepFace = brep.Faces.ToList();

            List<int> indexAdjecentFaces = (brepFace[bottomFaceIndex].AdjacentFaces()).ToList(); // get faces adjacent to bottom face
            List<int> indexAdjecentEdges = (brepFace[bottomFaceIndex].AdjacentEdges()).ToList(); // get edges adjacent to bottom face
            indexAdjecentFaces.Add(bottomFaceIndex); // add bottom face to list
            for (int i = 0; i < brepFace.Count; i++)
            {
                if (!indexAdjecentFaces.Contains(brepFace.IndexOf(brepFace[i])))
                {
                    indexAdjecentEdges.AddRange(brepFace[i].AdjacentEdges()); //top face edges
                    break;
                }
            }

            // Find rails
            List<BrepEdge> brepEdges = brep.Edges.ToList();
            List<Curve> rails = new List<Curve>(brepEdges);
            foreach (int index in indexAdjecentEdges) { rails.Remove(brepEdges[index]); }
            return rails;
        }
        private DataTree<Point3d> DivideRailIntoWPoints(List<Curve> rails, BrepFace brepBottomFace, int w)
        {
            DataTree<Point3d> railPoints = new DataTree<Point3d>();

            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(w, true, out Point3d[] pt);
                List<Point3d> point = pt.ToList();
                
                brepBottomFace.ClosestPoint(point[0], out double PointOnCurveU, out double PointOnCurveV);
                Point3d testPoint = brepBottomFace.PointAt(PointOnCurveU, PointOnCurveV);
                Vector3d distanceToFace = testPoint - point[0];
                if (distanceToFace.Length > 0.001) { point.Reverse(); }

                for (int j = 0; j < point.Count; j++)
                {
                    railPoints.Add(point[j], new GH_Path(j)); //tree with w points on each rail. Branch: floor
                }
            }

            // Check if the rails must be re-oredered to generate elements with nodes counting ccw
            if (railPoints.Branch(0).Count > 1)
            {
                Curve testCurve = Curve.CreateControlPointCurve(railPoints.Branch(0));
                Vector3d direction = railPoints.Branch(w)[0] - railPoints.Branch(0)[0];
                string curveOrientation = testCurve.ClosedCurveOrientation(direction).ToString();
                if (curveOrientation == "Clockwise")
                {
                    for (int i = 0; i < railPoints.BranchCount; i++)
                    {
                        railPoints.Branch(i).Reverse();
                    }
                }
            }
            else if (railPoints.Branch(0).Count == 1)
            {
                // Create rails
                DataTree<Point3d> railPointsTemp = new DataTree<Point3d>();

                Point3d pt1 = brepBottomFace.PointAt(0,0);
                Point3d pt4 = brepBottomFace.PointAt(1, 0);
                Point3d pt3 = brepBottomFace.PointAt(1, 1);
                Point3d pt2 = brepBottomFace.PointAt(0, 1);
                List<Point3d> newBasePt = new List<Point3d>() { pt1, pt2, pt3, pt4 };
                for (int i = 0; i < railPoints.BranchCount; i++)
                {
                    Point3d pt = railPoints.Branch(i)[0];
                    foreach (Point3d basePt in newBasePt)
                    {
                        Vector3d move = basePt - pt;
                        Point3d newPt = new Point3d(pt.X + move.X, pt.Y + move.Y, pt.Z);
                        railPointsTemp.Add(newPt, new GH_Path(i));
                    }
                }
                railPoints = railPointsTemp;
            }
            return railPoints;
        }
        private List<Plane> GetPlanes(DataTree<Point3d> railPoints)
        {
            List<Plane> planes = new List<Plane>();
            if (railPoints == null) { return null; }

            for (int i = 0; i < railPoints.BranchCount; i++)
            {
                Vector3d vec1 = railPoints.Branch(i)[1] - railPoints.Branch(i)[0];
                Vector3d vec2 = railPoints.Branch(i)[3] - railPoints.Branch(i)[0];
                Plane plane = new Plane(railPoints.Branch(i)[0], vec1, vec2);
                planes.Add(plane);
            }
            return planes;
        }
        private List<Node> CreateNodes(SmartMesh refMesh, List<Plane> planes)
        {
            List<Node> nodes = new List<Node>();
            List<Node> nodesTest = new List<Node>();
            int numNodesInPlane = refMesh.Nodes.Count;
            int w = planes.Count - 1;

            Plane basePlane = planes[0];
            SmartMesh meshToTransform = refMesh;
            List<Node> nodetest2 = meshToTransform.Nodes;
            for (int i = 0; i < planes.Count; i++)
            {
                if (i > 0) { basePlane = planes[i - 1]; }
                Transform tranformation = Transform.PlaneToPlane(basePlane, planes[i]);

                // Creating global nodes
                for (int j = 0; j < nodetest2.Count; j++)
                {
                    Node nodeToTransform = nodetest2[j];
                    Point3d pointToTransform = nodeToTransform.Coordinate;
                    pointToTransform.Transform(tranformation);
                    int idJump = numNodesInPlane;
                    if (i == 0 | i == w) { nodeToTransform.BC_W = true; } else { nodeToTransform.BC_W = false; } // assign BCW
                    if (i == 0) { idJump = 0; }
                    Node n = new Node(nodeToTransform.GlobalId + idJump, pointToTransform, nodeToTransform.BC_U, nodeToTransform.BC_V, nodeToTransform.BC_W);
                    nodes.Add(n);
                }
                nodetest2 = new List<Node>(nodes);
                nodesTest.AddRange(nodes);
                nodes.Clear();
            }
            return nodesTest;
        }
        private List<Element> CreateElements(SmartMesh refMesh, int w, List<Node> nodes)
        {
            List<Element> elements = new List<Element>();
            int numElementsInPlane = refMesh.Elements.Count;
            int numNodesInPlane = refMesh.Nodes.Count;

            int elemId = 0;
            for (int i = 0; i < w; i++)  // loop levels
            {
                for (int j = 0; j < numElementsInPlane; j++) // loop elements in a level
                {
                    List<Node> elementNodes = new List<Node>();
                    List<int> connectivity = new List<int>();
                    Element refElement = refMesh.Elements[j];

                    if (refMesh.Elements[0].Type == "Triangle")
                    {
                        connectivity.Add(refElement.Connectivity[0] + numNodesInPlane * i);
                        connectivity.Add(refElement.Connectivity[1] + numNodesInPlane * i);
                        connectivity.Add(refElement.Connectivity[2] + numNodesInPlane * i);
                        connectivity.Add(refElement.Connectivity[0] + numNodesInPlane * (i + 1));
                        connectivity.Add(refElement.Connectivity[1] + numNodesInPlane * (i + 1));
                        connectivity.Add(refElement.Connectivity[2] + numNodesInPlane * (i + 1));
                    }
                    else 
                    {
                        connectivity.Add(refElement.Connectivity[0] + numNodesInPlane * i);
                        connectivity.Add(refElement.Connectivity[1] + numNodesInPlane * i);
                        connectivity.Add(refElement.Connectivity[2] + numNodesInPlane * i);
                        connectivity.Add(refElement.Connectivity[3] + numNodesInPlane * i);
                        connectivity.Add(refElement.Connectivity[0] + numNodesInPlane * (i + 1));
                        connectivity.Add(refElement.Connectivity[1] + numNodesInPlane * (i + 1));
                        connectivity.Add(refElement.Connectivity[2] + numNodesInPlane * (i + 1));
                        connectivity.Add(refElement.Connectivity[3] + numNodesInPlane * (i + 1));

                    }

                    foreach (int id in connectivity)
                    {
                        elementNodes.Add(nodes[id]);
                    }

                    Element element = new Element(elemId, elementNodes, connectivity);

                    //add element and mesh to element list
                    elements.Add(element);
                    elemId++;
                }
            }
            return elements;
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
                return Properties.Resources.Icon_SweepSolid;
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
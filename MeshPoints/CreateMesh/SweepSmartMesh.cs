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
    public class SweepSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the SweepSmartMesh class.
        /// </summary>
        public SweepSmartMesh()
          : base("Sweep SmartMesh", "ss",
              "Sweep a referance surface SmartMesh",
              "SmartMesh", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "Brep to mesh with sweeping", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("BottomFace", "index", "Index of bottomFace", GH_ParamAccess.item);
            pManager.AddIntegerParameter("w", "w","Number element in w direction", GH_ParamAccess.item, 4);
            pManager.AddGenericParameter("SmartMesh", "mb", "Reference SmartMesh to sweep", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "sm", "SmartMesh from sweeping", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "m", "Mesh from sweeping", GH_ParamAccess.item);
            pManager.AddGenericParameter("test list", "sm", "SmartMesh from sweeping", GH_ParamAccess.list);

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

            // 6. Create global mesh
            Mesh globalMesh = CreateGlobalMeshWithWeld(elements);

            // 7. Create SmartMesh
            SmartMesh solidMesh = new SmartMesh(nodes, elements, globalMesh, "Solid");

            DA.SetData(0, solidMesh);
            DA.SetData(1, solidMesh.Mesh);
            DA.SetDataList(2, elements); // temp
        }

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
            Curve testCurve = Curve.CreateControlPointCurve(railPoints.Branch(0), 1);
            Vector3d direction = railPoints.Branch(w)[0] - railPoints.Branch(0)[0];
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

            //Divide each rail into w points.
            for (int i = 0; i < rails.Count; i++)
            {
                rails[i].DivideByCount(w, true, out Point3d[] wPt);  //divide each rail in w number of points
                List<Point3d> wPoints = wPt.ToList();
                for (int j = 0; j < wPoints.Count; j++)
                {
                    railPoints.Add(wPoints[j], new GH_Path(j)); //tree with w points on each rail. Branch: floor
                }
            }

            // Check if the rails must be re-oredered to generate elements with nodes counting ccw
            Curve testCurve = Curve.CreateControlPointCurve(railPoints.Branch(0), 1);
            Vector3d direction = railPoints.Branch(w)[0] - railPoints.Branch(0)[0];
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
                    Node nodeToTransform = nodetest2[j];//meshToTransform.Nodes[j];
                    Point3d pointToTransform = nodeToTransform.Coordinate;
                    pointToTransform.Transform(tranformation);

                    if (i == 0 | i == w) { nodeToTransform.BC_W = true; } else { nodeToTransform.BC_W = false; } // assign BCW
                    Node n = new Node(nodeToTransform.GlobalId + numNodesInPlane * i, pointToTransform, nodeToTransform.BC_U, nodeToTransform.BC_V, nodeToTransform.BC_W);
                    nodes.Add(n);
                }
                nodetest2 = new List<Node>(nodes);
                nodesTest.AddRange(nodes);
                nodes.Clear();
            }
            return nodes;
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

                    connectivity.Add(refElement.Connectivity[0] + numNodesInPlane * i);
                    connectivity.Add(refElement.Connectivity[1] + numNodesInPlane * i);
                    connectivity.Add(refElement.Connectivity[2] + numNodesInPlane * i);
                    connectivity.Add(refElement.Connectivity[3] + numNodesInPlane * i);
                    connectivity.Add(refElement.Connectivity[0] + numNodesInPlane * (i+1));
                    connectivity.Add(refElement.Connectivity[1] + numNodesInPlane * (i+1));
                    connectivity.Add(refElement.Connectivity[2] + numNodesInPlane * (i+1));
                    connectivity.Add(refElement.Connectivity[3] + numNodesInPlane * (i+1));

                    foreach (int id in connectivity)
                    {
                        elementNodes.Add(nodes[id]);
                    }

                    Element element = new Element(elemId, elementNodes, connectivity);

                    // create local mesh
                    Mesh localMesh = new Mesh();
                    foreach (Node node in elementNodes)
                    {
                        localMesh.Vertices.Add(node.Coordinate); 
                    }
                    localMesh.Faces.AddFace(0, 1, 5, 4);
                    localMesh.Faces.AddFace(1, 2, 6, 5);
                    localMesh.Faces.AddFace(2, 3, 7, 6);
                    localMesh.Faces.AddFace(3, 0, 4, 7);
                    localMesh.Faces.AddFace(0, 1, 2, 3);
                    localMesh.Faces.AddFace(4, 5, 6, 7);

                    localMesh.Normals.ComputeNormals();
                    localMesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
                    localMesh.Compact(); // to ensure that it calculate
                    element.Mesh = localMesh;

                    //add element and mesh to element list
                    elements.Add(element);
                    elemId++;
                }
            }
            return elements;
        }
        private Mesh CreateGlobalMeshWithWeld(List<Element> elements)
        {
            Mesh globalMesh = new Mesh();
            foreach (Element element in elements)
            {
                globalMesh.Append(element.Mesh);
            }
            globalMesh.Weld(0.01);

            globalMesh.Normals.ComputeNormals();
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
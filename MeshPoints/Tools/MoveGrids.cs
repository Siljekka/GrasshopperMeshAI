using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;

namespace MeshPoints.Tools
{
    public class MoveGrids : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MoveGrids class.
        /// </summary>
        public MoveGrids()
          : base("Move Grids", "mg",
              "Move mesh grids",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "sm", "Input a SmartMesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Gene pool for translation in u direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("v genes", "qp", "Gene pool for translation in v direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("w genes", "qp", "Gene pool for translation in w direction", GH_ParamAccess.list);
            pManager[3].Optional = true; // if solid
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "sm", "Updated mesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "m", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            SmartMesh oldMesh = new SmartMesh();
            List<double> genesU = new List<double>();
            List<double> genesV = new List<double>();
            List<double> genesW = new List<double>();

            DA.GetData(0, ref oldMesh);
            DA.GetDataList(1, genesU);
            DA.GetDataList(2, genesV);
            DA.GetDataList(3, genesW);

            // Variables
            SmartMesh newMesh = new SmartMesh();

            // 1. Write error if wrong input
            if (!DA.GetData(0, ref oldMesh)) return;

            if (oldMesh.Type == "Solid" & !DA.GetDataList(3, genesW)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "For solid elements, must have input GenesW."); return; }
            if (genesU.Count < (oldMesh.nu - 2) * oldMesh.nw) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase u genes."); return; }
            if (genesV.Count < (oldMesh.nv - 2) * oldMesh.nw) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase v genes."); return; }
            if (oldMesh.Type == "Solid" & (genesW.Count < (oldMesh.nw - 2) * (oldMesh.nu - 2) * (oldMesh.nv - 2) )) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase w genes."); return; }
            if (oldMesh.nu == 0 | oldMesh.nv == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Do not support SmartMesh made as unstructured."); return; }

            // 2. Inherit properties from old mesh
            newMesh.nu = oldMesh.nu;
            newMesh.nv = oldMesh.nv;
            newMesh.nw = oldMesh.nw;
            newMesh.Type = oldMesh.Type;
            newMesh.Geometry = oldMesh.Geometry;
            Brep brep = oldMesh.Geometry.Brep;

            // 3. Create new nodes
            List<Point3d> newPoints = new List<Point3d>();
            double genU = 0;
            double genV = 0;
            double genW = 0;

            for (int k = 0; k < oldMesh.nw; k++)
            {
                if (k == 0 | k == oldMesh.nu - 1){ genW = 0;} // if edge grids in w dir
                else { genW = genesW[k-1]; }

                for (int j = 0; j < oldMesh.nv; j++) 
                {
                    if (j == 0 | j == oldMesh.nv - 1) { genV = 0; } // if edge grids in v dir
                    else { genV = genesV[j-1]; }

                    for (int i = 0; i < oldMesh.nu; i++)
                    {
                        if (i == 0 | i == oldMesh.nu - 1) { genU = 0; } // if edge grids in u dir
                        else { genU = genesU[i - 1]; }

                        int nodeIndex = i + j * oldMesh.nv + k * oldMesh.nu * oldMesh.nv;

                        Tuple<bool, BrepFace> pointFace = PointOnFace(oldMesh.Nodes[nodeIndex], brep); // Item1: IsOnFace, Item2: face. Silje: flytte dette inn i Node klasse? Og kall på fra GetNewCoord
                        Tuple<bool, BrepEdge> pointEdge = PointOnEdge(oldMesh.Nodes[nodeIndex], brep); // Item1: IsOnEdge, Item2: edge. Silje: flytte dette inn i Node klasse? Og kall på fra GetNewCoord
                        Point3d newPoint = GetNewCoordinateOfNode(nodeIndex, pointFace, pointEdge, oldMesh, genU, genV, genW);
                        newPoints.Add(newPoint);
                    }
                }
            }

            // 4. Set new nodes and elements
            newMesh.CreateNodes(newPoints, oldMesh.nu-1, oldMesh.nv-1, oldMesh.nw-1);
            if (newMesh.Type == "Surface")
            {
                newMesh.CreateQuadElements();
            }
            else
            {
                newMesh.CreateHexElements();
            }

            //4. Set new mesh 
            newMesh.CreateMesh();

            // Output
            DA.SetData(0, newMesh);
            DA.SetData(1, newMesh.Mesh);
        }
        #region Methods

        /// <summary>
        /// Check if point is on face.
        /// </summary>
        /// <returns> Returns true and face if point is on face. False and null face if point is not on face.</returns>
        private Tuple<bool, BrepFace> PointOnFace(Node node, Brep brep)
        {
            bool IsOnFace = false;
            BrepFace face = null;
            BrepFaceList brepFace = brep.Faces;

            foreach (BrepFace bFace in brepFace) // check if node is on face
            {
                if (node.Type == "Face")
                {
                    IsOnFace = node.IsOnFace(bFace);
                    if (IsOnFace)
                    {
                        face = bFace;
                        return new Tuple<bool, BrepFace>(IsOnFace, face);
                    }
                }
            }
            return new Tuple<bool, BrepFace>(IsOnFace, face);
        }

        /// <summary>
        /// Check if point is on edge.
        /// </summary>
        /// <returns> Returns true and edge if point is on edge. False and null edge if point is not on edge.</returns>
        private Tuple<bool, BrepEdge> PointOnEdge(Node node, Brep brep)
        {
            bool IsOnEdge = false;
            BrepEdge edge = null;
            BrepEdgeList brepEdge = brep.Edges;

            foreach (BrepEdge bEdge in brepEdge) // check if node is on edge
            {
                if (node.Type == "Edge")
                {
                    IsOnEdge = node.IsOnEdge(bEdge);
                    if (IsOnEdge)
                    {
                        edge = bEdge;
                        return new Tuple<bool, BrepEdge>(IsOnEdge, edge);
                    }
                }
            }
            return new Tuple<bool, BrepEdge>(IsOnEdge, edge);
        }

        /// <summary>
        /// Move the old node in allowable directions.
        /// </summary>
        /// <returns> Returns coordinates of moved node.</returns>
        private Point3d GetNewCoordinateOfNode(int i, Tuple<bool, BrepFace> pointFace, Tuple<bool, BrepEdge> pointEdge, SmartMesh mesh, double genU, double genV, double genW)
        {
            Point3d movedNode = new Point3d();
            bool IsOnEdge = pointEdge.Item1;
            bool IsOnFace = pointFace.Item1;
            BrepEdge edge = pointEdge.Item2;
            BrepFace face = pointFace.Item2;

            Vector3d translationVectorU = Vector3d.Zero;
            Vector3d translationVectorV = Vector3d.Zero;
            Vector3d translationVectorW = Vector3d.Zero;

            // Translation in x direction
            // 1. if: Node not restrained in x direction and gen positive.
            // 2. if: Node not restrained in x direction and gen negative.
            // 3. if: Node restrained in x direction.
            // Note: if point is on edge not restrained in x direction - meshPoint is made

            if (genU > 0 & !mesh.Nodes[i].BC_U) // 1. if
            {
                translationVectorU = 0.5 * (mesh.Nodes[i + 1].Coordinate - mesh.Nodes[i].Coordinate) * genU; // make vector translating node in U-direction
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genU, i, i + 1); return movedNode; } // make meshPoint
            }
            else if (genU < 0 & !mesh.Nodes[i].BC_U)  // 2. if
            {
                translationVectorU = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - 1].Coordinate) * genU;
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genU, i, i - 1); return movedNode; } // make meshPoint
            }
            else { translationVectorU = translationVectorU * 0; }  // 3. if


            if (genV > 0 & !mesh.Nodes[i].BC_V) // 1. if
            {
                translationVectorV = 0.5 * (mesh.Nodes[i + mesh.nu].Coordinate - mesh.Nodes[i].Coordinate) * genV;
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genV, i, i + mesh.nu); return movedNode; } // make meshPoint
            }
            else if (genV < 0 & !mesh.Nodes[i].BC_V) // 2. if
            {
                translationVectorV = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - mesh.nu].Coordinate) * genV;
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genV, i, i - mesh.nu); return movedNode; } // make meshPoint
            }
            else { translationVectorV = translationVectorV * 0; } // 3. if


            if (mesh.Type == "Solid")
            {

                if (genW > 0 & !mesh.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (mesh.Nodes[i + (mesh.nu) * (mesh.nv)].Coordinate - mesh.Nodes[i].Coordinate) * genW;
                    if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genW, i, i + (mesh.nu) * (mesh.nv)); return movedNode; } // make meshPoint
                }
                else if (genW < 0 & !mesh.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - (mesh.nu) * (mesh.nv)].Coordinate) * genW;
                    if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genW, i, i - (mesh.nu) * (mesh.nv)); return movedNode; } // make meshPoint
                }
                else { translationVectorW = translationVectorW * 0; } // 3. if                            
            }

            // 4. if: Make movedNode if node is on face or inside brep (if on edge, movedNode already made).
            double overlapTolerance = 0.99; // ensure no collision of vertices, reduce number to avoid "the look of triangles".
            movedNode = new Point3d
                (
                mesh.Nodes[i].Coordinate.X + (translationVectorU.X + translationVectorV.X + translationVectorW.X) * overlapTolerance,
                mesh.Nodes[i].Coordinate.Y + (translationVectorU.Y + translationVectorV.Y + translationVectorW.Y) * overlapTolerance,
                mesh.Nodes[i].Coordinate.Z + (translationVectorU.Z + translationVectorV.Z + translationVectorW.Z) * overlapTolerance
                );

            if (IsOnFace) // If node is on face: ensure it stays on face
            {
                Brep srf = face.DuplicateFace(false);
                movedNode = srf.ClosestPoint(movedNode); // "Project" meshPoint to surface.
            }
            return movedNode;
        }


        /// <summary>
        /// Make new node if point is on edge.
        /// </summary>
        /// <returns> Returns coordinates of moved node on edge.</returns>
        private Point3d EdgeNode(BrepEdge edge, SmartMesh mesh, double genes, int start, int stop)
        {
            Point3d movedNode = new Point3d();
            Curve edgeCurve1;
            Curve edgeCurve2;
            edgeCurve1 = edge.DuplicateCurve();
            edgeCurve2 = edge.DuplicateCurve();

            edgeCurve1.SetStartPoint(mesh.Nodes[start].Coordinate); //forces start point of edgeCurve
            edgeCurve1.SetEndPoint(mesh.Nodes[stop].Coordinate); //forces end point of edgeCurve

            edgeCurve2.SetStartPoint(mesh.Nodes[stop].Coordinate); //forces start point of edgeCurve
            edgeCurve2.SetEndPoint(mesh.Nodes[start].Coordinate); //forces end point of edgeCurve

            bool dummyCrit = edgeCurve2.GetLength() > 0.001; // dummy criteror for curve2 when construction failes

            if (genes > 0)
            {
                if (edgeCurve1.GetLength() > edgeCurve2.GetLength() & dummyCrit)
                {
                    edgeCurve2.Reverse();
                    movedNode = edgeCurve2.PointAtNormalizedLength((0.49 * genes));
                }
                else { movedNode = edgeCurve1.PointAtNormalizedLength((0.49 * genes)); } // move node along edgeCurve
            }
            else if (genes < 0)
            {
                if (edgeCurve1.GetLength() > edgeCurve2.GetLength() & dummyCrit)
                {
                    edgeCurve2.Reverse();
                    movedNode = edgeCurve2.PointAtNormalizedLength(-(0.49 * genes));
                }
                else { movedNode = edgeCurve1.PointAtNormalizedLength((-0.49 * genes)); } // move node along edgeCurve
            }

            return movedNode;
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
            get { return new Guid("ca71dc64-8038-4382-b77c-62e4c4951128"); }
        }
    }
}
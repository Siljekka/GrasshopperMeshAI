using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;


namespace MeshPoints.MoveNodes
{
    public class MoveNodes : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MoveMesh3DVertices class.
        /// </summary>
        public MoveNodes()
          : base("Move Nodes", "mn",
              "Move nodes of a SmartMesh with gene pools",
              "MyPlugIn", "Modify Mesh")
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
            List<Node> newNodes = new List<Node>();

            // 1. Write error if wrong input
            if (!DA.GetData(0, ref oldMesh)) return;

            if (oldMesh.Type == "Solid" & !DA.GetDataList(3, genesW)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "For solid elements, must have input GenesW."); return; }
            if ((genesU.Count < oldMesh.Nodes.Count) | (genesV.Count < oldMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; }
            if (oldMesh.Type == "Solid" & (genesW.Count < oldMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; }

            // 2. Inherit properties from old mesh
            newMesh.nu = oldMesh.nu;
            newMesh.nv = oldMesh.nv;
            newMesh.nw = oldMesh.nw;
            newMesh.Type = oldMesh.Type;
            newMesh.Geometry = oldMesh.Geometry;
            Brep brep = oldMesh.Geometry.Brep;


            // 3. Create new nodes
            for (int i = 0; i < oldMesh.Nodes.Count; i++)
            {
                // a. Check if node is on face or edge.
                Tuple<bool, BrepFace> pointFace = PointOnFace(oldMesh.Nodes[i], brep); // Item1: IsOnFace, Item2: face. Silje: flytte dette inn i Node klasse?
                Tuple<bool, BrepEdge> pointEdge = PointOnEdge(oldMesh.Nodes[i], brep); // Item1: IsOnEdge, Item2: edge. Silje: flytte dette inn i Node klasse?

                // b. Get coordinates of the moved node.
                Point3d meshPoint = GetNewCoordinateOfNode(i, pointFace, pointEdge, oldMesh, genesU, genesV, genesW);

                // c. Make new node from moved node.
                Node node = new Node(i, meshPoint, oldMesh.Nodes[i].BC_U, oldMesh.Nodes[i].BC_V, oldMesh.Nodes[i].BC_W);
                newNodes.Add(node);
            }

            // 4. Set new nodes and elements
            newMesh.Nodes = newNodes;
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
        private Point3d GetNewCoordinateOfNode(int i, Tuple<bool, BrepFace> pointFace, Tuple<bool, BrepEdge> pointEdge, SmartMesh mesh, List<double> genesU, List<double> genesV, List<double> genesW)
        {
            Point3d movedNode = new Point3d();
            bool IsOnEdge = pointEdge.Item1;
            bool IsOnFace = pointFace.Item1;
            BrepEdge edge = pointEdge.Item2;
            BrepFace face = pointFace.Item2;

            Vector3d translationVectorU = Vector3d.Zero;
            Vector3d translationVectorV = Vector3d.Zero;
            Vector3d translationVectorW = Vector3d.Zero;

            // Translation in U direction
            // 1. if: Node not restrained in U direction and gen positive.
            // 2. if: Node not restrained in U direction and gen negative.
            // 3. if: Node restrained in U direction.
            // Note: if point is on edge not restrained in U direction - meshPoint is made

            if (genesU[i] > 0 & !mesh.Nodes[i].BC_U) // 1. if
            {
                translationVectorU = 0.5 * (mesh.Nodes[i + 1].Coordinate - mesh.Nodes[i].Coordinate) * genesU[i]; // make vector translating node in U-direction
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genesU[i], i, i + 1); return movedNode; } // make meshPoint
            }
            else if (genesU[i] < 0 & !mesh.Nodes[i].BC_U)  // 2. if
            {
                translationVectorU = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - 1].Coordinate) * genesU[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genesU[i], i, i - 1); return movedNode; } // make meshPoint
            }
            else { translationVectorU = translationVectorU * 0; }  // 3. if

            // Translation in V direction
            // 1. if: Node not restrained in V direction and gen positive.
            // 2. if: Node not restrained in V direction and gen negative.
            // 3. if: Node restrained in V direction.
            // Note: if point is on edge not restrained in V direction - meshPoint is made
            if (genesV[i] > 0 & !mesh.Nodes[i].BC_V) // 1. if
            {
                translationVectorV = 0.5 * (mesh.Nodes[i + mesh.nu].Coordinate - mesh.Nodes[i].Coordinate) * genesV[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genesV[i], i, i + mesh.nu); return movedNode; } // make meshPoint
            }
            else if (genesV[i] < 0 & !mesh.Nodes[i].BC_V) // 2. if
            {
                translationVectorV = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - mesh.nu ].Coordinate) * genesV[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genesV[i], i, i - mesh.nu); return movedNode; } // make meshPoint
            }
            else { translationVectorV = translationVectorV * 0; } // 3. if

            // Translation in W direction
            // 1. if: Node not restrained in W direction and gen positive.
            // 2. if: Node not restrained in W direction and gen negative.
            // 3. if: Node restrained in W direction.
            // Note: if point is on edge not restrained in W direction - meshPoint is made

            if (mesh.Type == "Solid")
            {

                if (genesW[i] > 0 & !mesh.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (mesh.Nodes[i + (mesh.nu) * (mesh.nv)].Coordinate - mesh.Nodes[i].Coordinate) * genesW[i];
                    if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genesW[i], i, i + (mesh.nu) * (mesh.nv)); return movedNode; } // make meshPoint
                }
                else if (genesW[i] < 0 & !mesh.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - (mesh.nu) * (mesh.nv)].Coordinate) * genesW[i];
                    if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genesW[i], i, i - (mesh.nu) * (mesh.nv)); return movedNode; } // make meshPoint
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
            {// to do: Hilde
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
            
            if (genes > 0)
            {
                if (edgeCurve1.GetLength() > edgeCurve2.GetLength()) 
                {
                    edgeCurve2.Reverse();
                    movedNode = edgeCurve2.PointAtNormalizedLength((0.49 * genes)); // to do: hilde, når gene = 0 starter den ved feil end?
                }
                else { movedNode = edgeCurve1.PointAtNormalizedLength((0.49 * genes)); } // move node along edgeCurve
            }
            else if (genes < 0)
            {
                if (edgeCurve1.GetLength() > edgeCurve2.GetLength())
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
                return Properties.Resources.Icon_MoveSolidMeshVertices;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6de864c8-742d-4bba-96f7-54f3577082c0"); }
        }
    }
}
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using MeshPoints.Classes;


namespace MeshPoints.Tools
{
    public class MoveNodes : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MoveMesh3DVertices class.
        /// </summary>
        public MoveNodes()
          : base("Move Nodes", "moveN",
              "Move nodes of a structured SmartMesh by translation vectors in range [-1,1].",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Translation vectors for u-direction.", GH_ParamAccess.list); 
            pManager.AddGenericParameter("v genes", "qp", "Translation vectors for v-direction.", GH_ParamAccess.list);
            pManager.AddGenericParameter("w genes", "qp", "Translation vectors for w-direction.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Symmetry Line", "sym", "Introduce symmetry by inserting symmetri line.", GH_ParamAccess.item);
            pManager[3].Optional = true; // if solid
            pManager[4].Optional = true; // if symmetry
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "Updated SmartMesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "mesh", "Mesh.", GH_ParamAccess.item);
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
            Curve symLine = null;

            DA.GetData(0, ref oldMesh);
            DA.GetDataList(1, genesU);
            DA.GetDataList(2, genesV);
            DA.GetDataList(3, genesW);
            DA.GetData(4, ref symLine);

            // Variables
            SmartMesh newMesh = new SmartMesh();
            List<Point3d> newPoints = new List<Point3d>();
            double genU = 0;
            double genV = 0;
            double genW = 0;
            double overlapTolerance = 0.95; // ensure no collision of vertices, reduce number to avoid "the look of triangles".


            // 1. Write error if wrong input
            if (!DA.GetData(0, ref oldMesh)) return;
            if (oldMesh.Type == "Solid" & !DA.GetDataList(3, genesW)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "For solid elements, must have input GenesW."); return; }
            
            if (symLine == null)
            {
                if ((genesU.Count < oldMesh.Nodes.Count) | (genesV.Count < oldMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; }
                if (oldMesh.Type == "Solid" & (genesW.Count < oldMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; }
                if (oldMesh.nu == 0 | oldMesh.nv == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Do not support SmartMesh made as unstructured."); return; }
            }
            else 
            {
                if ((genesU.Count < oldMesh.Nodes.Count/2) | (genesV.Count < oldMesh.Nodes.Count/2)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; }
                if (oldMesh.Type == "Solid" & (genesW.Count < oldMesh.Nodes.Count/2)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; }
                if (oldMesh.nu == 0 | oldMesh.nv == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Do not support SmartMesh made as unstructured."); return; }
            }

            // 2. Inherit properties from old mesh
            newMesh.nu = oldMesh.nu;
            newMesh.nv = oldMesh.nv;
            newMesh.nw = oldMesh.nw;
            newMesh.Type = oldMesh.Type;
            newMesh.Geometry = oldMesh.Geometry;
            Brep brep = oldMesh.Geometry.Brep;

            // 2. If sym-case, find the symmetry
            if (symLine != null)
            {
                // Get nodes on symmetry line
                List<int> nodeIdOnSymEdge = GetNodesOnSymmetryLine(oldMesh, symLine);

                // Find direction of symmetry line
                if (nodeIdOnSymEdge.Count == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Did not find symmetry line."); return; }
                int indexDiff = nodeIdOnSymEdge[1] - nodeIdOnSymEdge[0];
                string symDirection = "";
                if (indexDiff == 1) { symDirection = "u"; }
                else if (indexDiff == oldMesh.nu) { symDirection = "v"; }
                else if (indexDiff == oldMesh.nu * oldMesh.nv) { symDirection = "w"; }

                // Get mirror connectivity
                List<List<int>> mirrorConnectivity = GetMirrorConnectivity(oldMesh, symDirection, nodeIdOnSymEdge);

                // Create new nodes
                int counter = 0;
                int uCounter = 0;
                int vCounter = 0;
                int wCounter = 0;
                Point3d[] newPointsArray = new Point3d[oldMesh.Nodes.Count];
                int indexToFix = -1;
                for (int i = 0; i < mirrorConnectivity.Count - 1 + mirrorConnectivity[0].Count; i++)
                {
                    // a. Get correct genes
                    if (i < mirrorConnectivity[0].Count) // if symmetry edge
                    {
                        if (symDirection == "u") { genU = genesU[i]; genV = 0;  genW = 0; uCounter++; }
                        else if (symDirection == "v") { genU = 0; genV = genesV[i]; genW = 0; vCounter++; }
                        else { genU = 0; genV = 0; genW = genesW[i]; wCounter++;}
                        indexToFix++;
                    }
                    else
                    {
                        genU = genesU[uCounter];
                        genV = genesV[vCounter];
                        genW = genesW[wCounter];
                        indexToFix = 0; counter++; uCounter++; vCounter++; wCounter++;
                    }

                    // b. Check if node is on face or edge.
                    Tuple<bool, BrepFace> pointFace = PointOnFace(oldMesh.Nodes[mirrorConnectivity[counter][indexToFix]], brep);  
                    Tuple<bool, BrepEdge> pointEdge = PointOnEdge(oldMesh.Nodes[mirrorConnectivity[counter][indexToFix]], brep); 

                    // c. Get coordinates of the moved node.
                    Point3d point = GetNewCoordinateOfNode(mirrorConnectivity[counter][indexToFix], pointFace, pointEdge, oldMesh, genU, genV, genW, overlapTolerance);

                    // d. Get coordinate of mirror node
                    if (i < mirrorConnectivity[0].Count) // if symmetry edge
                    {
                        newPointsArray[mirrorConnectivity[0][i]] = point;
                    }
                    else 
                    {
                        // Transform
                        Plane symPlane = new Plane(symLine.PointAtStart, symLine.PointAtStart - symLine.PointAtEnd, Vector3d.ZAxis);
                        Transform mirrorMatrix = Transform.Mirror(symPlane);
                        Point3d mirrorPoint = point;
                        mirrorPoint.Transform(mirrorMatrix);

                        newPointsArray[mirrorConnectivity[counter][0]] = point;
                        newPointsArray[mirrorConnectivity[counter][1]] = mirrorPoint;
                    }
                }             
                foreach (Point3d pt in newPointsArray) { newPoints.Add(pt); } // from array to list of points
            }
            else
            {
                // 3. Create new nodes
                for (int i = 0; i < oldMesh.Nodes.Count; i++)
                {
                    // a. Check if node is on face or edge.
                    Tuple<bool, BrepFace> pointFace = PointOnFace(oldMesh.Nodes[i], brep); // Item1: IsOnFace, Item2: face. Silje: flytte dette inn i Node klasse? Og kall på fra GetNewCoord
                    Tuple<bool, BrepEdge> pointEdge = PointOnEdge(oldMesh.Nodes[i], brep); // Item1: IsOnEdge, Item2: edge. Silje: flytte dette inn i Node klasse? Og kall på fra GetNewCoord
                  
                    if (oldMesh.Type == "Solid")
                    {
                        genW = genesW[i];
                    }

                    // b. Get coordinates of the moved node.
                    Point3d point = GetNewCoordinateOfNode(i, pointFace, pointEdge, oldMesh, genesU[i], genesV[i], genW, overlapTolerance);
                    newPoints.Add(point);
                }
            }     

            // 4. Set new nodes and elements
            newMesh.CreateNodes(newPoints, newMesh.nu-1, newMesh.nv-1, newMesh.nw-1);
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
        private Point3d GetNewCoordinateOfNode(int i, Tuple<bool, BrepFace> pointFace, Tuple<bool, BrepEdge> pointEdge, SmartMesh mesh, double genU, double genV, double genW, double overlapTolerance)
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
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genU, i, i + 1, overlapTolerance); return movedNode; } // make meshPoint
            }
            else if (genU < 0 & !mesh.Nodes[i].BC_U)  // 2. if
            {
                translationVectorU = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - 1].Coordinate) * genU;
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genU, i, i - 1, overlapTolerance); return movedNode; } // make meshPoint
            }
            else { translationVectorU = translationVectorU * 0; }  // 3. if

           
            if (genV > 0 & !mesh.Nodes[i].BC_V) // 1. if
            {
                translationVectorV = 0.5 * (mesh.Nodes[i + mesh.nu].Coordinate - mesh.Nodes[i].Coordinate) * genV;
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genV, i, i + mesh.nu, overlapTolerance); return movedNode; } // make meshPoint
            }
            else if (genV < 0 & !mesh.Nodes[i].BC_V) // 2. if
            {
                translationVectorV = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - mesh.nu ].Coordinate) * genV;
                if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genV, i, i - mesh.nu, overlapTolerance); return movedNode; } // make meshPoint
            }
            else { translationVectorV = translationVectorV * 0; } // 3. if


            if (mesh.Type == "Solid")
            {

                if (genW > 0 & !mesh.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (mesh.Nodes[i + (mesh.nu) * (mesh.nv)].Coordinate - mesh.Nodes[i].Coordinate) * genW;
                    if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genW, i, i + (mesh.nu) * (mesh.nv), overlapTolerance); return movedNode; } // make meshPoint
                }
                else if (genW < 0 & !mesh.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (mesh.Nodes[i].Coordinate - mesh.Nodes[i - (mesh.nu) * (mesh.nv)].Coordinate) * genW;
                    if (IsOnEdge) { movedNode = EdgeNode(edge, mesh, genW, i, i - (mesh.nu) * (mesh.nv), overlapTolerance); return movedNode; } // make meshPoint
                }
                else { translationVectorW = translationVectorW * 0; } // 3. if                            
            }

            // 4. if: Make movedNode if node is on face or inside brep (if on edge, movedNode already made).
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
        private Point3d EdgeNode(BrepEdge edge, SmartMesh mesh, double genes, int start, int stop, double overlapTolerance)
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
                    movedNode = edgeCurve2.PointAtNormalizedLength(0.5 * overlapTolerance * genes); 
                }
                else { movedNode = edgeCurve1.PointAtNormalizedLength(0.5 * overlapTolerance * genes); } // move node along edgeCurve
            }
            else if (genes < 0)
            {
                if (edgeCurve1.GetLength() > edgeCurve2.GetLength() & dummyCrit)
                {
                    edgeCurve2.Reverse();
                    movedNode = edgeCurve2.PointAtNormalizedLength(-(0.5 * overlapTolerance * genes));
                }
                else { movedNode = edgeCurve1.PointAtNormalizedLength((-0.5 * overlapTolerance * genes)); } // move node along edgeCurve
            }

            return movedNode;
        }

        private List<int> GetNodesOnSymmetryLine(SmartMesh oldMesh, Curve symLine)
        {
            List<int> nodeIdOnSymEdge = new List<int>();
            foreach (Node node in oldMesh.Nodes)
            {
                Point3d point = node.Coordinate;
                symLine.ClosestPoint(point, out double PointOnCurve);
                Point3d testPoint = symLine.PointAt(PointOnCurve);  // make test point 
                double distanceToEdge = (testPoint - point).Length; // calculate distance between testPoint and node
                if (distanceToEdge <= 0.0001 & distanceToEdge >= -0.0001) // if distance = 0: node is on edge
                {
                    nodeIdOnSymEdge.Add(node.GlobalId);
                }
            }
            return nodeIdOnSymEdge;
        }
        private List<List<int>> GetMirrorConnectivity(SmartMesh oldMesh, string symDirection, List<int> nodeIdOnSymEdge)
        {
            // Find mirror connectivity
            int nu = oldMesh.nu;
            int nv = oldMesh.nv;
            int nw = oldMesh.nw;

            List<List<int>> mirrorConnectivity = new List<List<int>>();
            int count = nodeIdOnSymEdge.Count;
            for (int i = 0; i < count; i++)
            {
                int id = nodeIdOnSymEdge[i];
                if (symDirection == "w")
                {
                    for (int j = 1; j < nv; j++)
                    {
                        nodeIdOnSymEdge.Add(oldMesh.Nodes[id + (nu) * j].GlobalId);
                    }
                }
            }
            nodeIdOnSymEdge.Sort();
            mirrorConnectivity.Add(nodeIdOnSymEdge); // first list is the id of nodes on symEdge 
            int idCounter = 0;

            if (symDirection == "u")
            {
                for (int k = 0; k < nw; k++)
                {
                    for (int j = 0; j < nv; j++)
                    {
                        for (int i = 0; i < nu; i++)
                        {
                            List<int> nodesToPare = new List<int>() { idCounter, ((nv - (j + 1)) * nu + i) + k * nu * nv };
                            mirrorConnectivity.Add(nodesToPare);
                            idCounter++;
                        }
                    }
                    idCounter = (k + 1) * nu * nv;
                }
            }
            else if (symDirection == "v")
            {
                for (int k = 0; k < nw; k++)
                {
                    for (int j = 0; j < nv; j++)
                    {
                        for (int i = 0; i < Math.Floor((double)nu / (double)2); i++)
                        {
                            List<int> nodesToPare = new List<int>() { idCounter, (nu - 1) - i + j * nv + k * nu * nv };
                            mirrorConnectivity.Add(nodesToPare);
                            idCounter++;
                        }
                        idCounter = (j + 1) * nu + k * nv * nu;
                    }
                    idCounter = (k + 1) * nu * nv;
                }
            }
            else if (symDirection == "w")
            {
                for (int k = 0; k < nw; k++)
                {
                    for (int j = 0; j < nv; j++)
                    {
                        for (int i = 0; i < Math.Floor((double)nu / (double)2); i++)
                        {
                            List<int> nodesToPare = new List<int>() { idCounter, (nu - 1) - i + j * nu + k * nu * nv };
                            mirrorConnectivity.Add(nodesToPare);
                            idCounter++;
                        }
                        idCounter = (j + 1) * nu + k * nu * nv;
                    }
                    idCounter = (k + 1) * nu * nv;
                }
            }

            return mirrorConnectivity;
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
                return Properties.Resources.Icon_MoveNodes;
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
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry.Collections;
using System.Drawing;
using MeshPoints.Classes;
using Rhino.Geometry.Intersect;


namespace MeshPoints.Galapagos
{
    public class GalapagosMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MoveMesh3DVertices class.
        /// </summary>
        public GalapagosMesh()
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
            pManager.AddGenericParameter("Geometry", "geo", "Input source geometry", GH_ParamAccess.item);
            pManager.AddGenericParameter("SmartMesh", "sm", "Input a SmartMesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Gene pool for translation in u direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("v genes", "qp", "Gene pool for translation in v direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("w genes", "qp", "Gene pool for translation in w direction", GH_ParamAccess.list);
            pManager[4].Optional = true; // if solid

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
            Mesh3D oldMesh = new Mesh3D();
            Brep brep = new Brep();
            List<double> genesU = new List<double>();
            List<double> genesV = new List<double>();
            List<double> genesW = new List<double>();

            DA.GetData(0, ref brep);
            DA.GetData(1, ref oldMesh);
            DA.GetDataList(2, genesU);
            DA.GetDataList(3, genesV);
            DA.GetDataList(4, genesW);

            // Variables
            Mesh3D newMesh = new Mesh3D();
            Mesh allMesh = new Mesh();
            Node n = new Node();
            List<Node> newNodes = new List<Node>();
            List<Element> elements = new List<Element>();

            newMesh.nu = oldMesh.nu;
            newMesh.nv = oldMesh.nv;
            newMesh.nw = oldMesh.nw;
            newMesh.Type = oldMesh.Type;

            // 1. Write error if wrong input
            if (!brep.IsValid) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Brep input is not valid."); return; }
            if (oldMesh == null) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "SmartMesh input is not valid."); return; }
            if ((genesU.Count < oldMesh.Nodes.Count) | (genesV.Count < oldMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; } 
            if (oldMesh.Type == "Solid" & (genesW.Count < oldMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; }

            // 2. Move and make new nodes
            for (int i = 0; i < oldMesh.Nodes.Count; i++)
            {
                // a. Check if node is on face or edge.
                //    if node is on face: true and face is output
                //    if node is on edge: true and edge is output
                Tuple<bool, BrepFace> pointFace = PointOnFace(i, oldMesh.Nodes, brep); // Item1: IsOnFace, Item2: face
                Tuple<bool, BrepEdge> pointEdge = PointOnEdge(i, oldMesh.Nodes, brep); // Item1: IsOnEdge, Item2: edge

                // b. Get coordinates of the moved node.
                Point3d meshPoint = GetMovedNode(i, pointFace, pointEdge, oldMesh, genesU, genesV, genesW);

                // c. Make new node from moved node.
                n = new Node(i, meshPoint, oldMesh.Nodes[i].BC_U, oldMesh.Nodes[i].BC_V, oldMesh.Nodes[i].BC_W); // todo: fix local id;
                newNodes.Add(n);
            }

            // 3. Make elements from moved nodes
            newMesh.Nodes = newNodes;
            if (newMesh.Type == "Surface")
            {
                newMesh.SetQuadElements();
            }
            else
            {
                newMesh.SetHexElements();
            }
            //elements = CreateNewElements(newNodes, oldMesh);


            //4. Create global mesh 
            allMesh = CreateGlobalMesh(newMesh.Elements); //todo: do this without using weld!
            newMesh.mesh = allMesh;
            //5. Add properties to updated SmartMesh
            
            /*Mesh3D newMesh;
            if (oldMesh.Type == "Surface")
            {
                 newMesh = new Mesh3D(oldMesh.nu, oldMesh.nv, newNodes, elements, allMesh);
            }
            else
            {
                newMesh = new Mesh3D(oldMesh.nu, oldMesh.nv, oldMesh.nw, newNodes, elements, allMesh);
            }
            */

            // Output
            DA.SetData(0, newMesh);
            DA.SetData(1, newMesh.mesh);
        }
        #region Methods
        /// <summary>Create hexahedral or quadlateral elements depending on a the input mesh.</summary>
        /// <param name="newNodes">List of global nodes</param>
        /// <param name="nu">Number of points in nu.</param>
        /// <param name="nv">Number of points in nv.</param>
        /// <param name="nw">Number of points in nw.</param>
        /// <returns>Return a list of elements.</returns>
        private List<Element> CreateNewElements(List<Node> newNodes, Mesh3D oldMesh) // to do: endre til ny metode
        {
            List<Element> elements = new List<Element>();

            if (oldMesh.Type == "Surface")
            {
                //elements = CreateQuadElements(newNodes, oldMesh.nu, oldMesh.nv);
            }
            else
            {
                //elements = CreateHexElementsOld(newNodes, oldMesh.nu, oldMesh.nv, oldMesh.nw);
            } 
            return elements;
        }

        /// <summary>Create global mesh.</summary>
        /// <param name="elements">List of elements</param>
        /// <returns>Return global mesh</returns>
        private Mesh CreateGlobalMesh(List<Element> elements) // to do: endre til ny metode
        {  
            Mesh allMesh = new Mesh();
            foreach (Element el in elements)
            {
                allMesh.Append(el.mesh);
            }
            allMesh.Weld(0.01);

            return allMesh;
        }
        
        /// <summary>Make mesh consistent</summary>
        /// <param name="mesh">Mesh to make consistent.</param>
        /// <returns></returns>
        private void MakeConsistent(Mesh mesh)
        {   //todo: code used before - remove?
            mesh.Normals.ComputeNormals();  // todo: control if needed
            mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            mesh.Compact(); // to ensure that it calculate
        }

        /// <summary>
        /// todo: write description of method
        /// </summary>
        /// <returns>todo: write description of output.</returns>
        private Mesh CreateGlobalMesh(Mesh m, int counter, int nu, int nv)
        {   //todo: code used before - remove?
            int meshPtsAtLevel = nu * nv;
            m.Faces.AddFace(counter, counter + 1, counter + meshPtsAtLevel + 1, counter + meshPtsAtLevel);
            m.Faces.AddFace(counter + 1, counter + nu + 1, counter + meshPtsAtLevel + nu + 1, counter + meshPtsAtLevel + 1);
            m.Faces.AddFace(counter + nu + 1, counter + nu, counter + meshPtsAtLevel + nu, counter + meshPtsAtLevel + nu + 1);
            m.Faces.AddFace(counter + nu, counter, counter + meshPtsAtLevel, counter + meshPtsAtLevel + nu);
            m.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);
            m.Faces.AddFace(counter + meshPtsAtLevel, counter + meshPtsAtLevel + 1, counter + meshPtsAtLevel + nu + 1, counter + meshPtsAtLevel + nu);
            return m;
        }

        /// <summary>
        /// Check if point is on face.
        /// </summary>
        /// <returns> Returns true and face if point is on face. False and null face if point is not on face.</returns>
        private Tuple<bool, BrepFace> PointOnFace(int i, List<Node> nodes, Brep brep)
        {
            bool IsOnFace = false;
            BrepFace face = null;
            BrepFaceList brepFace = brep.Faces;

            foreach (BrepFace bFace in brepFace) // check if node is on edge
            {
                bFace.ClosestPoint(nodes[i].Coordinate, out double PointOnCurveU, out double PointOnCurveV);
                Point3d testPoint = bFace.PointAt(PointOnCurveU, PointOnCurveV);  // make test point 
                double distanceToFace = testPoint.DistanceTo(nodes[i].Coordinate); // calculate distance between testPoint and node
                if (distanceToFace <= 0.0001 & distanceToFace >= -0.0001) // if distance = 0: node is on edge
                {
                    if (nodes[i].BC_U & nodes[i].BC_V & nodes[i].BC_W) { IsOnFace = false; } // cornerpoints
                    else if ((!nodes[i].BC_U & !nodes[i].BC_V) | (!nodes[i].BC_U & !nodes[i].BC_W) | (!nodes[i].BC_V & !nodes[i].BC_W))
                    {
                        IsOnFace = true;
                        face = bFace;
                    }
                }
            }

            return new Tuple<bool, BrepFace>(IsOnFace, face);
        }

        /// <summary>
        /// Check if point is on edge.
        /// </summary>
        /// <returns> Returns true and edge if point is on edge. False and null edge if point is not on edge.</returns>
        private Tuple<bool, BrepEdge> PointOnEdge(int i, List<Node> nodes, Brep brep)
        {
            bool IsOnEdge = false;
            BrepEdge edge = null;
            BrepEdgeList brepEdge = brep.Edges;

            IsOnEdge = false;
            foreach (BrepEdge bEdge in brepEdge) // check if node is on edge
            {
                bEdge.ClosestPoint(nodes[i].Coordinate, out double PointOnCurve);
                Point3d testPoint = bEdge.PointAt(PointOnCurve);  // make test point 
                double distanceToEdge = testPoint.DistanceTo(nodes[i].Coordinate); // calculate distance between testPoint and node
                if (distanceToEdge <= 0.0001 & distanceToEdge >= -0.0001) // if distance = 0: node is on edge
                {
                    if (nodes[i].BC_U & nodes[i].BC_V & nodes[i].BC_W) { IsOnEdge = false; } // cornerpoints: IsOnCurve must be false
                    else if ((nodes[i].BC_U & nodes[i].BC_V) | (nodes[i].BC_U & nodes[i].BC_W) | (nodes[i].BC_V & nodes[i].BC_W))
                    {
                        IsOnEdge = true;
                        edge = bEdge;
                    }
                }
            }
            return new Tuple<bool, BrepEdge>(IsOnEdge, edge);
        }

        /// <summary>
        /// Move the old node in allowable directions.
        /// </summary>
        /// <returns> Returns coordinates of moved node.</returns>
        private Point3d GetMovedNode(int i, Tuple<bool, BrepFace> pointFace, Tuple<bool, BrepEdge> pointEdge, Mesh3D m, List<double> genesU, List<double> genesV, List<double> genesW)
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
            // 1. if: Node not restrained in U direction and genes positive.
            // 2. if: Node not restrained in U direction and genes negative.
            // 3. if: Node restrained in U direction.
            // Note: if point is on edge not restrained in U direction - meshPoint is made
            if (genesU[i] >= 0 & !m.Nodes[i].BC_U) // 1. if
            {
                translationVectorU = 0.5 * (m.Nodes[i + 1].Coordinate - m.Nodes[i].Coordinate) * genesU[i]; // make vector translating node in U-direction
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesU[i], i, i + 1); return movedNode; } // make meshPoint
            }
            else if (genesU[i] <= 0 & !m.Nodes[i].BC_U)  // 2. if
            {
                translationVectorU = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - 1].Coordinate) * genesU[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesU[i], i, i - 1); return movedNode; } // make meshPoint
            }
            else { translationVectorU = translationVectorU * 0; }  // 3. if

            // Translation in V direction
            // 1. if: Node not restrained in V direction and genes positive.
            // 2. if: Node not restrained in V direction and genes negative.
            // 3. if: Node restrained in V direction.
            // Note: if point is on edge not restrained in V direction - meshPoint is made
            if (genesV[i] >= 0 & !m.Nodes[i].BC_V) // 1. if
            {
                translationVectorV = 0.5 * (m.Nodes[i + m.nu].Coordinate - m.Nodes[i].Coordinate) * genesV[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesV[i], i, i + m.nu); return movedNode; } // make meshPoint
            }
            else if (genesV[i] <= 0 & !m.Nodes[i].BC_V) // 2. if
            {
                translationVectorV = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - m.nu ].Coordinate) * genesV[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesV[i], i, i - m.nu); return movedNode; } // make meshPoint
            }
            else { translationVectorV = translationVectorV * 0; } // 3. if

            // Translation in W direction
            // 1. if: Node not restrained in W direction and genes positive.
            // 2. if: Node not restrained in W direction and genes negative.
            // 3. if: Node restrained in W direction.
            // Note: if point is on edge not restrained in W direction - meshPoint is made
            if (m.Type == "Solid")
            {
                if (genesW[i] >= 0 & !m.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (m.Nodes[i + (m.nu) * (m.nv)].Coordinate - m.Nodes[i].Coordinate) * genesW[i];
                    if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesW[i], i, i + (m.nu) * (m.nv)); return movedNode; } // make meshPoint
                }
                else if (genesW[i] <= 0 & !m.Nodes[i].BC_W) // 1. if
                {
                    translationVectorW = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - (m.nu) * (m.nv)].Coordinate) * genesW[i];
                    if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesW[i], i, i - (m.nu) * (m.nv)); return movedNode; } // make meshPoint
                }
                else { translationVectorW = translationVectorW * 0; } // 3. if
            }

            // 4. if: Make movedNode if node is on face or inside brep (if on edge, movedNode already made).
            double overlapTolerance = 0.99; // ensure no collision of vertices, reduce number to avoid "the look of triangles".
            movedNode = new Point3d
                (
                m.Nodes[i].Coordinate.X + (translationVectorU.X + translationVectorV.X + translationVectorW.X) * overlapTolerance,
                m.Nodes[i].Coordinate.Y + (translationVectorU.Y + translationVectorV.Y + translationVectorW.Y) * overlapTolerance,
                m.Nodes[i].Coordinate.Z + (translationVectorU.Z + translationVectorV.Z + translationVectorW.Z) * overlapTolerance
                );
                
            if (IsOnFace) // If node is on face: ensure it stays on face
            {// to do: slett?
                Brep srf = face.DuplicateFace(false);
                movedNode = srf.ClosestPoint(movedNode); // "Project" meshPoint to surface.
            }
            return movedNode;

        }

        /// <summary>
        /// Make new node if point is on edge.
        /// </summary>
        /// <returns> Returns coordinates of moved node on edge.</returns>
        private Point3d EdgeNode(BrepEdge edge, Mesh3D m, double genes, int start, int stop)
        {
            Point3d movedNode = new Point3d();
            Curve edgeCurve1;
            Curve edgeCurve2;
            edgeCurve1 = edge.DuplicateCurve();
            edgeCurve2 = edge.DuplicateCurve();

            edgeCurve1.SetStartPoint(m.Nodes[start].Coordinate); //forces start point of edgeCurve
            edgeCurve1.SetEndPoint(m.Nodes[stop].Coordinate); //forces end point of edgeCurve

            edgeCurve2.SetStartPoint(m.Nodes[stop].Coordinate); //forces start point of edgeCurve
            edgeCurve2.SetEndPoint(m.Nodes[start].Coordinate); //forces end point of edgeCurve
            
            if (genes >= 0)
            {
                if (edgeCurve1.GetLength() > edgeCurve2.GetLength()) 
                {
                    edgeCurve2.Reverse();
                    movedNode = edgeCurve2.PointAtNormalizedLength((0.49 * genes));
                }
                else { movedNode = edgeCurve1.PointAtNormalizedLength((0.49 * genes)); } // move node along edgeCurve}
            }
            else if (genes <= 0)
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

        private List<Element> CreateQuadElementsOld(List<Node> nodes, int nu, int nv) // to do: erstatt med Mesh3D metode
        {
            List<Element> elements = new List<Element>();
            int uSequence = 0;
            int counter = 0;

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

        private List<Element> CreateHexElementsOld(List<Node> nodes, int nu, int nv, int nw)// to do: erstatt medh Mesh3d
        {
            List<Element> elements = new List<Element>();
            int elemId = 0;

            for (int i = 0; i < nw - 1; i++)  // loop levels
            {
                int sequence = 0;
                int counter = (nu * nv) * i;

                for (int j = 0; j < (nu * nv) - nu - 1; j++) // loop elements in a level
                {
                    List<Node> elementNodes = new List<Node>();
                    List<int> connectivity = new List<int>();

                    if (sequence < nu - 1)
                    {
                        connectivity.Add(counter);
                        connectivity.Add(counter + 1);
                        connectivity.Add(counter + nu + 1);
                        connectivity.Add(counter + nu);
                        connectivity.Add(counter + nu * nv);
                        connectivity.Add(counter + 1 + nu * nv);
                        connectivity.Add(counter + nu + 1 + nu * nv);
                        connectivity.Add(counter + nu + nu * nv);

                        foreach (int id in connectivity)
                        {
                            elementNodes.Add(nodes[id]);
                        }

                        Element element = new Element(elemId, elementNodes, connectivity);

                        // create local mesh
                        Mesh localMesh = new Mesh();
                        foreach (Node node in elementNodes)
                        {
                            localMesh.Vertices.Add(node.Coordinate); //0
                        }
                        localMesh.Faces.AddFace(0, 1, 5, 4);
                        localMesh.Faces.AddFace(1, 2, 6, 5);
                        localMesh.Faces.AddFace(2, 3, 7, 6);
                        localMesh.Faces.AddFace(3, 0, 4, 7);
                        localMesh.Faces.AddFace(0, 1, 2, 3);
                        localMesh.Faces.AddFace(4, 5, 6, 7);

                        localMesh.Normals.ComputeNormals();  //Control if needed
                        localMesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                        localMesh.Compact(); //to ensure that it calculate
                        element.mesh = localMesh;

                        //add element and mesh to element list
                        elements.Add(element);

                        sequence++;
                        elemId++;
                        counter++;
                    }
                    else { sequence = 0; counter++; }
                }
            }
            return elements;
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
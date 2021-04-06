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
          : base("Move SolidMesh Vertices", "m3dv",
              "Move mesh vertices with gene pools",
              "MyPlugIn", "Modify Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Brep", "bp", "Input source brep", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh3D", "m2d", "Input Mesh2D", GH_ParamAccess.item);
            pManager.AddGenericParameter("u genes ", "qp", "Gene pool for translation in u direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("v genes", "qp", "Gene pool for translation in v direction", GH_ParamAccess.list);
            pManager.AddGenericParameter("w genes", "qp", "Gene pool for translation in w direction", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh3D", "m3D", "Updated solid mesh", GH_ParamAccess.item);
            //pManager.AddGenericParameter("test", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Input
            // Input
            Mesh3D inputMesh = new Mesh3D();
            Brep brep = new Brep();
            List<double> genesU = new List<double>();
            List<double> genesV = new List<double>();
            List<double> genesW = new List<double>();

            DA.GetData(0, ref brep);
            DA.GetData(1, ref inputMesh);
            DA.GetDataList(2, genesU);
            DA.GetDataList(3, genesV);
            DA.GetDataList(4, genesW);
            #endregion

            #region Variables
            // Variables
            Mesh3D solidMesh = new Mesh3D();
            Mesh allMesh = new Mesh();
            Node n = new Node();
            List<Node> nodes = new List<Node>();
            List<Element> elements = new List<Element>();
            //int newRow = 0; old variable
            //int counter = 0; old variable
            //Element e = new Element(); old variable
            //Mesh mesh = new Mesh(); old variable
            //Mesh globalMesh = new Mesh(); old variable
            //Mesh3D meshUpdated = new Mesh3D(); old variable
            #endregion

            // 1. Write warning and error if wrong input
            if (!brep.IsValid) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Brep input is not valid."); return; }
            if (inputMesh == null) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh3D input is not valid."); return; }
            if ((genesU.Count < inputMesh.Nodes.Count) | (genesV.Count < inputMesh.Nodes.Count) | (genesW.Count < inputMesh.Nodes.Count)) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Must increase genes."); return; } //todo: add warning message


            // 2. Move and make new nodes
            for (int i = 0; i < inputMesh.Nodes.Count; i++)
            {
                // a. Check if node is on face or edge.
                //    if node is on face: true and face is output
                //    if node is on edge: true and edge is output
                Tuple<bool, BrepFace> pointFace = PointOnFace(i, inputMesh.Nodes, brep); // Item1: IsOnFace, Item2: face
                Tuple<bool, BrepEdge> pointEdge = PointOnEdge(i, inputMesh.Nodes, brep); // Item1: IsOnEdge, Item2: edge

                // b. Get coordinates of the moved node.
                Point3d meshPoint = GetMovedNode(i, pointFace, pointEdge, inputMesh, genesU, genesV, genesW);

                // c. Make new node from moved node.
                n = new Node(i, meshPoint, inputMesh.Nodes[i].BC_U, inputMesh.Nodes[i].BC_V, inputMesh.Nodes[i].BC_W); // todo: fix local id;
                nodes.Add(n);
                //globalMesh.Vertices.Add(meshPoint); old line
            }

            // 3. Make elements from moved nodes
            elements = CreateHexElements(nodes, inputMesh.nu, inputMesh.nv, inputMesh.nw);

            //4. Create global mesh 
            allMesh = CreateGlobalMesh(elements); //todo: do this without using weld!

            //5. Add properties to Mesh3D
            solidMesh.Nodes = nodes;
            solidMesh.Elements = elements;
            solidMesh.mesh = allMesh;

            #region Old code: Make Element and mesh3D
            /*m.nu = m.nu + 1;
            m.nv = m.nv + 1;
            counter = 0;
            for (int j = 0; j < m.nw; j++) // loop trough levels
            {
                counter = m.nu * m.nv * j;
                for (int i = 0; i < (m.nu - 1) * (m.nv - 1); i++) // mesh a level
                {
                    int id = i * (j + 1); // element id
                    e = CreateElement(id, nodes, counter, m.nu, m.nv);
                    elements.Add(e); // add element and mesh to element list
                    globalMesh = CreateGlobalMesh(globalMesh, counter, m.nu, m.nv);

                    // clear
                    e = new Element();
                    mesh = new Mesh();

                    // element counter
                    counter++;
                    newRow++; ;
                    if (newRow == (m.nu - 1)) //new row
                    {
                        counter++;
                        newRow = 0;
                    }
                }
            }*/
            #endregion

            // Output
            DA.SetData(0, solidMesh);
        }
        #region Methods
        /// <summary>
        /// Create Elements: assign ElementId, ElementMesh and Nodes incl. Coordiantes, GlobalId, LocalId and Boundary Conditions), elementId, elementMesh.
        /// </summary>
        /// <returns>List with elements incl properties</returns>
        private List<Element> CreateHexElements(List<Node> nodes, int nu, int nv, int nw)
        {
            Element e = new Element();
            Mesh mesh = new Mesh();
            List<Element> elements = new List<Element>();
            int elemId = 0;

            nu = nu + 1; //input nu = nu - 1. Exs: nu = 3, total points in u-direction is 4;
            nv = nv + 1; //input nv = nv - 1. Exs: nv = 3, total points in v-direction is 4;
            nw = nw + 1; //input nv = nv - 1. Exs: nv = 3, total points in v-direction is 4;

            for (int i = 0; i < nw - 1; i++)  // loop levels
            {
                int count2 = 0;
                int counter = (nu * nv) * i;

                for (int j = 0; j < (nu * nv) - nu - 1; j++) // loop elements in a level
                {
                    e.Id = elemId;
                    e.IsCube = true;
                    if (count2 < nu - 1)
                    {
                        Node n1 = new Node(1, nodes[counter].GlobalId, nodes[counter].Coordinate, nodes[counter].BC_U, nodes[counter].BC_V, nodes[counter].BC_W);
                        e.Node1 = n1;

                        Node n2 = new Node(2, nodes[counter + 1].GlobalId, nodes[counter + 1].Coordinate, nodes[counter + 1].BC_U, nodes[counter + 1].BC_V, nodes[counter + 1].BC_W);
                        e.Node2 = n2;

                        Node n3 = new Node(3, nodes[counter + nu + 1].GlobalId, nodes[counter + nu + 1].Coordinate, nodes[counter + nu + 1].BC_U, nodes[counter + nu + 1].BC_V, nodes[counter + nu + 1].BC_W);
                        e.Node3 = n3;

                        Node n4 = new Node(4, nodes[counter + nu].GlobalId, nodes[counter + nu].Coordinate, nodes[counter + nu].BC_U, nodes[counter + nu].BC_V, nodes[counter + nu].BC_W);
                        e.Node4 = n4;

                        Node n5 = new Node(5, nodes[counter + nu * nv].GlobalId, nodes[counter + nu*nv].Coordinate, nodes[counter + nu * nv].BC_U, nodes[counter + nu * nv].BC_V, nodes[counter + nu * nv].BC_W);
                        e.Node5 = n5;

                        Node n6 = new Node(6, nodes[counter + 1 + nu * nv].GlobalId, nodes[counter + 1 + nu*nv].Coordinate, nodes[counter + 1 + nu * nv].BC_U, nodes[counter + 1 + nu * nv].BC_V, nodes[counter + 1 + nu * nv].BC_W);
                        e.Node6 = n6;

                        Node n7 = new Node(7, nodes[counter + nu + 1 + nu * nv].GlobalId, nodes[counter + nu + 1 + nu*nv].Coordinate, nodes[counter + nu + 1 + nu * nv].BC_U, nodes[counter + nu + 1 + nu * nv].BC_V, nodes[counter + nu + 1 + nu * nv].BC_W);
                        e.Node7 = n7;

                        Node n8 = new Node(8, nodes[counter + nu + nu * nv].GlobalId, nodes[counter + nu + nu*nv].Coordinate, nodes[counter + nu + nu * nv].BC_U, nodes[counter + nu + nu * nv].BC_V, nodes[counter + nu + nu * nv].BC_W);
                        e.Node8 = n8;

                        mesh.Vertices.Add(e.Node1.Coordinate); //0
                        mesh.Vertices.Add(e.Node2.Coordinate); //1
                        mesh.Vertices.Add(e.Node3.Coordinate); //2
                        mesh.Vertices.Add(e.Node4.Coordinate); //3
                        mesh.Vertices.Add(e.Node5.Coordinate); //4
                        mesh.Vertices.Add(e.Node6.Coordinate); //5
                        mesh.Vertices.Add(e.Node7.Coordinate); //6
                        mesh.Vertices.Add(e.Node8.Coordinate); //7

                        mesh.Faces.AddFace(0, 1, 5, 4);
                        mesh.Faces.AddFace(1, 2, 6, 5);
                        mesh.Faces.AddFace(2, 3, 7, 6);
                        mesh.Faces.AddFace(3, 0, 4, 7);
                        mesh.Faces.AddFace(0, 1, 2, 3);
                        mesh.Faces.AddFace(4, 5, 6, 7);

                        mesh.Normals.ComputeNormals();  //Control if needed
                        mesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
                        mesh.Compact(); //to ensure that it calculate
                        e.mesh = mesh;

                        //add element and mesh to element list
                        elements.Add(e);

                        //clear
                        e = new Element();
                        mesh = new Mesh();

                        count2++;
                        elemId++;
                        counter++;
                    }
                    else { count2 = 0; counter++; }
                }
            }
            return elements;
        }

        /// <summary>
        /// Create Global mesh. todo: make without weld
        /// </summary>
        /// <returns>Global mesh</returns>
        private Mesh CreateGlobalMesh(List<Element> elements)
        {   //todo: make this without weld
            Mesh allMesh = new Mesh();
            foreach (Element el in elements)
            {
                allMesh.Append(el.mesh);
            }
            allMesh.Weld(0.01);

            return allMesh;
        }
        /// <summary>
        /// todo: write description of method
        /// </summary>
        /// <returns> todo: write description of output.</returns>
        private Mesh MakeConsistent(Mesh mesh)
        {   //todo: code used before - remove?
            mesh.Normals.ComputeNormals();  // todo: control if needed
            mesh.FaceNormals.ComputeFaceNormals();  // want a consistant mesh
            mesh.Compact(); // to ensure that it calculate
            return mesh;
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
        /// todo: write description of method
        /// </summary>
        /// <returns> todo: write description of output.</returns>
        private Element CreateElement(int id, List<Node> nodes, int counter, int nu, int nv)
        {   //todo: code used before - remove?
            Element e = new Element();
            int meshPtsAtLevel = nu * nv;

            e.Id = id;
            e.Node1 = nodes[counter];
            e.Node1.LocalId = 1;
            e.Node2 = nodes[counter + 1];
            e.Node2.LocalId = 2;
            e.Node3 = nodes[counter + nu + 1];
            e.Node3.LocalId = 3;
            e.Node4 = nodes[counter + nu];
            e.Node4.LocalId = 4;
            e.Node5 = nodes[counter + meshPtsAtLevel];
            e.Node5.LocalId = 5;
            e.Node6 = nodes[counter + meshPtsAtLevel + 1];
            e.Node6.LocalId = 6;
            e.Node7 = nodes[counter + meshPtsAtLevel + nu + 1];
            e.Node7.LocalId = 7;
            e.Node8 = nodes[counter + meshPtsAtLevel + nu];
            e.Node8.LocalId = 8;

            Mesh mesh = new Mesh();
            mesh.Vertices.Add(e.Node1.Coordinate); //0
            mesh.Vertices.Add(e.Node2.Coordinate); //1
            mesh.Vertices.Add(e.Node3.Coordinate); //2
            mesh.Vertices.Add(e.Node4.Coordinate); //3
            mesh.Vertices.Add(e.Node5.Coordinate); //4
            mesh.Vertices.Add(e.Node6.Coordinate); //5
            mesh.Vertices.Add(e.Node7.Coordinate); //6
            mesh.Vertices.Add(e.Node8.Coordinate); //7

            mesh.Faces.AddFace(0, 1, 5, 4);
            mesh.Faces.AddFace(1, 2, 6, 5);
            mesh.Faces.AddFace(2, 3, 7, 6);
            mesh.Faces.AddFace(3, 0, 4, 7);
            mesh.Faces.AddFace(0, 1, 2, 3);
            mesh.Faces.AddFace(4, 5, 6, 7);
            mesh = MakeConsistent(mesh);
            e.mesh = mesh;
            return e;
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
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesU[i], i, i + 1); } // make meshPoint
            }
            else if (genesU[i] <= 0 & !m.Nodes[i].BC_U)  // 2. if
            {
                translationVectorU = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - 1].Coordinate) * genesU[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesU[i], i, i - 1); } // make meshPoint
            }
            else { translationVectorU = translationVectorU * 0; }  // 3. if

            // Translation in V direction
            // 1. if: Node not restrained in V direction and genes positive.
            // 2. if: Node not restrained in V direction and genes negative.
            // 3. if: Node restrained in V direction.
            // Note: if point is on edge not restrained in V direction - meshPoint is made
            if (genesV[i] >= 0 & !m.Nodes[i].BC_V) // 1. if
            {
                translationVectorV = 0.5 * (m.Nodes[i + (m.nu + 1)].Coordinate - m.Nodes[i].Coordinate) * genesV[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesV[i], i, i + (m.nu + 1)); } // make meshPoint
            }
            else if (genesV[i] <= 0 & !m.Nodes[i].BC_V) // 2. if
            {
                translationVectorV = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - (m.nu + 1)].Coordinate) * genesV[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesV[i], i, i - (m.nu + 1)); } // make meshPoint
            }
            else { translationVectorV = translationVectorV * 0; } // 3. if

            // Translation in W direction
            // 1. if: Node not restrained in W direction and genes positive.
            // 2. if: Node not restrained in W direction and genes negative.
            // 3. if: Node restrained in W direction.
            // Note: if point is on edge not restrained in W direction - meshPoint is made
            if (genesW[i] >= 0 & !m.Nodes[i].BC_W) // 1. if
            {
                translationVectorW = 0.5 * (m.Nodes[i + (m.nu + 1) * (m.nv + 1)].Coordinate - m.Nodes[i].Coordinate) * genesW[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesW[i], i, i + (m.nu + 1) * (m.nv + 1)); } // make meshPoint
            }
            else if (genesW[i] <= 0 & !m.Nodes[i].BC_W) // 1. if
            {
                translationVectorW = 0.5 * (m.Nodes[i].Coordinate - m.Nodes[i - (m.nu + 1) * (m.nv + 1)].Coordinate ) * genesW[i];
                if (IsOnEdge) { movedNode = EdgeNode(edge, m, genesW[i], i, i - (m.nu + 1) * (m.nv + 1)); } // make meshPoint
            }
            else { translationVectorW = translationVectorW * 0; } // 3. if

            // 4. if: Make movedNode if node is on face or inside brep (if on edge, movedNode already made).
            if (IsOnFace | !IsOnEdge)
            {
                double overlapTolerance = 0.99; // ensure no collision of vertices, reduce number to avoid "the look of triangles".
                movedNode = new Point3d(m.Nodes[i].Coordinate.X + (translationVectorU.X + translationVectorV.X + translationVectorW.X) * overlapTolerance,
                    m.Nodes[i].Coordinate.Y + (translationVectorU.Y + translationVectorV.Y + translationVectorW.Y) * overlapTolerance,
                    m.Nodes[i].Coordinate.Z + (translationVectorU.Z + translationVectorV.Z + translationVectorW.Z) * overlapTolerance);
                
                if (IsOnFace) // If node is on face: ensure it stays on face
                {
                    Brep srf = face.DuplicateFace(false);
                    movedNode = srf.ClosestPoint(movedNode); // "Project" meshPoint to surface.
                }
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
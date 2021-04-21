using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;

namespace MeshPoints.Classes
{
    class Mesh3D // to do: change to shell and solid
    {
        public List<Element> Elements { get; set; } //list of elements
        public List<Node> Nodes { get; set; } //list of nodes
        public Mesh mesh { get; set; } //mesh
        public int nu { get; set; } //number of nodes in x-dir
        public int nv { get; set; } //number of nodes in y-dir
        public int nw { get; set; } //number of nodes in z-dir
        public bool inp { get; set; }
        public string Type { get; set; } // to do: inplementer
        public Geometry Geometry { get; set; } // to do: temporary
        public Mesh3D()
        {
            //Empty constructor
        }

        public Mesh3D(int _nu, int _nv, List<Node> _nodes, List<Element> _elements, Mesh _mesh) // for shell mesh
        {
            nu = _nu;
            nv = _nv;
            nw = 1;
            Nodes = _nodes;
            Elements = _elements;
            mesh = _mesh;
            Type = "Surface";
        }

        public Mesh3D(int _nu, int _nv, int _nw, List<Node> _nodes, List<Element> _elements, Mesh _mesh) // for solid mesh
        {
            nu = _nu;
            nv = _nv;
            nw = _nw;
            Nodes = _nodes;
            Elements = _elements;
            mesh = _mesh;
            Type = "Solid";
        }

        public void SetQuadElements() 
        {
            List<Node> nodes = this.Nodes;
            int nu = this.nu;
            int nv = this.nv;

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
            this.Elements = elements;
        }

        public void SetHexElements()
        {
            int nu = this.nu;
            int nv = this.nv;
            int nw = this.nw;
            List<Node> nodes = this.Nodes;
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
            this.Elements = elements;
        }

        public void SetMesh()
        {
            int nu = this.nu;
            int nv = this.nv;
            int nw = this.nw;
            List<Node> nodes = this.Nodes;
            
            Mesh mesh = new Mesh();
            int counter;
            int jump;

            // Assign mesh vertices from node coordinates
            foreach (Node node in nodes)
            {
                mesh.Vertices.Add(node.Coordinate);
            }

            // Mesh planes in u - dir
            counter = 0;
            jump = 0; // new w sequence
            for (int i = 0; i < nu; i++)
            {
                for (int j = 0; j < (nv - 1) * (nw - 1); j++)
                {
                    if (jump < nv - 1)
                    {

                        mesh.Faces.AddFace(counter, counter + (nu * nv), counter + (nu * nv) + nu, counter + nu);

                        counter += nu;
                        jump++;
                    }
                    else { counter += nu; jump = 0; j--; }
                }
                counter = (i + 1);
                jump = 0;
            }


            // Mesh planes in v - dir
            counter = 0;
            int counterU = 0;
            jump = 0; // new u sequence
            for (int i = 0; i < nv; i++)
            {
                for (int j = 0; j < (nu - 1) * (nw - 1); j++)
                {
                    if (jump < nw - 1)
                    {
                        mesh.Faces.AddFace(counter, counter + 1, counter + (nu * nv) + 1, counter + nu * nv);

                        counter += nu * nv;
                        jump++;
                    }
                    else { counterU++; counter = counterU; jump = 0; j--; }
                }
                counter = (i + 1) * nu;
                jump = 0;
                counterU = counter;
            }

            // Mesh planes in w - dir
            counter = 0;
            jump = 0; // new v sequence
            for (int i = 0; i < nw; i++)
            {
                for (int j = 0; j < (nu - 1) * (nv - 1); j++)
                {
                    if (jump < nu - 1)
                    {
                        mesh.Faces.AddFace(counter, counter + 1, counter + nu + 1, counter + nu);

                        counter++;
                        jump++;
                    }
                    else { counter++; jump = 0; j--; }
                }
                counter = (i + 1) * nu * nv;
                jump = 0;
            }
            mesh.Normals.ComputeNormals();  //Control if needed
            mesh.FaceNormals.ComputeFaceNormals();  //want a consistant mesh
            mesh.Compact(); //to ensure that it calculate

            this.mesh = mesh; ;
        }
    }
}

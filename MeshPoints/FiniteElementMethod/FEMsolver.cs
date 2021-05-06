using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Drawing;
using Rhino;

// Csparse
using LA = MathNet.Numerics.LinearAlgebra;
using CSparse;
using CSD = CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;


namespace MeshPoints.FiniteElementMethod
{
    public class FEMsolver : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMsolver class.
        /// </summary>
        public FEMsolver()
          : base("FEM solver", "FEM",
              "Finite element method solver with quad 4 and hex 8 elements.",
              "MyPlugIn", "FEM")
        { 
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "smartMesh", "Input a SmartMesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Loads", "loads", "Input a load vector", GH_ParamAccess.list);
            pManager.AddGenericParameter("Boundary conditions", "BC", "Input a boundary condition vector", GH_ParamAccess.list);
            pManager.AddGenericParameter("Material", "material", "Input a list of material sorted: Young modulus, Poisson Ratio", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("u1", "disp", "Displacement of nodes in u1 dir", GH_ParamAccess.list);
            pManager.AddGenericParameter("u2", "disp", "Displacement of nodes in u2 dir", GH_ParamAccess.list);
            pManager.AddGenericParameter("u3", "disp", "Displacement of node in u3 dir", GH_ParamAccess.list);
            pManager.AddGenericParameter("Nodal stress", "node stress", "Stress at nodes", GH_ParamAccess.list);
            pManager.AddGenericParameter("Mises stress", "mises", "Calculate mises stress at nodes", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Input
            SmartMesh mesh = new SmartMesh(); // to do: change to MeshGeometry elns
            List<double> loads = new List<double>();
            List<List<int>> boundaryConditions = new List<List<int>>();
            Material material = new Material();

            DA.GetData(0, ref mesh);
            DA.GetDataList(1, loads);
            DA.GetDataList(2, boundaryConditions);
            DA.GetData(3, ref material);
            #endregion

            #region Code
            List<Node> nodes = mesh.Nodes;
            List<Element> elements = mesh.Elements;
            int numNodes = nodes.Count;
            int nodeDOFS = 0;

            // 1. Check if mesh is Surface or Solid
            if (String.Equals(mesh.Type, "Surface"))  { nodeDOFS = 2;}
            else if (String.Equals( mesh.Type,"Solid")) { nodeDOFS = 3;}
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid mesh: Need to spesify if mesh is surface or solid."); }

            // 2. Get global stiffness matrix
            LA.Matrix<double> K_global = CalculateGlobalStiffnessMatrix(elements, numNodes, nodeDOFS, material);

            // 3. Get load vector
            LA.Matrix<double> R = LA.Double.DenseMatrix.Build.Dense(numNodes * nodeDOFS, 1);
            for (int i = 0; i < loads.Count; i++)
            {
                R[i, 0] = loads[i];
            }

            // Fic BoundaryConditions:
            if (boundaryConditions.Count > mesh.Nodes.Count) { boundaryConditions = FixBoundaryConditions(boundaryConditions, mesh.Nodes.Count); }
            
            // 4. Calculate displacement 
            LA.Matrix<double> u = CalculateDisplacement(K_global, R, boundaryConditions); 

            var stress = CalculateGlobalStress(elements, u, material, nodeDOFS);
            LA.Matrix<double> globalStress = stress.Item1;
            LA.Vector<double> mises = stress.Item2;

            ColorMeshAfterStress(mesh, mises, material);

            // prepare output
            double[] nodalDeformation = u.Column(0).ToArray();

            List<double> u1 = new List<double>();
            List<double> u2 = new List<double>();
            List<double> u3 = new List<double>();
            List<double> nodalMises = new List<double>();

            for (int i = 0; i < mesh.Nodes.Count; i++)
            {
                u1.Add(u[i * nodeDOFS, 0]);
                u2.Add(u[i * nodeDOFS + 1, 0]);
                u3.Add(u[i * nodeDOFS + 2, 0]);

                nodalMises.Add(mises[i]);
            }

            List<double[]> nodalStress = new List<double[]>();
            for (int i = 0; i < globalStress.ColumnCount; i++)
            {
               nodalStress.Add(globalStress.Column(i).ToArray());
            }
            #endregion

            #region Output
            DA.SetDataList(0, u1);
            DA.SetDataList(1, u2);
            DA.SetDataList(2, u3);
            DA.SetDataList(3, nodalStress);
            DA.SetDataList(4, nodalMises);
            #endregion 
        }

        #region Methods
        private List<List<int>> FixBoundaryConditions(List<List<int>> boundaryConditions, int numNodes)
        {
            List<List<int>> totalBC = new List<List<int>>();
            for (int i = 0; i < numNodes; i++) // loop number nodes
            {
                List<int> dofList = new List<int>(boundaryConditions[i]);  // get dofList of first input list of BC
                for (int j = 0; j < dofList.Count; j++) // loop dofs
                {
                    for (int k = 1; k < boundaryConditions.Count / numNodes; k++) // loop the remaining inout list of BC 
                    {
                        dofList[j] = dofList[j] + boundaryConditions[i + k * numNodes][j];
                    }
                }
                totalBC.Add(dofList);
            }
            return totalBC;
        }

        private Tuple<LA.Matrix<double>, List<LA.Matrix<double>>> Synne(List<Node> nodeList, Material material)
        {
            LA.Matrix<double> Ke = LA.Matrix<double>.Build.Dense(24, 24);
            List<LA.Matrix<double>> Be = new List<LA.Matrix<double>>();

            //3D Constitutive LA.Matrix: C
            double E = material.YoungModulus;
            double nu = material.PossionRatio;
            double value = (double)E / ((1 + nu) * (1 - 2 * nu));
            LA.Matrix<double> C = LA.Double.DenseMatrix.OfArray(new double[,]
            {
                {1-nu, nu, nu, 0, 0, 0},
                {nu, 1-nu, nu, 0, 0, 0},
                {nu, nu, 1-nu, 0, 0, 0},
                {0, 0, 0, (1-2*nu)/2, 0, 0},
                {0, 0, 0, 0, (1-2*nu)/2, 0},
                {0, 0, 0, 0, 0, (1-2*nu)/2},
            });

            C = C.Multiply(value); //Constitutive LA.Matrix

            //Gauss points
            LA.Vector<double> gaussPoints = DenseVector.OfArray(new double[] { -1 / Math.Sqrt(3), 1 / Math.Sqrt(3) }); //Gauss points

            Point3d point = new Point3d(0, 0, 0);
            List<Point3d> pNatural = new List<Point3d>();

            for (int i = 0; i < nodeList.Count; i++)
            {
                pNatural.Add( nodeList[i].Coordinate);
            }

            LA.Matrix<double> coordinates = LA.Double.DenseMatrix.OfArray(new double[,]
            {
                {pNatural[0].X,pNatural[0].Y , pNatural[0].Z},
                {pNatural[1].X,pNatural[1].Y , pNatural[1].Z},
                {pNatural[2].X,pNatural[2].Y , pNatural[2].Z},
                {pNatural[3].X,pNatural[3].Y , pNatural[3].Z},
                {pNatural[4].X,pNatural[4].Y , pNatural[4].Z},
                {pNatural[5].X,pNatural[5].Y , pNatural[5].Z},
                {pNatural[6].X,pNatural[6].Y , pNatural[6].Z},
                {pNatural[7].X,pNatural[7].Y , pNatural[7].Z},
            });

            //Numerical integration
            foreach (double g3 in gaussPoints)
            {
                foreach (double g2 in gaussPoints)
                {
                    foreach (double g1 in gaussPoints)
                    {
                        //Shape functions
                        LA.Matrix<double> shapeF = LA.Double.DenseMatrix.OfArray(new double[,]
                       {
                            {-(1-g2)*(1-g3), (1-g2)*(1-g3), (1+g2)*(1-g3),-(1+g2)*(1-g3),-(1-g2)*(1+g3),(1-g2)*(1+g3),(1+g2)*(1+g3),-(1+g2)*(1+g3)},
                            {-(1-g1)*(1-g3), -(1+g1)*(1-g3), (1+g1)*(1-g3),(1-g1)*(1-g3),-(1-g1)*(1+g3),-(1+g1)*(1+g3),(1+g1)*(1+g3),(1-g1)*(1+g3)},
                            {-(1-g1)*(1-g2), -(1+g1)*(1-g2), -(1+g1)*(1+g2),-(1-g1)*(1+g2),(1-g1)*(1-g2),(1+g1)*(1-g2),(1+g1)*(1+g2),(1-g1)*(1+g2)},

                       });

                        shapeF = shapeF.Divide(8); //Divided by 8

                        //Jacobi Matrix

                        LA.Matrix<double> JacobiMatrix = LA.Matrix<double>.Build.Dense(3, 3);

                        JacobiMatrix = shapeF.Multiply(coordinates);

                        // Auxiliar LA.Matrix for assemblinng of B-matrix 
                        LA.Matrix<double> auxiliar = LA.Matrix<double>.Build.Dense(3, 8);

                        auxiliar = JacobiMatrix.Inverse().Multiply(shapeF);

                        // B matrix
                        LA.Matrix<double> B = LA.Matrix<double>.Build.Dense(6, 24);

                        //First three rows
                        for (int i = 0; i < 3; i++)
                        {
                            for (int j = 0; j <= 7; j++)
                            {
                                B[i, 3 * j + 1 + (i - 1)] = auxiliar[i, j];
                            }
                        }

                        //Fourth row
                        for (int j = 0; j <= 7; j++)
                        {
                            B[3, 3 * j] = auxiliar[1, j];
                        }

                        for (int j = 0; j <= 7; j++)
                        {
                            B[3, 3 * j + 1] = auxiliar[0, j];
                        }

                        //Fifth row
                        for (int j = 0; j <= 7; j++)
                        {
                            B[4, 3 * j + 2] = auxiliar[1, j];
                        }

                        for (int j = 0; j <= 7; j++)
                        {
                            B[4, 3 * j + 1] = auxiliar[2, j];
                        }

                        //Sixth row
                        for (int j = 0; j <= 7; j++)
                        {
                            B[5, 3 * j] = auxiliar[2, j];
                        }

                        for (int j = 0; j <= 7; j++)
                        {
                            B[5, 3 * j + 2] = auxiliar[0, j];
                        }

                        Be.Add(B);

                        //Adding the stiffness matrix. Ke = Ke + B'*C*B*Det(JacobiMatrix)
                        Ke = Ke.Add(B.Transpose().Multiply(C).Multiply(B).Multiply(JacobiMatrix.Determinant()));
                    }
                }
            }

            //Changing order of Be to fit the global numbering
            LA.Matrix<double> B_2 = Be[2];
            Be[2] = Be[3];
            Be[3] = B_2;
            LA.Matrix<double> B_6 = Be[6];
            Be[6] = Be[7];
            Be[7] = B_6;

            bool sym = Ke.IsSymmetric();
            return Tuple.Create(Ke, Be);
        } // to do: slett

        private Tuple<LA.Matrix<double>, List<LA.Matrix<double>>> CalculateElementMatrices(Element element, Material material, int nodeDOFS)
        {
            // summary: calculate local K and B matrix

            // material
            LA.Matrix<double> C = GetMaterialConstant(material.YoungModulus, material.PossionRatio, nodeDOFS);

            // shapefunction
            FEM _FEM = new FEM();

            // create local stiffness matrix
            int numElementNodes = element.Nodes.Count;
            LA.Matrix<double> K_local = LA.Matrix<double>.Build.Dense(nodeDOFS * numElementNodes, nodeDOFS * numElementNodes);

            // create local deformation matrix
            List<LA.Matrix<double>> B_local = new List<LA.Matrix<double>>();

            // Global coordinates of the corner nodes of the actual element
            LA.Matrix<double> globalCoordinates = LA.Matrix<double>.Build.Dense(numElementNodes, nodeDOFS);
            for (int i = 0; i < numElementNodes; i++)
            {
                globalCoordinates[i, 0] = element.Nodes[i].Coordinate.X; // column of x coordinates
                globalCoordinates[i, 1] = element.Nodes[i].Coordinate.Y; // column of y coordinates
                if (nodeDOFS == 3)
                {
                    globalCoordinates[i, 2] = element.Nodes[i].Coordinate.Z; // colum of z coordinates
                }
            }

            //Numerical integration
           LA.Matrix<double> gaussNodes = _FEM.GetGaussPoints((double)Math.Sqrt((double)1 / (double)3), nodeDOFS);

           for (int n = 0; n < gaussNodes.RowCount; n++)  // loop gauss nodes
            {
                // Substitute the natural coordinates into the symbolic expression
                var r = gaussNodes.Row(n)[0];
               var s = gaussNodes.Row(n)[1];
               double t = 0;
               if (nodeDOFS == 3) { t = gaussNodes.Row(n)[2]; }

               // Partial derivatives of the shape functions
               LA.Matrix<double> shapeFunctionsDerivatedNatural = _FEM.DerivateWithNatrualCoordinates(r, s, t, nodeDOFS); 

              // Calculate Jacobian matrix
              LA.Matrix<double> jacobianMatrix = shapeFunctionsDerivatedNatural.Multiply(globalCoordinates);

               // Calculate B - LA.Matrix
               LA.Matrix<double> shapeFuncDerivatedCartesian = jacobianMatrix.Inverse().Multiply(shapeFunctionsDerivatedNatural);

                double checkDet = jacobianMatrix.Determinant();
                if (checkDet < 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Negativ jac det"); }
               int dimRowB = 0;
               if (nodeDOFS == 2) { dimRowB = 3; }
               else { dimRowB = 6; }

               LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense( dimRowB , nodeDOFS*numElementNodes);

                for (int i = 0; i < numElementNodes; i++)
                {
                    for (int j = 0; j < nodeDOFS; j++)
                    {
                        if (nodeDOFS == 2) // surface
                        {
                            if (j == 0)
                            {
                                B_i[0, 2 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                                B_i[2, 2 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                            }
                            else if (j == 1)
                            {
                                B_i[1, j + 2 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                                B_i[2, j + 2 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                            }
                        }
                        else // solid
                        {
                            if (j == 0)
                            {
                                B_i[0, 3 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                                B_i[4, 3 * i] = shapeFuncDerivatedCartesian.Row(2)[i];
                                B_i[5, 3 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                            }
                            else if (j == 1)
                            {
                                B_i[1, j + 3 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                                B_i[3, j + 3 * i] = shapeFuncDerivatedCartesian.Row(2)[i];
                                B_i[5, j + 3 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                            }
                            else if (j == 2)
                            {
                                B_i[2, j + 3 * i] = shapeFuncDerivatedCartesian.Row(2)[i];
                                B_i[3, j + 3 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                                B_i[4, j + 3 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                            }
                        }
                    }
                }
                

            B_local.Add(B_i);
            K_local = K_local + B_i.Transpose().Multiply(C).Multiply(B_i).Multiply(jacobianMatrix.Determinant());
           }

            // to do: check
            
            // Matrix<double> B_2 = B_local[2];
            //B_local[2] = B_local[3];
           // B_local[3] = B_2;
           // Matrix<double> B_6 = B_local[6];
            //B_local[6] = B_local[7];
           // B_local[7] = B_6;

            return Tuple.Create(K_local, B_local);
        }

        private Tuple<LA.Matrix<double>, List<LA.Matrix<double>>> Magnus(Element e, Material material)
        {
            // material
            LA.Matrix<double> C = GetMaterialConstant(material.YoungModulus, material.PossionRatio, 3);

            // shapefunction
            FEM _FEM = new FEM();
            int nodeDOFS = 3;
            // create local stiffness matrix
            int numElementNodes = e.Nodes.Count;
            LA.Matrix<double> K_local = LA.Matrix<double>.Build.Dense(3 * numElementNodes, 3 * numElementNodes);

            // create local deformation matrix
            List<LA.Matrix<double>> B_local = new List<LA.Matrix<double>>();

            // Global X, Y, and Z-coordinates of the corner nodes of the actual element
            List<double> gX = new List<double>();
            List<double> gY = new List<double>();
            List<double> gZ = new List<double>();

            for (int i = 0; i < 8; i++)
            {
                gX.Add(e.Nodes[i].Coordinate.X);
                gY.Add(e.Nodes[i].Coordinate.Y);
                gZ.Add(e.Nodes[i].Coordinate.Z);

            }


            //Numerical integration
            //Matrix<double> gaussNodes = _FEM.GetGaussPoints((double)Math.Sqrt((double)1 / (double)3), nodeDOFS);

            // Make to gauss points...
            // Corner nodes of the isoparametric element
            var naturalNodes = new List<List<Double>>
            {
                new List<double> { -1, -1, -1 }, new List<double> { 1, -1, -1}, new List<double> { 1, 1, -1 }, new List<double> { -1, 1, -1 },
                new List<double> { -1, -1, 1 }, new List<double> { 1, -1, 1 }, new List<double> { 1, 1, 1 }, new List<double> { -1, 1, 1 }
            };


            List<double> jacobiansOfElement = new List<double>();
            foreach (List<Double> node in naturalNodes)
            {
                // Substitute the natural coordinates into the symbolic expression
                var r = node[0] * Math.Sqrt((double)1 / (double)3);
                var s = node[1] * Math.Sqrt((double)1 / (double)3);
                var t = node[2] * Math.Sqrt((double)1 / (double)3);

                // Partial derivatives of the shape functions
                var N1Dr = -0.125 * (s - 1) * (t - 1);
                var N1Ds = -0.125 * (r - 1) * (t - 1);
                var N1Dt = -0.125 * (r - 1) * (s - 1);
                var N2Dr = 0.125 * (s - 1) * (t - 1);
                var N2Ds = 0.125 * (r + 1) * (t - 1);
                var N2Dt = 0.125 * (r + 1) * (s - 1);
                var N3Dr = -0.125 * (s + 1) * (t - 1);
                var N3Ds = -0.125 * (r + 1) * (t - 1);
                var N3Dt = -0.125 * (r + 1) * (s + 1);
                var N4Dr = 0.125 * (s + 1) * (t - 1);
                var N4Ds = 0.125 * (r - 1) * (t - 1);
                var N4Dt = 0.125 * (r - 1) * (s + 1);
                var N5Dr = 0.125 * (s - 1) * (t + 1);
                var N5Ds = 0.125 * (r - 1) * (t + 1);
                var N5Dt = 0.125 * (r - 1) * (s - 1);
                var N6Dr = -0.125 * (s - 1) * (t + 1);
                var N6Ds = -0.125 * (r + 1) * (t + 1);
                var N6Dt = -0.125 * (r + 1) * (s - 1);
                var N7Dr = 0.125 * (s + 1) * (t + 1);
                var N7Ds = 0.125 * (r + 1) * (t + 1);
                var N7Dt = 0.125 * (r + 1) * (s + 1);
                var N8Dr = -0.125 * (s + 1) * (t + 1);
                var N8Ds = -0.125 * (r - 1) * (t + 1);
                var N8Dt = -0.125 * (r - 1) * (s + 1);

                var sfDr = new List<double>
                    {
                        N1Dr, N2Dr, N3Dr, N4Dr, N5Dr, N6Dr, N7Dr, N8Dr
                    };
                var sfDs = new List<double>
                    {
                        N1Ds, N2Ds, N3Ds, N4Ds, N5Ds, N6Ds, N7Ds, N8Ds
                    };
                var sfDt = new List<double>
                    {
                        N1Dt, N2Dt, N3Dt, N4Dt, N5Dt, N6Dt, N7Dt, N8Dt
                    };

                // Evaluates each partial derivative in the isoparametric node
                var calcDerivs = new List<Double>
                    {
                        MultiplyLists(gX, sfDr),
                        MultiplyLists(gX, sfDs),
                        MultiplyLists(gX, sfDt),

                        MultiplyLists(gY, sfDr),
                        MultiplyLists(gY, sfDs),
                        MultiplyLists(gY, sfDt),

                        MultiplyLists(gZ, sfDr),
                        MultiplyLists(gZ, sfDs),
                        MultiplyLists(gZ, sfDt)
                    };

                // Helper function to piecewise multiply elements of two lists of length 8
                double MultiplyLists(List<double> a, List<double> b)
                {
                    double sum = 0.0;
                    for (int i = 0; i < 8; i++)
                    {
                        sum += (a[i] * b[i]);
                    }
                    return sum;
                }

                // Structure data in the form of a Jacobian matrix
                LA.Matrix<double> jacobianMatrix = LA.DenseMatrixModule.ofArray2(new double[,]
                {
                        {calcDerivs[0], calcDerivs[3], calcDerivs[6] },
                        {calcDerivs[1], calcDerivs[4], calcDerivs[7] },
                        {calcDerivs[2], calcDerivs[5], calcDerivs[8] },
                });

                var jacobianDeterminant = jacobianMatrix.Determinant();
                jacobiansOfElement.Add(Math.Abs(jacobianDeterminant));

                // Partial derivatives of the shape functions
                LA.Matrix<double> shapeFunctionsDerivatedNatural = _FEM.DerivateWithNatrualCoordinates(r, s, t, 3);

                // Calculate Jacobian matrix
                //Matrix<double> jacobianMatrix = shapeFunctionsDerivatedNatural.Multiply(globalCoordinates);

                // Structure data in the form of a Jacobian matrix
                //Matrix<double> jacobianMatrix = calcDerivs.Transpose(); // to do: riktig ?? slik magnus hadde det

                // Calculate B - matrix
                LA.Matrix<double> shapeFuncDerivatedCartesian = jacobianMatrix.Inverse().Multiply(shapeFunctionsDerivatedNatural);

                int dimRowB = 0;
                if (nodeDOFS == 2) { dimRowB = 3; }
                else { dimRowB = 6; }

                LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense(dimRowB, nodeDOFS * numElementNodes);

                for (int i = 0; i < numElementNodes; i++)
                {
                    for (int j = 0; j < nodeDOFS; j++)
                    {
                        if (nodeDOFS == 2) // surface
                        {
                            if (j == 0)
                            {
                                B_i[0, 2 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                                B_i[2, 2 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                            }
                            else if (j == 1)
                            {
                                B_i[1, j + 2 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                                B_i[2, j + 2 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                            }
                        }
                        else // solid
                        {
                            if (j == 0)
                            {
                                B_i[0, 3 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                                B_i[4, 3 * i] = shapeFuncDerivatedCartesian.Row(2)[i];
                                B_i[5, 3 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                            }
                            else if (j == 1)
                            {
                                B_i[1, j + 3 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                                B_i[3, j + 3 * i] = shapeFuncDerivatedCartesian.Row(2)[i];
                                B_i[5, j + 3 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                            }
                            else if (j == 2)
                            {
                                B_i[2, j + 3 * i] = shapeFuncDerivatedCartesian.Row(2)[i];
                                B_i[3, j + 3 * i] = shapeFuncDerivatedCartesian.Row(1)[i];
                                B_i[4, j + 3 * i] = shapeFuncDerivatedCartesian.Row(0)[i];
                            }
                        }
                    }
                }
                    B_local.Add(B_i);
                    K_local = K_local + B_i.Transpose().Multiply(C).Multiply(B_i).Multiply(jacobianMatrix.Determinant());
            }
            return Tuple.Create(K_local, B_local);
        }

        private LA.Matrix<double> CalculateGlobalStiffnessMatrix(List<Element> elements, int numNode, int nodeDOFS, Material material)
        {
            // create stiffness matrix
            LA.Matrix<double> K_global = LA.Matrix<double>.Build.Dense(numNode * nodeDOFS, numNode * nodeDOFS);
            foreach (Element element in elements)
            {
                List<int> con = element.Connectivity;
                LA.Matrix<double> K_local = CalculateElementMatrices(element, material, nodeDOFS).Item1;
                //LA.Matrix<double> K_local = Synne(element.Nodes, material).Item1;
                //LA.Matrix<double> K_local = Magnus(element, material).Item1;

                // loop nodes of elements
                for (int i = 0; i < con.Count; i++)
                {
                    for (int j = 0; j < con.Count; j++)
                    {
                        // loop relevant local stiffness contribution
                        for (int dofRow = 0; dofRow < nodeDOFS; dofRow++)
                        {
                            for (int dofCol = 0; dofCol < nodeDOFS; dofCol++) 
                            {
                                K_global[nodeDOFS * con[i] + dofRow, nodeDOFS * con[j] + dofCol] = K_global[nodeDOFS * con[i] + dofRow, nodeDOFS * con[j] + dofCol] + K_local[nodeDOFS * i + dofRow , nodeDOFS * j + dofCol];
                            }
                        }
                    }
                }
       
            }
            return K_global;
        }

        private LA.Matrix<double> CalculateDisplacement(LA.Matrix<double> K_global, LA.Matrix<double> R, List<List<int>> applyBCToDOF)
        {
            // summary: include boundary condistions and calculate global displacement
            
            // Make list of boundary condistions
            List<int> BCList = new List<int>();
            for (int i = 0; i < applyBCToDOF.Count; i++)
            {
                for (int j = 0; j < applyBCToDOF[0].Count; j++)
                {
                    BCList.Add(applyBCToDOF[i][j]);
                }
            }

            for (int i = 0; i < BCList.Count; i++)
            {
                for (int j = 0; j < BCList.Count; j++)
                {

                    if (BCList[i] == 1)
                    {
                        if (i != j)
                        {
                            K_global[i, j] = 0;
                        }
                        else
                        {
                            K_global[i, j] = 1;
                            R[i, 0] = 0;
                        }
                    }
                }
            }

            // Time recorder
            var sw0 = new System.Diagnostics.Stopwatch();
            var sw1 = new System.Diagnostics.Stopwatch();

            // Mathnet.Numerics to CSparse

            var b = R.Column(0);
            sw0.Start();

            var CMA = K_global.Storage.ToColumnMajorArray();
            CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_global.RowCount, K_global.ColumnCount, CMA);

            SparseLU CS_K_global = SparseLU.Create(CCS, ColumnOrdering.MinimumDegreeAtPlusA, 0.0);
            double[] CS_u = CSD.Vector.Create(K_global.RowCount * 1, 0.0);
            double[] CS_R = R.Column(0).ToArray();


            sw0.Stop();
            Rhino.RhinoApp.WriteLine($"Elapsed [msec] = " + sw0.Elapsed.TotalMilliseconds);
            sw0.Start();

            sw1.Restart();
            CS_K_global.Solve(CS_R, CS_u);
            sw1.Stop();
            sw0.Stop();
            Rhino.RhinoApp.WriteLine($"### {K_global.RowCount} x {K_global.ColumnCount} Matrix. CSparse Elapsed [msec] = {sw1.Elapsed.TotalMilliseconds}");

            // CSparse to Mathnet.Numerics
            LA.Matrix<double> u = LA.Double.DenseMatrix.OfColumnArrays(CS_u); 
            return u;  
        }

        private Tuple<LA.Matrix<double>, LA.Matrix<double>> CalculateElementStrainStress(Element element, LA.Matrix<double> u, Material material, int nodeDOFS)
        {
            // summary: calculate a list of strain and stress vectors for each node in a element.
            LA.Matrix<double> C = GetMaterialConstant(material.YoungModulus, material.PossionRatio, nodeDOFS);

            FEM _FEM = new FEM();
            List<LA.Matrix<double>> B_local = CalculateElementMatrices(element, material, nodeDOFS).Item2;
            LA.Matrix<double> elementGaussStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementGaussStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> localDeformation = LA.Double.DenseMatrix.Build.Dense(nodeDOFS * B_local.Count,1);
            
            // get deformation of nodes connected to element
            for (int i = 0; i < element.Connectivity.Count; i++)
            {
                localDeformation[nodeDOFS * i, 0] = u[nodeDOFS * element.Connectivity[i],0];
                localDeformation[nodeDOFS * i + 1, 0] = u[nodeDOFS * element.Connectivity[i] + 1, 0];
                localDeformation[nodeDOFS * i + 2, 0] = u[nodeDOFS * element.Connectivity[i] + 2, 0];
            }
            // get gauss strain and stress
            for (int n = 0; n < B_local.Count; n++)
            {
                // B-matrix is calculated from gauss points
                LA.Matrix<double> gaussStrain = B_local[n].Multiply(localDeformation);
                LA.Matrix<double> gaussStress = C.Multiply(B_local[n]).Multiply(localDeformation);

                for (int i = 0; i < B_local[0].RowCount; i++)
                {
                    elementGaussStrain[i, n] = gaussStrain[i,0];
                    elementGaussStress[i, n] = gaussStress[i,0];
                }
            }

            // get node strain and stress by extrapolation
            LA.Matrix<double> extrapolationNodes = _FEM.GetGaussPoints(Math.Sqrt(3), nodeDOFS);

            for (int n = 0; n < B_local.Count; n++)
            { 
                // get stress and strain in nodes
                var r = extrapolationNodes.Row(n)[0];
                var s = extrapolationNodes.Row(n)[1];
                double t = 0;
                if (nodeDOFS == 3) { t = extrapolationNodes.Row(n)[2]; }

                LA.Vector<double> shapefunctionValuesInNode = _FEM.GetShapeFunctions(r, s, t, nodeDOFS);
                LA.Vector<double> nodeStrain = elementGaussStrain.Multiply(shapefunctionValuesInNode);
                LA.Vector<double> nodeStress = elementGaussStress.Multiply(shapefunctionValuesInNode);
                for (int i = 0; i < B_local[0].RowCount; i++)
                {
                    elementStrain[i, n] = nodeStrain[i];
                    elementStress[i, n] = nodeStress[i];
                }
            }
            return Tuple.Create(elementStrain, elementStress);
        }

        private Tuple<LA.Matrix<double>, LA.Vector<double>> CalculateGlobalStress(List<Element> elements, LA.Matrix<double> u, Material material, int nodeDOFS)
        {
            int numNodes =  u.RowCount / 3;
            int stressRowDim = 4;
            if (nodeDOFS == 3) { stressRowDim = 6; }
            LA.Matrix<double> globalStress = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);
            LA.Matrix<double> counter = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);


            foreach (Element element in elements)
            {
                LA.Matrix<double> elementStress = CalculateElementStrainStress(element, u, material, nodeDOFS).Item2;

                List<int> connectivity = element.Connectivity;

                for (int i = 0; i < elementStress.RowCount; i++) // loop the stress
                {
                    for (int j = 0; j < elementStress.ColumnCount; j++) // loop the element nodes
                    {
                        globalStress[i, connectivity[j]] = globalStress[i, connectivity[j]] + elementStress[i, j];
                        counter[i, connectivity[j]]++;
                    }
                }
            }

            // get average
            for (int i = 0; i < globalStress.RowCount; i++) // loop the stress
            {
                for (int j = 0; j < globalStress.ColumnCount; j++) // loop the element nodes
                {
                    if (counter[i, j] > 1)
                    {
                        globalStress[i, j] = globalStress[i, j] / (double)counter[i, j];
                        counter[i, j] = 0;
                    }
                }
            }

            // Mises
            LA.Vector<double> mises = DenseVector.Build.Dense(numNodes);
            for (int i = 0; i < numNodes; i++)
            {
                if (nodeDOFS == 2)
                {
                    LA.Vector<double> nodeStress = globalStress.Column(i);
                    double Sxx = nodeStress[0];
                    double Syy = nodeStress[1];
                    double Sxy = nodeStress[2];
                    mises[i] = Math.Sqrt( Math.Pow(Sxx, 2) - Sxx * Syy + Math.Pow(Syy, 2) + 3 * Math.Pow(Sxy, 2));
                }
                else
                {
                    LA.Vector<double> nodeStress = globalStress.Column(i);
                    double Sxx = nodeStress[0];
                    double Syy = nodeStress[1];
                    double Szz = nodeStress[2];
                    double Sxy = nodeStress[3];
                    double Sxz = nodeStress[4];
                    double Syz = nodeStress[5];
                    mises[i] = Math.Sqrt(0.5 * (Math.Pow(Sxx - Syy, 2) + Math.Pow(Syy - Szz, 2) + Math.Pow(Szz - Sxx, 2)) + 3 * (Math.Pow(Sxy, 2) + Math.Pow(Sxz, 2) + Math.Pow(Syz, 2)));
                }
            }
            return Tuple.Create(globalStress, mises);
        }

        private LA.Matrix<double> GetMaterialConstant(double youngModulus, double possionRatio, int nodeDOFS )
        {
            if (nodeDOFS == 2)
            {
                LA.Matrix<double> C = LA.Double.DenseMatrix.OfArray(new double[,]
                {
                    {1, possionRatio, 0},
                    {possionRatio, 1, 0},
                    {0, 0, (1-possionRatio)/2}
                });
                C.Multiply((double) youngModulus / (1 - Math.Pow( possionRatio, 2)));
                return C;
            }
            else
            {
                LA.Matrix<double> C = LA.Double.DenseMatrix.OfArray(new double[,]
                {
                    {1-possionRatio, possionRatio, possionRatio, 0, 0, 0},
                    {possionRatio, 1-possionRatio, possionRatio, 0, 0, 0},
                    {possionRatio, possionRatio, 1- possionRatio, 0, 0, 0},
                    {0, 0, 0, (1-2*possionRatio)/(double)2, 0, 0},
                    {0, 0, 0, 0, (1-2*possionRatio)/(double)2, 0},
                    {0, 0, 0, 0, 0, (1-2*possionRatio)/(double)2},
                });
                C = C.Multiply((double)youngModulus / (double)((1 + possionRatio) * (1 - 2 * possionRatio)));
                return C;
            }
        }

        private void ColorMeshAfterStress(SmartMesh mesh, LA.Vector<double> mises, Material material)
        {
            double maxValue = material.YieldingStress;
            double minValue = 0;
            Color color = Color.White;

            double range = (maxValue - minValue) / (double) 13;
            for (int i = 0; i < mesh.Nodes.Count; i++)
            {
                // to do: change, for likt synne
                if (mises[i] < minValue + range) color = Color.Blue;
                else if (mises[i] < minValue + 2 * range) color = Color.RoyalBlue;
                else if (mises[i] < minValue + 3 * range) color = Color.DeepSkyBlue;
                else if (mises[i] < minValue + 4 * range) color = Color.Cyan;
                else if (mises[i] < minValue + 5 * range) color = Color.PaleGreen;
                else if (mises[i] < minValue + 6 * range) color = Color.LimeGreen;
                else if (mises[i] < minValue + 7 * range) color = Color.Lime;
                else if (mises[i] < minValue + 8 * range) color = Color.Lime;
                else if (mises[i] < minValue + 9 * range) color = Color.GreenYellow;
                else if (mises[i] < minValue + 10 * range) color = Color.Yellow;
                else if (mises[i] < minValue + 11 * range) color = Color.Orange;
                else if (mises[i] < minValue + 12 * range) color = Color.OrangeRed;
                else color = Color.Red;

                mesh.Mesh.VertexColors.Add(color);
            }
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
            get { return new Guid("a82cb774-ef88-487c-bbe2-a283b76cc7bc"); }
        }
    }
}
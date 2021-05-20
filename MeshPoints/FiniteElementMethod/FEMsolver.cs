using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Drawing;

// Csparse
using LA = MathNet.Numerics.LinearAlgebra;
using CSparse;
using CSD = CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;


namespace MeshPoints.FiniteElementMethod
{
    public class FEMSolver : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMsolver class.
        /// </summary>
        public FEMSolver()
          : base("FEM solver", "FEM",
              "Finite element method solver with quad 4 and hex 8 elements.",
              "SmartMesh", "FEM")
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
            pManager.AddGenericParameter("Element stress", "element stress", "Stress in elements", GH_ParamAccess.list);
            pManager.AddGenericParameter("Mises stress", "mises", "Calculate mises stress at nodes", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            SmartMesh smartMesh = new SmartMesh();
            List<double> loads = new List<double>();
            List<List<int>> boundaryConditions = new List<List<int>>();
            Material material = new Material();

            DA.GetData(0, ref smartMesh);
            DA.GetDataList(1, loads);
            DA.GetDataList(2, boundaryConditions);
            DA.GetData(3, ref material);
           


            // Code

            List<Node> nodes = smartMesh.Nodes;
            List<Element> elements = smartMesh.Elements;
            int numNodes = nodes.Count;
            int nodeDOFS = 0;

            // 1. Check if mesh is Surface or Solid
            if (String.Equals(smartMesh.Type, "Surface"))  { nodeDOFS = 2;}
            else if (String.Equals( smartMesh.Type,"Solid")) { nodeDOFS = 3;}
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid mesh: Need to spesify if mesh is surface or solid."); }

            // 2. Get global stiffness matrix
            LA.Matrix<double> K_global = CalculateGlobalStiffnessMatrix(elements, numNodes, nodeDOFS, material);

            // 3. Get load vector
            LA.Matrix<double> R = LA.Double.DenseMatrix.Build.Dense(numNodes * nodeDOFS, 1);
            for (int i = 0; i < loads.Count; i++)
            {
                R[i, 0] = loads[i];
            }

            // 5. Fix BoundaryConditions
            boundaryConditions = FixBoundaryConditions(boundaryConditions, smartMesh.Nodes.Count);
            //if (boundaryConditions.Count > mesh.Nodes.Count) { boundaryConditions = FixBoundaryConditions(boundaryConditions, mesh.Nodes.Count); } Hilde?? 
            
            // 6. Calculate displacement 
            LA.Matrix<double> u = CalculateDisplacement(K_global, R, boundaryConditions); 

            // 7. Calculate stress
            var stress = CalculateGlobalStress(elements, u, material, nodeDOFS);
            LA.Matrix<double> globalStress = stress.Item1;
            LA.Vector<double> mises = stress.Item2;
            LA.Vector<double> misesElement = stress.Item3;
            ColorMeshAfterStress(smartMesh, mises, material);

            // 8. Prepare output
            List<double> u1 = new List<double>();
            List<double> u2 = new List<double>();
            List<double> u3 = new List<double>();
            List<double> nodalMises = new List<double>();
            List<double> elementMises = new List<double>();

            for (int i = 0; i < smartMesh.Elements.Count; i++)
            {
                elementMises.Add(misesElement[i]);
            }

            for (int i = 0; i < smartMesh.Nodes.Count; i++)
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



            // Output
            DA.SetDataList(0, u1);
            DA.SetDataList(1, u2);
            DA.SetDataList(2, u3);
            DA.SetDataList(3, nodalStress);
            DA.SetDataList(4, elementMises);
            DA.SetDataList(5, nodalMises);

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
                    BCList.Add(applyBCToDOF[i][j]); // list of 0 and 1 values for boundary condition for dof; true = 1, false = 0
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

            // Reduce K_global and R
            var reducedData = ReduceKandR(K_global, R, BCList);
            LA.Matrix<double> K_global_red = reducedData.Item1;
            LA.Matrix<double> R_red = reducedData.Item2;

            // Time recorder
            var sw0 = new System.Diagnostics.Stopwatch();

            // Mathnet.Numerics to CSparse
            var CMA = K_global_red.Storage.ToColumnMajorArray();
            CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_global_red.RowCount, K_global_red.ColumnCount, CMA);

            SparseLU CS_K_global = SparseLU.Create(CCS, ColumnOrdering.MinimumDegreeAtPlusA, 0.0);
            double[] CS_u = CSD.Vector.Create(K_global_red.RowCount * 1, 0.0);
            double[] CS_R = R_red.Column(0).ToArray();
;
            sw0.Start();
            CS_K_global.Solve(CS_R, CS_u);
            sw0.Stop();
            Rhino.RhinoApp.WriteLine($"### {K_global_red.RowCount} x {K_global_red.ColumnCount} Matrix. CSparse Elapsed [msec] = {sw0.Elapsed.TotalMilliseconds}");

            // CSparse to Mathnet.Numerics
            LA.Matrix<double> u = LA.Double.DenseMatrix.OfColumnArrays(CS_u);

            // Get total displacement
            
            LA.Vector<double> insertVec = DenseVector.Build.Dense(1);

            for (int i = 0; i < BCList.Count; i++)
            {
                if (BCList[i] == 1)
                {
                    u = u.InsertRow(i, insertVec);
                }
            }

            return u;  
        }

        private Tuple<LA.Matrix<double>, LA.Matrix<double>> ReduceKandR(LA.Matrix<double> K_global, LA.Matrix<double> R, List<int> BC)
        {
            int removeIndex = 0;
            for (int i = 0; i < K_global.RowCount; i++)
            {
                if (BC[i] == 1)
                {
                    K_global = K_global.RemoveRow(removeIndex); 
                    K_global = K_global.RemoveColumn(removeIndex);
                    R = R.RemoveRow(removeIndex);
                    removeIndex--;
                }
                removeIndex++;
            }
            return Tuple.Create(K_global, R);
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

        private Tuple<LA.Matrix<double>, LA.Vector<double>, LA.Vector<double>> CalculateGlobalStress(List<Element> elements, LA.Matrix<double> u, Material material, int nodeDOFS)
        {
            int numNodes =  u.RowCount / 3;
            int stressRowDim = 4;
            if (nodeDOFS == 3) { stressRowDim = 6; }
            LA.Matrix<double> globalStress = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);
            LA.Matrix<double> counter = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);
            List<LA.Matrix<double>> elementStressList = new List<LA.Matrix<double>>(); 
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
                elementStressList.Add(elementStress);
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

            // Nodal Mises
            LA.Vector<double> mises = DenseVector.Build.Dense(numNodes);
            for (int i = 0; i < numNodes; i++)
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

            // Element mises
            LA.Vector<double> elementMises = DenseVector.Build.Dense(numNodes);
            for (int i = 0; i < elementStressList.Count; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    LA.Vector<double> nodeStress = elementStressList[i].Column(j);
                    double Sxx = nodeStress[0];
                    double Syy = nodeStress[1];
                    double Szz = nodeStress[2];
                    double Sxy = nodeStress[3];
                    double Sxz = nodeStress[4];
                    double Syz = nodeStress[5];
                    elementMises[i] = elementMises[i] + Math.Sqrt(0.5 * (Math.Pow(Sxx - Syy, 2) + Math.Pow(Syy - Szz, 2) + Math.Pow(Szz - Sxx, 2)) + 3 * (Math.Pow(Sxy, 2) + Math.Pow(Sxz, 2) + Math.Pow(Syz, 2)));
                }
                elementMises[i] = elementMises[i] / (double)8; // get average of nodal mises
            }


            return Tuple.Create(globalStress, mises, elementMises);
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
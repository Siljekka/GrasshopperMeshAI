using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;


namespace MeshPoints
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
            pManager.AddGenericParameter("Mesh", "solidmesh", "Input a SolidMesh", GH_ParamAccess.item); // to do: change name
            pManager.AddGenericParameter("Loads", "loads", "Input a load vector", GH_ParamAccess.item);
            pManager.AddGenericParameter("Boundary conditions", "BC", "Input a boundary condition vector", GH_ParamAccess.item);
            pManager.AddGenericParameter("Material", "material", "Input a list of material sorted: Young modulus, Poisson Ratio", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodal displacement", "disp", "Displacement of nodes", GH_ParamAccess.list);
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
            Mesh3D mesh = new Mesh3D(); // to do: change to MeshGeometry elns
            List<string> loads = new List<string>();
            List<string> boundaryConditions = new List<string>();
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

            // shell or solid
            if (mesh.Type == "shell")  { nodeDOFS = 2;}
            else if (mesh.Type == "solid") { nodeDOFS = 3;}
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid mesh: Need to spesify if mesh is shell or solid."); }

            // get global stiffness matrix
            Matrix<double> K_global = CalculateGlobalStiffnessMatrix(elements, numNodes, nodeDOFS, material);

            // get boundary conditions
            // get list of index of dofs that are fixed
            //List<List<int>> uIsZero =  new List<List<int>>(); // 0: false, 1: true

            List<int> uIsZero = new List<int>(); 

            // get load 
            Matrix<double> R = DenseMatrix.Build.Dense(numNodes * nodeDOFS, 1);

            // calculate displacement 
            Matrix<double> u = CalculateDisplacement(K_global, R, uIsZero); // fix: uIszero

            var stress = CalculateGlobalStress(elements, u, material, nodeDOFS);
            Matrix<double> globalStress = stress.Item1;
            Vector<double> mises = stress.Item2;

            // prepare output
            double[] nodalDeformation = u.Column(0).ToArray();
            double[] nodalMises = mises.ToArray();
            List<double[]> nodalStress = new List<double[]>();
            for (int i = 0; i < globalStress.ColumnCount; i++)
            {
                nodalStress.Add(globalStress.Column(i).ToArray());
            }
            #endregion 

            #region Output
            DA.SetDataList(0, nodalDeformation);
            DA.SetDataList(1, nodalStress);
            DA.SetDataList(2, nodalMises);
            #endregion 
        }

        #region Methods

        private Tuple<Matrix<double>, List<Matrix<double>>> CalculateElementMatrices(Element element, Material material, int nodeDOFS)
        {
            // summary: calculate local K and B matrix

            // material
            Matrix<double> C = GetMaterialConstant(material.YoungModulus, material.PossionRatio, nodeDOFS);

            // shapefunction
            FEM _FEM = new FEM();

            // create local stiffness matrix
            int numElementNodes = element.Nodes.Count;
            Matrix<double> K_local = Matrix<double>.Build.Dense(nodeDOFS * numElementNodes, nodeDOFS * numElementNodes);

            // create local deformation matrix
            List<Matrix<double>> B_local = new List<Matrix<double>>();

            // Gauss points
            Matrix<double> gaussNodes = _FEM.GetGaussPoints(Math.Sqrt(1 / 3), nodeDOFS);

            // Global coordinates of the corner nodes of the actual element
            Matrix<double> globalCoordinates = Matrix<double>.Build.Dense(numElementNodes, nodeDOFS);
            for (int i = 0; i < numElementNodes; i++)
            {
                globalCoordinates.Row(i)[0] = element.Nodes[i].Coordinate.X; // column of x coordinates
                globalCoordinates.Row(i)[1] = element.Nodes[i].Coordinate.Y; // column of y coordinates
                if (nodeDOFS == 3)
                {
                    globalCoordinates.Row(i)[2] = element.Nodes[i].Coordinate.Y; // colum of z coordinates
                }
            }

            // loop gauss nodes
            for (int n = 0; n < gaussNodes.RowCount; n++)
            {
                // Substitute the natural coordinates into the symbolic expression
                var r = gaussNodes.Row(n)[0];
                var s = gaussNodes.Row(n)[1];
                double t = 0;
                if (nodeDOFS == 3) { t = gaussNodes.Row(n)[2]; }

                // Partial derivatives of the shape functions
                Matrix<double> shapeFunctionsDerivatedNatural = _FEM.DerivateWithNatrualCoordinates(r, s, t, nodeDOFS); // to do: sjekk class om riktig

                // Calculate Jacobian matrix
                Matrix<double> jacobianMatrix = shapeFunctionsDerivatedNatural.Multiply(globalCoordinates);

                // Structure data in the form of a Jacobian matrix
                //Matrix<double> jacobianMatrix = calcDerivs.Transpose(); // to do: riktig ?? slik magnus hadde det

                // Calculate B - matrix
                Matrix<double> shapeFuncDerivatedCartesian = jacobianMatrix.Inverse().Multiply(shapeFunctionsDerivatedNatural); // to do: sjekk om riktig
                
                int dimRowB = 0;
                if (nodeDOFS == 2) { dimRowB = 3; }
                else { dimRowB = 6; }

                Matrix<double> B_i = DenseMatrix.Build.Dense( dimRowB , nodeDOFS*numElementNodes);

                for (int i = 0; i < numElementNodes; i++)
                { 
                    for (int j = 0; j < nodeDOFS; j++)
                    {

                        if (nodeDOFS == 2) // shell
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
                K_local.Add(B_i.Transpose().Multiply(C).Multiply(jacobianMatrix.Determinant()));
            }
            return Tuple.Create(K_local, B_local);
        }

        private Matrix<double> CalculateGlobalStiffnessMatrix(List<Element> elements, int numNode, int nodeDOFS, Material material)
        {
            // create stiffness matrix
            Matrix<double> K_global = Matrix<double>.Build.Dense(numNode * nodeDOFS, numNode * nodeDOFS);
            foreach (Element element in elements)
            {
                List<int> connectivity = element.Connectivity;
                Matrix<double> K_local = CalculateElementMatrices(element, material, nodeDOFS).Item1;

                // isoparametric mapping // to do: endre, for likt synne
                // loop nodes of elements
                for (int i = 0; i < connectivity.Count; i++)
                {
                    for (int j = 0; j < connectivity.Count; j++)
                    {
                        // loop relevant local stiffness contribution
                        for (int row = 0; row < nodeDOFS; row++)
                        {
                            for (int column = 0; column < nodeDOFS; column++) // to do: sjekk
                            {
                                K_global[nodeDOFS * connectivity[i] + row, nodeDOFS * connectivity[j] + column] = K_local[nodeDOFS * row + i, nodeDOFS * column + j];
                            }
                        }
                    }
                }
            }
            return K_global;
        }

        private Matrix<double> CalculateDisplacement(Matrix<double> K_global, Matrix<double> R, List<int> uIsZero)
        {
            foreach (int i in uIsZero)
            {
                foreach (int j in uIsZero)
                {
                    if (i != j)
                    {
                        K_global[i, j] = 0;
                    }
                    else 
                    {
                        K_global[i, j] = 1;
                        R[i, 1] = 0;
                    }
                }
            }
            Matrix<double> u = K_global.Inverse().Multiply(R);
            return u;  
        }

        private Tuple<Matrix<double>, Matrix<double>> CalculateElementStrainStress(Element element, Matrix<double> u, Material material, int nodeDOFS)
        {
            // summary: calculate a list of strain and stress vectors for each node in a element.
            Matrix<double> C = GetMaterialConstant(material.YoungModulus, material.PossionRatio, nodeDOFS);

            FEM _FEM = new FEM();
            List<Matrix<double>> B_local = CalculateElementMatrices(element, material, nodeDOFS).Item2;
            Matrix<double> elementGaussStrain = DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            Matrix<double> elementGaussStress = DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            Matrix<double> elementStrain = DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            Matrix<double> elementStress = DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);

            // get gauss strain and stress
            for (int n = 0; n < B_local.Count; n++)
            {
                // B-matrix is calculated from gauss points
                Matrix<double> gaussStrain = B_local[n].Multiply(u);
                Matrix<double> gaussStress = C.Multiply(B_local[n]).Multiply(u);

                for (int i = 0; i < B_local[0].RowCount; i++)
                {
                    elementGaussStrain[n, i] = gaussStrain[0,i];
                    elementGaussStress[n, i] = gaussStress[0, i];
                }
            }

            // get node strain and stress by interpolation
            Matrix<double> interpolationNodes = _FEM.GetGaussPoints(1 / Math.Sqrt(3), nodeDOFS); // marcin: ok constant?

            for (int n = 0; n < B_local.Count; n++)
            { 
                // get stress and strain in nodes
                var r = interpolationNodes.Row(n)[0];
                var s = interpolationNodes.Row(n)[1];
                double t = 0;
                if (nodeDOFS == 3) { t = interpolationNodes.Row(n)[2]; }

                Vector<double> shapefunctionValuesInNode = _FEM.GetShapeFunctions(r, s, t, nodeDOFS); // marcin: correct interpolation?
                Vector<double> nodeStrain = elementGaussStrain.Multiply(shapefunctionValuesInNode);
                Vector<double> nodeStress = elementGaussStress.Multiply(shapefunctionValuesInNode);
                for (int i = 0; i < B_local[0].RowCount; i++)
                {
                    elementStrain[n,i] = nodeStrain[i];
                    elementStress[n,i] = nodeStress[i];
                }
            }
            return Tuple.Create(elementStrain, elementStress);
        }

        private Tuple<Matrix<double>, Vector<double>> CalculateGlobalStress(List<Element> elements, Matrix<double> u, Material material, int nodeDOFS)
        {
            int numNodes = u.RowCount;
            int stressRowDim = 4;
            if (nodeDOFS == 3) { stressRowDim = 6; }
            Matrix<double> globalStress = DenseMatrix.Build.Dense(stressRowDim, numNodes);
            foreach (Element element in elements)
            {
                var elementStressStrain = CalculateElementStrainStress(element, u, material, nodeDOFS);
                Matrix<double> elementStrain = elementStressStrain.Item1;
                Matrix<double> elementStress = elementStressStrain.Item2;

                List<int> connectivity = element.Connectivity;

                for (int i = 0; i < elementStress.RowCount; i++) // loop the stress
                {
                    for (int j = 0; j < connectivity.Count; j++) // loop the element nodes
                    {
                        globalStress[i, connectivity[j]] = globalStress[i, connectivity[j]] + elementStress[i, j]; 
                    }
                    // to do: find the avarage if more contributions
                }
            }

            // Mises
            Vector<double> mises = DenseVector.Build.Dense(numNodes);
            for (int i = 0; i < numNodes; i++)
            {
                if (nodeDOFS == 2)
                {
                    Vector<double> nodeStress = globalStress.Column(i);
                    double Sxx = nodeStress[0];
                    double Syy = nodeStress[1];
                    double Sxy = nodeStress[2];
                    mises[i] = Math.Sqrt( Math.Pow(Sxx, 2) - Sxx * Syy + Math.Pow(Syy, 2) + 3 * Math.Pow(Sxy, 2));
                }
                else
                {
                    Vector<double> nodeStress = globalStress.Column(i);
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

        private Matrix<double> GetMaterialConstant(double youngModulus, double possionRatio, int nodeDOFS )
        {
            if (nodeDOFS == 2)
            {
                Matrix<double> C = DenseMatrix.OfArray(new double[,]
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
                Matrix<double> C = DenseMatrix.OfArray(new double[,]
                {
                    {1-possionRatio, possionRatio, possionRatio, 0, 0, 0},
                    {possionRatio, 1-possionRatio, possionRatio, 0, 0, 0},
                    {possionRatio, possionRatio, 1- possionRatio, 0, 0, 0},
                    {0, 0, 0, (1-2*possionRatio)/2, 0, 0},
                    {0, 0, 0, 0, (1-2*possionRatio)/2, 0},
                    {0, 0, 0, 0, 0, (1-2*possionRatio)/2},
                });
                C.Multiply((double)youngModulus / ((1 + possionRatio) * (1 - 2 * possionRatio)));
                return C;
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
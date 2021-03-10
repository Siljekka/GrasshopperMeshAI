using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints
{
    public class GenerateINPfile : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GenerateINPfile class.
        /// </summary>
        public GenerateINPfile()
          : base("Generate inp-file", "inp",
              "Generate inp-file for 3D-models. Default material is steel with E=210000, nu=0.3. Made for linear elastic analysis.",
              "MyPlugIn", "inp-file")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SolidMesh", "solid", "Solid mesh", GH_ParamAccess.item);
            pManager.AddTextParameter("ElementType", "element", "String with element type from Abaqus (IMPORTANT: must be written exactly as given in Abaqus). Default is: C3D8I", GH_ParamAccess.item, "C3D8I");
            pManager.AddNumberParameter("Young modulus", "E", "Value of Young modulus [MPa]. Default value is 210000 MPa", GH_ParamAccess.item, 210000);
            pManager.AddNumberParameter("Poisson Ratio", "nu", "Value of poisson ratio [-]. Default value is 0.3", GH_ParamAccess.item, 0.3);
            //pManager.AddGenericParameter("Material (string)", "material", "String with material from Abaqus", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("inp-file", "inp", "String containing inp-file", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Input
            Mesh3D solidMesh = new Mesh3D();
            string elementType = "C3D8I";
            double Emodul = 210000;
            double nu = 0.3;
            DA.GetData(0, ref solidMesh);
            DA.GetData(1, ref elementType);
            DA.GetData(2, ref Emodul);
            DA.GetData(3, ref nu);
            #endregion

            List<string> inpText = new List<string>();
            List<Node> nodes = solidMesh.Nodes;
            List<Element> elements = solidMesh.Elements;

            string partName = "Geometry"; //todo: fix name
            string sectionName = "Section"; //todo: fix name
            string materialName = "Steel";
            
            
            if (Emodul != 210000 | nu != 0.3) { materialName = "custom material"; } //todo: egendefinert på engelsk


            inpText.Add("*Heading");
            inpText.Add("**< text describing the problem being simulated.>"); //todo: fix description
            inpText.Add("**SI Units");
            inpText.Add("**x1=x, x2=y, x3=z");
            inpText.Add("*Preprint, echo = YES, model = YES, history = YES"); //recomended. Gives printout of the input file and of the model and history definition data
            
            // Start of part
            inpText.Add("**");
            inpText.Add("PARTS");
            inpText.Add("**");
            inpText.Add(String.Format("*Part, name={0}", partName)); //clean to have parts included, but not needed.
            inpText.Add("**");
             
            // Nodes
            inpText.Add("*Node");
            foreach (Node node in nodes)
            {
                int globalId = node.GlobalId;
                double nodeX = node.Coordinate.X;
                double nodeY = node.Coordinate.Y;
                double nodeZ = node.Coordinate.Z;
                inpText.Add(String.Format("{0}, {1}, {2}, {3}", globalId, nodeX, nodeY, nodeZ)); //GlobalId, x-coord, y-coord, z-coord
            }
            inpText.Add("**");

            // Elements
            inpText.Add(String.Format("*Element, type={0}", elementType));
            foreach (Element e in elements)
            {
                int elementId = e.Id;
                int n1 = e.Node1.GlobalId;
                int n2 = e.Node2.GlobalId;
                int n3 = e.Node3.GlobalId;
                int n4 = e.Node4.GlobalId;
                int n5 = e.Node5.GlobalId;
                int n6 = e.Node6.GlobalId;
                int n7 = e.Node7.GlobalId;
                int n8 = e.Node8.GlobalId;
                inpText.Add(String.Format("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}", elementId, n1, n2, n3, n4, n5, n6, n7, n8)); //ElementId, n1, n2, n3, n4, n5, n6, n7, n8
            }
            inpText.Add("**");

            // optional
            // *Nset, nset =/*OPTIONAL: name*/, generate //Optional: assigns nodes to a node set.
            // <first - node, last - node, Increment.>
            // *Elset, elset = Fixed, generate
            // <first - element, last - element, Increment.>
            // have as many nset and elset one want

            inpText.Add("*Nset, nset=Whole, generate");
            inpText.Add(String.Format("1, {0}, 1 >", nodes.Count));
            inpText.Add("*Elset, elset=Whole, generate");
            inpText.Add(String.Format("1, {0}, 1 >", elements.Count));
            inpText.Add("**");

            // Section
            inpText.Add(String.Format("**Section: {0}", sectionName));
            inpText.Add(String.Format("*Solid Section, elset=Whole, material={0}", materialName));//todo: fix materialinput;
            inpText.Add(",");
            inpText.Add("*End Part");

            // Assembly
            inpText.Add("**");
            inpText.Add("**");
            inpText.Add("ASSEMBLY");
            inpText.Add("**");
            inpText.Add("*Assembly, name=Assembly");
            inpText.Add("**");
            inpText.Add(String.Format("*Instance, name=instance, part={0}", partName));
            inpText.Add("*End Instance");
            inpText.Add("**");
            inpText.Add("*End Assembly");
            //todo: check extra elset, nset

            // Material
            inpText.Add("**");
            inpText.Add("** MATERIALS");
            inpText.Add("**");
            inpText.Add(String.Format("*Material, name={0}", materialName));
            inpText.Add("*Elastic");
            inpText.Add(String.Format("{0}, {1}", Emodul, nu));

            // Boundary Condition
            inpText.Add("**");
            inpText.Add("**BOUNDARY CONDITIONS");
            inpText.Add("**");
            inpText.Add("** Name: BCType Type: Displacement / Rotation"); //fix name
            inpText.Add("*Boundary");
            inpText.Add("Input: Fill out relevant BC"); //todo: fill in BC

            // STEP
            inpText.Add("**----------------------------------------------------------------");
            inpText.Add("**");
            inpText.Add("STEP: Load");
            inpText.Add("**");
            inpText.Add("**Step, name=Load, nlgeom=NO");
            inpText.Add("*Static");
            inpText.Add("0.1, 1., 1e-05, 1.");

            // Load
            inpText.Add("**LOADS");
            inpText.Add("**");
            inpText.Add("**Name: LoadType   Type: <Input: Type of load>"); //todo: fix 
            inpText.Add("<Input: Fill out relevant loads>");
            inpText.Add("<Input: Fill out relevant loads>");
            inpText.Add("**");

            // OUTPUT
            //todo: fix output
            inpText.Add("**OUTPUT REQUESTS");
            inpText.Add("**");
            inpText.Add("*Restart, write, frequency = 0");
            inpText.Add("**");
            inpText.Add("**FIELD OUTPUT: F - Output - 1");
            inpText.Add("**");
            inpText.Add("*Output, field");
            inpText.Add("*Node Output");
            inpText.Add("CF, RF, U");
            inpText.Add("*Element Output, directions = YES");
            inpText.Add("E, MISES, S");
            inpText.Add("**");
            inpText.Add("**HISTORY OUTPUT: H - Output - 1");
            inpText.Add("**");
            inpText.Add("*Output, history");
            inpText.Add("*Energy Output, elset = BeamInstance.Fixed, variable = PRESELECT");
            inpText.Add("*End Step");


            DA.GetDataList(0, inpText);
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
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("1e9c1cd4-77ed-49a2-bd0c-34fb0871f992"); }
        }
    }
}
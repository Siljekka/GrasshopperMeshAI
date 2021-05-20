using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;
using System.IO;

namespace MeshPoints
{
    public class GenerateINPfile : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GenerateINPfile class.
        /// </summary>
        public GenerateINPfile()
          : base("Generate inp-file", "inp",
              "Generate inp-file for 3D analysis. Solid elements (C3D8) and shell elements (S4) are made for SmartMesh dependent on mesh type. Default material is steel with E=210000, nu=0.3. Made for linear elastic analysis.",
              "SmartMesh", "Tools")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "mesh", "SmartMesh Class. Geometry must have been modelled in mm.", GH_ParamAccess.item);
            pManager.AddGenericParameter("FilePath", "path", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("Save", "save", "", GH_ParamAccess.item);
            pManager.AddNumberParameter("Shell thickness", "t", "Value of shell thickness [mm]. Only for SurfaceMesh. Default value is 10 mm", GH_ParamAccess.item);
            pManager.AddNumberParameter("Young modulus", "E", "Value of Young modulus [MPa]. Default value is 210000 MPa", GH_ParamAccess.item);
            pManager.AddNumberParameter("Poisson Ratio", "nu", "Value of poisson ratio [-]. Default value is 0.3", GH_ParamAccess.item);
           

            pManager[1].Optional = true; // Write to file is optional
            pManager[2].Optional = true; // Write to file is optional
            pManager[3].Optional = true; // Shell thickness is optional
            pManager[4].Optional = true; // Young modulus is optional
            pManager[5].Optional = true; // Poisson Ratio is optional
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("inp", "inp", "String containing inp-file", GH_ParamAccess.list); //todo: change name.
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            SmartMesh smartMesh = new SmartMesh();
            string filePath = "Empty";
            bool writeFile = false;
            double sectionThickness = 10;
            double Emodul = 210000;
            double nu = 0.3;
            
            DA.GetData(0, ref smartMesh);
            DA.GetData(1, ref filePath);
            DA.GetData(2, ref writeFile);
            DA.GetData(3, ref sectionThickness);
            DA.GetData(4, ref Emodul);
            DA.GetData(5, ref nu);
            

            // Variables
            List<string> inpText = new List<string>();
            string elementType = "Empty";
            string partName = "Geometry"; //todo: fix name
            string sectionName = "Section"; //todo: fix name
            string materialName = "Steel";


            // 0. Check input.
            if (!DA.GetData(0, ref smartMesh)) return;
            

            // 1. Set material properties.
            if (Emodul != 210000 | nu != 0.3) { materialName = "custom material"; } //todo: egendefinert på engelsk

            // 2. Set element type dependent on solid or surface mesh.
            if (smartMesh.Type == "Solid") { elementType = "C3D8"; }
            else if (smartMesh.Type == "Surface") { elementType = "S4"; }
            else { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No element type has been assigned."); }
            
            // 3. Generate inp-file.
            if (smartMesh.Type == "Solid")
            {
                inpText = GenerateSolidfile(smartMesh, elementType, Emodul, nu, partName, sectionName, materialName);
            }
            else if (smartMesh.Type == "Surface")
            {
                inpText = GenerateSurfacefile(smartMesh, elementType, sectionThickness, Emodul, nu, partName, sectionName, materialName);
            }

            // 4. Write to .txt-file
            if (DA.GetData(1, ref filePath) & writeFile)
            {
                var file = @filePath;
                File.WriteAllLines(file, inpText.ToArray());
            }
            
            // Output
            DA.SetDataList(0, inpText);
        }

        #region Methods

        /// <summary> Generate inp file for solid models. </summary>
        /// <param name="solidMesh"> SmartMesh Class. </param>
        /// <param name="elementType"> Element type used for the solid elements. </param>
        /// <param name="Emodul"> Modulus of Elasticity. </param>
        /// <param name="nu"> Poisson Ratio. </param>
        /// <param name="partName"> Name of the Part. </param>
        /// <param name="sectionName"> Name of the Section. </param>
        /// <param name="materialName"> Name of the Material. </param>
        /// <returns> .txt-file that can be converted to an .inp-file </returns>
        private List<string> GenerateSolidfile(SmartMesh solidMesh, string elementType, double Emodul, double nu, string partName, string sectionName, string materialName)
        {
            //Transform.ChangeBasis(Vector3d.XAxis, Vector3d.YAxis, Vector3d.ZAxis, Vector3d.ZAxis, Vector3d.XAxis, Vector3d.YAxis);
            List<string> inpText = new List<string>();
            List<Node> nodes = solidMesh.Nodes;
            List<Element> elements = solidMesh.Elements;

            inpText.Add("*Heading");
            inpText.Add("**<Input: text describing the problem being simulated.>"); //todo: fix description
            inpText.Add("**SI Units");
            inpText.Add("**x1=x, x2=y, x3=z");
            inpText.Add("*Preprint, echo=YES, model=YES, history=YES"); //recomended. Gives printout of the input file and of the model and history definition data

            // Start of part
            inpText.Add("**");
            inpText.Add("** PARTS");
            inpText.Add("**");
            inpText.Add(String.Format("*Part, name={0}", partName)); //clean to have parts included, but not needed.
            
            // Nodes
            inpText.Add("*Node");
            foreach (Node node in nodes)
            {
                int globalId = node.GlobalId+1;
                double nodeX = node.Coordinate.X;
                double nodeY = node.Coordinate.Y;
                double nodeZ = node.Coordinate.Z;
                inpText.Add(String.Format("{0}, {1}, {2}, {3}", globalId, nodeX, nodeZ, nodeY)); //GlobalId, x-coord, y-coord, z-coord
            }

            // Elements
            inpText.Add(String.Format("*Element, type={0}", elementType));
            foreach (Element e in elements)
            {
                int elementId = e.Id+1;
                int n1 = e.Nodes[0].GlobalId + 1;
                int n2 = e.Nodes[1].GlobalId + 1;
                int n3 = e.Nodes[2].GlobalId + 1;
                int n4 = e.Nodes[3].GlobalId + 1;
                int n5 = e.Nodes[4].GlobalId + 1;
                int n6 = e.Nodes[5].GlobalId + 1;
                int n7 = e.Nodes[6].GlobalId + 1;
                int n8 = e.Nodes[7].GlobalId + 1;
                inpText.Add(String.Format("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}", elementId, n4, n3, n2, n1, n8, n7, n6, n5)); //ElementId, nodes. Count cw because of transformation of coordsyst. between rhino and abaqus.
            }
            inpText.Add("**");

            // optional
            // *Nset, nset =/*OPTIONAL: name*/, generate //Optional: assigns nodes to a node set.
            // <first - node, last - node, Increment.>
            // *Elset, elset = Fixed, generate
            // <first - element, last - element, Increment.>
            // have as many nset and elset one want

            inpText.Add("*Nset, nset=Whole, generate");
            inpText.Add(String.Format("1, {0}, 1", nodes.Count));
            inpText.Add("*Elset, elset=Whole, generate");
            inpText.Add(String.Format("1, {0}, 1", elements.Count));
            inpText.Add("**");

            // Section
            inpText.Add(String.Format("** Section: {0}", sectionName));
            inpText.Add(String.Format("*Solid Section, elset=Whole, material={0}", materialName));//todo: fix materialinput;
            inpText.Add(",");
            inpText.Add("*End Part");

            // Assembly
            inpText.Add("**");
            inpText.Add("**");
            inpText.Add("** ASSEMBLY");
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
            inpText.Add("** BOUNDARY CONDITIONS");
            inpText.Add("**");
            inpText.Add("** Name: BCType1 Type: Displacement/Rotation"); //fix name
            inpText.Add("*Boundary");
            inpText.Add("<Input: Fill out relevant BC>"); //todo: fill in BC
            inpText.Add("** Name: BCType2 Type: Displacement/Rotation"); //fix name
            inpText.Add("**Boundary");
            inpText.Add("**<Input: Fill out relevant BC>"); //todo: fill in BC

            // STEP
            inpText.Add("**----------------------------------------------------------------");
            inpText.Add("**");
            inpText.Add("** STEP: Load");
            inpText.Add("**");
            inpText.Add("*Step, name=Load, nlgeom=NO");
            inpText.Add("<Input: description of load applid in step>");
            inpText.Add("*Static");
            inpText.Add("0.1, 1., 1e-05, 1.");

            // Load
            inpText.Add("** LOADS");
            inpText.Add("**");
            inpText.Add("** Name: LoadType   Type: <Input: Type of load>"); //todo: fix 
            inpText.Add("<Input: Fill out relevant loads>");
            inpText.Add("<Input: Fill out relevant loads>");
            inpText.Add("**");

            // OUTPUT
            //todo: fix output
            inpText.Add("**OUTPUT REQUESTS");
            inpText.Add("**");
            inpText.Add("*Restart, write, frequency=0"); //fix for solid
            inpText.Add("**");

            inpText.Add("**FIELD OUTPUT: F-Output-1");
            inpText.Add("**");
            inpText.Add("*Output, field, variable=PRESELECT");
            inpText.Add("**");

            inpText.Add("**HISTORY OUTPUT: H-Output-1");
            inpText.Add("**");
            inpText.Add("*Output, history, variable=PRESELECT");

            inpText.Add("*End Step");

            return inpText;
        }

        /// <summary> Generate inp file for surface models. </summary>
        /// <param name="solidMesh"> SmartMesh Class. </param>
        /// <param name="elementType"> Element type used for the shell elements. </param>
        /// <param name="sectionThickness"> Thickness of the shell element. </param>
        /// <param name="Emodul"> Modulus of Elasticity. </param>
        /// <param name="nu"> Poisson Ratio. </param>
        /// <param name="partName"> Name of the Part. </param>
        /// <param name="sectionName"> Name of the Section. </param>
        /// <param name="materialName"> Name of the Material. </param>
        /// <returns> .txt-file that can be converted to an .inp-file </returns>
        private List<string> GenerateSurfacefile(SmartMesh surfaceMesh, string elementType, double sectionThickness, double Emodul, double nu, string partName, string sectionName, string materialName)
        {
            List<string> inpText = new List<string>();
            List<Node> nodes = surfaceMesh.Nodes;
            List<Element> elements = surfaceMesh.Elements;
            inpText.Add("*Heading");
            inpText.Add("**< text describing the problem being simulated.>"); //todo: fix description
            inpText.Add("**SI Units");
            inpText.Add("**x1=x, x2=y, x3=z");
            inpText.Add("*Preprint, echo=YES, model=YES, history=YES"); //recomended. Gives printout of the input file and of the model and history definition data

            // Start of part
            inpText.Add("**");
            inpText.Add("** PARTS");
            inpText.Add("**");
            inpText.Add(String.Format("*Part, name={0}", partName)); //clean to have parts included, but not needed.

            // Nodes
            inpText.Add("*Node");
            foreach (Node node in nodes)
            {
                int globalId = node.GlobalId + 1;
                double nodeX = node.Coordinate.X;
                double nodeY = node.Coordinate.Y;
                double nodeZ = node.Coordinate.Z;
                inpText.Add(String.Format("{0}, {1}, {2}, {3}", globalId, nodeX, nodeY, nodeZ)); //GlobalId, x-coord, y-coord, z-coord
            }

            // Elements
            inpText.Add(String.Format("*Element, type={0}", elementType));
            foreach (Element e in elements)
            {
                int elementId = e.Id + 1;
                int n1 = e.Nodes[0].GlobalId + 1;
                int n2 = e.Nodes[1].GlobalId + 1;
                int n3 = e.Nodes[2].GlobalId + 1;
                int n4 = e.Nodes[3].GlobalId + 1;
                inpText.Add(String.Format("{0}, {1}, {2}, {3}, {4}", elementId, n1, n2, n3, n4)); //ElementId, n1, n2, n3, n4
            }
            inpText.Add("**");

            // optional
            // *Nset, nset =/*OPTIONAL: name*/, generate //Optional: assigns nodes to a node set.
            // <first - node, last - node, Increment.>
            // *Elset, elset = Fixed, generate
            // <first - element, last - element, Increment.>
            // have as many nset and elset one want

            inpText.Add("*Nset, nset=Whole, generate");
            inpText.Add(String.Format("1, {0}, 1", nodes.Count));
            inpText.Add("*Elset, elset=Whole, generate");
            inpText.Add(String.Format("1, {0}, 1", elements.Count));
            inpText.Add("**");

            // Section
            inpText.Add(String.Format("**Section: {0}", sectionName));
            inpText.Add(String.Format("*Shell General Section, elset=Whole, material={0}", materialName));//todo: fix materialinput
            inpText.Add(String.Format("{0},", sectionThickness));
            inpText.Add("*End Part");

            // Assembly
            inpText.Add("**");
            inpText.Add("**");
            inpText.Add("** ASSEMBLY");
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
            inpText.Add("** Name: BCType1 Type: Displacement/Rotation"); //fix name
            inpText.Add("*Boundary");
            inpText.Add("<Input: Fill out relevant BC>");
            inpText.Add("** Name: BCType2 Type: Displacement/Rotation"); //fix name
            inpText.Add("**Boundary");
            inpText.Add("**<Input: Fill out relevant BC>"); 

            // STEP
            inpText.Add("**----------------------------------------------------------------");
            inpText.Add("**");
            inpText.Add("** STEP: Load");
            inpText.Add("**");
            inpText.Add("*Step, name=Load, nlgeom=NO");
            inpText.Add("<Input: description of load applid in step>"); 
            inpText.Add("*Static");
            inpText.Add("0.1, 1., 1e-05, 1.");

            // Load
            inpText.Add("**");
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
            inpText.Add("*Restart, write, frequency=0"); //fix for solid
            inpText.Add("**");

            inpText.Add("**FIELD OUTPUT: F-Output-1");
            inpText.Add("**");
            inpText.Add("*Output, field, variable=PRESELECT");
            inpText.Add("**");

            inpText.Add("**HISTORY OUTPUT: H-Output-1");
            inpText.Add("**");
            inpText.Add("*Output, history, variable=PRESELECT");
            
            inpText.Add("*End Step");

            return inpText;
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
                return Properties.Resources.Icon_INP;
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
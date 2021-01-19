using Grasshopper.Kernel;
using System;
using System.Drawing;

namespace MeshPoints
{
    public class MeshPointsInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "MeshPoints";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("bfc8b8d2-140f-4b0d-9a57-af823cfa2185");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "HP Inc.";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "";
            }
        }
    }
}

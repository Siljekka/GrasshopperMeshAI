using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using Grasshopper.Kernel;


namespace MeshPoints.Classes
{
    class Quality
    {
        //Properties
        public double AspectRatio { get; set; }
        public double Skewness { get; set; }
        public Element element { get; set; }


        //Constructor
        public Quality()
            {
                //Empty constructor
            }

        public Quality(Element _elem, double _aspectRatio, double _skewness)
        {
            AspectRatio = _aspectRatio;
            Skewness = _skewness;
            element = _elem;
        }

        //Methods


    }
}

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
        
        public Element Elem { get; set; }


        //Constructor
        public Quality()
            {
                //Empty constructor
            }

        public Quality(double _aspectRatio, double _skewness, Element _elem)
        {
            AspectRatio = _aspectRatio;
            Skewness = _skewness;
            Elem = _elem;
        }

        //Methods


    }
}

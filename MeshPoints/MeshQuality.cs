using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;


namespace MeshPoints
{
    class MeshQuality
    {
        //Properties
        public List<double> aspectRatio { get; set; }
        public List<double> angleRatio { get; set; }
        public MeshFace meshFace { get; set; }

        //Constructor
        public MeshQuality()
        {
            //Empty constructor
        }
        public MeshQuality(List<double> _aspectRatio, List<double> _angleRatio)
        {
            aspectRatio = _aspectRatio;
            angleRatio = _angleRatio;
        }

        //Methods
    }
}

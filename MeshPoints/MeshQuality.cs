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
        public List<double> AspectRatio { get; set; }
        public List<double> AngleRatio { get; set; }
        public MeshFace MeshFace { get; set; }

        //Constructor
        public MeshQuality()
        {
            //Empty constructor
        }
        public MeshQuality(List<double> _aspectRatio, List<double> _angleRatio)
        {
            AspectRatio = _aspectRatio;
            AngleRatio = _angleRatio;
        }

        public MeshQuality(List<double> _aspectRatio)
        {
            AspectRatio = _aspectRatio;
        }

        //Methods
    }
}

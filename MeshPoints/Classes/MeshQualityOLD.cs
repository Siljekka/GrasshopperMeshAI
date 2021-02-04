using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

//_______OLD VERSION_________

namespace MeshPoints.Classes
{
    class MeshQualityOLD
    {
        //Properties
        public double AspectRatio { get; set; }
        public double Skewness { get; set; }
        public MeshFace MeshFace { get; set; }

        //add vertices

        //Constructor
        public MeshQualityOLD()
        {
            //Empty constructor
        }
        public MeshQualityOLD(MeshFace _meshFace, double _aspectRatio, double _skewness)
        {
            MeshFace = _meshFace;
            AspectRatio = _aspectRatio;
            Skewness = _skewness;
        }
        //Methods
    }
}

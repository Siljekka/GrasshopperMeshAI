using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeshPoints.Classes
{
    class qConstants
    {
        public qConstants()
        {
            //Empty constructor
        }

        public double GetThetaTolerance()
        {
            return 0.16667 * Math.PI;
        }
        public double GetThetaToleranceForClosing()
        {
            return 1.5 * 0.16667 * Math.PI; // 1.5 * ThetaTolerance
        }
    }
}

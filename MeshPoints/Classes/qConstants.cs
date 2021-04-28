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
            return Math.PI / (double)6;
        }
        public double GetThetaToleranceForClosing()
        {
            return 1.5 * Math.PI / (double)6; // 1.5 * ThetaTolerance
        }

        public double GetTransitionTolerance()
        {
            return 2.5;
        }

    }
}

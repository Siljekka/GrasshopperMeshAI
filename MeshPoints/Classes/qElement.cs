using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    class qElement
    {
        public List<qEdge> EdgeList { get; set; }
        public List<double> AngleList { get; } // todo: when angle is larger than pi it does not work..
        public bool IsQuad { get; }


        // Constructer

        public qElement()
        {
            // empty constructor
        }

        public qElement(List<qEdge> _edgeList)
        {
            EdgeList = _edgeList;
            AngleList = CalculateAngles(_edgeList);

            if (_edgeList.Count == 4) { IsQuad = true; }
            else { IsQuad = false; }
        }

        // Methods
        
        private List<double> CalculateAngles(List<qEdge> _edgeList)
        {
            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;
            double ang = 0;
            List<double> angList = new List<double>();
            List<qEdge> edgeListCopy = new List<qEdge>(_edgeList);
            edgeListCopy.Add(edgeListCopy[0]);

            // todo: change calculation of angles AnalyzeTriangle (do not work for chevron)
            for (int i = 0; i < edgeListCopy.Count-1; i++)
            {
                Point3d start1 = edgeListCopy[i].StartNode.Coordinate;
                Point3d end1 = edgeListCopy[i].EndNode.Coordinate;
                Point3d start2 = edgeListCopy[i + 1].StartNode.Coordinate;
                Point3d end2 = edgeListCopy[i + 1].EndNode.Coordinate;

                vec1 = start1 - end1;
                vec2 = end2 - start2;
                
                ang = Vector3d.VectorAngle(vec1, vec2); // radian
                angList.Add(ang);
            }
            return angList;
        
        }
        
    }
}

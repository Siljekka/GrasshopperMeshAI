using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;


namespace MeshPoints.Classes
{
    class qEdge
    {
        public qNode StartNode { get; set; }
        public qNode EndNode { get; set; }
        public double Length { get; set; }
        public Line EdgeLine { get; set; } // for visualization
        public qElement Element1 { get; set; }
        public qElement Element2 { get; set; }
        public qEdge LeftFrontNeighbor { get; set; }
        public qEdge RightFrontNeighbor { get; set; }
        public int Level { get; set; } // level if front edge
        public bool IsQuadSideEdge { get; set; }




        public qEdge()
        {
            // empty constructor
        }

        public qEdge(qNode _startNode, qNode _endNode)
        {
            StartNode = _startNode;
            EndNode = _endNode;
            Length = CalculateLength(_startNode, _endNode);
            EdgeLine = VisualizeLine(_startNode, _endNode);
            Level = 0; // check if needed
            IsQuadSideEdge = false;
        }

        public double CalculateLength(qNode _startNode, qNode _endNode)
        {
            return _startNode.Coordinate.DistanceTo(_endNode.Coordinate);
        }
        public Line VisualizeLine(qNode _startNode, qNode _endNode)
        {
            return new Line(_startNode.Coordinate, _endNode.Coordinate);
        }


    }
}

﻿using System;
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
        public int Index { get; set; }
        public qNode StartNode { get; set; }
        public qNode EndNode { get; set; }
        public double Length { get; set; }
        public Line EdgeLine { get; set; } // makes the user see the edge, might delete?
        public qElement Element1 { get; set; }
        public qElement Element2 { get; set; }
        public qEdge LeftFrontNeighbor { get; set; }
        public qEdge RightFrontNeighbor { get; set; }

        

    public qEdge()
        {
            // empty constructor
        }

        public qEdge(int _index, qNode _startNode, qNode _endNode)
        {
            Index = _index;
            StartNode = _startNode;
            EndNode = _endNode;
            Length = CalculateLength(_startNode, _endNode);
            EdgeLine = VisualizeLine(_startNode, _endNode);
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

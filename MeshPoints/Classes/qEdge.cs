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

        public qNode GetOppositeNode(qNode node)
        {
            qNode oppositeNode = new qNode();
            if (node == this.StartNode) { oppositeNode = this.EndNode; }
            else if (node == this.EndNode) { oppositeNode = this.StartNode; }
            return oppositeNode;
        }
        public bool IsFrontEdge()
        {
            // summary: check if an edge is a front edge
            bool check = false;
            List<qElement> connectedElements = this.GetConnectedElements();
            if (connectedElements.Count == 1)
            {
                if (!this.Element1.IsQuad)
                {
                    check = true;
                }
            }
            else if (connectedElements.Count != 1)
            {
                if (!this.Element1.IsQuad & this.Element2.IsQuad) { check = true; }
                else if (this.Element1.IsQuad & !this.Element2.IsQuad) { check = true; }
            }

            return check;
        }

        public List<qElement> GetConnectedElements()
        {
            // summary: get conneccted elements to an edge. Assume edge has updated elements element 1 and/or element 2.
            List<qElement> connectedElements = new List<qElement>();
            connectedElements.Add(this.Element1);
            if (this.Element2 != null)
            {
                connectedElements.Add(this.Element2);
            }
            return connectedElements;
        } 



    }
}

﻿using System;
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
        public List<double> AngleList { get; set; } // todo: when angle is larger than pi it does not work..
        public List<Line> Contour { get; set; }
        public bool IsQuad { get; }


        // Constructors:
        public qElement()
        {
            // empty constructor
        }
        public qElement(List<qEdge> _edgeList)
        {
            EdgeList = _edgeList;
            //FixEdgeOrder();
            AngleList = CalculateAngles(_edgeList);
            Contour = GetContourOfElement(_edgeList);

            if (_edgeList.Count == 4) { IsQuad = true; }
            else { IsQuad = false; }
        }

        // Methods
        public List<double> CalculateAngles(List<qEdge> _edgeList)
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
        public List<Line> GetContourOfElement(List<qEdge> _edgeList)
        {
            List<Line> contour = new List<Line>();
            foreach (qEdge edge in _edgeList)
            {
                Line contourLine = new Line(edge.StartNode.Coordinate, edge.EndNode.Coordinate);
                contour.Add(contourLine);
            }
            return contour;
        }
        public Point3d GetElementCenter()
        {
            // summary: get center of an element
            double sx = 0;
            double sy = 0;
            double sz = 0;

            List<qEdge> edgeList = this.EdgeList;
            foreach (qEdge edge in edgeList)
            {
                Point3d startPoint = edge.StartNode.Coordinate;
                Point3d endPoint = edge.EndNode.Coordinate;
                List<Point3d> pts = new List<Point3d>() { startPoint, endPoint };
                foreach (Point3d pt in pts)
                {
                    sx = sx + pt.X;
                    sy = sy + pt.Y;
                    sz = sz + pt.Z;
                }
            }
            int n = edgeList.Count * 2;
            Point3d centerPt = new Point3d(sx / n, sy / n, sz / n);
            return centerPt;
        } 
        public List<qNode> GetNodesOfElement()
        {
            // summary: get nodes of an element
            List<qNode> nodeList = new List<qNode>();
            if (!this.IsQuad)
            {
                foreach (qEdge edge in this.EdgeList)
                {
                    if (!nodeList.Contains(edge.StartNode))
                    {
                        nodeList.Add(edge.StartNode);
                    }

                    if (!nodeList.Contains(edge.EndNode))
                    {
                        nodeList.Add(edge.EndNode);
                    }
                }
            }
            else if (this.IsQuad)
            {
                List<qEdge> quadEdges = this.EdgeList;

                qEdge baseEdge = quadEdges[0];
                qEdge rightEdge = quadEdges[1];
                qEdge leftEdge = quadEdges[2];
                qEdge topEdge = quadEdges[3];

                qNode node1 = new qNode();
                qNode node2 = new qNode();
                qNode node3 = new qNode();
                qNode node4 = new qNode();

                if (baseEdge.StartNode == leftEdge.StartNode | baseEdge.StartNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.StartNode;
                    node2 = baseEdge.EndNode;
                }
                else if (baseEdge.EndNode == leftEdge.StartNode | baseEdge.EndNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.EndNode;
                    node2 = baseEdge.StartNode;
                }
                else if (baseEdge.StartNode == rightEdge.StartNode | baseEdge.StartNode == rightEdge.EndNode)
                {
                    node2 = baseEdge.StartNode;
                    node1 = baseEdge.EndNode;
                }
                else if (baseEdge.EndNode == rightEdge.StartNode | baseEdge.EndNode == rightEdge.EndNode)
                {
                    node2 = baseEdge.EndNode;
                    node1 = baseEdge.StartNode;
                }

                if (topEdge.StartNode == leftEdge.StartNode | topEdge.StartNode == leftEdge.EndNode)
                {
                    node3 = topEdge.EndNode;
                    node4 = topEdge.StartNode;
                }
                else if (topEdge.EndNode == leftEdge.StartNode | topEdge.EndNode == leftEdge.EndNode)
                {
                    node3 = topEdge.StartNode;
                    node4 = topEdge.EndNode;
                }
                else if (topEdge.StartNode == rightEdge.StartNode | topEdge.StartNode == rightEdge.EndNode)
                {
                    node4 = topEdge.EndNode;
                    node3 = topEdge.StartNode;
                }
                else if (topEdge.EndNode == rightEdge.StartNode | topEdge.EndNode == rightEdge.EndNode)
                {
                    node4 = topEdge.StartNode;
                    node3 = topEdge.EndNode;
                }

                nodeList = new List<qNode> { node1, node2, node3, node4 }; // n1: bottom left, n2: bottom right, n3: top right, n4: top left
            }
            return nodeList;
        }

        public void FixEdgeOrder()
        {
            // summary: fix edge order of triangle elements

            qEdge edge = this.EdgeList[0];
            qEdge edgeConnectedToStartNode = new qEdge();
            qEdge edgeConnectedToEndNode = new qEdge();
            qEdge edgeNotConnected = new qEdge();

            for (int i = 1; i < this.EdgeList.Count; i ++)
            {
                if (edge.StartNode == this.EdgeList[i].StartNode | edge.StartNode == this.EdgeList[i].EndNode)
                {
                    edgeConnectedToStartNode = this.EdgeList[i];
                }
                else if (edge.EndNode == this.EdgeList[i].StartNode | edge.EndNode == this.EdgeList[i].EndNode)
                {
                    edgeConnectedToEndNode = this.EdgeList[i];
                }
                else 
                {
                    edgeNotConnected = this.EdgeList[i]; 
                }
            }

            Point3d midPointEdg = 0.5 * (edge.StartNode.Coordinate + edge.EndNode.Coordinate); // mid point of edge
            Point3d centerPoint = GetElementCenter();
            Vector3d centerToMidVector = midPointEdg - centerPoint;
            Vector3d centerToEndNodeVector = edge.EndNode.Coordinate - centerPoint;
            Vector3d centerToStartNodeVector = edge.StartNode.Coordinate - centerPoint;
            double startAngle = Vector3d.VectorAngle(centerToMidVector, centerToStartNodeVector, Vector3d.ZAxis); // todo: make normal more general
            double endAngle = Vector3d.VectorAngle(centerToMidVector, centerToEndNodeVector, Vector3d.ZAxis); // todo: make normal more general

            if (endAngle < startAngle)
            {
                if (!this.IsQuad)
                {
                    this.EdgeList = new List<qEdge>() { edge, edgeConnectedToEndNode, edgeConnectedToStartNode };
                }
                else 
                {
                    this.EdgeList = new List<qEdge>() { edge, edgeConnectedToEndNode, edgeConnectedToStartNode, edgeNotConnected };
                }
            }
            else
            {
                if (!this.IsQuad)
                {
                    this.EdgeList = new List<qEdge>() { edge, edgeConnectedToStartNode, edgeConnectedToEndNode };
                }
                else
                {
                    this.EdgeList = new List<qEdge>() { edge, edgeConnectedToStartNode, edgeConnectedToEndNode, edgeNotConnected };
                }
            }
        }

        public bool IsChevron()
        {
            if (this.AngleList.Max() > (double)200 / (double)180 * Math.PI) { return true; }
            else { return false; }
        }
        public bool IsTriangleInverted()
        {
            // summary: check if a triangle element is inverted
            bool isInverted = false;

            Point3d A = this.EdgeList[0].GetSharedNode(this.EdgeList[1]).Coordinate;
            Point3d B = this.EdgeList[1].GetSharedNode(this.EdgeList[2]).Coordinate;
            Point3d C = this.EdgeList[2].GetSharedNode(this.EdgeList[0]).Coordinate;

            // check area
            double area = 0.5 * (A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y));
            if (area <= 0) { isInverted = true; }
            return isInverted;
        }
    }
}

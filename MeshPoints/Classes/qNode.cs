﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace MeshPoints.Classes
{
    class qNode
    {
        public Point3d Coordinate { get; set; }
        public bool BoundaryNode { get; set; }

        // Constructors:
        public qNode()
        {
            // empty constructor
        }

        public qNode(Point3d _coordinate, bool _boundaryNode)
        {
            Coordinate = _coordinate;
            BoundaryNode = _boundaryNode;
        }

        // Methods:
        public List<qEdge> GetConnectedEdges(List<qEdge> globalEdgeList)
        {
            qNode node = this;
            // summary: get connected edges to a node
            List<qEdge> connectedEdges = new List<qEdge>();
            foreach (qEdge edge in globalEdgeList)
            {
                if ((edge.StartNode == node) | (edge.EndNode == node))
                {
                    connectedEdges.Add(edge);
                }
            }
            //if (connectedEdges.Count == 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "GetConnectedEdges: connectedEdges is empty"); }
            return connectedEdges;
        }
        public bool IsFrontNode(List<qEdge> frontEdges)
        {
            qNode node = this;
            // summary: check if node is a front node
            bool isFrontNode = false;
            foreach (qEdge front in frontEdges)
            {
                if (node == front.StartNode | node == front.EndNode)
                {
                    isFrontNode = true;
                    break;
                }
            }
            return isFrontNode;
        } // class: kan bli implementert i: qNode
        public List<qElement> GetQuadsConnectedToNode(List<qEdge> globalEdgeList)
        {
            qNode node = this;
            // summary: get all quad elements connected to a node
            List<qElement> quadElements = new List<qElement>();
            var connectedEdges = node.GetConnectedEdges(globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                List<qElement> connectedElements = edge.GetConnectedElements();
                foreach (qElement element in connectedElements)
                {
                    if (element.IsQuad) { quadElements.Add(element); }
                }
            }

            // delete dublicates
            List<qElement> quadElementsNoDublicates = new List<qElement>();
            if (quadElements.Count > 0)
            {
                foreach (qElement element in quadElements)
                {
                    if (!quadElementsNoDublicates.Contains(element)) { quadElementsNoDublicates.Add(element); }
                }
            }

            return quadElementsNoDublicates;
        } // class: kan bli implementert i: qNode
        public List<qElement> GetTrianglesConnectedToNode(List<qEdge> globalEdgeList)
        {
            qNode node = this;
            // summary: get all quad elements connected to a node
            List<qElement> triangleElements = new List<qElement>();
            var connectedEdges = node.GetConnectedEdges(globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                List<qElement> connectedElements = edge.GetConnectedElements();
                foreach (qElement element in connectedElements)
                {
                    if (!element.IsQuad) { triangleElements.Add(element); }
                }
            }

            // delete dublicates
            List<qElement> triangleElementsNoDublicates = new List<qElement>();
            if (triangleElements.Count > 0)
            {
                foreach (qElement element in triangleElements)
                {
                    if (!triangleElementsNoDublicates.Contains(element)) { triangleElementsNoDublicates.Add(element); }
                }
            }

            return triangleElementsNoDublicates;
        } // class: kan bli implementert i: qNode
    }
}

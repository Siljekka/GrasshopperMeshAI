﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeshPoints
{
    class Element
    {
        public Node Node1 { get; set; }
        public Node Node2 { get; set; }
        public Node Node3 { get; set; }
        public Node Node4 { get; set; }

        public Node Node5 { get; set; }
        public Node Node6 { get; set; }
        public Node Node7 { get; set; }
        public Node Node8 { get; set; }

        public Quality QualityElement { get; set; }
        public int Id { get; set; }
        public bool IsCube { get; set; } //if 3D element

        //Constructer

        //_empty
        public Element()
        {
            //Empty constructor
        }

        //_for 2D
        public Element(int _id, Node _node1, Node _node2, Node _node3, Node _node4)
        {
            Id = _id;
            Node1 = _node1;
            Node2 = _node2;
            Node3 = _node3;
            Node4 = _node4;
        }

        //_for 3D
        public Element(int _id, Node _node1, Node _node2, Node _node3, Node _node4, Node _node5, Node _node6, Node _node7, Node _node8)
        {
            Id = _id;
            Node1 = _node1;
            Node2 = _node2;
            Node3 = _node3;
            Node4 = _node4;
            Node5 = _node5;
            Node6 = _node6;
            Node7 = _node7;
            Node8 = _node8;
        }

    }
}

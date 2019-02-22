package com.uwaterloo;

public class Edge {
    Vertex startVertice;
    Vertex endVertice;


    public Edge(Vertex startVertice, Vertex endVertice) {
        this.startVertice = startVertice;
        this.endVertice = endVertice;
    }

    public Vertex getStartVertice() {
        return startVertice;
    }

    public Vertex getEndVertice() {
        return endVertice;
    }
}

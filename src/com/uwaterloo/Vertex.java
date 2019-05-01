package com.uwaterloo;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Vertex {
    int pos;
    char AA;
    //MutationsPattern pattern;
    List<Edge> inEdges;
    List<Edge> outEdges;

    public Vertex(int pos, char AA) {
        this.pos = pos;
        this.AA = AA;
        inEdges = new ArrayList<>();
        outEdges = new ArrayList<>();
    }


    /*
    public Vertex(int pos, MutationsPattern pattern) {
        this.pos = pos;
        this.pattern = pattern;
    }
    */

    public int getPos() {
        return pos;
    }

    public char getAA() {
        return AA;
    }

    /*
    public MutationsPattern getPattern() {
        return this.pattern;
    }

 */
    public List<Edge> getInEdges() {
        return inEdges;
    }


    public List<Edge> getOutEdges() {
        return outEdges;
    }

    @Override
    public String toString() {
        return pos + " " + AA;
    }

    @Override
    public boolean equals(Object o) {
        if (this.getPos() != ((Vertex) o).getPos()) {
            return false;
        }
        if (this.getAA() == ((Vertex) o).getAA()) {
            return true;
        }

        return false;
    }

    @Override
    public int hashCode() {
        int result = 7;
        result = 31 * result + pos;
        result = 31 * result + AA;
        return result;
    }
}

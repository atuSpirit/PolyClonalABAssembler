package com.uwaterloo;

public class Vertex {
    int pos;
    //char AA;
    MutationsPattern pattern;

    public Vertex(int pos, MutationsPattern pattern) {
        this.pos = pos;
        this.pattern = pattern;
    }

    public int getPos() {
        return pos;
    }

    public MutationsPattern getPattern() {
        return this.pattern;
    }
}

package com.uwaterloo.DenovoAssembler;

import com.uwaterloo.Utils.Vertex;

import java.util.ArrayList;
import java.util.List;

public class WeightedVertex extends Vertex {
    int weight; //Storing the best weight of path ended in this Vertex.  Used in dn assembly
    List<WeightedVertex> maxPreVertexList;    //The preVertex list has the largest weight among InEdges. The list's weight is same
    boolean visited = false;

    public WeightedVertex(int pos, char AA) {
        super(pos, AA);
        this.weight = 0;
        this.maxPreVertexList = new ArrayList<>();
    }

    public WeightedVertex(int pos, char AA, int weight) {
        super(pos, AA);
        this.weight = weight;
        this.maxPreVertexList = new ArrayList<>();
    }

    public int getWeight() {
        return weight;
    }

    public void setWeight(int weight) {
        this.weight = weight;
    }

    public void setVisited() {
        this.visited = true;
    }

    public boolean isVisited() {
        return this.visited;
    }

    /**
     * If a new maxPrevertex come in, if the weight is same with exisint maxPretex,
     * add it to the existing list, otherwise, discard the existing one, add it
     * as a new list.
     * @param vertex
     */
    public void setMaxPreVertexList(WeightedVertex vertex) {
        if ((this.maxPreVertexList.size() == 0) ||
                (vertex.getWeight() == this.maxPreVertexList.get(0).getWeight()) ) {
            this.maxPreVertexList.add(vertex);
        } else {
            this.maxPreVertexList = new ArrayList<>();
            this.maxPreVertexList.add(vertex);
        }

    }

    public List<WeightedVertex> getMaxPreVertexList() {
        return this.maxPreVertexList;
    }
}

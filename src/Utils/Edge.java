package Utils;

public class Edge {
    Vertex startVertice;
    Vertex endVertice;
    int weight;

    public Edge(Vertex startVertice, Vertex endVertice) {
        this.startVertice = startVertice;
        this.endVertice = endVertice;
    }

    public Edge(Vertex startVertice, Vertex endVertice, int weight) {
        this.startVertice = startVertice;
        this.endVertice = endVertice;
        this.weight = weight;
    }

    public Vertex getStartVertice() {
        return startVertice;
    }

    public Vertex getEndVertice() {
        return endVertice;
    }

    public int getWeight() {
        return weight;
    }

    public void setWeight(int weight) {
        this.weight = weight;
    }

    @Override
    public boolean equals(Object o) {
        if (!this.getStartVertice().equals(((Edge) o).getStartVertice())) {
            return false;
        }
        if (this.getEndVertice().equals(((Edge) o).getEndVertice())) {
            return true;
        }

        return false;
    }

    @Override
    public int hashCode() {
        int result = 7;
        result = 31 * result + startVertice.toString().hashCode();
        result = 31 * result + endVertice.toString().hashCode();
        return result;
    }

    @Override
    public String toString() {
        return startVertice.toString() + " " + endVertice.toString() + " " + weight;
    }

}

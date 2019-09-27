package com.uwaterloo.DenovoAssembler;

import com.uwaterloo.ScanTemplateMapper.PSMAligned;
import com.uwaterloo.Utils.Edge;

import java.util.*;

public class UncertainRegionGraphAssembler {
    HashMap<String, DenovoOnly> scanDnMap;
    //TreeMap<Integer, List<WeightedVertex>> AAsGraph;

    public UncertainRegionGraphAssembler(HashMap<String, DenovoOnly> scanDnMap) {
        this.scanDnMap = scanDnMap;
        //this.AAsGraph = new TreeMap<>();
    }

    /*
    public TreeMap<Integer, List<WeightedVertex>> getAAsGraph() {
        return this.AAsGraph;
    }
    */

    /* Simple case that does not deal with insertion and deletion in gap region */
    public List<Contig> assembleOneRegionByGraph(UncertainRegion regionToAssemble) {
        System.out.println("region starting from " + regionToAssemble.getStartPos() + " to " + regionToAssemble.getEndPos());

        List<DenovoAligned> dnAlignToRightList = new ArrayList<>(regionToAssemble.getDnAlignToRightSet());
        List<DenovoAligned> dnAlignToLeftList = new ArrayList<>(regionToAssemble.getDnAlignToLeftSet());
        List<PSMAligned> psmAlignedList = new ArrayList<>(regionToAssemble.getPsmAlignedSet());

        TreeMap<Integer, List<WeightedVertex>> assemblyGraph = new TreeMap<>();

        /*
        System.out.println("Dn to right:");
        for (DenovoAligned dnAligned : dnAlignToRightList) {
            System.out.print(dnAligned.gettStart() + " " + new String(scanDnMap.get(dnAligned.getDnScan()).getAAs()) +
                    " " + dnAligned.getScore() + " " + dnAligned.getDnScan());
            for (short score : scanDnMap.get(dnAligned.getDnScan()).getConfScores()) {
                System.out.print(" " + score);
            }
            System.out.println();
        }
        /**/

        buildDnGraph(assemblyGraph, dnAlignToRightList, true);

        /*
        System.out.println("dn to left");
        for (DenovoAligned dnAligned : dnAlignToLeftList) {
            System.out.print(dnAligned.gettEnd() + " " + new String(scanDnMap.get(dnAligned.getDnScan()).getAAs()) +
                    " " +  dnAligned.getScore() + " " + dnAligned.getDnScan());
            for (short score : scanDnMap.get(dnAligned.getDnScan()).getConfScores()) {
                System.out.print(" " + score);
            }
            System.out.println();
        }
         /**/
        buildDnGraph(assemblyGraph, dnAlignToLeftList, false);

        /*
        System.out.println("psm");
        for (PSMAligned psm : psmAlignedList) {
            System.out.print(psm.getStart() + " " + new String(psm.getAAs()));
            for (int score : psm.getIonScores()) {
                System.out.print(" " + score);
            }
            System.out.println();
        }

        /* */

        buildPSMGraph(assemblyGraph, psmAlignedList);

        System.out.println("assembled contigs");
        //Extract the longest pathes which are the assembly with best score
        List<Contig> contigs = findLongestPath(assemblyGraph);

        //Debug
        for (Contig contig : contigs) {
            System.out.println(contig.toString());
        }

        return contigs;
    }

    /* TODO: a assembleOneRegionByGraph that will consider the insertion or deletion */
    public List<Contig> complicated_assembleOneRegionByGraph(UncertainRegion regionToAssemble) {
        System.out.println("region starting from " + regionToAssemble.getStartPos() + " to " + regionToAssemble.getEndPos());

        List<DenovoAligned> dnAlignToRightList = new ArrayList<>(regionToAssemble.getDnAlignToRightSet());
        List<DenovoAligned> dnAlignToLeftList = new ArrayList<>(regionToAssemble.getDnAlignToLeftSet());
        List<PSMAligned> psmAlignedList = new ArrayList<>(regionToAssemble.getPsmAlignedSet());

        //Sort dnToRightList by tStart, if same tStart, sort by score descending.
        Collections.sort(dnAlignToRightList, DenovoAligned.cmpReverseScore());
        Collections.sort(dnAlignToRightList, DenovoAligned.cmpTStart());

        //Sort dnToRightList by tStart, if same tStart, sort by score descending.
        Collections.sort(dnAlignToLeftList, DenovoAligned.cmpReverseScore());
        Collections.sort(dnAlignToLeftList, DenovoAligned.cmpReverseTEnd());

        //Sort psmList by tStart
        Collections.sort(psmAlignedList, PSMAligned.cmpStart());


        // Build assembly graph based on dn and psms for to right region.
        TreeMap<Integer, List<WeightedVertex>> assemblyToRightGraph = new TreeMap<>();
        buildDnGraph(assemblyToRightGraph, dnAlignToRightList, true);
        buildPSMGraph(assemblyToRightGraph, psmAlignedList);

        // Build the assembly graph based on dn and psms for to left region.
        TreeMap<Integer, List<WeightedVertex>> assemblyToLeftGraph = new TreeMap<>();
        buildDnGraph(assemblyToLeftGraph, dnAlignToLeftList, false);
        buildPSMGraph(assemblyToLeftGraph, psmAlignedList);

        //Extract the longest paths which is the assembly with best score
        return null;
    }

    /**
     * Add edges and vertices in dnAlignList to existing assemblyGraph.
     * @param assemblyGraph
     * @param dnAlignList
     */
    void buildDnGraph(TreeMap<Integer, List<WeightedVertex>> assemblyGraph,
                      List<DenovoAligned> dnAlignList,
                      boolean toRight) {
        for (DenovoAligned dnAligned : dnAlignList) {
            char[] AAs = scanDnMap.get(dnAligned.getDnScan()).getAAs();
            short[] confs = scanDnMap.get(dnAligned.getDnScan()).getConfScores();
            if (toRight) {
                processOnePeptide(assemblyGraph, dnAligned.gettStart(), AAs, confs);
            } else {
                processOnePeptide(assemblyGraph, dnAligned.gettEnd() + 1 - AAs.length, AAs, confs);
            }
        }
    }

    /**
     * Add edges and vertices in PSMs including (db result and spider result) to existing
     * assemblyGraph.
     * @param assemblyGraph
     * @param psmAlignedList
     */
    void buildPSMGraph(TreeMap<Integer, List<WeightedVertex>> assemblyGraph, List<PSMAligned> psmAlignedList) {
        for (PSMAligned psmAligned : psmAlignedList) {
            char[] AAs = psmAligned.getAAs();
            short[] confs = psmAligned.getIonScores();
            processOnePeptide(assemblyGraph, psmAligned.getStart(), AAs, confs);
        }
    }

    /**
     * According to the AA and confs of a peptide, add it to the graph.
     * The peptide can be a denovo only, a db or a spider peptide.
     */
    public void processOnePeptide(TreeMap<Integer, List<WeightedVertex>> assemblyGraph, int posStart, char[] AAs, short[] confs) {
        int length = AAs.length;
        List<WeightedVertex> vertexList;

        WeightedVertex sourceWeightedVertex = new WeightedVertex(posStart - 1, '*');

        /* Locate the sourceWeightedVertex in the graph, if not exists, add it. */
        if (!assemblyGraph.containsKey(posStart - 1)) {
            vertexList = new ArrayList<>();
            vertexList.add(sourceWeightedVertex);
            assemblyGraph.put(posStart - 1, vertexList);
        } else {
            vertexList = assemblyGraph.get(posStart - 1);
            WeightedVertex existingWeightedVertex = findWeightedVertex(vertexList, sourceWeightedVertex);
            if (existingWeightedVertex == null) {
                vertexList.add(sourceWeightedVertex);
            } else {
                sourceWeightedVertex = existingWeightedVertex;
            }
        }

        WeightedVertex preWeightedVertex = sourceWeightedVertex;

        /* Add all AAs in the peptide to the graph */
        for (int i = 0; i < length; i++) {
            int currPos = posStart + i;
            WeightedVertex currWeightedVertex = new WeightedVertex(currPos, AAs[i]);

            if (!assemblyGraph.containsKey(currPos)) {
                vertexList = new ArrayList<>();
                vertexList.add(currWeightedVertex);
                assemblyGraph.put(currPos, vertexList);
            } else {
                vertexList = assemblyGraph.get(currPos);
                WeightedVertex existingWeightedVertex = findWeightedVertex(vertexList, currWeightedVertex);
                if (existingWeightedVertex == null) {
                    vertexList.add(currWeightedVertex);
                } else {
                    currWeightedVertex = existingWeightedVertex;
                }
            }
            //Add edge to the two vertices. If the edge exists, update the weight
            Edge newEdge = new Edge(preWeightedVertex, currWeightedVertex, confs[i]);
            Edge edge = findEdge(currWeightedVertex.getInEdges(), newEdge);
            if (edge == null) {
                currWeightedVertex.getInEdges().add(newEdge);
                preWeightedVertex.getOutEdges().add(newEdge);
            } else {
                edge.setWeight(edge.getWeight() + newEdge.getWeight());
            }

            preWeightedVertex = currWeightedVertex;
        }
    }

    /**
     * Given a vertex(pos, AA), find whether it exists in vertexList
     * @param vertexList
     * @param vertex
     * @return null if not find, the vertex in the vertexList if exists.
     */
    private WeightedVertex findWeightedVertex(List<WeightedVertex> vertexList, WeightedVertex vertex) {
        if (vertexList == null) return null;
        for (WeightedVertex v : vertexList) {
            if (v.equals(vertex)) {
                return v;
            }
        }
        return null;
    }

    /**
     * Looking for an edge in a edge list. If not exist, return null, otherwise, return
     * the edge in the edgeList.
     * @param edgeList
     * @param edge
     * @return
     */
    public Edge findEdge(List<Edge> edgeList, Edge edge) {
        if (edgeList == null) return null;
        for (Edge e : edgeList) {
            if (e.equals(edge)) {
                return e;
            }
        }
        return null;
    }

    /**
     * Find the path with max weight(longest path). They are the assembled contigs
     * with bestScore.
     * @param assemblyGraph
     * @return
     */
    List<Contig> findLongestPath(TreeMap<Integer, List<WeightedVertex>> assemblyGraph) {
        List<WeightedVertex> vertexWithMaxWeightList = updateVertexWeight(assemblyGraph);
        List<Contig> longestPathContigs = new ArrayList<>();
        for (WeightedVertex vertexWithMaxWeight : vertexWithMaxWeightList) {
            longestPathContigs.addAll(backtrack(vertexWithMaxWeight));
        }

        return longestPathContigs;
    }

    /**
     * Tranverse all vertices from left to right. Update their weight according to
     * the weight of their precedent and inEdge weight. Keep a list of vertexWithMaxWeight.
     * For each connected graph, there is only one vertexWtihMaxWeight.
     * The list of vertexWithMaxWeight is for where the gap cannot be bridged.
     * @param assemblyGraph
     * @return
     */
    List<WeightedVertex> updateVertexWeight(TreeMap<Integer, List<WeightedVertex>> assemblyGraph) {
        List<WeightedVertex> vertexWithMaxWeightList = new ArrayList<>();
        WeightedVertex vertexWithMaxWeight = null;
        int maxWeight = 0;
        int prePos = -10;
        for (int pos : assemblyGraph.keySet()) {
            //Update the weight of all vertices in this pos
            List<WeightedVertex> vertexList = assemblyGraph.get(pos);
            if (prePos != pos - 1) {
                if (vertexWithMaxWeight != null) {
                    vertexWithMaxWeightList.add(vertexWithMaxWeight);
                }
                maxWeight = 0;
                vertexWithMaxWeight = vertexList.get(0);
            }

            for (WeightedVertex currentVertex : vertexList) {
                List<Edge> inEdges = currentVertex.getInEdges();
                //Find the max weight of currentVertex derived from inEdges
                for (Edge inEdge : inEdges) {
                    WeightedVertex preVertex = (WeightedVertex) inEdge.getStartVertice();
                    int inEdgeWeight = preVertex.getWeight() + inEdge.getWeight();
                    if (inEdgeWeight >= currentVertex.getWeight()) {
                        currentVertex.setMaxPreVertexList(preVertex);
                        currentVertex.setWeight(inEdgeWeight);
                    }
                }
                //Update the maxVertex if needed
                if (currentVertex.getWeight() > maxWeight) {
                    maxWeight = currentVertex.getWeight();
                    vertexWithMaxWeight = currentVertex;
                }
            }
            prePos = pos;
        }
        if (vertexWithMaxWeight != null) {
            vertexWithMaxWeightList.add(vertexWithMaxWeight);
        }

        return vertexWithMaxWeightList;

    }

    /**
     * Depth-First transverse to find the pathes with maxWeight.
     * @param vertexWithMaxWeight
     * @return
     */
    List<Contig> backtrack(WeightedVertex vertexWithMaxWeight) {
        List<Contig> contigList = new ArrayList<>();
        Stack<WeightedVertex> vertexStack = new Stack<>();
        vertexStack.push(vertexWithMaxWeight);

        List<Character> AAList = new ArrayList<>();
        List<Integer> confsList = new ArrayList<>();

        while (!vertexStack.empty()) {
            WeightedVertex vertex = vertexStack.pop();
            //If reached the head of contig, store the contig in contigList
            if (vertex.getAA() == '*') {
                int length = AAList.size();
                char[] AAs = new char[length];
                int[] confs = new int[length];
                int score = 0;
                for (int i = 0; i < length; i++) {
                    AAs[i] = AAList.get(length - 1 - i);
                }
                confs[0] = confsList.get(length - 1);
                score = confs[0];
                for (int i = 1; i < length; i++) {
                    confs[i] = confsList.get(length - 1 - i) - confsList.get(length - i);
                    score += confs[i];
                }

                Contig contig = new Contig(vertex.getPos() + 1, vertex.getPos() + length,
                        AAs, confs, score);
                contigList.add(contig);
            } else {
                AAList.add(vertex.getAA());
                confsList.add(vertex.getWeight());
                vertex.setVisited();

                //Add all preVertexWithMaxWeight in the stack
                for (WeightedVertex preVertex : vertex.getMaxPreVertexList()) {
                    if (!preVertex.isVisited()) {
                        vertexStack.push(preVertex);
                    }
                }
            }
        }
        return contigList;
    }
}

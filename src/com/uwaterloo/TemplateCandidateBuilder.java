package com.uwaterloo;

import Utils.CharEqual;
import Utils.Edge;
import Utils.Vertex;
import com.uwaterloo.ScanTemplateMapper.PSMAligned;
import com.uwaterloo.ScanTemplateMapper.TemplateHooked;
import com.uwaterloo.SignificantMutationsFinder.MutationsPattern;

import java.util.*;

public class TemplateCandidateBuilder {
    List<HashMap<List<Integer>, List<MutationsPattern>>> mutationsList;

    public TemplateCandidateBuilder() {

    }

    public TemplateCandidateBuilder(List<HashMap<List<Integer>, List<MutationsPattern>>> mutationsList) {
        this.mutationsList = mutationsList;
    }

    /* Get the mutation positions set from mutationList*/
    private List<Integer> getPosSet() {
        Set<Integer> posSet = new HashSet<>();
        for (int i = 0; i < mutationsList.size(); i++) {
            HashMap<List<Integer>, List<MutationsPattern>> mutations = mutationsList.get(i);
            for (List<Integer> posList : mutations.keySet()) {
                for (int pos : posList) {
                    posSet.add(pos);
                }
            }
        }
        List<Integer> posList = new ArrayList<>(posSet);
        Collections.sort(posList);

        return posList;
    }

    private TreeMap<Integer, List<MutationsPattern>> sortMutationsAccordingPos() {
        TreeMap<Integer, List<MutationsPattern>> mutationsSorted = new TreeMap<>();
        for (int i = 0; i < mutationsList.size(); i++) {
            HashMap<List<Integer>, List<MutationsPattern>> mutations = mutationsList.get(i);
            for (List<Integer> posList : mutations.keySet()) {
                List<MutationsPattern> patternList = mutations.get(posList);
                int pos = posList.get(0);
                if (!mutationsSorted.containsKey(pos)) {
                    mutationsSorted.put(pos, patternList);
                } else {
                    mutationsSorted.get(pos).addAll(patternList);
                }
            }
        }
        return mutationsSorted;
    }


    private void printMutationsAlongPos(TreeMap<Integer, List<MutationsPattern>> posMutationMap) {
        for (int pos : posMutationMap.keySet()) {
            System.out.println(pos);
            List<MutationsPattern> patterns = posMutationMap.get(pos);
            for (MutationsPattern pattern : patterns) {
                System.out.println(pattern.toString());
            }
        }
    }



    private int indexOfMaxFreq(Set<Integer> indexList, List<MutationsPattern> mutationCongigList) {
        int maxIndex = 0;
        int maxFreq = 0;
        for (int i : indexList) {
            int freq = mutationCongigList.get(i).getFreq();
            if (freq > maxFreq) {
                maxIndex = i;
                maxFreq = freq;
            }
        }
        return maxIndex;
    }


    /**
     * Generate new template sequence according to old template sequence and mutated patterns.
     * @param templateHooked
     * @param topMutationContigList
     * @return A new template sequence
     */
    private char[] getCandidateTemplate(TemplateHooked templateHooked, List<MutationsPattern> topMutationContigList) {
        char[] candidateTemplate = templateHooked.getSeq().clone();
        if (topMutationContigList.size() == 0) {
            return candidateTemplate;
        }

        for (MutationsPattern mutationsPattern: topMutationContigList) {
            List<Integer> posList = mutationsPattern.getPosList();
            String pattern = mutationsPattern.getAAs();
            //DEBUG
            //System.out.println("Debug: " + pattern + " " + posList.toString() + " " + mutationsPattern.getFreq());

            for (int i = 0; i < posList.size(); i++) {
                candidateTemplate[posList.get(i)] = pattern.charAt(i);
            }
        }

        //System.out.println(new String(candidateTemplate));
        return candidateTemplate;
    }



    private List<MutationsPattern> filterMutationContigs(List<MutationsPattern> assembledMutationContigs, int confThresh) {
        List<MutationsPattern> filteredMutationContigs = new ArrayList<>();
        for (MutationsPattern mutationsPattern : assembledMutationContigs) {
            if (mutationsPattern.getFreq() >= confThresh) {
                filteredMutationContigs.add(mutationsPattern);
            }
        }
        return filteredMutationContigs;
    }

    /**
     * For each pos in PosArray, choose the significant AAs and build a vertice for each of them.
     * @return
     */
    private TreeMap<Integer, Map<Character, Vertex>> buildVertice(TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        TreeMap<Integer, Map<Character, Vertex>> vertexesList = new TreeMap<>();
        for (int pos : significantAAsPerPos.keySet()) {
            List<MutationsPattern> mutations = significantAAsPerPos.get(pos);
            Map<Character, Vertex> vertices = new HashMap<>();
            for (MutationsPattern mutation : mutations) {
                char AA = mutation.getAAs().charAt(0);
                Vertex v = new Vertex(pos, AA);
                vertices.put(AA, v);
            }

            vertexesList.put(pos, vertices);
        }
        return vertexesList;
    }

    /**
     * Add the template pattern in vertices position to extendedMutationPatterns.
     * To reduce the combination number of candidate templates.
     * If one of templateAA does not appears in vertices, split the templatePattern into
     * multiple one excluding those position whose AA is not significant.
     * @param extendedMutations
     * @param templateHooked
     * @param verticesMap
     */
    private void addTemplatePattern(TreeMap<Integer, List<MutationsPattern>> extendedMutations,
                                    TemplateHooked templateHooked,
                                    TreeMap<Integer, Map<Character, Vertex>> verticesMap) {
        List<Integer> posList = new ArrayList<>();
        String AAs = "";

        for (int pos : verticesMap.keySet()) {
            char templateAA = templateHooked.getSeq()[pos];

            if (verticesMap.get(pos).containsKey(templateAA)) {
                posList.add(pos);
                AAs += templateAA;
            } else {
                MutationsPattern templatePattern = new MutationsPattern(posList, AAs, 1, 1);
                if (posList.size() < 2) {
                    continue;
                }

                if (extendedMutations.containsKey(posList.get(0))) {
                    extendedMutations.get(posList.get(0)).add(templatePattern);
                } else {
                    List<MutationsPattern> list = new ArrayList<>();
                    list.add(templatePattern);
                    extendedMutations.put(posList.get(0), list);
                }

                posList = new ArrayList<>();
                AAs = "";
            }
        }
        if (!AAs.equals("")) {
            MutationsPattern templatePattern = new MutationsPattern(posList, AAs, 1, 1);
            int pos = posList.get(0);
            if (extendedMutations.containsKey(pos)) {
                extendedMutations.get(pos).add(templatePattern);
            } else {
                List<MutationsPattern> templatePatternList = new ArrayList<>();
                templatePatternList.add(templatePattern);
                extendedMutations.put(posList.get(0), templatePatternList);
            }
        }
    }



    /**
     * Judge whether an extentedMutation contains only significant AA in verticesMap or not.
     * If the mutation contains position not in verticeMap, if the AA is same as that
     * on templateHooked. If the mutation AA does not appear in verticeMap, it is a insignificant
     * extendedMutation.
     * @param extendedMutation
     * @param verticesMap
     * @param templateHooked
     */
    private boolean isSignificantPattern(MutationsPattern extendedMutation,
                                           TreeMap<Integer, Map<Character, Vertex>> verticesMap,
                                            TemplateHooked templateHooked) {

        List<Integer> posList = extendedMutation.getPosList();
        if (posList.size() == 1) {
            //If pattern contain one position, there will be no edge
            return false;
        }
        //Check whether all AAs in the extendedMutation are all significantAA. If not, the mutationPatter will not be used.
        String AAString = extendedMutation.getAAs();

        boolean isSignificantAA = true;
        boolean containSignificantPos = false;

        for (int i = 0; i < posList.size(); i++) {
            int pos = posList.get(i);
            char mutationAA = AAString.charAt(i);

            //If pos is not a significant pos
            if (!verticesMap.containsKey(pos)) {
                //if the AA in extendedMutation equals that one on template,
                //keep going, otherwise, discard this mutation pattern
                if (mutationAA == templateHooked.getSeq()[pos]) {
                    continue;
                } else {
                    isSignificantAA = false;
                    break;
                }
            }

            containSignificantPos = true;
            Map<Character, Vertex> vertexMap = verticesMap.get(pos);
            //If the currAA does not in vertex, it is not significant
            if (!vertexMap.containsKey(mutationAA)) {
                isSignificantAA = false;
                break;
            }
        }
        return isSignificantAA && containSignificantPos;
    }

    /**
     * Add an edge to vertexEdgeList.  If the edge exists, only update the weight.
     * @param edge
     * @param edgeWeight
     * @param vertexEdgeList
     */
    private void addEdgeToVertexEdgeList(Edge edge, int edgeWeight, List<Edge> vertexEdgeList) {
        boolean edgeExist = false;

        for (int i = 0; i < vertexEdgeList.size(); i++) {
            Edge existEdge = vertexEdgeList.get(i);
            if (edge.equals(existEdge)) {
                edgeExist = true;
                vertexEdgeList.get(i).setWeight(existEdge.getWeight() + edgeWeight);
            }
        }
        if (!edgeExist) {
            edge.setWeight(edgeWeight);
            vertexEdgeList.add(edge);
        }
    }

    /**
     * Build edges between adjacent two vertices according to significant extendedMutationPattern.
     * @param extendedMutationsMap The extendedMutationPatterns
     * @param verticesMap a sorted map <position, verticeMap at this position>, verticeMap is
     *                    <AA, vertex> map.
     * @param templateHooked
     * @return The verticesMap is modified. Each vertex is attached with inEdges and outEdges
     *          with proper accumulated weight
     */
    public List<List<Vertex>> buildEdges(TreeMap<Integer, List<MutationsPattern>> extendedMutationsMap,
                                   TreeMap<Integer, Map<Character, Vertex>> verticesMap,
                                    TemplateHooked templateHooked) {
        for (int startPos : extendedMutationsMap.keySet()) {
            List<MutationsPattern> extendedMutations = extendedMutationsMap.get(startPos);
            for (MutationsPattern extendedMutation : extendedMutations) {
                //Skip those extendMutationPatterns containing insignificant AA
                if (!isSignificantPattern(extendedMutation, verticesMap, templateHooked)) continue;

                System.out.println("Debug " + extendedMutation.toString());
                List<Integer> posList = extendedMutation.getPosList();
                String AAString = extendedMutation.getAAs();

                int index = 0;
                int currPos = posList.get(index);
                //find the currPos which is significant
                while (!verticesMap.containsKey(currPos)) {
                    index += 1;
                    currPos = posList.get(index);
                }
                char currAA = AAString.charAt(index);

                while (index < (posList.size() - 1)) {
                    int nextPos = posList.get(index + 1);
                    //find next significant pos
                    if (!verticesMap.containsKey(nextPos)) {
                        index++;
                        continue;
                    }
                    char nextAA = AAString.charAt(index + 1);

                    Vertex vertex1 = verticesMap.get(currPos).get(currAA);
                    Vertex vertex2 = verticesMap.get(nextPos).get(nextAA);
                    Edge edge = new Edge(vertex1, vertex2);
                    Edge edgeCopy = new Edge(vertex1, vertex2);
                    addEdgeToVertexEdgeList(edge, extendedMutation.getScore(), vertex1.getOutEdges());
                    addEdgeToVertexEdgeList(edgeCopy, extendedMutation.getScore(), vertex2.getInEdges());

                    currPos = nextPos;
                    currAA = nextAA;
                    index++;
                }
            }
        }

        //Store verticesMap to a list of vertexlist at each position.
        List<List<Vertex>> verticesList = new ArrayList<>();
        for (int pos : verticesMap.keySet()) {
            List<Vertex> vertices = new ArrayList<>();
            Map<Character, Vertex> vertexMap = verticesMap.get(pos);
            for (Vertex vertex : vertexMap.values()) {
                vertices.add(vertex);
            }
            verticesList.add(vertices);
        }
        return verticesList;
    }

    public void printEdges(TreeMap<Integer, Map<Character, Vertex>> verticesMap) {
        for (int pos : verticesMap.keySet()) {
            System.out.println("pos " + pos + ":");
            Map<Character, Vertex> vertices = verticesMap.get(pos);
            for (Vertex vertex : vertices.values()) {
                /*
                System.out.println("Inner edges");
                for (Edge inEdge : vertex.getInEdges()) {
                    System.out.println(inEdge.toString());
                }
                */
                System.out.println("Outer edges");
                for (Edge outEdge : vertex.getOutEdges()) {
                    System.out.println(outEdge.toString());
                }

            }

        }
    }

    /**
     * For those vertices has no incoming edge, add edge between them and
     * all of the vertices at previous position.
     * @param verticesList
     */
    public void addConnection(List<List<Vertex>> verticesList) {
        int index = 0;
        while (index < (verticesList.size() - 1)) {
            List<Vertex> preVertices = verticesList.get(index);
            List<Vertex> nextVertices = verticesList.get(index + 1);
            Set<Edge> edgesToBeAdded = new HashSet<>();
            for (Vertex preVertex : preVertices) {
                if (preVertex.getOutEdges().size() == 0) {
                    for (Vertex nextVertex : nextVertices) {
                        Edge edge = new Edge(preVertex, nextVertex);
                        addEdgeToVertexEdgeList(edge, 0, edge.getStartVertice().getOutEdges());
                        addEdgeToVertexEdgeList(edge, 0, edge.getEndVertice().getInEdges());
                        //edgesToBeAdded.add(edge);
                    }
                }
            }
            
            for (Vertex nextVertex : nextVertices) {
                if (nextVertex.getInEdges().size() == 0) {
                    for (Vertex preVertex : preVertices) {
                        Edge edge = new Edge(preVertex, nextVertex);
                        addEdgeToVertexEdgeList(edge, 0, edge.getStartVertice().getOutEdges());
                        addEdgeToVertexEdgeList(edge, 0, edge.getEndVertice().getInEdges());
                        //edgesToBeAdded.add(edge);
                    }
                }
            }

            /*
            for (Edge edge : edgesToBeAdded) {
                addEdgeToVertexEdgeList(edge, 0, edge.getStartVertice().getOutEdges());
                addEdgeToVertexEdgeList(edge, 0, edge.getEndVertice().getInEdges());
            }
            */

            index++;
        }
    }

    public List<MutationsPattern> generatePathCombination(List<List<Vertex>> verticesList) {
        List<MutationsPattern> pathCombinations = new ArrayList<>();

        //Put all edges from first position
        int index = 0;
        for (Vertex vertex : verticesList.get(0)) {
            for (Edge edge : vertex.getOutEdges()) {
                List<Integer> posList = new ArrayList<>();
                posList.add(edge.getStartVertice().getPos());
                posList.add(edge.getEndVertice().getPos());
                MutationsPattern prefixPattern = new MutationsPattern(posList,
                        String.valueOf(edge.getStartVertice().getAA()) + String.valueOf(edge.getEndVertice().getAA()),
                                                0, edge.getWeight());
                pathCombinations.add(prefixPattern);
            }
        }
        //Append AA from each position to pathCombination
        index = 1;
        while (index < (verticesList.size() - 1)) {
            List<MutationsPattern> newAddedPrefixPattern = new ArrayList<>();
            for (Vertex vertex1 : verticesList.get(index)) {
                for (MutationsPattern prefixPattern : pathCombinations) {
                    //Find the prefix whose suffix AA is same as vertex1.AA to extend
                    if (prefixPattern.getAAs().endsWith(String.valueOf(vertex1.getAA()))) {
                        int edgeIndex = 0;
                        //If vertex1.getOutEdges contain more than one edge, need to clone new prefixPattern
                        while (edgeIndex < (vertex1.getOutEdges().size() - 1)) {
                            Edge edge = vertex1.getOutEdges().get(edgeIndex);
                            Vertex vertex2 = edge.getEndVertice();
                            List<Integer> newPosList = new ArrayList<>();
                            for (int pos : prefixPattern.getPosList()) {
                                newPosList.add(pos);
                            }
                            newPosList.add(vertex2.getPos());
                            MutationsPattern prefixCopy = new MutationsPattern(newPosList,
                                    prefixPattern.getAAs() + String.valueOf(vertex2.getAA()),
                                    prefixPattern.getFreq(),
                                    prefixPattern.getScore() + edge.getWeight());
                            newAddedPrefixPattern.add(prefixCopy);
                            edgeIndex++;
                        }

                        Edge edge = vertex1.getOutEdges().get(edgeIndex);
                        Vertex vertex2 = edge.getEndVertice();
                        prefixPattern.getPosList().add(vertex2.getPos());
                        prefixPattern.setAAs(prefixPattern.getAAs() + vertex2.getAA());
                        prefixPattern.setScore(prefixPattern.getScore() + edge.getWeight());
                    }
                }
            }
            pathCombinations.addAll(newAddedPrefixPattern);
            /*
            System.out.println("Prefix Path at index: " +index);
            for (MutationsPattern path : pathCombinations) {
                System.out.println("Debug: " + path.toString());
            }
            */

            index++;
        }
        return pathCombinations;
    }

    /**
     * Merge nested patterns with same starting position into a larger one. For example,
     * Pos 66
     * MSSR at [66, 70, 74, 76] freq: 7 score: 2650
     * MSLY at [66, 70, 74, 76] freq: 1 score: 400
     * MSR at [66, 70, 76] freq: 7 score: 1950
     * M at [66] freq: 15 score: 1400
     * MSR is a nested pattern of MSSR, because all the reads of SR already checked in SSR, so
     * its pattern belongs to SSR too.  So delete SR pattern, only remain SSR and SLY. Current method is to keep only
     * those patterns with max length
     * @param mutationSorted a list of MutationPatterns sorted by ascending position. For each position, the patterns are
     *                       sorted descending by the length of the pattern
     * @param posArray
     * @return
     */
    private TreeMap<Integer, List<MutationsPattern>> mergeNestedPatterns(TreeMap<Integer, List<MutationsPattern>> mutationSorted, List<Integer> posArray) {
        TreeMap<Integer, List<MutationsPattern>> mergedMutationPatterns = new TreeMap<>();
        for (Integer pos : posArray) {
            List<MutationsPattern> patterns = mutationSorted.get(pos);
            if (patterns == null) {
                continue;
            }
            int maxLen = patterns.get(0).getPosList().size();
            List<MutationsPattern> patternsWithMaxLen = new ArrayList<>();
            for (MutationsPattern pattern : patterns) {
                if (pattern.getPosList().size() == maxLen) {
                    patternsWithMaxLen.add(pattern);
                } else {
                    boolean nested = false;

                    for (MutationsPattern longPattern : patternsWithMaxLen) {
                        List<Integer> posList = pattern.getPosList();
                        List<Integer> longPosList = longPattern.getPosList();
                        int j = 0;
                        boolean matched = false;
                        for (int i = 0; i < longPosList.size(); i++) {
                            if (posList.get(j) == longPosList.get(i)) {
                                if (pattern.getAAs().charAt(j) == longPattern.getAAs().charAt(i)) {
                                    matched = true;
                                } else {
                                    matched = false;
                                }
                                j++;
                                if (j >= posList.size()) {
                                    break;
                                }
                            }
                        }
                        nested = matched | nested;
                    }
                    if (!nested) {
                        System.err.println("Pattern " + pattern + " is not nested to longer patter at position " + pos);
                    }
                }
            }
            mergedMutationPatterns.put(pos, patternsWithMaxLen);
        }
        return mergedMutationPatterns;
    }


    /**
     * Examine each position in posArray.  For each position, extract all scans covering the position.
     * For each unexamined scan, extract the pattern in posArray to extendPatterns.
     * Since all scans covering the position are those scans containing confident score. Scans without
     * significant mutation has been removed in MutationValidator. Since each mutation pattern
     * is extracted from a scan based on posArray, it cannot be extended any more.
     */
    public TreeMap<Integer, List<MutationsPattern>> extendPatterns(List<Integer> posArray, TemplateHooked templateHooked,
                                                                    HashMap<String, PSMAligned> scanPSMMap) {
        Map<MutationsPattern, MutationsPattern> extendedPatterns = new HashMap<>();
        Set<String> checkScanSet = new HashSet<>();
        for (Integer pos : posArray) {
            //Extract all scans covering this position
            List<String> scanList = templateHooked.getMappedScanList().get(pos);

            for (String scan : scanList) {
                if (checkScanSet.contains(scan)) {
                     continue;
                } else {
                    checkScanSet.add(scan);
                }
                PSMAligned psmAligned = scanPSMMap.get(scan);
                int start = psmAligned.getStart();
                int end = psmAligned.getEnd();

                //Get the pos List of posArray in the range of this peptide
                List<Integer> posList = getPosListInRange(posArray, start, end);

                MutationsPattern extendedPattern =  extractAAs(posList, psmAligned);

                if (extendedPatterns.containsKey(extendedPattern)) {
                    //Update the freq and score
                    extendedPatterns.get(extendedPattern).setFreq(extendedPatterns.get(extendedPattern).getFreq() + 1);
                    extendedPatterns.get(extendedPattern).setScore(extendedPattern.getScore() +
                            extendedPatterns.get(extendedPattern).getScore());
                    //The extendedPattern only contain one intensity in its intensitySet
                    long intensity = extendedPattern.getIntensitySet().iterator().next();
                    //If the intensity has not appeared in extendedPatterns, add it to the intensity set
                    if (!extendedPatterns.get(extendedPattern).getIntensitySet().contains(intensity)) {
                        extendedPatterns.get(extendedPattern).getIntensitySet().add(intensity);
                    }
                } else {
                    extendedPatterns.put(extendedPattern, extendedPattern);
                }
            }
        }
        //Storing the extendedPatterns according to their first position
        TreeMap<Integer, List<MutationsPattern>> extendedPatternsMap = new TreeMap<>();
        for (MutationsPattern extendedPattern : extendedPatterns.values()) {
            int pos = extendedPattern.getPosList().get(0);
            if (extendedPatternsMap.containsKey(pos)) {
                extendedPatternsMap.get(pos).add(extendedPattern);
            } else {
                List<MutationsPattern> extendedPatternsList = new ArrayList<>();
                extendedPatternsList.add(extendedPattern);
                extendedPatternsMap.put(pos, extendedPatternsList);
            }
        }
        return extendedPatternsMap;
    }

    /**
     * Based on the extended mutation patterns, extract frequency and scores of AAs in each position.
     * The frequency is the sum of the AA appears in extended patterns of this position.
     * The score is the sum of the score of the extended patterns containing AA divided by the length of the pattern.
     * @param extendedMutations the mutation patterns extended according to peptides
     * @return mutationPattern list containing variations per pos
     */
    private TreeMap<Integer, List<MutationsPattern>> getVariationsPerPos(TreeMap<Integer, List<MutationsPattern>> extendedMutations) {
        TreeMap<Integer, List<MutationsPattern>> variationsPerPos = new TreeMap<>();
        for (int pos : extendedMutations.keySet()) {
            List<MutationsPattern> mutationsPatternList = extendedMutations.get(pos);
            for (MutationsPattern pattern : mutationsPatternList) {
                List<Integer> posList = pattern.getPosList();
                //For a given pattern, extract the AA appearing in the pattern and put in corresponding mutation pattern
                for (int i = 0; i < posList.size(); i++) {
                    int AApos = posList.get(i);
                    char AA = pattern.getAAs().charAt(i);
                    List<Integer> aPosList = new ArrayList<>();
                    aPosList.add(AApos);
                    String str = String.valueOf(AA);

                    if (!variationsPerPos.containsKey(AApos)) {
                        List<MutationsPattern> patternList = new ArrayList<>();
                        MutationsPattern AAPattern = new MutationsPattern(aPosList, str, pattern.getFreq(),
                                pattern.getScore() / pattern.getAAs().length(), pattern.getIntensitySet());
                        patternList.add(AAPattern);
                        variationsPerPos.put(AApos, patternList);
                    } else {
                        List<MutationsPattern> patternList = variationsPerPos.get(AApos);
                        boolean AAexist = false;
                        for (MutationsPattern AAPattern : patternList) {
                            //if (AAPattern.getAAs().equals(str)) { Change it to the following to allow D/N and Q/E cases
                            if (CharEqual.charEqual(AAPattern.getAAs().charAt(0), AA)) {
                                AAPattern.setFreq(AAPattern.getFreq() + pattern.getFreq());
                                AAPattern.setScore(AAPattern.getScore() + pattern.getScore() / pattern.getAAs().length());
                                AAPattern.getIntensitySet().addAll(pattern.getIntensitySet());
                                AAexist = true;
                            }
                        }
                        if (!AAexist) {
                            patternList.add(new MutationsPattern(aPosList, str, pattern.getFreq(),
                                    pattern.getScore() / pattern.getAAs().length(), pattern.getIntensitySet()));
                        }
                        variationsPerPos.put(AApos, patternList);
                    }
                }
            }
        }
        return variationsPerPos;
    }

    /**
     *
     * Based on the extended mutation patterns, extract frequency and scores of AAs in each position.
     * The frequency is the sum of the AA appears in extended patterns of this position.
     * The score is the sum of the score of the extended patterns containing AA divided by the length of the pattern.
     * @param extendedMutations the mutation patterns extended according to peptides
     * @return mutationPattern list containing variations per pos
     */
    private TreeMap<Integer, List<MutationsPattern>> getVariationsPerPos(TreeMap<Integer, List<MutationsPattern>> extendedMutations,
                                                                         HashMap<String, PSMAligned> scanPSMMap,
                                                                         TemplateHooked templateHooked) {
        TreeMap<Integer, List<MutationsPattern>> variationsPerPos = new TreeMap<>();
    //    List<Integer> mutationPosList = getPosListFromPatterns(extendedMutations);
        for (int pos : extendedMutations.keySet()) {
            List<MutationsPattern> mutationsPatternList = extendedMutations.get(pos);
            for (MutationsPattern pattern : mutationsPatternList) {
                List<Integer> posList = pattern.getPosList();
                //For a given pattern, extract the AA appearing in the pattern and put in corresponding mutation pattern
                for (int i = 0; i < posList.size(); i++) {
                    int AApos = posList.get(i);
                    char AA = pattern.getAAs().charAt(i);
                    List<Integer> aPosList = new ArrayList<>();
                    aPosList.add(AApos);
                    String str = String.valueOf(AA);

                    if (!variationsPerPos.containsKey(AApos)) {
                        List<MutationsPattern> patternList = new ArrayList<>();
                        MutationsPattern AAPattern = new MutationsPattern(aPosList, str, pattern.getFreq(),
                                pattern.getScore() / pattern.getAAs().length(), pattern.getIntensitySet());
                        patternList.add(AAPattern);
                        variationsPerPos.put(AApos, patternList);
                    } else {
                        List<MutationsPattern> patternList = variationsPerPos.get(AApos);
                        boolean AAexist = false;
                        for (MutationsPattern AAPattern : patternList) {
                            if (AAPattern.getAAs().equals(str)) {
                                AAPattern.setFreq(AAPattern.getFreq() + pattern.getFreq());
                                AAPattern.setScore(AAPattern.getScore() + pattern.getScore() / pattern.getAAs().length());
                                AAexist = true;
                            }
                        }
                        if (!AAexist) {
                            patternList.add(new MutationsPattern(aPosList, str, pattern.getFreq(),
                                    pattern.getScore() / pattern.getAAs().length(), pattern.getIntensitySet()));
                        }
                        variationsPerPos.put(AApos, patternList);
                    }
                }
            }
        }
        return variationsPerPos;
    }


    /* Filter the AA under ratio_thresh. All pattern in sigAAsPerPos are a single AA */
    private TreeMap<Integer, List<MutationsPattern>> getSigVariationsPerPos(TreeMap<Integer, List<MutationsPattern>> AAsPerPos,
                                                                            char[] templateAAs, double ratio_thresh, int minFreq) {
        TreeMap<Integer, List<MutationsPattern>> sigAAsPerPos = new TreeMap<>();
        for (int pos : AAsPerPos.keySet()) {
            List<MutationsPattern> patterns = AAsPerPos.get(pos);
            List<MutationsPattern> sigPatterns = new ArrayList<>();
            double scoreSum = 0.0;
            for (MutationsPattern pattern : patterns) {
                scoreSum += pattern.getScore();
            }

            for (MutationsPattern pattern : patterns) {
                if (((pattern.getScore() / scoreSum) >= ratio_thresh) && (pattern.getFreq() >= minFreq)) {
                    sigPatterns.add(pattern);
                }
            }
            if ((sigPatterns.size() == 0)) {
                continue;
            }
            //If only the mutation is same as template AA
            if ((sigPatterns.size() == 1) && (sigPatterns.get(0).getAAs().charAt(0) == templateAAs[pos])) {
                continue;
            }
            sigAAsPerPos.put(pos, sigPatterns);
        }
        return sigAAsPerPos;
    }

    /* Helper function to extract AAs, score, intensity information in subPosList from AAs. */
    private MutationsPattern extractAAs(List<Integer> posList, PSMAligned psmAligned) {
        char[] AAs = psmAligned.getAAs();
        int start = psmAligned.getStart();
        short[] scores = psmAligned.getIonScores();
        if (scores == null) {
            //For the case Peaks does not have fragment ions, set all score to zero.
            scores = new short[AAs.length];
            for (int i = 0; i < AAs.length; i++) {
                scores[i] = 0;
            }
        }
        String extractedAAs = "";
        int score = 0;
        for (int pos : posList) {
            extractedAAs += AAs[pos - start];
            score += scores[pos - start];
        }

        Set<Long> intensitySet = new HashSet<>();
        intensitySet.add((long) psmAligned.getIntensity());
        MutationsPattern extendedPattern = new MutationsPattern(posList, extractedAAs, 1, score, intensitySet);
        return extendedPattern;
    }

    /* Helper function to get pos in posArray in range of [start, end] */
    private List<Integer> getPosListInRange(List<Integer> posArray, int start, int end) {
        List<Integer> posListInRange = new ArrayList<>();
        for (int pos : posArray) {
            if ((pos >= start) && (pos <= end)) {
                posListInRange.add(pos);
            }
        }
        return posListInRange;
    }

    /* Extract the AA on template of positions in posArray. If pos contains only one
       significant AA, the position of template should be changed to this AA.
    */
    private MutationsPattern getAAsOnTemplate(TemplateHooked templateHooked, List<Integer> posArray,
                                             TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        String AAsOnTemplate = "";
        List<Integer> newPosArray = new ArrayList<>();
        for (int pos : posArray) {
            List<MutationsPattern> significantAAs = significantAAsPerPos.get(pos);
            if (significantAAs == null) {
                continue;
            }
            //If pos contains only one significant AA, the position of template should be changed to this AA.
            if (significantAAs.size() == 1) {
                AAsOnTemplate += significantAAs.get(0).getAAs().charAt(0);
            } else {
                AAsOnTemplate += templateHooked.getSeq()[pos];
            }
            newPosArray.add(pos);
        }
        MutationsPattern templatePattern = new MutationsPattern(newPosArray, AAsOnTemplate, 1, 1);
        return templatePattern;
    }



    public List<char[]> buildCandidateTemplate(TemplateHooked templateHooked, HashMap<String, PSMAligned> scanPSMMap,
                                               double significantThreshold, int minFreq) {
        System.out.println(templateHooked.getTemplateAccession());
        List<Integer> posArray = getPosSet();
        System.out.println(posArray.toString());

        //Build a map that the patterns covering a position will be linked to this position. So a pattern might be mapped to more than one position
        TreeMap<Integer, List<MutationsPattern>> mutationSorted = sortMutationsAccordingPos();

        printMutationsAlongPos(mutationSorted);

        System.out.println("Extend mutation patterns...");
        //TreeMap<Integer, List<MutationsPattern>> extendedMutations = extendPatterns(processedMutations, posArray,
        TreeMap<Integer, List<MutationsPattern>> extendedMutations = extendPatterns(posArray, templateHooked, scanPSMMap);
                //extendPatterns(mutationSorted, posArray,
                 //                                                                   templateHooked, scanPSMMap);
        printMutationsAlongPos(extendedMutations);

        System.out.println("Get variations per pos...");
        //Get AA variation for each pos. The mutationPattern contains only one AA
        TreeMap<Integer, List<MutationsPattern>> AAsPerPos = getVariationsPerPos(extendedMutations);
        printMutationsAlongPos(AAsPerPos);

        System.out.println("Filter significant AA per pos...");
        TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos = getSigVariationsPerPos(AAsPerPos,
                                                                    templateHooked.getSeq(), significantThreshold, minFreq);

        printMutationsAlongPos(significantAAsPerPos);


        /* Sum up the intensity for each variation of each position. If an intensity appears multiple times,
            they might be from the same feature.  Then it should be counted only once.
         */
        sumIntensitySet(significantAAsPerPos);
 //       printIntensitySet(significantAAsPerPos);


        /* Get the AAs on template of positions in posArray. If significantAAsPerPos has only one variation, change
            the AA of this position to it.
         */
        MutationsPattern templatePattern = getAAsOnTemplate(templateHooked, posArray, significantAAsPerPos);
        posArray = templatePattern.getPosList();
        System.out.println(posArray.toString());
        //System.out.println(templatePattern.toString());

        System.out.println("Generating candidate templates...");
        List<char[]> candidateTemplates = new ArrayList<>();
        //If there is no significant mutations, return the original templateHooked.
        if (significantAAsPerPos.size() == 0) {
            candidateTemplates.add(templateHooked.getSeq());
            return candidateTemplates;
        }


        boolean polyClonal = true;
        //For Polyclonal, each template might diverge into three.
        boolean homo = true;
        boolean graphAssembly = true;
        boolean useTopScoreNotIntensity = false;    //true use score, false use intensity

        if (!polyClonal) {
            //For mAB, no need second candidate, pick the candidate with highest score
            char[] candidateTemplate1 = templateHooked.getSeq().clone();

            if (useTopScoreNotIntensity) {
                //Choose mutation according to higher score
                changeCandidateTemplateAccordMaxScore(candidateTemplate1, significantAAsPerPos);
            } else {
                //Choose mutation according to higher intensity
                changeCandidateTemplateAccordMaxIntensity(candidateTemplate1, significantAAsPerPos);
            }

            candidateTemplates.add(candidateTemplate1);
        } else {
            if (homo) {
                //Apply homogeneous mutations to template to generate one candidate template
                List<MutationsPattern> homogeneousMutations = getHomogeneousMutations(significantAAsPerPos, templatePattern.getAAs());
                //Debug
                System.out.println("homogenous mutations:");
                for (MutationsPattern pattern : homogeneousMutations) {
                    System.out.println(pattern.toString());
                }
                char[] homoCandidateTemplate = getCandidateTemplate(templateHooked, homogeneousMutations);
                candidateTemplates.add(homoCandidateTemplate);
            } else if (!graphAssembly) {
                //Add the template with homogeneous position altered as one of the template
                List<MutationsPattern> patternList = new ArrayList<>();
                patternList.add(templatePattern);
                char[] templateCandidate = getCandidateTemplate(templateHooked, patternList);
                //TODO: not to include the template one when using intensity
                //candidateTemplates.add(templateCandidate);

                char[] candidateTemplate1 = templateCandidate.clone();
                char[] candidateTemplate2 = templateCandidate.clone();

                if (useTopScoreNotIntensity) {
                    //Choose mutation according to higher score
                    changeCandidateTemplateAccordMaxScore(candidateTemplate1, significantAAsPerPos);
                    changeCandidateTemplateAccordMaxScore(candidateTemplate2, significantAAsPerPos);
                } else {
                    //Choose mutation according to higher intensity
                    changeCandidateTemplateAccordMaxIntensity(candidateTemplate1, significantAAsPerPos);
                    changeCandidateTemplateAccordMaxIntensity(candidateTemplate2, significantAAsPerPos);
                }

                if (!java.util.Arrays.equals(templateCandidate, candidateTemplate1)) {
                    candidateTemplates.add(candidateTemplate1);
                }
                if (!java.util.Arrays.equals(templateCandidate, candidateTemplate2)) {
                    candidateTemplates.add(candidateTemplate2);
                }
            } else {
                /* assemble according to graph */

                System.out.println("Building vertices...");
                TreeMap<Integer, Map<Character, Vertex>> verticesMap = buildVertice(significantAAsPerPos);

                //Add the template pattern in extended pattern to reduce the complexity
                addTemplatePattern(extendedMutations, templateHooked, verticesMap);

                System.out.println("Building Edges...");
                List<List<Vertex>> verticesList = buildEdges(extendedMutations, verticesMap, templateHooked);
                addConnection(verticesList);
                printEdges(verticesMap);
                List<MutationsPattern> pathCombination = generatePathCombination(verticesList);

                Collections.sort(pathCombination, MutationsPattern.cmpReverseScore());
                System.out.println(templateHooked.getTemplateAccession());
                System.out.println("candidate number: " + pathCombination.size());
                for (MutationsPattern path : pathCombination) {
                    //System.out.println(path.toString());
                    List<MutationsPattern> pathList = new ArrayList<>();
                    pathList.add(path);
                    char[] candidateTemplate = getCandidateTemplate(templateHooked, pathList);
                    String infoString = ">" + path.toString() +">";
                    char[] candidateTemplateWithInfo = (infoString + new String(candidateTemplate)).toCharArray();
                    //candidateTemplates.add(candidateTemplate);
                    candidateTemplates.add(candidateTemplateWithInfo);
                }

            }

        }

/*
        List<MutationsPattern> topMutationContigList1 = pickMuationContigsWithTopFreq(filteredExtendedMutations, posMutationsMap);
        char[] candidateTemplate1 = getCandidateTemplate(templateHooked, topMutationContigList1);

        List<MutationsPattern> topMutationContigList2 = pickMuationContigsWithTopFreq(filteredExtendedMutations, posMutationsMap);
        char[] candidateTemplate2 = getCandidateTemplate(templateHooked, topMutationContigList2);
        */

        return candidateTemplates;

    }


    private void sumIntensitySet(TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        for (int pos : significantAAsPerPos.keySet()) {
            for (MutationsPattern pattern : significantAAsPerPos.get(pos)) {
                Set<Long> intensitySet = pattern.getIntensitySet();
                long sumIntensity = 0;
                for (long intensity : intensitySet) {
                    sumIntensity += intensity;
                }
                intensitySet = new HashSet<Long>();
                intensitySet.add(sumIntensity);

                pattern.setIntensitySet(intensitySet);
                //System.out.println(pos + "," + pattern.getAAs() + "," + sumIntensity);
            }
        }
    }

    private void printIntensitySet(TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        for (int pos : significantAAsPerPos.keySet()) {
            for (MutationsPattern pattern : significantAAsPerPos.get(pos)) {
                Set<Long> intensitySet = pattern.getIntensitySet();
                /*
                long sumIntensity = 0;
                for (int intensity : intensitySet) {
                    sumIntensity += intensity;
                    System.out.println(pos + "," + pattern.getAAs() + "," + intensity);
                }
                */
                if (intensitySet.size() > 1) {
                    System.err.println("intensity set size larger than 1!");
                }
                System.out.println(pos + "," + pattern.getAAs() + "," + intensitySet.iterator().next());
            }
        }
    }

    private void changeCandidateTemplateAccordMaxScore(char[] candidateTemplate, TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        for (int pos : significantAAsPerPos.keySet()) {
            List<MutationsPattern> AAsList = significantAAsPerPos.get(pos);

            MutationsPattern patternWithMaxScore = AAsList.get(0);
            int maxScore = patternWithMaxScore.getScore();
            for (MutationsPattern AA : AAsList) {
                if (AA.getScore() > maxScore) {
                    patternWithMaxScore = AA;
                    maxScore = AA.getScore();
                }
            }
            candidateTemplate[pos] = patternWithMaxScore.getAAs().charAt(0);
            patternWithMaxScore.setScore(1);
        }
    }

    private void changeCandidateTemplateAccordMaxIntensity(char[] candidateTemplate, TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        for (int pos : significantAAsPerPos.keySet()) {
            List<MutationsPattern> AAsList = significantAAsPerPos.get(pos);

            MutationsPattern patternWithMaxIntensity = AAsList.get(0);
            long maxIntensity = patternWithMaxIntensity.getIntensitySet().iterator().next();
            for (MutationsPattern AA : AAsList) {
                if (AA.getIntensitySet().iterator().next() > maxIntensity) {
                    patternWithMaxIntensity = AA;
                    maxIntensity = AA.getIntensitySet().iterator().next();
                }
            }
            candidateTemplate[pos] = patternWithMaxIntensity.getAAs().charAt(0);
            Set<Long> intensitySet = new HashSet<>();
            intensitySet.add((long) 1);
            patternWithMaxIntensity.setIntensitySet(intensitySet);

        }
    }

    private List<MutationsPattern> filterSignificantMutations(TreeMap<Integer, List<MutationsPattern>> extendedMutations,
                                                                                TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        List<MutationsPattern> significantExtendedMutations = new ArrayList<>();
        for (List<MutationsPattern> mutationsPatternList : extendedMutations.values()) {
            for (MutationsPattern pattern : mutationsPatternList) {
                boolean significant = true;
                List<Integer> posList = pattern.getPosList();
                for (int i = 0; i < posList.size(); i++) {
                    boolean found = false;
                    List<MutationsPattern> significantAAList = significantAAsPerPos.get(posList.get(i));
                    for (MutationsPattern AA : significantAAList) {
                        //If the ith AA in pattern is found in significantAAList
                        if (pattern.getAAs().charAt(i) == AA.getAAs().charAt(0)) {
                            found = true;
                            break;
                        }
                    }
                    if (found == false) {
                        significant = false;
                        break;
                    }
                }
                if (significant == true) {
                    significantExtendedMutations.add(pattern);
                }
            }
        }
        return significantExtendedMutations;

    }

    /* Find those homogeneous mutations in significant mutations */
    private List<MutationsPattern> getHomogeneousMutations(TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos,
                                                           String templateAAs) {
        List<MutationsPattern> homogeneousMutations = new ArrayList<>();

        for (int pos : significantAAsPerPos.keySet()) {
            if (significantAAsPerPos.get(pos).size() == 1) {
                homogeneousMutations.addAll(significantAAsPerPos.get(pos));
            }
        }
        return homogeneousMutations;
    }


}

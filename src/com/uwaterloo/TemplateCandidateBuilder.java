package com.uwaterloo;

import java.util.*;

public class TemplateCandidateBuilder {
    List<HashMap<List<Integer>, List<MutationsPattern>>> mutationsList;

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

    /* For some pattern, they does not contain the consecutive position from the posList,
        add the AA on template to them to make them consecutive for later procession.
        For exmaple, pos List contain {18, 19, 27}, but one pattern only has {18, 27}.
        Then the AA at 19 on template will be added to the middle.
        This is wrong, because the middle AA could not same with that on template.  It should
        be same as the AA on reads.
     */
    private TreeMap<Integer, List<MutationsPattern>> makePatternsConsecutive(TreeMap<Integer, List<MutationsPattern>> mutationSorted,
                                                           List<Integer> posList, TemplateHooked templateHooked) {
        int size = posList.size();
        for (int startIndex = 0; startIndex < size; startIndex++) {
            if (!mutationSorted.containsKey(posList.get(startIndex))) {
                continue;
            }
            List<MutationsPattern> mutationsPatterns = mutationSorted.get(posList.get(startIndex));
            List<MutationsPattern> processedMutationPatterns = new ArrayList<>();

            /* Preprocess mutationPattern, insert template AA if needed */
            for (MutationsPattern mutationsPattern : mutationsPatterns) {
                List<Integer> subPosList = mutationsPattern.getPosList();
                String pattern = "";
                List<Integer> newSubPosList = new ArrayList<>();
                int subPosListIndex = 0;
                int posListIndex = startIndex + subPosListIndex;
                while (subPosListIndex < subPosList.size()) {
                    while (subPosList.get(subPosListIndex) > posList.get(posListIndex)) {
                        newSubPosList.add(posList.get(posListIndex));
                        pattern += templateHooked.getSeq()[posList.get(posListIndex)];
                        posListIndex++;
                    }
                    pattern += mutationsPattern.AAs.charAt(subPosListIndex);
                    newSubPosList.add(subPosList.get(subPosListIndex));
                    subPosListIndex++;
                    posListIndex++;
                }

                MutationsPattern newMutationPattern = new MutationsPattern(newSubPosList, pattern,
                        mutationsPattern.getFreq(), mutationsPattern.getScore());
                processedMutationPatterns.add(newMutationPattern);

            }
            mutationSorted.put(posList.get(startIndex), processedMutationPatterns);

            /* Testing
            System.out.println(posList.get(startIndex));
            for (MutationsPattern mutationsPattern : processedMutationPatterns) {
                System.out.println(mutationsPattern.toString());
            }
            */
        }
        return mutationSorted;
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

    /**
     * Assemble mutation patterns according to overlaps on patterns.  However, it is not correct,
     * because the adjacent pattern does not contain AA extracted, it might be different. For example,
     * pattern1 and pattern2 overlap at pos i.  pattern1 [..., i-1, i], pattern2 [i, i+1...]. The AA[i-1] of
     * pattern1 is not equals to AA[i-1] of reads containing pattern2.
     *
     * @param posArray
     * @param mutationSorted
     * @return
     */
    private List<MutationsPattern> assembleMutationPatterns(List<Integer> posArray,
                                                            TreeMap<Integer, List<MutationsPattern>> mutationSorted) {
        List<MutationsPattern> mutationContigList = new ArrayList<>();

        int size = posArray.size();
        int endPos = -1;
        int contigNum = 0;
        for (int i = 0; i < size; i++) {
            int pos = posArray.get(i);
            List<MutationsPattern> mutationsPatterns = mutationSorted.get(pos);

            if (mutationsPatterns == null) {
                continue;
            }

            if (pos > endPos) {
                //There will be no overlap, add all new patterns starting from pos
                for (MutationsPattern mutationsPattern : mutationsPatterns) {
                    /*
                    if (mutationsPattern.getPosList().size() == 1) {
                        continue;
                    }
                    */
                    mutationContigList.add(mutationsPattern);
                    List<Integer> posList = mutationsPattern.getPosList();
                    int lastPos = posList.get(posList.size() - 1);
                    if (lastPos > endPos) {
                        endPos = lastPos;
                    }
                }
            } else {
                //There will be overlap position. assemble patterns according to overlap
                //Append overlapped patterns to mutationsPatternstoExtend, update the endPos
                endPos = assemblePatternsByOverlapAA(mutationContigList, endPos,
                        mutationSorted.get(pos));
            }
        }

        //DEBUG
        for (MutationsPattern mutationsPattern : mutationContigList) {
            System.out.println(mutationsPattern.toString());
        }
        return mutationContigList;

    }

    /* need to be rewritten using reads covering */
    private List<MutationsPattern> assembleMutationPatterns(List<Integer> posArray, TemplateHooked templateHooked,
                                                            TreeMap<Integer, List<MutationsPattern>> mutationSorted) {
        TreeMap<Integer, List<Vertex>> vertexesList = buildVertice(posArray, mutationSorted);



        //TODO
        return null;
    }

    /**
     * According to the assembled mutation Contigs list, attach the indexes of mutationContigs covering pos to pos.
     * @param mutationContigList The list of assembled mutation contigs
     * @return a map of <pos, List<index of contigs covering the pos>
     */
    private TreeMap<Integer, Set<Integer>> buildPosAssembledMutationsMap(List<MutationsPattern> mutationContigList) {
        /* Attach the index of mutationContigs to each position the pattern covers */
        TreeMap<Integer, Set<Integer>> posMutationsMap = new TreeMap<>();
        for (int index = 0; index < mutationContigList.size(); index++) {
            List<Integer> posList = mutationContigList.get(index).getPosList();
            for (int pos : posList) {
                if (posMutationsMap.containsKey(pos)) {
                    posMutationsMap.get(pos).add(index);
                } else {
                    Set<Integer> indexSet = new HashSet<Integer>();
                    indexSet.add(index);
                    posMutationsMap.put(pos, indexSet);
                }
            }
        }

        return posMutationsMap;
    }

    /**
     * Pick the mutation contigs with highest frequency on each position. After adding the pattern with
     * maximum freq to a list of contigs with highest freq, remove the pattern from posMutationMap. In this way,
     * second call of this function will result a list of contigs with second highest freq.
     * @param mutationContigList
     * @param posMutationsMap
     * @return a list of mutationContig whose freq are the most in the position.
     *          The mutationContig with max freq is removed from posMutations Map.
     */
    private List<MutationsPattern> pickMuationContigsWithTopFreq(List<MutationsPattern> mutationContigList,
                                                                 TreeMap<Integer, Set<Integer>> posMutationsMap) {
        List<MutationsPattern> topMutationContigList = new ArrayList<>();
        int preIndex = -1;
        for (int pos : posMutationsMap.keySet()) {
            Set<Integer> indexList = posMutationsMap.get(pos);
            if (indexList.size() == 0) {
                continue;
            }
            int maxIndex = indexOfMaxFreq(indexList, mutationContigList);

            if (maxIndex != preIndex) {
                topMutationContigList.add(mutationContigList.get(maxIndex));
                preIndex = maxIndex;
            }
            posMutationsMap.get(pos).remove(maxIndex);
        }

        return topMutationContigList;
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
     * Append patterns in mutationsPatterns to mutationsPatternsToExtend according to overlap.
     * @param mutationContigList current mutation contig list
     * @param mutationsPatterns The mutation patterns to be appended to existing contigs in mutationContigList
     *                          or to be added as a new contig if no overlapping
     * @return The last position in succMutPatList. The mutationsPatternsToExtend are updated
     * to the assembled pattern list
     */
    private int assemblePatternsByOverlapAA(List<MutationsPattern> mutationContigList, int endPos,
                                                               List<MutationsPattern> mutationsPatterns) {
        //Get the latest contig to extend
        List<MutationsPattern> mutationsPatternsToExtend = mutationContigList;

        for (MutationsPattern succMutPat : mutationsPatterns) {
            List<Integer> succPosList = succMutPat.getPosList();
            //Does not consider pattern of size one
            if (succPosList.size() == 1) {
                continue;
            }

            if (succPosList.get(succPosList.size() - 1) > endPos) {
                endPos = succPosList.get(succPosList.size() - 1);
            }
            boolean isOverlap = false;
            for (MutationsPattern preMutPat : mutationsPatternsToExtend) {
                int index = overlapStart(preMutPat, succMutPat);
                //If overlapped, append the succMutPat to existing contigs
                if (index != -1) {
                    isOverlap = true;
                    if (succMutPat.getFreq() > preMutPat.getFreq()) {
                        preMutPat.setFreq(succMutPat.getFreq());
                    }

                    List<Integer> prePosList = preMutPat.getPosList();
                    //If succMatPattern is a substring of preMutPat, no need to update preMutPattern
                    if (prePosList.get(prePosList.size() - 1) > succPosList.get(succPosList.size() - 1)) {
                        continue;
                    }

                    preMutPat.setAAs(preMutPat.getAAs().substring(0, index) + succMutPat.getAAs());
                    List<Integer> newPosList = new ArrayList<>();
                    newPosList.addAll(prePosList.subList(0, index));
                    newPosList.addAll(succPosList);
                    preMutPat.setPosList(newPosList);

//                    prePosList = prePosList.subList(0, index);
//                    prePosList.addAll(succPosList);
                }
            }
            //If succMat does not overlap with any patterns in mutationsPatternstoExtend,
            //it should be added to contigList as a new contig
            if (!isOverlap) {
                mutationContigList.add(succMutPat);
            }
        }

        return endPos;

    }

    /**
     * If preMutPat and succMutPat is overlapped, return the first index of overlapping,
     * otherwise return -1
     * @param preMutPat
     * @param succMutPat
     * @return
     */
    private int overlapStart(MutationsPattern preMutPat, MutationsPattern succMutPat) {
        List<Integer> prePosList = preMutPat.getPosList();
        List<Integer> succPosList = succMutPat.getPosList();
        int overLapStart = -1;
        boolean isOverlapped = false;

       // System.out.println("Debug preMutPat " + preMutPat.getAAs() + " succMutPat: " + succMutPat.getAAs());
        if (prePosList.get(prePosList.size()  - 1) < succPosList.get(0)) {
            return -1;
        }
        int preIndex = 0;
        int succIndex = 0;

        while (preIndex < prePosList.size()) {
            //There is a pos in prePosList == the first pos in succPosList
            if (prePosList.get(preIndex) == succPosList.get(succIndex)) {
                break;
            }
            preIndex++;
        }
        overLapStart = preIndex;
        while ((preIndex < prePosList.size()) && (succIndex < succPosList.size())) {
            if (preMutPat.getAAs().charAt(preIndex) != succMutPat.getAAs().charAt(succIndex)) {
                return -1;
            }
            preIndex++;
            succIndex++;
        }


        return overLapStart;
    }

    /**
     * Generate new template sequence according to old template sequence and mutated patterns.
     * @param templateHooked
     * @param topMutationContigList
     * @return A new template sequence
     */
    private char[] getCandidateTemplate(TemplateHooked templateHooked, List<MutationsPattern> topMutationContigList) {
        char[] candidateTemplate = templateHooked.getSeq().clone();

        for (MutationsPattern mutationsPattern: topMutationContigList) {
            List<Integer> posList = mutationsPattern.getPosList();
            String pattern = mutationsPattern.getAAs();
            //DEBUG
            System.out.println("Debug: " + pattern + " " + posList.toString() + " " + mutationsPattern.getFreq());

            for (int i = 0; i < posList.size(); i++) {
                candidateTemplate[posList.get(i)] = pattern.charAt(i);
            }
        }
        System.out.println(new String(candidateTemplate));
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
    private TreeMap<Integer, List<Vertex>> buildVertice(List<Integer> posArray, TreeMap<Integer, List<MutationsPattern>> mutationSorted) {
        TreeMap<Integer, List<Vertex>> vertexesList = new TreeMap<>();
        for (int pos : posArray) {
            List<MutationsPattern> mutations = mutationSorted.get(pos);
            List<Vertex> vertexs = new ArrayList<>();
            for (MutationsPattern mutation : mutations) {
                Vertex v = new Vertex(pos, mutation);
                vertexs.add(v);
            }

            vertexesList.put(pos, vertexs);
        }
        return vertexesList;
    }

    /**
     * Build edges between adjacent two vertexes.
     * @param templateHooked
     * @param scanPSMMap
     * @param posArray
     * @param vertexes
     * @return a list of edgeList.  The list has the same length with posArray. list[i] contains all the edges starting
     * from pos[i]
     */
    private List<Set<Edge>> buildEdges(TemplateHooked templateHooked, HashMap<String, PSMAligned> scanPSMMap,
                                        List<Integer> posArray, List<Map<Character, Vertex>> vertexes) {
        if (posArray.size() < 2) {
            return null;
        }

        List<Set<Edge>> edgesList = new ArrayList<>();
        for (int i = 0; i < (posArray.size() - 1); i++) {
            Set<Edge> edges = new HashSet<>();
            edgesList.add(edges);
        }

        int index = 0;
        while (index < (posArray.size() - 1)) {
            int currPos = posArray.get(index);
            int nextPos = posArray.get(index + 1);

            List<String> scans = templateHooked.getMappedScanList().get(currPos);
            for (String scan : scans) {
                PSMAligned psmAligned = scanPSMMap.get(scan);
                int pos2 = nextPos - psmAligned.getStart();
                if (pos2 > psmAligned.getAAs().length) {
                    continue;
                }
                int pos1 = currPos - psmAligned.getStart();
                char AA1 = psmAligned.getAAs()[pos1];
                Vertex vertex1 = vertexes.get(currPos).get(AA1);
                if (vertex1 == null) {
                    System.err.println("Error: null vertex");
                }

                char AA2 = psmAligned.getAAs()[pos2];
                Vertex vertex2 = vertexes.get(nextPos).get(AA2);
                if (vertex2 == null) {
                    System.err.println("Error: null vertex");
                }
                Edge edge = new Edge(vertex1, vertex2);
                edgesList.get(index).add(edge);
            }

            index += 1;
        }
        return edgesList;
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
     *  According to spectra in templateHooked and posArray, extend the patterns to include adjacent positions in posArray.
     *  For example,
     *  18
     *  DG at [18, 27] freq: 20 score: 3433 will become [18, 19, 27]
     *  19
     * VGS at [19, 27, 34] freq: 10 score: 2391 -> DVGS at [18, 19, 27, 34]
     */
    private TreeMap<Integer, List<MutationsPattern>> extendPatterns(TreeMap<Integer, List<MutationsPattern>> processedMutations,
                                                                    List<Integer> posArray, TemplateHooked templateHooked,
                                                                    HashMap<String, PSMAligned> scanPSMMap) {
        Map<MutationsPattern, MutationsPattern> extendedPatterns = new HashMap<>();
        for (Integer pos : posArray) {
            List<MutationsPattern> mutationsPatterns = processedMutations.get(pos);
            List<String> scanList = templateHooked.getMappedScanList().get(pos);
            for (MutationsPattern pattern : mutationsPatterns) {
                List<Integer> subPosList = pattern.getPosList();
                for (String scan : scanList) {
                    PSMAligned psmAligned = scanPSMMap.get(scan);
                    int start = psmAligned.getStart();
                    int end = psmAligned.getEnd();
                    if (end > subPosList.get(subPosList.size() - 1)) {
                        continue;
                    }
                    //Get the pos List of posArray in the range of this peptide
                    List<Integer> posList = getPosListInRange(posArray, start, end);
                    String AAs = extractAAs(subPosList, psmAligned).getAAs();
                    if (AAs.equals(pattern.getAAs())) {
                        MutationsPattern extendedPattern =  extractAAs(posList, psmAligned);;
                        if (extendedPatterns.containsKey(extendedPattern)) {
                            //Update the freq and score
                            extendedPatterns.get(extendedPattern).setFreq(extendedPatterns.get(extendedPattern).getFreq() + 1);
                            extendedPatterns.get(extendedPattern).setScore(extendedPattern.getScore() +
                                        extendedPatterns.get(extendedPattern).getScore());
                        } else {
                            extendedPatterns.put(extendedPattern, extendedPattern);
                        }
                    }
                }

            }
        }
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

    /* Helper function to extract AAs in subPosList from AAs. */
    private MutationsPattern extractAAs(List<Integer> posList, PSMAligned psmAligned) {
        char[] AAs = psmAligned.getAAs();
        int start = psmAligned.getStart();
        short[] scores = psmAligned.getIonScores();
        String extractedAAs = "";
        int score = 0;
        for (int pos : posList) {
            extractedAAs += AAs[pos - start];
            score += scores[pos - start];
        }
        MutationsPattern extendedPattern = new MutationsPattern(posList, extractedAAs, 1, score);
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

    public List<char[]> buildCandidateTemplate(TemplateHooked templateHooked, HashMap<String, PSMAligned> scanPSMMap) {
        List<Integer> posArray = getPosSet();
        System.out.println(posArray.toString());

        TreeMap<Integer, List<MutationsPattern>> mutationSorted = sortMutationsAccordingPos();
        /* New method to adopt the method of generating a connecting graph and find all the connected path
        List<Map<Character, Vertex>> vertexes = buildVertice(posArray, mutationSorted);
        //There is problem of generating edges in this way
        List<Set<Edge>> edges = buildEdges(templateHooked, scanPSMMap, posArray, vertexes);
        List<MutationsPattern> mutationContigs = findConnectedPaths(vertexes, edges, posArray);
*/

        /*
        System.out.println("Process all patterns consecutive according to posArray...");
        TreeMap<Integer, List<MutationsPattern>> processedMutations = makePatternsConsecutive(mutationSorted, posArray, templateHooked);
        printMutationsAlongPos(mutationSorted);
*/
        printMutationsAlongPos(mutationSorted);
        System.out.println("Merge nested mutation patterns...");
        TreeMap<Integer, List<MutationsPattern>> processedMutations = mergeNestedPatterns(mutationSorted, posArray);
        printMutationsAlongPos(processedMutations);

        TreeMap<Integer, List<MutationsPattern>> extendedMutations = extendPatterns(processedMutations, posArray,
                                                                                    templateHooked, scanPSMMap);
        printMutationsAlongPos(extendedMutations);

        /* Do not assemble mutation patterns */
        System.out.println("Assembling mutations patterns...");
       // List<MutationsPattern> assembledMutationContigs = assembleMutationPatterns(posArray, mutationSorted);
        List<MutationsPattern> assembledMutationContigs = assembleMutationPatterns(posArray, templateHooked, extendedMutations);

       /* List<MutationsPattern> assembledMutationContigs = new ArrayList<>();
        for (List<MutationsPattern> patternList : mutationSorted.values()) {
            //for (List<MutationsPattern> patternList : processedMutations.values()) {  Wrong one
            assembledMutationContigs.addAll(patternList);
        }
        */

        int confThresh = 1; //Currently use 3 as minimal frequency of the pattern to be considered.
        List<MutationsPattern> filteredMutationContigs = filterMutationContigs(assembledMutationContigs, confThresh);
        TreeMap<Integer, Set<Integer>> posMutationsMap = buildPosAssembledMutationsMap(filteredMutationContigs);

        /* Get the top freq to form template 1, then second top freq to form template 2.  It might contain problem that
            the second one is a mosiac one. In fact, the second one should have part of first one.  Need more judgement
            whether the frequency in one template are similar.
         */
        System.out.println("Generating candidate templates...");
        List<MutationsPattern> topMutationContigList1 = pickMuationContigsWithTopFreq(filteredMutationContigs, posMutationsMap);
        char[] candidateTemplate1 = getCandidateTemplate(templateHooked, topMutationContigList1);

        List<MutationsPattern> topMutationContigList2 = pickMuationContigsWithTopFreq(filteredMutationContigs, posMutationsMap);
        char[] candidateTemplate2 = getCandidateTemplate(templateHooked, topMutationContigList2);

        List<char[]> top2CandidateTemplates = new ArrayList<>();
        top2CandidateTemplates.add(candidateTemplate1);
        top2CandidateTemplates.add(candidateTemplate2);

        return top2CandidateTemplates;

    }




}

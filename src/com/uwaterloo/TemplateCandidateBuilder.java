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
    private TreeMap<Integer, List<Vertex>> buildVertice(List<Integer> posArray, TreeMap<Integer, List<MutationsPattern>> significantAAsPerPos) {
        TreeMap<Integer, List<Vertex>> vertexesList = new TreeMap<>();
        for (int pos : posArray) {
            List<MutationsPattern> mutations = significantAAsPerPos.get(pos);
            List<Vertex> vertices = new ArrayList<>();
            for (MutationsPattern mutation : mutations) {
                Vertex v = new Vertex(pos, mutation);
                vertices.add(v);
            }

            vertexesList.put(pos, vertices);
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
     * Will not used any more.
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
            if (pos == 47) {
                System.out.println("Debug in extendPatterns");
            }
            //Extract all scans covering this position
            List<String> scanList = templateHooked.getMappedScanList().get(pos);

            //Extract all mutation patterns starting from this position
            List<MutationsPattern> mutationsPatterns = processedMutations.get(pos);
            /* If the mutation Patterns is null for pos, the position contains mutation, however,
               no mutation patterns starting from this position.  We should extract a single length AA pattern of this position
               to form a new mutation pattern from the scans covering this position */
            if (mutationsPatterns == null) {
                List<Integer> AAPosList = new ArrayList<>();
                AAPosList.add(pos);
                for (String scan : scanList) {
                    PSMAligned psmAligned = scanPSMMap.get(scan);
                    MutationsPattern extendedPattern = extractAAs(AAPosList, psmAligned);

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
                continue;
            }

            //For each mutation pattern, extend it based on the scans covering this position
            for (MutationsPattern pattern : mutationsPatterns) {
                List<Integer> subPosList = pattern.getPosList();

                for (String scan : scanList) {
                    PSMAligned psmAligned = scanPSMMap.get(scan);
                    int start = psmAligned.getStart();
                    int end = psmAligned.getEnd();
                    if (end < subPosList.get(subPosList.size() - 1)) {
                        continue;
                    }
                    //Get the pos List of posArray in the range of this peptide
                    List<Integer> posList = getPosListInRange(posArray, start, end);
                    String AAs = extractAAs(subPosList, psmAligned).getAAs();
                    if (AAs.equals(pattern.getAAs())) {
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
     * Rewrite it to not to scan based on mutationPatterns.
     * Examine each position in posArray.  For each position, extract all scans covering the position.
     * For each unexamined scan, extract the pattern in posArray to extendPatterns.
     * Since all scans covering the position are those scans containing confident score. Scans without
     * significant mutation has been removed in MutationValidator. Since each mutation pattern
     * is extracted from a scan based on posArray, it cannot be extended any more.
     */
    private TreeMap<Integer, List<MutationsPattern>> extendPatterns(List<Integer> posArray, TemplateHooked templateHooked,
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
     *
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
                                                                            char[] templateAAs, double ratio_thresh) {
        TreeMap<Integer, List<MutationsPattern>> sigAAsPerPos = new TreeMap<>();
        for (int pos : AAsPerPos.keySet()) {
            List<MutationsPattern> patterns = AAsPerPos.get(pos);
            List<MutationsPattern> sigPatterns = new ArrayList<>();
            double scoreSum = 0.0;
            for (MutationsPattern pattern : patterns) {
                scoreSum += pattern.getScore();
            }

            for (MutationsPattern pattern : patterns) {
                if ((pattern.getScore() / scoreSum) >= ratio_thresh) {
                    sigPatterns.add(pattern);
                }
            }
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

    public List<char[]> buildCandidateTemplate(TemplateHooked templateHooked, HashMap<String, PSMAligned> scanPSMMap, double significantThreshold) {
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
                                                                    templateHooked.getSeq(), significantThreshold);

        printMutationsAlongPos(significantAAsPerPos);

        /* Sum up the intensity for each variation of each position. If an intensity appears multiple times,
            they might be from the same feature.  Then it should be counted only once.
         */
        sumIntensitySet(significantAAsPerPos);
 //       printIntensitySet(significantAAsPerPos);

        /*
        List<MutationsPattern> filteredExtendedMutations = filterSignificantMutations(extendedMutations,
                                                                                significantAAsPerPos);

        TreeMap<Integer, Set<Integer>> posMutationsMap = buildPosAssembledMutationsMap(filteredExtendedMutations);
        */
        /* Do not assemble mutation patterns */
        /*
        System.out.println("Assembling mutations patterns...");
        //List<MutationsPattern> assembledMutationContigs = assembleMutationPatterns(posArray, extendedMutations);
       List<MutationsPattern> assembledMutationContigs = assembleMutationPatterns(posArray, templateHooked,
               extendedMutations, significantAAsPerPos);
*/
       /* List<MutationsPattern> assembledMutationContigs = new ArrayList<>();
        for (List<MutationsPattern> patternList : mutationSorted.values()) {
            //for (List<MutationsPattern> patternList : processedMutations.values()) {  Wrong one
            assembledMutationContigs.addAll(patternList);
        }
        */
/*
        int confThresh = 1; //Currently use 3 as minimal frequency of the pattern to be considered.
        List<MutationsPattern> filteredMutationContigs = filterMutationContigs(assembledMutationContigs, confThresh);
        TreeMap<Integer, Set<Integer>> posMutationsMap = buildPosAssembledMutationsMap(filteredMutationContigs);
*/
        /* Get the top freq to form template 1, then second top freq to form template 2.  It might contain problem that
            the second one is a mosiac one. In fact, the second one should have part of first one.  Need more judgement
            whether the frequency in one template are similar.
         */

        /* Get the AAs on template of positions in posArray. If significantAAsPerPos has only one variation, change
            the AA of this position to it.
         */
        MutationsPattern templatePattern = getAAsOnTemplate(templateHooked, posArray, significantAAsPerPos);
        posArray = templatePattern.getPosList();
        System.out.println(posArray.toString());
        System.out.println(templatePattern.toString());

        System.out.println("Generating candidate templates...");
        List<char[]> candidateTemplates = new ArrayList<>();

        boolean homo = false;

        if (homo) {
            //Apply homogeneous mutations to template to generate one candidate template
            List<MutationsPattern> homogeneousMutations = getHomogeneousMutations(significantAAsPerPos, templatePattern.getAAs());
            char[] homoCandidateTemplate = getCandidateTemplate(templateHooked, homogeneousMutations);
            candidateTemplates.add(homoCandidateTemplate);
        } else {
            //Add the template with homogeneous position altered as one of the template
            List<MutationsPattern> patternList = new ArrayList<>();
            patternList.add(templatePattern);
            char[] templateCandidate = getCandidateTemplate(templateHooked, patternList);
            candidateTemplates.add(templateCandidate);


            char[] candidateTemplate1 = templateCandidate.clone();
            char[] candidateTemplate2 = templateCandidate.clone();

            boolean useTopScoreNotIntensity = false;    //true use score, false use intensity
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
        int i = 0;
        for (int pos : significantAAsPerPos.keySet()) {
            if (significantAAsPerPos.get(pos).size() == 1) {
                homogeneousMutations.addAll(significantAAsPerPos.get(pos));
            }
            i++;
        }
        return homogeneousMutations;
    }


}

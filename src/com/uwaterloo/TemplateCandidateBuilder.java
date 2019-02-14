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

    public List<char[]> buildCandidateTemplate(TemplateHooked templateHooked) {
        List<Integer> posArray = getPosSet();
        System.out.println(posArray.toString());

        TreeMap<Integer, List<MutationsPattern>> mutationSorted = sortMutationsAccordingPos();

        System.out.println("Process all patterns consecutive according to posArray...");
        TreeMap<Integer, List<MutationsPattern>> processedMutations = makePatternsConsecutive(mutationSorted, posArray, templateHooked);
        printMutationsAlongPos(mutationSorted);

        System.out.println("Assembling mutations patterns...");
        List<MutationsPattern> assembledMutationContigs = assembleMutationPatterns(posArray, processedMutations);

        int confThresh = 3; //Currently use 3 as minimal frequency of the pattern to be considered.
        List<MutationsPattern> filteredMutationContigs = filterMutationContigs(assembledMutationContigs, confThresh);
        TreeMap<Integer, Set<Integer>> posMutationsMap = buildPosAssembledMutationsMap(filteredMutationContigs);

        /* Get the top freq to form template 1, then second top freq to form template 2.  It might contain problem that
            the second one is a mosiac one. In fact, the second one should have part of first one.  Need more judgement
            whether the frequency in one template are similar.
         */
        System.out.println("Generating candidate templates...");
        List<MutationsPattern> topMutationContigList1 = pickMuationContigsWithTopFreq(assembledMutationContigs, posMutationsMap);
        char[] candidateTemplate1 = getCandidateTemplate(templateHooked, topMutationContigList1);

        List<MutationsPattern> topMutationContigList2 = pickMuationContigsWithTopFreq(assembledMutationContigs, posMutationsMap);
        char[] candidateTemplate2 = getCandidateTemplate(templateHooked, topMutationContigList2);

        List<char[]> top2CandidateTemplates = new ArrayList<>();
        top2CandidateTemplates.add(candidateTemplate1);
        top2CandidateTemplates.add(candidateTemplate2);

        return top2CandidateTemplates;

    }




}

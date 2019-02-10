package com.uwaterloo;

import java.util.*;

public class TemplateCandidateBuilder {
    List<HashMap<List<Integer>, List<String>>> mutationsList;

    public TemplateCandidateBuilder(List<HashMap<List<Integer>, List<String>>> mutationsList) {
        this.mutationsList = mutationsList;
    }

    /* Get the mutation positions set from mutationList*/
    private List<Integer> getPosSet() {
        Set<Integer> posSet = new HashSet<>();
        for (int i = 0; i < mutationsList.size(); i++) {
            HashMap<List<Integer>, List<String>> mutations = mutationsList.get(i);
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
            HashMap<List<Integer>, List<String>> mutations = mutationsList.get(i);
            for (List<Integer> posList : mutations.keySet()) {
                List<String> patterns = mutations.get(posList);
                int pos = posList.get(0);
                for (String pattern : patterns) {
                    String[] items = pattern.split("_");
                    MutationsPattern mutationsPattern = new MutationsPattern(posList, items[0],
                            Integer.valueOf(items[1]));
                    if (!mutationsSorted.containsKey(pos)) {
                        List<MutationsPattern> newPatternList = new ArrayList<>();
                        mutationsSorted.put(pos, newPatternList);
                    }
                    mutationsSorted.get(pos).add(mutationsPattern);
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

                MutationsPattern newMutationPattern = new MutationsPattern(newSubPosList, pattern, mutationsPattern.getFreq());
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

    private TreeMap<Integer, List<MutationsPattern>> assembleMutationPatterns(List<Integer> posArray,
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

        for (MutationsPattern mutationsPattern : mutationContigList) {
            System.out.println(mutationsPattern.toString());

        }
        return null;

    }

    /* Given a list of patterns with same start but different pattern or different length,
        Return a list of merged pattern
     */
    private void processPatternsStartingAtPos(List<MutationsPattern> mutationsPatterns,
                                              TemplateHooked templateHooked,
                                              List<Integer> posList, int startIndex) {
        List<MutationsPattern> processedMutationPatterns = new ArrayList<>();


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
        int preIndex = -1;
        int succIndex = 0;
        while (preIndex < prePosList.size()) {
            preIndex++;
            //There is a pos in prePosList == the first pos in succPosList
            if (prePosList.get(preIndex) == succPosList.get(succIndex)) {
                break;
            }
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



    public void buildCandidateTemplate(TemplateHooked templateHooked) {
        List<Integer> posArray = getPosSet();
        System.out.println(posArray.toString());

        TreeMap<Integer, List<MutationsPattern>> mutationSorted = sortMutationsAccordingPos();

        System.out.println("Process all patterns consecutive according to posArray");
        TreeMap<Integer, List<MutationsPattern>> processedMutations = makePatternsConsecutive(mutationSorted, posArray, templateHooked);
        printMutationsAlongPos(mutationSorted);

        System.out.println("Assembling mutations patterns");
        TreeMap<Integer, List<MutationsPattern>> assembledMutations = assembleMutationPatterns(posArray, processedMutations);


        //makePatternsConsecutive(templateHooked, posArray);
    }




}

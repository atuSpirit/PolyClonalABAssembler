package com.uwaterloo;

import java.util.*;

/* Validate mutations in the spectrum to decide whether to discard the psm
    with mutation or change the amino acid on the template.
  */
public class MutationValidator {
    /**
     * Find the maximum number of mutations each PSMs contains.
     * @param templateHooked
     * @return
     */
    int findMaxMutationNumPerPSM(TemplateHooked templateHooked) {
        int maxMutationNum = 0;
        ArrayList<ArrayList<PSMAligned>> listOfSpiderList = templateHooked.getSpiderList();
        for (ArrayList<PSMAligned> psmAlignedList : listOfSpiderList) {
            for (PSMAligned psmAligned : psmAlignedList) {
                if (psmAligned.getPositionOfVariations() == null) {
                    continue;
                }
                int mutationNum = psmAligned.getPositionOfVariations().size();
                if (maxMutationNum < mutationNum) {
                    maxMutationNum = mutationNum;
                }
            }
        }

        return maxMutationNum;
    }

    /* Extract positions where containing psm with mutationNum mutations */
    public ArrayList<Integer> extractPositionWithMutationNum(TemplateHooked templateHooked, int mutationNum) {
        ArrayList<Integer> posWithMutationNumList = new ArrayList<>();
        ArrayList<ArrayList<PSMAligned>> listOfSpiderList = templateHooked.getSpiderList();
        int size = listOfSpiderList.size();
        for (int pos = 0; pos < size; pos++) {
            ArrayList<PSMAligned> psmAlignedList = listOfSpiderList.get(pos);
            for (PSMAligned psmAligned : psmAlignedList) {
                if (psmAligned.getPositionOfVariations() == null) {
                    continue;
                }
                if (psmAligned.getPositionOfVariations().size() == mutationNum) {
                    posWithMutationNumList.add(pos);
                    break;
                }
            }
        }
        return posWithMutationNumList;
    }

    /**
     * Validate mutation combinations given mutated number and start position
     * @param templateHooked    The template to be validated
     * @param pos   The start position of the template containing the spectrum with mutated
     * @param mutatedNum
     * @return the first PSMAligned which start at pos and have mutatedNum mutations, null if not found.
     */
    public List<PSMAligned> getListOfPSMAlignedGivenPosAndMutatedNum(TemplateHooked templateHooked, int pos, int mutatedNum) {
        List<PSMAligned> psmAlignedList = templateHooked.getSpiderList().get(pos);
        List<PSMAligned> psmListGivenPosAndMutatedNum = new LinkedList<>();
        for (PSMAligned psmAligned : psmAlignedList) {
            if (psmAligned.getPositionOfVariations().size() == mutatedNum) {
                psmListGivenPosAndMutatedNum.add(psmAligned);
            }
        }
        return psmListGivenPosAndMutatedNum;
    }

    /**
     * Decide whether the PSM with mutation should be kept or discarded on the templateHooked.
     * Get all scans crossing the mutation position list.
     * Return a list of significant mutated list on template
     * @param templateHooked
     * @param psmAligned
     * @param scanPSMMap
     * @return <list<position>, mutationPattern>
     */
    private HashMap<List<Integer>, List<String>> validateOnePSMAligned(TemplateHooked templateHooked, PSMAligned psmAligned,
                                                                 HashMap<String, PSMAligned> scanPSMMap) {
        double threshold = 0.1;

        /* validate the psms having same list of variation positions */
        HashMap<String, String> extractedScanAAPatternMap = buildScanAAPatternMapForOnePSM(templateHooked,
                psmAligned, scanPSMMap);
        //printExtractedPSMs(extractedPSMsMap);
        String templatePattern = extractedScanAAPatternMap.get("Template");
        extractedScanAAPatternMap.remove("Template");
        HashMap<String, Integer> patternFreqTable = countPatternFreqTable(extractedScanAAPatternMap);
        List<String> significantMutatedPatternsOnTemplate = getSignificantMutatedList(patternFreqTable,
                threshold);

        /* Remove scans and PSMAlign containing unsignificant mutation patterns */
        HashSet<String> significantPatterns = new HashSet<>(significantMutatedPatternsOnTemplate);
        removeUnsignificantMutatedScans(templateHooked, psmAligned.getStart(),
                extractedScanAAPatternMap, significantPatterns);

        /* If significant patterns contain one common AA which appears in template,
            this mutated pattern is not of length mutatedNum. return null
            mutatedAAsOnTemplate. The mutated position will be considered in
             later cycle with less mutatedNum.
         */
        boolean isShrinked = checkPatternShrink(templatePattern, significantPatterns);
        if (isShrinked) {
            return null;
        }
        /* The frequency of significant patterns are extracted from PatternFreqTable and
            attached to pattern. */
        List<String> sigMutationsWithFreqOnTemplate = attachFreqToSignificantPattern(significantMutatedPatternsOnTemplate,
                patternFreqTable);

        HashMap<List<Integer>, List<String>> mutatedAAsOnTemplate = new HashMap<>();
        /* heterogous AA on template */
//        if (sigMutationsWithFreqOnTemplate.size() > 1) {
        int start = psmAligned.getStart();
        List<Integer> posList = new ArrayList<>();
        for (int i : psmAligned.getPositionOfVariations()) {
            posList.add(i + start);
        }

        mutatedAAsOnTemplate.put(posList, sigMutationsWithFreqOnTemplate);
//      }
        return mutatedAAsOnTemplate;
    }

    private List<String> attachFreqToSignificantPattern(List<String> significantMutatedPatternsOnTemplate,
                                                HashMap<String, Integer> patternFreqTable) {
        List<String> patternWithFreqList = new ArrayList<>();
        for (String pattern : significantMutatedPatternsOnTemplate) {
            String patternFreq = pattern + "_" + patternFreqTable.get(pattern);
            patternWithFreqList.add(patternFreq);
        }
        return patternWithFreqList;

    }

    /* If all ith position of patterns equals to the AA on templatePattern, this position is
        not a mutation. The pattern should be shrink to a pattern with less mutationNum.
     */
    private boolean checkPatternShrink(String templatePattern, HashSet<String> significantPatterns) {
        int size = templatePattern.length();
        for (int i = 0; i < size; i++) {
            char AA = templatePattern.charAt(i);
            boolean equalToTemplate = true;
            for (String pattern : significantPatterns) {
                if (pattern.charAt(i) != AA) {
                    equalToTemplate = false;
                }
            }
            if (equalToTemplate) {
                return true;
            }
        }
        return false;
    }

    private void removeUnsignificantMutatedScans(TemplateHooked templateHooked, int pos,
                                                 HashMap<String, String> scanAAPatternMap,
                                                 Set<String> significantPatterns) {
        Set<String> unsignificantScanSet = new HashSet<>();
        for (String scan : scanAAPatternMap.keySet()) {
            if (significantPatterns.contains(scanAAPatternMap.get(scan)) == false) {
                unsignificantScanSet.add(scan);
            }
        }

        /* remove from spiderList of position pos on templateHooked */
        ArrayList<PSMAligned> spiderList = templateHooked.getSpiderList().get(pos);
        int spiderListSize = spiderList.size();
        ArrayList<PSMAligned> newSpiderList = new ArrayList<>();
        for (int index = 0; index < spiderListSize; index++) {
            PSMAligned psmAligned = spiderList.get(index);
            if (!unsignificantScanSet.contains(psmAligned.getScan())) {
                newSpiderList.add(psmAligned);
                spiderList.set(index, null);
            }
        }
        templateHooked.getSpiderList().set(pos, newSpiderList);

        /* remove from scanList of position pos on templateHooked */
        ArrayList<String> scanList = templateHooked.getMappedScanList().get(pos);
        int scanListSize = scanList.size();
        ArrayList<String> newScanList = new ArrayList<>();
        for (int index = 0; index < scanListSize; index++) {
            String scan = scanList.get(index);
            if (!unsignificantScanSet.contains(scan)) {
                newScanList.add(scan);
                scanList.set(index, null);
            }
        }
        templateHooked.getMappedScanList().set(pos, newScanList);
    }


    /* return those patterns whose freq greater than threshold in the top 3 patterns */
    private List<String> getSignificantMutatedList(HashMap<String, Integer> patternFreqTable,
                                                                     double threshold) {
        List<String> significantMutatedPattern = new ArrayList<>();
        /* Get the sum of scans */
        double total = 0.0;
        for (int value : patternFreqTable.values()) {
            total += value;
        }

        /* Get the 3 patterns with top 3 frequency */
        String[] top3Patterns = getPatternsWithTop3Freq(patternFreqTable);

        /* Print significant freq*/
        for (String pattern : top3Patterns) {
            //Skip the case that there is less than 3 patterns in total
            if (pattern == null) {
                break;
            }
            int freq = patternFreqTable.get(pattern);
            double ratio = freq / total;

            if (ratio > threshold) {
                significantMutatedPattern.add(pattern);
//                System.out.printf("%s %d %.4f\n", pattern, freq, freq / total);
            }
        }
        return significantMutatedPattern;
    }

    /**
     * Decide whether the PSM with mutation should be kept or discarded on the templateHooked.
     * Get all scans crossing the mutation position list.
     * For psmAligned, extract those psms having same list of variation positions.
     * For each scan, extract the AA patterns formed by mutation position list
     * @param templateHooked
     * @param psmAligned
     * @param scanPSMMap
     * @return a Map between scan and its AA pattern
     */
    private HashMap<String, String> buildScanAAPatternMapForOnePSM(TemplateHooked templateHooked, PSMAligned psmAligned,
                                                                   HashMap<String, PSMAligned> scanPSMMap) {
        int start = psmAligned.getStart();
        List<Integer> posList = new ArrayList<>();
        for (int i : psmAligned.getPositionOfVariations()) {
            posList.add(start + i);
        }

        int mutationNum = posList.size();
        char[] templateAAPattern = new char[mutationNum];
        Set<String> scanSet = new HashSet<>();

        for (int i = 0; i < mutationNum; i++) {
            int mutationPos = posList.get(i);

            scanSet.addAll(templateHooked.getMappedScanList().get(mutationPos));
            templateAAPattern[i] = templateHooked.getSeq()[mutationPos];
        }

//        System.out.println("Template AA Pattern: " + new String(templateAAPattern));

        HashMap<String, String> scanAAPatternMap = new HashMap<>();
        /* Add the template pattern in scanAAPatternMap */
        scanAAPatternMap.put("Template", new String(templateAAPattern));

        for (String scan : scanSet) {
            PSMAligned correspondingPsmAligned = scanPSMMap.get(scan);
            String AAPattern = extractAAPattern(correspondingPsmAligned, posList, templateAAPattern);
            //If correspondingPSMAligned has shorter posList, null will be returned.
            if (AAPattern != null) {
                scanAAPatternMap.put(scan, AAPattern);
                //DEBUG
                if (AAPattern.equals("LNS")) {
                    //Got result, this happens only on template
                    System.out.println("Debug: LNS scan: " + scan);
                }
            }

        }

        return scanAAPatternMap;
    }


    /**
     * Count the frequency of each pattern
     * @param scanAAPatternMap a map between <scan, pattern>
     * @return a lookup table <pattern, its frequency>
     */
    private HashMap<String, Integer> countPatternFreqTable(HashMap<String, String> scanAAPatternMap) {
        HashMap<String, Integer> patternFreqTable = new HashMap<>();

        for (String AAPattern : scanAAPatternMap.values()) {
            if (patternFreqTable.containsKey(AAPattern)) {
                patternFreqTable.put(AAPattern, patternFreqTable.get(AAPattern) + 1);
            } else {
                patternFreqTable.put(AAPattern, 1);
            }
        }

        return patternFreqTable;
    }

     /**
     * Extract the amino acids from positions in the posList from one PSM
     * (ins) won't cause problem.  (del) will cause problem, need to rethink.
      * If this psmAligned's peptide does not contains full list of posList,
      * null will be returned. This psmAligned will be considered in another
      * cycle with less mutationNum.
     * @param psmAligned
     * @param posList
     * @param templateAAPattern
     * @return  The AA combination of posList of this psm
     */
    private String extractAAPattern(PSMAligned psmAligned, List<Integer> posList,
                                      char[] templateAAPattern) {
        int start = psmAligned.getStart();
        int end = psmAligned.getEnd() - start;

        if (psmAligned.getPeptide().contains("del")) {
            System.err.println("Error in extractAAComb() : scan " + psmAligned.getScan() +
                    " contains del " + psmAligned.getPeptide());
        }

        int size = posList.size();
        char[] AAComb = new char[size];
        for (int i = 0; i < posList.size(); i++) {
            int pos = posList.get(i) - start;
            if ((pos < 0) || (pos > end)) {
                return null;
                //AAComb[i] = templateAAPattern[i];
                //AAComb[i] = '*';
            } else {
                AAComb[i] = psmAligned.getAAs()[pos];
            }

        }
        String AAPattern = new String(AAComb);

        return AAPattern;
    }

    /**
     * Get patterns with top 3 frequency
     * @param patternFreqTable
     * @return No more than 3 patterns with top frequency
     */
    private String[] getPatternsWithTop3Freq(HashMap<String, Integer> patternFreqTable) {
        String[] top3Patterns = new String[3];
        int[] top3Freqs = {-1, -1, -1};
        for (Map.Entry<String, Integer> entry : patternFreqTable.entrySet()) {
            int freq = entry.getValue();
            if (freq >= top3Freqs[0]) {
                top3Freqs[2] = top3Freqs[1];
                top3Freqs[1] = top3Freqs[0];
                top3Freqs[0] = freq;
                top3Patterns[2] = top3Patterns[1];
                top3Patterns[1] = top3Patterns[0];
                top3Patterns[0] = entry.getKey();
            } else if (freq >= top3Freqs[1]) {
                top3Freqs[2] = top3Freqs[1];
                top3Patterns[2] = top3Patterns[1];
                top3Freqs[1] = freq;
                top3Patterns[1] = entry.getKey();
            } else if (freq >= top3Freqs[0]) {
                top3Freqs[2] = freq;
                top3Patterns[2] = entry.getKey();
            }
        }
        return top3Patterns;
    }

    /* Test function to validate the patterns extracted */
    private void printExtractedPSMs(HashMap<String,char[]> extractedPSMsMap) {
        for (Map.Entry<String, char[]> entry : extractedPSMsMap.entrySet()) {
            String scan = entry.getKey();
            char[] AAComb = entry.getValue();
            System.out.println(scan + " : " + new String(AAComb));
        }

    }


    /**
     * Get one psm from psmAlignedList whose variation positions haven't appear in
     * variationPosChecked set.  Add its variation positions to the set.
     * @param variationPosChecked  a Set containing all checked variation positions
     * @param psmAlignedList    a list of psmAligned
     * @return a PSM whose variation positions haven't appear in
     *      variationPosChecked set.  variationPosChecked set is updated by adding
     *      those new variation positions.
     */
    private PSMAligned getOnePSMWithNewVariationPosition(Set<Integer> variationPosChecked, List<PSMAligned> psmAlignedList) {
        /* Pick one psm with new list of variation positions */
        boolean containNewVaritions = false;
        PSMAligned psmAlignedChosen = null;

        for (PSMAligned psmAligned : psmAlignedList) {
            List<Integer> posOfVariations = psmAligned.getPositionOfVariations();
            for (int i : posOfVariations) {
                if (!variationPosChecked.contains(i)) {
                    variationPosChecked.add(i);
                    containNewVaritions = true;
                }
            }
            if (containNewVaritions == true) {
                psmAlignedChosen = psmAligned;
                break;
            }
        }
        return psmAlignedChosen;
    }


    public List<HashMap<List<Integer>, List<String>>> validateMutations(TemplateHooked templateHooked,
                                                                        HashMap<String, PSMAligned> scanPSMMap) {
        int maxMutationNum = findMaxMutationNumPerPSM(templateHooked);
        System.out.println("Max mutation num: " + maxMutationNum);
        //Store mutationsOnTemplate according to the mutationNum separately.
        List<HashMap<List<Integer>, List<String>>> mutationsOnTemplateList = new ArrayList<>();

        for (int mutationNum = maxMutationNum; mutationNum > 0; mutationNum--) {
            HashMap<List<Integer>, List<String>> mutationsOnTemplate = new HashMap<>();
            List<Integer> posWithMaxMutationNumList = extractPositionWithMutationNum(templateHooked, mutationNum);

            //int pos = posWithMaxMutationNumList.get(0);
            for (int pos : posWithMaxMutationNumList) {
                //System.out.println("pos: " + pos);

                Set<Integer> variationPosChecked = new HashSet<>();
                while (true) {
                    List<PSMAligned> psmAlignedList = getListOfPSMAlignedGivenPosAndMutatedNum(templateHooked, pos,
                                                                                                mutationNum);
                    if (psmAlignedList.size() == 0) {
                        break;
                    }

                    PSMAligned psmAligned = getOnePSMWithNewVariationPosition(variationPosChecked, psmAlignedList);
                    if (psmAligned == null) {
                        break;
                    }


                    HashMap<List<Integer>, List<String>> partMutatedAAsOnTemplate = validateOnePSMAligned(templateHooked,
                            psmAligned, scanPSMMap);
                    if (partMutatedAAsOnTemplate != null) {
                        for (List<Integer> posList : partMutatedAAsOnTemplate.keySet()) {
                            mutationsOnTemplate.put(posList, partMutatedAAsOnTemplate.get(posList));
                        }
                    }
                }
            }
            mutationsOnTemplateList.add(mutationsOnTemplate);
        }
        return mutationsOnTemplateList;
    }


}

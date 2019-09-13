package com.uwaterloo.SignificantMutationsFinder;

import com.uwaterloo.ScanTemplateMapper.PSMAligned;
import com.uwaterloo.ScanTemplateMapper.TemplateHooked;

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
    public int findMaxMutationNumPerPSM(TemplateHooked templateHooked) {
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

    /* Extract positions where containing psm with mutationNum of mutations */
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
     * @param significantThreshold
     * @return <list<position>, mutationPattern>
     */
    private HashMap<List<Integer>, List<MutationsPattern>> validateOnePSMAligned(TemplateHooked templateHooked, PSMAligned psmAligned,
                                                                                 HashMap<String, PSMAligned> scanPSMMap,
                                                                                 double significantThreshold) {
        /* validate the psms having same list of variation positions */
        HashMap<String, MutationsPattern> extractedScanAAPatternMap = buildScanAAPatternMapForOnePSM(templateHooked,
                psmAligned, scanPSMMap);
        //printExtractedPSMs(extractedPSMsMap);
        String templatePattern = extractedScanAAPatternMap.get("Template").getAAs();
        extractedScanAAPatternMap.remove("Template");
        /* Use frequency as pattern score
        double threshold = 0.1;
        HashMap<String, Integer> patternFreqTable = countPatternFreqTable(extractedScanAAPatternMap);
        List<String> significantMutatedPatternsOnTemplate = getSignificantFreqMutatedList(patternFreqTable,
                threshold);
        */

        if (extractedScanAAPatternMap.size() == 0) {
            return null;
        }
        /* Use ion score as pattern score and compute significant using baysian probability */
        HashMap<MutationsPattern, Integer> patternScoreTable = countPatternScoreTable(extractedScanAAPatternMap);
        /* Debug
        for (MutationsPattern pattern : patternScoreTable.keySet()) {
            System.out.println(pattern.toString() + " " + patternScoreTable.get(pattern));
        }
        */

        List<MutationsPattern> significantMutatedPatternsOnTemplate = getSignificantScoreMutatedList(patternScoreTable,
                                                                                            significantThreshold);

        /* Remove scans and PSMAlign containing unsignificant mutation patterns */
        HashSet<MutationsPattern> significantPatterns = new HashSet<>(significantMutatedPatternsOnTemplate);
        removeUnsignificantMutatedScans(templateHooked, psmAligned.getStart(),
                extractedScanAAPatternMap, significantPatterns);

        /* If only one pattern existing in significantMutatedPatternsOnTemplate and the pattern is the same
            with template, there is no need to store this pattern
         */
        if ((significantMutatedPatternsOnTemplate.size() == 1) &&
                (significantMutatedPatternsOnTemplate.get(0).getAAs().equals(templatePattern))) {
            return null;
        }

        /* If significant patterns contain one common AA which appears in template,
            this mutated pattern is not of length mutatedNum. return null
            mutatedAAsOnTemplate. The mutated position will be considered in
             later cycle with less mutatedNum.
         */
        /* The shrink caused some problems later TODO check whether shrink cause problem if not assemble mutationPattern*/
        boolean isShrinked = checkPatternShrink(templatePattern, significantPatterns);
        if (isShrinked) {
            return null;
        }


        /* The following code are used for method using frequency
        // The frequency of significant patterns are extracted from PatternFreqTable and attached to pattern.
        List<String> sigMutationsWithFreqOnTemplate = attachFreqToSignificantPattern(significantMutatedPatternsOnTemplate,
                patternFreqTable);

        HashMap<List<Integer>, List<String>> mutatedAAsOnTemplate = new HashMap<>();
        int start = psmAligned.getStart();
        List<Integer> posList = new ArrayList<>();
        for (int i : psmAligned.getPositionOfVariations()) {
            posList.add(i + start);
        }

        mutatedAAsOnTemplate.put(posList, sigMutationsWithFreqOnTemplate);
//      }
*/
        HashMap<List<Integer>, List<MutationsPattern>> mutatedAAsOnTemplate = new HashMap<>();

        List<Integer> posList = significantMutatedPatternsOnTemplate.get(0).getPosList();
        //The posList in significantMutatedPatternsOnTemplate should be the same.
        mutatedAAsOnTemplate.put(posList, significantMutatedPatternsOnTemplate);

        return mutatedAAsOnTemplate;
    }




    /* return those patterns whose freq greater than threshold in the top 3 patterns */
    private List<String> getSignificantFreqMutatedList(HashMap<String, Integer> patternFreqTable,
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
     * Compute significant mutationPatterns according to top 3 pattern scores.
     * @param patternScoreTable
     * @param threshold  The ratio threshold as significant
     * @return a list of significant mutationPatterns
     */
    private List<MutationsPattern> getSignificantScoreMutatedList(HashMap<MutationsPattern, Integer> patternScoreTable,
                                                                  double threshold) {
        List<MutationsPattern> significantMutatedPattern = new ArrayList<>();

        /* Get the sum of scores */
        double total = 0.0;
        for (int value : patternScoreTable.values()) {
            total += value;
        }

        /* Get the 3 patterns with top 3 frequency */
        MutationsPattern[] top3Patterns = getPatternsWithTop3Score(patternScoreTable);

        /* Print significant freq*/
        for (MutationsPattern pattern : top3Patterns) {
            //Skip the case that there is less than 3 patterns in total
            if (pattern == null) {
                break;
            }
            int score = patternScoreTable.get(pattern);
            double ratio = score / total;

            if (ratio > threshold) {
                significantMutatedPattern.add(pattern);
                //Debug
                //System.out.printf("Debug %s %.4f\n", pattern, score / total);
            } else {
                //Debug
                //System.out.printf("Debug unsignificant %s %.4f\n", pattern, score / total);
            }
        }

        return significantMutatedPattern;
    }



    /* There is a bug that the unsignificant pattern might in dblist other than spiderList.
    * TODO */
    private void removeUnsignificantMutatedScans(TemplateHooked templateHooked, int pos,
                                                 HashMap<String, MutationsPattern> scanAAPatternMap,
                                                 Set<MutationsPattern> significantPatterns) {
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
            } else {
                //Debug
                //System.out.println("Debug spider remove scan " + psmAligned.getScan());
            }
        }
        templateHooked.getSpiderList().set(pos, newSpiderList);

        /* remove from dbList of position pos on templateHooked .   Consider later. Currently, for
        current mutationNum only remove insignificant spider. the db list will be kept for consideration of other
           smaller mutationNum
        ArrayList<PSMAligned> dbList = templateHooked.getDbList().get(pos);
        int dbListSize = dbList.size();
        ArrayList<PSMAligned> newDbList = new ArrayList<>();
        for (int index = 0; index < dbListSize; index++) {
            PSMAligned psmAligned = dbList.get(index);
            if (psmAligned != null && !unsignificantScanSet.contains(psmAligned.getScan())) {
                newDbList.add(psmAligned);
                dbList.set(index, null);
            } else if (psmAligned != null && unsignificantScanSet.contains(psmAligned.getScan())) {
                //Debug
                System.out.println("db remove scan " + psmAligned.getScan());
            }
        }
        templateHooked.getDbList().set(pos, newSpiderList);
*/

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


    /* If all ith position of patterns equals to the AA on templatePattern, this position is
            not a mutation. The pattern should be shrink to a pattern with less mutationNum.
         */
    private boolean checkPatternShrink(String templatePattern, HashSet<MutationsPattern> significantPatterns) {
        int size = templatePattern.length();
        for (int i = 0; i < size; i++) {
            char AA = templatePattern.charAt(i);
            boolean equalToTemplate = true;
            for (MutationsPattern mutationsPattern : significantPatterns) {
                String pattern = mutationsPattern.getAAs();
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

    private List<String> attachFreqToSignificantPattern(List<String> significantMutatedPatternsOnTemplate,
                                                HashMap<String, Integer> patternFreqTable) {
        List<String> patternWithFreqList = new ArrayList<>();
        for (String pattern : significantMutatedPatternsOnTemplate) {
            String patternFreq = pattern + "_" + patternFreqTable.get(pattern);
            patternWithFreqList.add(patternFreq);
        }
        return patternWithFreqList;

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
    private HashMap<String, MutationsPattern> buildScanAAPatternMapForOnePSM(TemplateHooked templateHooked, PSMAligned psmAligned,
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

        HashMap<String, MutationsPattern> scanAAPatternMap = new HashMap<>();
        /* Add the template pattern in scanAAPatternMap */
        scanAAPatternMap.put("Template", new MutationsPattern(posList, new String(templateAAPattern), 0, 0));

        for (String scan : scanSet) {
            PSMAligned correspondingPsmAligned = scanPSMMap.get(scan);

            MutationsPattern AAPattern = extractAAPattern(correspondingPsmAligned, posList, templateAAPattern);
            //If correspondingPSMAligned has shorter posList, null will be returned.
            if (AAPattern != null) {
                scanAAPatternMap.put(scan, AAPattern);

                /* //DEBUG
                if (scan.equals("F3:9267")) {
                    System.out.println("Debug: " + AAPattern + " Peptide: " +
                    correspondingPsmAligned.getPeptide() + " pos List" + posList.toString());

                }

                if (AAPattern.equals("SY")) {
                    //Got result, this happens only on template
                    System.out.println("Debug: SY scan: " + scan);
                }
                */
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
     * Compute the sum of ions scores of each pattern
     * @param scanAAPatternMap  a map between <scan, MutationsPattren>
     * @return a lookup table <pattern, the sum of scores of this pattern>
     */
    private HashMap<MutationsPattern, Integer> countPatternScoreTable(HashMap<String, MutationsPattern> scanAAPatternMap) {
        //TODO deal with patterns with wildcards.  Maybe add all possible head AAs for the pattern to form multiple new pattern
        //Compute the sum of scores for each AA patterns
        HashMap<String, MutationsPattern> sumMutationPatternScores = new HashMap<>();
        for (MutationsPattern mutationsPattern : scanAAPatternMap.values()) {
            String AAs = mutationsPattern.getAAs();
            if (sumMutationPatternScores.containsKey(AAs)) {
                sumMutationPatternScores.get(AAs).setFreq(sumMutationPatternScores.get(AAs).getFreq() + mutationsPattern.getFreq());
                sumMutationPatternScores.get(AAs).setScore(sumMutationPatternScores.get(AAs).getScore() + mutationsPattern.getScore());
            } else {
                MutationsPattern pattern = new MutationsPattern(mutationsPattern.getPosList(),
                        mutationsPattern.getAAs(), mutationsPattern.getFreq(), mutationsPattern.getScore());
                sumMutationPatternScores.put(AAs, pattern);
            }
        }

        //Generate a the hashMap of mutationPattern and its summed score.
        HashMap<MutationsPattern, Integer> patternScoreTable = new HashMap<>();
        for (MutationsPattern mutationsPattern : sumMutationPatternScores.values()) {
            patternScoreTable.put(mutationsPattern, mutationsPattern.getScore());
        }

        return patternScoreTable;
    }

     /**
     * Extract the amino acids from positions in the posList from one PSM
     * (ins) won't cause problem.  (del) will cause problem, need to rethink.
      * If this psmAligned's peptide does not contains full list of posList,
      * null will be returned. This psmAligned will be considered in another
      * cycle with less mutationNum. If the size of posList is less than
      * the variationNum of psmAligned, its pattern should not be extracted,
      * because it will be the subset of longer mutation pattern already extracted.
      * TODO extract all scan covering all posList. If mutations equals or the header or end is one less, them should be considered.
      *
     * @param psmAligned
     * @param posList
     * @param templateAAPattern
     * @return  A mutationPattern, containing the AA combination of posList of this psm, the score of the
      * mutationPattern are the sum of ionsScore. freq = 1
     */
    private MutationsPattern extractAAPattern(PSMAligned psmAligned, List<Integer> posList,
                                      char[] templateAAPattern) {
        if ((psmAligned.getPositionOfVariations()) != null && (posList.size() < psmAligned.getPositionOfVariations().size())) {
            return null;
        }
        int start = psmAligned.getStart();
        int end = psmAligned.getEnd() - start;

        /*
        if (psmAligned.getPeptide().contains("del")) {
            System.err.println("Error in extractAAComb() : scan " + psmAligned.getScan() +
                    " contains del " + psmAligned.getPeptide());
            //TODO: to deal with del part
            return null;
        } */

        int size = posList.size();
        char[] AAComb = new char[size];
        int ionScoreSum = 0;
        for (int i = 0; i < posList.size(); i++) {
            int pos = posList.get(i) - start;
            if ((pos < 0) || (pos > end)) {
                return null;
                //AAComb[i] = templateAAPattern[i];
                //AAComb[i] = '*'; //TODO reconsider
            } else {
                AAComb[i] = psmAligned.getAAs()[pos];
                //If this scan does not have ions information, return null. It is probably due to PEAKS bug.
                if (psmAligned.getIonScores() == null) {
                    return null;
                }
                ionScoreSum += psmAligned.getIonScores()[pos];
            }

        }
        MutationsPattern AAPattern = new MutationsPattern(posList, new String(AAComb), 1, ionScoreSum);

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


    /**
     * Get mutationPatterns with top 3 score
     * @param patternScoreTable
     * @return No more than 3 patterns with top frequency
     */
    private MutationsPattern[] getPatternsWithTop3Score(HashMap<MutationsPattern, Integer> patternScoreTable) {
        MutationsPattern[] top3Patterns = new MutationsPattern[3];
        int[] top3Scores = {-1, -1, -1};
        for (Map.Entry<MutationsPattern, Integer> entry : patternScoreTable.entrySet()) {
            int score = entry.getValue();
            if (score >= top3Scores[0]) {
                top3Scores[2] = top3Scores[1];
                top3Scores[1] = top3Scores[0];
                top3Scores[0] = score;
                top3Patterns[2] = top3Patterns[1];
                top3Patterns[1] = top3Patterns[0];
                top3Patterns[0] = entry.getKey();
            } else if (score >= top3Scores[1]) {
                top3Scores[2] = top3Scores[1];
                top3Patterns[2] = top3Patterns[1];
                top3Scores[1] = score;
                top3Patterns[1] = entry.getKey();
            } else if (score >= top3Scores[0]) {
                top3Scores[2] = score;
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
    private PSMAligned getOnePSMWithNewVariationPosition(Set<String> variationPosChecked, List<PSMAligned> psmAlignedList) {
        /* Pick one psm with new list of variation positions */
        for (PSMAligned psmAligned : psmAlignedList) {
            int start = psmAligned.getStart();
            List<Integer> posOfVariations = psmAligned.getPositionOfVariations();
            List<Integer> posList = new ArrayList<>();
            for (int i = 0; i < posOfVariations.size(); i++) {
                posList.add(posOfVariations.get(i) + start);
            }
            String posListString = posList.toString();
            if (variationPosChecked.contains(posListString)) {
                continue;
            } else {
                variationPosChecked.add(posListString);
                return psmAligned;
            }
        }
        return null;
    }

    private List<Integer> extractPositionWithMutation(TemplateHooked templateHooked) {
        List<Integer> posList = new ArrayList<>();
        ArrayList<ArrayList<PSMAligned>> listOfSpiderList = templateHooked.getSpiderList();
        int size = listOfSpiderList.size();
        int prePos = -1;
        for (int pos = 0; pos < size; pos++) {
            ArrayList<PSMAligned> psmAlignedList = listOfSpiderList.get(pos);
            for (PSMAligned psmAligned : psmAlignedList) {
                if (psmAligned.getPositionOfVariations() == null) {
                    continue;
                }
                if (pos != prePos) {
                    posList.add(pos);
                    prePos = pos;
                }
            }
        }
        return posList;

    }

    /**
     * For the given templateHooked, find a list of significant mutation patterns.
     * @param templateHooked
     * @param scanPSMMap
     * @param significantThreshold
     * @return a list of mutation patterns on the template, denoted by HashMap<posList, list of significant patterns.
     */
    public List<HashMap<List<Integer>, List<MutationsPattern>>> findSignificantMutations(TemplateHooked templateHooked,
                                                                                         HashMap<String, PSMAligned> scanPSMMap,
                                                                                         double significantThreshold) {
        int maxMutationNum = findMaxMutationNumPerPSM(templateHooked);
        System.out.println("Max mutation num: " + maxMutationNum);
        //Store mutationsOnTemplate according to the mutationNum separately.
        List<HashMap<List<Integer>, List<MutationsPattern>>> mutationsOnTemplateList = new ArrayList<>();
        //Try new method to extract all pos containing mutations
        //List<Integer> posWithMaxMutationNumList = extractPositionWithMutation(templateHooked);

        //for (int mutationNum = 0; mutationNum < maxMutationNum; mutationNum++) {
        for (int mutationNum = maxMutationNum; mutationNum > 0; mutationNum--) {
            HashMap<List<Integer>, List<MutationsPattern>> mutationsOnTemplate = new HashMap<>();
            //Old way extract only those spectrum with given number of mutations.
            List<Integer> posWithMaxMutationNumList = extractPositionWithMutationNum(templateHooked, mutationNum);

            //int pos = posWithMaxMutationNumList.get(0);
            Set<String> variationPosChecked = new HashSet<>();
            for (int pos : posWithMaxMutationNumList) {
                //System.out.println("pos: " + pos);
                while (true) {
                    List<PSMAligned> psmAlignedList = getListOfPSMAlignedGivenPosAndMutatedNum(templateHooked, pos,
                            mutationNum);
                    //List<PSMAligned> psmAlignedList = templateHooked.getSpiderList().get(pos);
                    if (psmAlignedList.size() == 0) {
                        break;
                    }

                    PSMAligned psmAligned = getOnePSMWithNewVariationPosition(variationPosChecked, psmAlignedList);
                    if (psmAligned == null) {
                        break;
                    }

                    HashMap<List<Integer>, List<MutationsPattern>> partMutatedAAsOnTemplate = validateOnePSMAligned(templateHooked,
                            psmAligned, scanPSMMap, significantThreshold);
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

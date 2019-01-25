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
        ArrayList<LinkedList<PSMAligned>> listOfSpiderList = templateHooked.getSpiderList();
        for (LinkedList<PSMAligned> psmAlignedList : listOfSpiderList) {
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
        ArrayList<LinkedList<PSMAligned>> listOfSpiderList = templateHooked.getSpiderList();
        int size = listOfSpiderList.size();
        for (int pos = 0; pos < size; pos++) {
            LinkedList<PSMAligned> psmAlignedList = listOfSpiderList.get(pos);
            for (PSMAligned psmAligned : psmAlignedList) {
                if (psmAligned.getPositionOfVariations() == null) {
                    continue;
                }
                if (psmAligned.getPositionOfVariations().size() == mutationNum) {
                    posWithMutationNumList.add(pos);
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
    public PSMAligned getFirstPSMAlignedGivenPosAndMutatedNum(TemplateHooked templateHooked, int pos, int mutatedNum) {
        List<PSMAligned> psmAlignedList = templateHooked.getSpiderList().get(pos);
        for (PSMAligned psmAligned : psmAlignedList) {
            if (psmAligned.getPositionOfVariations().size() == mutatedNum) {
                return psmAligned;
            }
        }
        return null;
    }

    /**
     * Decide whether the PSM with mutation should be kept or discarded on the templateHooked.
     * Get all scans crossing the mutation position list.
     * @param templateHooked
     * @param psmAligned
     * @param scanPSMMap
     */
    private HashMap<String, char[]> validateOnePSMAligned(TemplateHooked templateHooked, PSMAligned psmAligned,
                                                    HashMap<String, PSMAligned> scanPSMMap) {
        int start = psmAligned.getStart();
        List<Integer> posList = new ArrayList<>();
        for (int i : psmAligned.getPositionOfVariations()) {
            posList.add(start + i);
        }

        int mutationNum = posList.size();
        char[] templateAAComb = new char[mutationNum];
        Set<String> scanSet = new HashSet<>();

        for (int i = 0; i < mutationNum; i++) {
            int mutationPos = posList.get(i);
            scanSet.addAll(templateHooked.getMappedScanList().get(mutationPos));
            templateAAComb[i] = templateHooked.getSeq()[mutationPos];
        }

        System.out.println("Template AA comb: " + new String(templateAAComb));

        HashMap<String, char[]> aminoCombination = new HashMap<>();
        for (String scan : scanSet) {
            PSMAligned correspondingPsmAligned = scanPSMMap.get(scan);
            char[] AAComb = extractAAComb(correspondingPsmAligned, posList, templateAAComb);
            aminoCombination.put(scan, AAComb);
        }

        return aminoCombination;


    }

    /**
     * Extract the amino acids in the posList.
     * (ins) won't cause problem.  (del) will cause problem, need to rethink.
     * @param psmAligned
     * @param posList
     * @param templateAAComb
     * @return
     */
    private char[] extractAAComb(PSMAligned psmAligned, List<Integer> posList,
                                      char[] templateAAComb) {
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
                AAComb[i] = templateAAComb[i];
            } else {
                AAComb[i] = psmAligned.getAAs()[pos];
            }

        }
        return AAComb;
    }

    private void printExtractedPSMs(HashMap<String,char[]> extractedPSMsMap) {
        for (Map.Entry<String, char[]> entry : extractedPSMsMap.entrySet()) {
            String scan = entry.getKey();
            char[] AAComb = entry.getValue();
            System.out.println(scan + " : " + new String(AAComb));
        }

    }

    public void validateMutations(TemplateHooked templateHooked, HashMap<String, PSMAligned> scanPSMMap) {
        int maxMutationNum = findMaxMutationNumPerPSM(templateHooked);
        System.out.println("Max mutation num: " + maxMutationNum);
        int mutationNum = 4;
        List<Integer> posWithMaxMutationNumList = extractPositionWithMutationNum(templateHooked, mutationNum);

        int pos = posWithMaxMutationNumList.get(0);
        PSMAligned psmAligned = getFirstPSMAlignedGivenPosAndMutatedNum(templateHooked, pos, mutationNum);

        System.out.println("The first psm: " + psmAligned.toString());

        HashMap<String, char[]> extractedPSMsMap = validateOnePSMAligned(templateHooked, psmAligned, scanPSMMap);
        printExtractedPSMs(extractedPSMsMap);
    }




}

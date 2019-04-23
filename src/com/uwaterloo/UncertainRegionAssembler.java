package com.uwaterloo;

import java.util.*;

public class UncertainRegionAssembler {

    public UncertainRegionAssembler() {

    }

    //Extract the Denovo only covering pos at templateHooked.
    private Set<DenovoAligned> extractDnToRightAtPos(TemplateHooked templateHooked, int pos) {
        List<DenovoAligned> denovoAlignedList = templateHooked.getDnToRightList().get(pos);
        Set<DenovoAligned> dnAlignSet = new HashSet<>();
        dnAlignSet.addAll(denovoAlignedList);

        return dnAlignSet;
    }

    //Extract the Denovo only covering pos at templateHooked.
    private Set<DenovoAligned> extractDnToLeftAtPos(TemplateHooked templateHooked, int pos) {
        List<DenovoAligned> denovoAlignedList = templateHooked.getDnToLeftList().get(pos);
        Set<DenovoAligned> dnAlignSet = new HashSet<>();
        dnAlignSet.addAll(denovoAlignedList);

        return dnAlignSet;
    }

    //Extract the db covering at pos at templateHooked
    private Set<PSMAligned> extractDbAtPos(TemplateHooked templateHooked, int pos,
                                           HashMap<String, PSMAligned> scanPSMMap) {
        Set<PSMAligned> dbSet = new HashSet<>();
        List<String> scanList = templateHooked.getMappedScanList().get(pos);
        for (String scan : scanList) {
            PSMAligned psmAligned = scanPSMMap.get(scan);
            //Skip spider psms
            if (psmAligned.getPositionOfVariations() != null) {
                continue;
            }
            dbSet.add(psmAligned);
        }
        return dbSet;
    }

    /**
     * Identify positions whose number of db / number of dn < dbDnRatioThresh to form a list of
     * uncertain regions. Each region is one position.
     * @param templateHooked Template with scanlist, dnToLeftList, dnToRightList
     * @param scanPSMMap <scan, PSMAlign> hashMap
     * @param dbDnRatioThresh The threshold of ratio of number of db / number of dn.
     * @return
     */
    private List<UncertainRegion> identifyUncertainRegions(TemplateHooked templateHooked,
                                          HashMap<String, PSMAligned> scanPSMMap,
                                          float dbDnRatioThresh) {
        int templateLength = templateHooked.getSeq().length;
        List<UncertainRegion> uncertainRegions = new ArrayList<>();
        for (int pos = 0; pos < templateLength; pos++) {
            List<DenovoAligned> denovoAlignedToRightList = templateHooked.getDnToRightList().get(pos);
            List<DenovoAligned> denovoAlignedToLeftList = templateHooked.getDnToLeftList().get(pos);
            int psmSize = templateHooked.getMappedScanList().get(pos).size();
            int dnSize = denovoAlignedToRightList.size() + denovoAlignedToLeftList.size();

            float ratio = (psmSize + 0.0f) / dnSize;
            if (ratio < dbDnRatioThresh) {
                Set<DenovoAligned> dnToRightSet = extractDnToRightAtPos(templateHooked, pos);
                Set<DenovoAligned> dnToLeftSet = extractDnToLeftAtPos(templateHooked, pos);
                Set<PSMAligned> dbSet = extractDbAtPos(templateHooked, pos, scanPSMMap);
                uncertainRegions.add(new UncertainRegion(pos, pos, dnToRightSet, dnToLeftSet, dbSet));
            }
        }
        return uncertainRegions;
    }

    /**
     * Merge adjacent uncertain positions into uncertain regions. The corresponding
     * dnToRightSet, dnToLeftSet, pmsSet are union.
     * @param uncertainRegionList
     * @param adjacentThresh The distance threshold between the end
     *                      of previous region and the start of the next region.
     * @return a list of merged uncertain regions
     */
    private List<UncertainRegion> mergeAdjacentRegion(List<UncertainRegion> uncertainRegionList,
                                                      int adjacentThresh) {
        List<UncertainRegion> mergedRegions = new ArrayList<>();
        int lastMergedIndex = 0;
        int preEndPos = 0;

        //Add the first uncertain region to merged regions
        mergedRegions.add(uncertainRegionList.get(0));

        for (int index = 1; index < uncertainRegionList.size(); index++) {
            UncertainRegion region = uncertainRegionList.get(index);
            int pos = region.getStartPos();
            if ((pos - preEndPos) <= adjacentThresh) {
                //Merge current region to previous merged region
                mergedRegions.get(lastMergedIndex).setEndPos(region.getEndPos());
                mergedRegions.get(lastMergedIndex).getDnAlignToRightSet().addAll(region.getDnAlignToRightSet());
                mergedRegions.get(lastMergedIndex).getDnAlignToLeftSet().addAll(region.getDnAlignToLeftSet());
                mergedRegions.get(lastMergedIndex).getPsmAlignedSet().addAll(region.getPsmAlignedSet());
            } else {
                mergedRegions.add(region);
                lastMergedIndex++;
            }
            preEndPos = region.getEndPos();
        }

       // System.out.println("merged region number: " + mergedRegions.size());
        return mergedRegions;
    }

    /**
     * Assemble the toRightDnSet, toLeftDnSet and psmSet. Find one contig for each set with highest
     * score.  Among these three, try to find a concensus with highest score. Return one contig if there
     * are adjacent, otherwise a list of contig.
     * @param regionToAssemble
     * @param scanDnMap
     * @return
     */
    private List<Contig> assembleOneRegion(UncertainRegion regionToAssemble, HashMap<String, DenovoOnly> scanDnMap) {
        System.out.println("region starting from " + regionToAssemble.getStartPos() + " to " + regionToAssemble.getEndPos());

        List<DenovoAligned> dnAlignToRightList = new ArrayList<>(regionToAssemble.getDnAlignToRightSet());
        List<DenovoAligned> dnAlignToLeftList = new ArrayList<>(regionToAssemble.getDnAlignToLeftSet());
        List<PSMAligned> psmAlignedList = new ArrayList<>(regionToAssemble.getPsmAlignedSet());

        //Sort dnToRightList by tStart
        Collections.sort(dnAlignToRightList, DenovoAligned.cmpTStart());
        //Sort dnToLeftList by reverse tEnd
        Collections.sort(dnAlignToLeftList, DenovoAligned.cmpReverseTEnd());
        //Sort psmList by tStart
        Collections.sort(psmAlignedList, PSMAligned.cmpStart());

        List<Contig> assembledContigs = new ArrayList<>();

/*
        System.out.println("dn to right");
        for (DenovoAligned dnAligned : dnAlignToRightList) {
            System.out.print(dnAligned.gettStart() + " " + new String(scanDnMap.get(dnAligned.getDnScan()).getAAs()) +
                    " " + dnAligned.getDnScan());
            for (short score : scanDnMap.get(dnAligned.getDnScan()).getConfScores()) {
                System.out.print(" " + score);
            }
            System.out.println();
        }
*/

        //Assemble the dnToRight into contigs, keep only the one with highest score
        List<Contig> toRightContigs = assembleDnToRightAlign(dnAlignToRightList, scanDnMap);
        if (toRightContigs != null) {
            //System.out.println("Contigs to right");
            /*
            for (Contig contig : toRightContigs) {
                System.out.println(contig.toString());
            }*/
            //Print only the contig with best score
            //System.out.println(toRightContigs.get(0).toString());
            assembledContigs.add(toRightContigs.get(0));
        }

/*
        System.out.println("dn to left");
        for (DenovoAligned dnAligned : dnAlignToLeftList) {
            System.out.println(dnAligned.gettStart() + " " + new String(scanDnMap.get(dnAligned.getDnScan()).getAAs()) +
                    " " + dnAligned.getDnScan());
            for (short score : scanDnMap.get(dnAligned.getDnScan()).getConfScores()) {
                System.out.print(" " + score);
            }
            System.out.println();
        }
*/

        //Assemble the dnToLeft into contigs, keep only the one with highest score
        List<Contig> toLeftContigs = assembleDnToLeftAlign(dnAlignToLeftList, scanDnMap);
        if (toLeftContigs != null) {
            //System.out.println("Contigs to Left");

            /*
            for (Contig contig : toLeftContigs) {
                System.out.println(contig.toString());
            }
            */
            //System.out.println(toLeftContigs.get(0).toString());
            assembledContigs.add(toLeftContigs.get(0));
        }

        //If toRightContig and toLeftContig contains overlap, merge them to a new contig
        if (toRightContigs != null && toLeftContigs != null) {
            List<Contig> bridgeContigs = findLeftRightOverlap(toRightContigs, toLeftContigs);
        }


        /*
        for (PSMAligned psm : psmAlignedList) {
            System.out.print(psm.getStart() + " " + new String(psm.getAAs()));
            for (int score : psm.getIonScores()) {
                System.out.print(" " + score);
            }
            System.out.println();

        }
        */


        //Assemble the PSMs into one contig, use only PSM with good fragmentation
        if (psmAlignedList != null) {
            //Use good quality peptide to form a contig
            List<Contig> psmContigs = assembleConfidentPSM(psmAlignedList);
            if (psmContigs.size() > 0) {
                /*System.out.println("db");
                for (Contig psmContig : psmContigs) {
                    System.out.println(psmContig.toString());
                }*/
                assembledContigs.addAll(psmContigs);
            }
        }

        return assembledContigs;

    }

    private List<Contig> assembleConfidentPSM(List<PSMAligned> psmAlignedList) {
        List<Contig> psmContigs = new ArrayList<>();
        int tStart = 5000;
        int tEnd = 0;
        char[] AAs = null;
        int scoreSum = 0;
        int seqLen = 0;

        //psmAlignedList is sorted by psm.tStart
        for (PSMAligned psm : psmAlignedList) {
            //Discard inconfident PSMs.
            if (isConfidentPSM(psm)) {
                if (psm.getStart() < tStart) {
                    tStart = psm.getStart();
                }
                //If the psmAlign is not overlapped with the previous contig, create a contig
                if ((psm.getStart() > tStart) && (psm.getStart() > tEnd)) {
                    psmContigs.add(new Contig(tStart, tEnd, AAs, scoreSum));
                    AAs = psm.getAAs();
                    tStart = psm.getStart();
                    tEnd = psm.getEnd();
                    scoreSum = 0;
                    for (short ionScore : psm.getIonScores()) {
                        scoreSum += ionScore;
                    }
                    seqLen = AAs.length;
                    continue;
                }
                if (psm.getEnd() > tEnd) {
                    tEnd = psm.getEnd();
                    seqLen = tEnd - tStart + 1;
                    System.out.println(seqLen);

                    if (AAs == null) {
                        AAs = psm.getAAs();
                    } else {
                        char[] newAAs = new char[seqLen];
                        for (int i = tStart; i < psm.getStart(); i++) {
                            newAAs[i - tStart] = AAs[i - tStart];
                        }
                        for (int i = psm.getStart(); i <= tEnd; i++) {
                            newAAs[i -tStart] = psm.getAAs()[i - psm.getStart()];
                        }
                        AAs = newAAs;
                    }
                }

                for (short ionScore : psm.getIonScores()) {
                    scoreSum += ionScore;
                }
            }
        }
        if (seqLen > 0) {
            psmContigs.add(new Contig(tStart, tEnd, AAs, scoreSum));
        }
        return psmContigs;
    }

    /* if the psm contain ionScore smaller than 33, there is no fragment for four AAs, it is a inconfident one.*/
    private boolean isConfidentPSM(PSMAligned psm) {
        int scoreThresh = 30;
        for (short ionScore : psm.getIonScores()) {
            if (ionScore < scoreThresh) {
                return false;
            }
        }
        return true;
    }

    private List<Contig> findLeftRightOverlap(List<Contig> toRightContigs, List<Contig> toLeftContigs) {
        List<Contig> bridgedContigs = new ArrayList<>();
        for (Contig toRightContig : toRightContigs) {
            char[] seq1 = toRightContig.getAAs();
            for (Contig toLeftContig : toLeftContigs) {
                char[] seq2 = toLeftContig.getAAs();
                int overlapIndex = findOverlapIndex(seq1, seq2);
                if (overlapIndex < seq1.length) {
                    System.out.println("bridge overlap" + new String(seq1) + " " + new String(seq2) + " " + overlapIndex) ;
                }
            }
        }
        return bridgedContigs;
    }

    /**
     * Find the index where the head of seq2 can overlap with the end of seq1
     * If no overlap, return seq1.length
     * @param seq1
     * @param seq2
     * @return
     */
    private int findOverlapIndex(char[] seq1, char[] seq2) {
        int i = 0;
        int j = 0;
        int start = -1;

        while (true) {
            while ((i < seq1.length) && (seq1[i] != seq2[0])) {
                i++;
            }
            if (i == seq1.length) {
                return i;
            }

            start = i;
            while ((i < seq1.length) && (j < seq2.length) && (seq1[i] == seq2[j])) {
                i++;
                j++;
            }
            if ((i == seq1.length) || (j == seq2.length)) {
                return start;
            } else {
                i = start + 1;
                j = 0;
            }
        }
    }


    private List<Contig> assembleDnToRightAlign(List<DenovoAligned> dnAlignToRightList,
                                                       HashMap<String, DenovoOnly> scanDnMap) {
        if (dnAlignToRightList.size() == 0) {
            return null;
        }
        List<Contig> contigs = new ArrayList<>();
        DenovoAligned firstDnAlign = dnAlignToRightList.get(0);
        Contig firstContig = new Contig(firstDnAlign.gettStart(), -1,
                scanDnMap.get(firstDnAlign.getDnScan()).getAAs(), firstDnAlign.getScore());
        contigs.add(firstContig);

        for (int i = 1; i < dnAlignToRightList.size(); i++) {
            DenovoAligned dnAligned = dnAlignToRightList.get(i);
            char[] dnAAs = scanDnMap.get(dnAligned.getDnScan()).getAAs();
            int dnTStart = dnAligned.gettStart();
            int overlappedContigIndex = findContigRightOverlap(contigs, dnAAs, dnTStart);


            if (overlappedContigIndex < contigs.size()) {
                Contig contig = contigs.get(overlappedContigIndex);
                char[] contigAAs = contig.getAAs();
                int startIndex = dnTStart - contig.gettStart();

                //If overlapped, add up their score
                contig.setScore(contig.getScore() + dnAligned.getScore());

                //Extend the AAs to right if there are more AAs
                int endIndex = dnAAs.length - (contigAAs.length - startIndex);
                if (endIndex > 0) {
                    char[] mergedAAs = new char[contigAAs.length + endIndex];
                    for (int j = 0; j < startIndex; j++) {
                        mergedAAs[j] = contigAAs[j];
                    }
                    for (int j = 0; j < dnAAs.length; j++) {
                        mergedAAs[j + startIndex] = dnAAs[j];
                    }
                    contig.setAAs(mergedAAs);
                }
            } else {
                //If not overlapped, add the dnAligned as a new contig
                Contig newContig = new Contig(dnAligned.gettStart(), -1,
                        scanDnMap.get(dnAligned.getDnScan()).getAAs(), dnAligned.getScore());
                contigs.add(newContig);
            }
        }

        //Sort the contigs according to their scores descending.
        Collections.sort(contigs, Contig.cmpReverseScore());

        for (Contig contig : contigs) {
            contig.settEnd(contig.gettStart() + contig.getAAs().length - 1);
        }

        return contigs;
    }

    /* Return the contig index whose right end overlap with the dnAA.
        If the index == contigs.size(), there is no overlapped contig.
     */
    private int findContigRightOverlap(List<Contig> contigs, char[] dnAAs, int dnTStart) {
        int contigIndex = 0;
        while (contigIndex < contigs.size()) {
            Contig contig = contigs.get(contigIndex);

            //If there is no position overlap
            if ((contig.getAAs().length + contig.gettStart()) < dnTStart) {
                contigIndex++;
                continue;
            }

            char[] contigAAs = contig.getAAs();


            boolean isOverlapped = true;
            int startIndex = dnTStart - contig.gettStart();
            int endIndex = ((contigAAs.length - startIndex) > dnAAs.length) ? (startIndex + dnAAs.length) : contigAAs.length;

            for (int j = startIndex; j < endIndex; j++) {
                if (contigAAs[j] != dnAAs[j - startIndex]) {
                    isOverlapped = false;
                    break;
                }
            }
            if (isOverlapped) {
                break;
            } else {
                contigIndex += 1;
            }
        }

        return contigIndex;
    }

    /**
     * assemble the denovo tags which is extending left to template.
     * @param dnAlignToLeftList The dnAlign list sorted by descending dn.tEnd
     * @param scanDnMap
     * @return
     */
    private List<Contig> assembleDnToLeftAlign(List<DenovoAligned> dnAlignToLeftList, HashMap<String, DenovoOnly> scanDnMap) {
        if (dnAlignToLeftList.size() == 0) {
            return null;
        }
        List<Contig> contigs = new ArrayList<>();
        DenovoAligned firstDnAlign = dnAlignToLeftList.get(0);
        Contig firstContig = new Contig(-1, firstDnAlign.gettEnd(),
                scanDnMap.get(firstDnAlign.getDnScan()).getAAs(), firstDnAlign.getScore());
        contigs.add(firstContig);

        for (int i = 1; i < dnAlignToLeftList.size(); i++) {
            DenovoAligned dnAligned = dnAlignToLeftList.get(i);
            char[] dnAAs = scanDnMap.get(dnAligned.getDnScan()).getAAs();
            int dnTEnd = dnAligned.gettEnd();
            int overlappedContigIndex = findContigLeftOverlap(contigs, dnAAs, dnTEnd);

            if (overlappedContigIndex < contigs.size()) {
                Contig contig = contigs.get(overlappedContigIndex);
                char[] contigAAs = contig.getAAs();
                int startIndex = contig.gettEnd() - dnTEnd;

                //If overlapped, add up their score
                contig.setScore(contig.getScore() + dnAligned.getScore());

                //Extend the AAs to right if there are more AAs
                int endIndex = dnAAs.length - (contigAAs.length - startIndex);
                if (endIndex > 0) {
                    int newLength = contigAAs.length + endIndex;
                    char[] mergedAAs = new char[newLength];
                    for (int j = 0; j < startIndex; j++) {
                        mergedAAs[newLength - 1 - j] = contigAAs[contigAAs.length - 1 - j];
                    }
                    for (int j = 0; j < dnAAs.length; j++) {
                        mergedAAs[newLength - 1 - j - startIndex] = dnAAs[dnAAs.length - 1 - j];
                    }
                    contig.setAAs(mergedAAs);
                }
            } else {
                //If not overlapped, add the dnAligned as a new contig
                Contig newContig = new Contig(-1, dnAligned.gettEnd(),
                        scanDnMap.get(dnAligned.getDnScan()).getAAs(), dnAligned.getScore());
                contigs.add(newContig);
            }
        }

        //Sort the contigs according to their scores descending.
        Collections.sort(contigs, Contig.cmpReverseScore());

        for (Contig contig : contigs) {
            contig.settStart(contig.gettEnd() - contig.getAAs().length + 1);
        }

        return contigs;
    }

    /**
     * From right to left, find the index of a contig which is overlapped with dnAA
     * @param contigs
     * @param dnAAs
     * @param dnTEnd
     * @return
     */
    private int findContigLeftOverlap(List<Contig> contigs, char[] dnAAs, int dnTEnd) {
        int contigIndex = 0;
        while (contigIndex < contigs.size()) {
            Contig contig = contigs.get(contigIndex);
            char[] contigAAs = contig.getAAs();
            boolean isOverlapped = true;

            //Start and end index from the end of contig
            int startIndex = contig.gettEnd() - dnTEnd;
            int endIndex = (dnAAs.length < (contigAAs.length - startIndex)) ? dnAAs.length : (contigAAs.length - startIndex);

            for (int j = startIndex; j < endIndex; j++) {
                if (contigAAs[contigAAs.length - 1 - j] != dnAAs[dnAAs.length - 1 - (j - startIndex)]) {
                    isOverlapped = false;
                    break;
                }
            }
            if (isOverlapped) {
                break;
            } else {
                contigIndex += 1;
            }
        }

        return contigIndex;
    }


    /**
     * Apply each contig one by one to the modified template.  For overlapped region, adopt the one with max score.
     * TODO
     * @param templateHooked
     * @param assembledContigs
     */
    private void generateCandidateTemplate(TemplateHooked templateHooked, List<Contig> assembledContigs) {
        List<Contig> mergedContigs = mergeContigs(assembledContigs);
        System.out.println("Merged");
        for (Contig mergedContig : mergedContigs) {
            System.out.println(mergedContig.toString());
        }

        char[] AAs = templateHooked.getSeq().clone();

        for (Contig contig : mergedContigs) {
            int tStart = contig.gettStart();
            int tEnd = contig.gettEnd();

            for (int j = tStart; j <= tEnd; j++) {
                AAs[j] = contig.getAAs()[j - tStart];
            }
        }
        List<char[]> candidateTemplateSeqs = new ArrayList<>();
        candidateTemplateSeqs.add(AAs);
        templateHooked.setModifiedTemplates(candidateTemplateSeqs);
    }

    /* Merge overlapped contigs, choose the AAs with higher score as consensus */
    private List<Contig> mergeContigs(List<Contig> assembledContigs) {
        Collections.sort(assembledContigs, Contig.cmpTStart());
        List<Contig> mergedContigs = new ArrayList<>();
        mergedContigs.add(assembledContigs.get(0));

        for (int i = 1; i < assembledContigs.size(); i++) {
            Contig contig = assembledContigs.get(i);
            boolean overlapped = false;
            for (Contig mergedContig : mergedContigs) {
                if (contig.gettStart() <= mergedContig.gettEnd()) {
                    if (contig.gettEnd() > mergedContig.gettEnd()) {
                        //If contig is longer than mergedContig. Merge parts of higher score to form a new mergedContig.
                        overlapped = true;
                        int middle = mergedContig.gettEnd();
                        int len = mergedContig.getAAs().length + (contig.gettEnd() - middle);
                        char[] AAs = new char[len];
                        for (int k = 0; k < (contig.gettStart() - mergedContig.gettStart()); k++) {
                            AAs[k] = mergedContig.getAAs()[k];
                        }
                        if (contig.getScore() > mergedContig.getScore()) {
                            for (int k = contig.gettStart(); k <= middle; k++) {
                                AAs[k - mergedContig.gettStart()] = contig.getAAs()[k - contig.gettStart()];
                            }
                        } else {
                            for (int k = contig.gettStart(); k <= middle; k++) {
                                AAs[k - mergedContig.gettStart()] = mergedContig.getAAs()[k - mergedContig.gettStart()];
                            }
                        }
                        for (int k = middle + 1; k <= contig.gettEnd(); k++) {
                            AAs[k - mergedContig.gettStart()] = contig.getAAs()[k - contig.gettStart()];
                        }
                        mergedContig.setAAs(AAs);
                    } else {
                        //If contig is nested in the mergedContig, only when it has higher score, it will be substituted.
                        if (contig.getScore() > mergedContig.getScore()) {
                            for (int k = contig.gettStart(); k <= contig.gettEnd(); k++) {
                                mergedContig.getAAs()[k - mergedContig.gettStart()] = contig.getAAs()[k - contig.gettStart()];
                            }
                        }
                    }
                }
            }
            if (!overlapped) {
                mergedContigs.add(contig);
            }
        }
        return mergedContigs;

    }

    public void assembleUncertainRegions(List<TemplateHooked> templateHookedList,
                                         List<HashMap<String, PSMAligned>> listOfscanPSMMap,
                                         HashMap<String, DenovoOnly> scanDnMap,
                                         float dbDnRatioThresh) {
        int adjacentThresh = 10;
        for (int i = 0; i < templateHookedList.size(); i++) {
            TemplateHooked templateHooked = templateHookedList.get(i);
            System.out.println("assemble uncertain region for " + templateHooked.getTemplateAccession());
            HashMap<String, PSMAligned> scanPSMMap = listOfscanPSMMap.get(i);
            //Find uncertain regions where db / dn ratio less than dbDnRatioThresh
            List<UncertainRegion> uncertainRegionList = identifyUncertainRegions(templateHooked, scanPSMMap, dbDnRatioThresh);
            if (uncertainRegionList.size() == 0) {
                continue;
            }
            //Merge adjacent uncertain region whose distance less than adjacentThresh
            List<UncertainRegion> mergedUncertainRegions = mergeAdjacentRegion(uncertainRegionList, adjacentThresh);

            //For each merged uncertain region, generate assembled contigs
            List<Contig> assembledContigs = new ArrayList<>();
            for (UncertainRegion regionToAssemble : mergedUncertainRegions) {
                assembledContigs.addAll(assembleOneRegion(regionToAssemble, scanDnMap));
            }

            System.out.println("Assembled contigs");
            for (Contig assembledContig : assembledContigs) {
                System.out.println(assembledContig.toString());
            }

            //Select contigs with max score and apply to the template to generate candidate template.
            generateCandidateTemplate(templateHooked, assembledContigs);

        }

    }




}

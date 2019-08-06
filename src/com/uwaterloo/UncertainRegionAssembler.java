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

        //Add the first uncertain region to merged regions
        mergedRegions.add(uncertainRegionList.get(0));
        int preEndPos = uncertainRegionList.get(0).getEndPos();

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

        //Sort dnToRightList by tStart, if same tStart, sort by score descending.
        Collections.sort(dnAlignToRightList, DenovoAligned.cmpReverseScore());
        Collections.sort(dnAlignToRightList, DenovoAligned.cmpTStart());

        //Sort dnToRightList by tStart, if same tStart, sort by score descending.
        Collections.sort(dnAlignToLeftList, DenovoAligned.cmpReverseScore());
        Collections.sort(dnAlignToLeftList, DenovoAligned.cmpReverseTEnd());

        //Sort psmList by tStart
        Collections.sort(psmAlignedList, PSMAligned.cmpStart());

        List<Contig> assembledContigs = new ArrayList<>();

/*
        System.out.println("dn to right");
        for (DenovoAligned dnAligned : dnAlignToRightList) {
            System.out.print(dnAligned.gettStart() + " " + new String(scanDnMap.get(dnAligned.getDnScan()).getAAs()) +
                    " " + dnAligned.getScore() + " " + dnAligned.getDnScan());
            for (short score : scanDnMap.get(dnAligned.getDnScan()).getConfScores()) {
                System.out.print(" " + score);
            }
            System.out.println();
        }
/**/

        //Assemble the dnToRight into contigs, keep only the one with highest score
        List<Contig> toRightContigs = assembleDnToRightAlign(dnAlignToRightList, scanDnMap);
        List<Contig> topToRightContigs = null;
        if (toRightContigs != null) {
            System.out.println("Contigs to right");
            /*
            for (Contig contig : toRightContigs) {
                System.out.println(contig.toString());
            }
            /**/

            //Pick the contig with best score for non overlap regions.
            topToRightContigs = pickTopContigs(toRightContigs);
            for (Contig topContig : topToRightContigs) {
                System.out.println(topContig.toString());
                assembledContigs.add(topContig);
            }
        }

/*
        System.out.println("dn to left");
        for (DenovoAligned dnAligned : dnAlignToLeftList) {
            System.out.print(dnAligned.gettStart() + " " + new String(scanDnMap.get(dnAligned.getDnScan()).getAAs()) +
                    " " +  dnAligned.getScore() + " " + dnAligned.getDnScan());
            for (short score : scanDnMap.get(dnAligned.getDnScan()).getConfScores()) {
                System.out.print(" " + score);
            }
            System.out.println();
        }
/**/

        //Assemble the dnToLeft into contigs, keep only the one with highest score
        List<Contig> toLeftContigs = assembleDnToLeftAlign(dnAlignToLeftList, scanDnMap);
        List<Contig> topToLeftContigs = null;
        if (toLeftContigs != null) {
            System.out.println("Contigs to Left");

            /*
            for (Contig contig : toLeftContigs) {
                System.out.println(contig.toString());
            }
            /**/
            //Pick the contig with best score for non overlap regions.
            topToLeftContigs = pickTopContigs(toLeftContigs);
            for (Contig topContig : topToLeftContigs) {
                System.out.println(topContig.toString());
                assembledContigs.add(topContig);
            }
        }
/*
        //If toRightContig and toLeftContig contains overlap, merge them to a new contig
        if (toRightContigs != null && toLeftContigs != null) {
            List<Contig> bridgeContigs = findLeftRightOverlap(topToRightContigs, topToLeftContigs);
            if (bridgeContigs.size() > 0) {
                System.out.println("Top Bridged contigs: ");
                /*
                for (Contig contig : bridgeContigs) {
                    System.out.println(contig.toString());
                }
                */
/*
                List<Contig> topContigs = pickTopContigs(bridgeContigs);
                for (Contig topContig : topContigs) {
                    System.out.println(topContig.toString());
                    //Todo need to add the overlapped contigs in
                    assembledContigs.add(topContig);
                }
            }
        }

*/
        /*
        for (PSMAligned psm : psmAlignedList) {
            System.out.print(psm.getStart() + " " + new String(psm.getAAs()));
            for (int score : psm.getIonScores()) {
                System.out.print(" " + score);
            }
            System.out.println();

        }
        */


        /* Only remain psm segment within range (startOfregion, endOfregion)
        Assemble the PSMs into one contig, use only PSM with good fragmentation
         */
        if (psmAlignedList != null) {
            //Use good quality peptide to form a contig
            int regionStart = regionToAssemble.getStartPos();
            int regionEnd = regionToAssemble.getEndPos();
            List<Contig> psmContigs = assembleConfidentPSM(psmAlignedList, regionStart, regionEnd);
            if (psmContigs.size() > 0) {
                System.out.println("db");
                printContigs(psmContigs);

                assembledContigs.addAll(psmContigs);
            }
        }

        return assembledContigs;

    }

    /**
     * Pick the contigs with highest score for each non-overlapped region.
     * For example, 45-90 one contig with highest score in this region, 102-140 anothr contig with highest score
     * @param contigs  contigs sorted by their scores descending.
     * @return
     */
    List<Contig> pickTopContigs(List<Contig> contigs) {
        List<Contig> topContigs = new ArrayList<>();
        topContigs.add(contigs.get(0));

        for (int i = 1; i < contigs.size(); i++) {
            boolean isOverlapped = false;
            for (Contig topContig : topContigs) {
                if (overlap(contigs.get(i), topContig)) {
                    isOverlapped = true;
                    break;
                }
            }
            if (!isOverlapped) {
                topContigs.add(contigs.get(i));
            }
        }
        return topContigs;
    }


    /**
     * Helper function to judge whether two contigs have overlapped position.
     * @param contig1
     * @param contig2
     * @return
     */
    boolean overlap(Contig contig1, Contig contig2) {
        return ((contig1.gettStart() >= contig2.gettStart()) && (contig1.gettStart() <= contig2.gettEnd())) ||
                ((contig1.gettEnd() >= contig2.gettStart()) && (contig1.gettEnd() <= contig2.gettEnd()));
    }

    private int[] copyConfs(short[] AAConfs) {
        int[] confs = new int[AAConfs.length];
        for (int i = 0; i < AAConfs.length; i++) {
            confs[i] = AAConfs[i];
        }
        return confs;
    }

    private List<Contig> assembleConfidentPSM(List<PSMAligned> psmAlignedList, int regionStart, int regionEnd) {
        List<Contig> psmContigs = new ArrayList<>();
        int tStart = 5000;
        int tEnd = 0;
        char[] AAs = null;
        int[] confs = null;
        int scoreSum = 0;
        int seqLen = 0;

        //psmAlignedList is sorted by psm.tStart
        for (PSMAligned psm : psmAlignedList) {
            //Discard inconfident PSMs.
            if (isConfidentPSM(psm)) {
                confs = copyConfs(psm.getIonScores());
                if (psm.getStart() < tStart) {
                    tStart = psm.getStart();
                }
                //If the psmAlign is not overlapped with the previous contig, create a contig
                if ((psm.getStart() > tStart) && (psm.getStart() > tEnd)) {
                    psmContigs.add(new Contig(tStart, tEnd, AAs, confs, scoreSum));
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
            psmContigs.add(new Contig(tStart, tEnd, AAs, confs, scoreSum));
        }
        return psmContigs;
    }

    /* if the psm contain ionScore smaller than 33, there is no fragment for four AAs, it is a inconfident one.*/
    private boolean isConfidentPSM(PSMAligned psm) {
        int scoreThresh = 30;

        if (psm.getIonScores() == null) {
            return false;
        }
        for (short ionScore : psm.getIonScores()) {
            if (ionScore < scoreThresh) {
                return false;
            }
        }
        return true;
    }

    /**
     * Merge the overlapped contigs
     * @param toRightContigs
     * @param toLeftContigs
     * @return
     */
    private List<Contig> findLeftRightOverlap(List<Contig> toRightContigs, List<Contig> toLeftContigs) {
        int minOverlapLen = 3;
        List<Contig> bridgedContigs = new ArrayList<>();
        for (Contig toRightContig : toRightContigs) {
            char[] seq1 = toRightContig.getAAs();
            for (Contig toLeftContig : toLeftContigs) {
                char[] seq2 = toLeftContig.getAAs();
                int overlapIndex = findOverlapIndex(seq1, seq2);
                if (overlapIndex < seq1.length && (seq1.length - overlapIndex) >= minOverlapLen) {
                    char[] mergedAA = new char[overlapIndex + 1 + seq2.length];
                    int[] confs = new int[overlapIndex + 1 + seq2.length];
                    for (int i = 0; i < overlapIndex; i++) {
                        mergedAA[i] = seq1[i];
                    }
                    for (int i = 0; i < seq1.length; i++) {
                        confs[i] = toRightContig.getConfs()[i];
                    }
                    for (int i = 0; i < seq2.length; i++) {
                        mergedAA[i + overlapIndex] = seq2[i];
                        confs[i + overlapIndex] += toLeftContig.getConfs()[i];
                    }

                    /* Create the merged contigs to have the score of the summation of the two contigs.
                       Although the mergedContigs might have insertion or deletion, in which case the
                       length of merged contig does not equal to the summation of the two lengths, still
                       use the start and end of the two merged contigs for later convenience to chain this
                       assembled contigs to correct positions on the template sequence.
                     */
                    Contig mergedContig = new Contig(toRightContig.gettStart(), toLeftContig.gettEnd(),
                                    mergedAA, confs, (toRightContig.getScore() + toLeftContig.getScore()));
                    bridgedContigs.add(mergedContig);

                    /*
                    System.out.println("bridge " + new String(seq1) + " " + toRightContig.getScore() + " "
                            + new String(seq2) + " " + toLeftContig.getScore() + " to " + mergedContig.toString());
                    if (mergedAA.length > (mergedContig.gettEnd() - mergedContig.gettStart())) {
                        System.out.println("Insertion");
                    }
                    */

                }
            }
        }

        Collections.sort(bridgedContigs, Contig.cmpReverseScore());
        int num = bridgedContigs.size() > 3 ? 3 : bridgedContigs.size();
        /*
        for (int i = 0; i < num; i++) {
            System.out.println(bridgedContigs.get(i).toString());
            if (bridgedContigs.get(i).getAAs().length > (bridgedContigs.get(i).gettEnd() - bridgedContigs.get(i).gettStart())) {
                System.out.println("Insertion");
            }
        }
        */
        return bridgedContigs;
    }

    /**
     * Find the index on seq1 starting from where the head of seq2 can overlap with the end of seq1
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
            while ((i < seq1.length) && !charEqual(seq1[i], seq2[0])) {
                i++;
            }
            if (i == seq1.length) {
                return i;
            }

            start = i;
            while ((i < seq1.length) && (j < seq2.length) && charEqual(seq1[i], seq2[j])) {
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

    private boolean charEqual(char c1, char c2) {
        if (c1 == c2) {
            return true;
        }
        if (((c1 == 'I') && (c2 == 'L'))
                || ((c1 == 'L') && (c2 == 'I'))
                || ((c1 == 'N') && (c2 == 'D'))
                || ((c1 == 'D') && (c2 == 'N'))) {
            return true;
        }
        return false;
    }
    /**
     * Assemble dnToRightAlign set into several contigs. If overlapped, merge into one contig.
     * @param dnAlignToRightList The list of dnAlign sorted by tStart ascending,
     *                           for same tStart, sort by score desceding.
     * @param scanDnMap
     * @return
     */
    private List<Contig> assembleDnToRightAlign(List<DenovoAligned> dnAlignToRightList,
                                                       HashMap<String, DenovoOnly> scanDnMap) {
        if (dnAlignToRightList.size() == 0) {
            return null;
        }
        List<Contig> contigs = new ArrayList<>();
        DenovoAligned firstDnAlign = dnAlignToRightList.get(0);
        Contig firstContig = new Contig(firstDnAlign.gettStart(), -1,
                scanDnMap.get(firstDnAlign.getDnScan()).getAAs(),
                copyConfs(scanDnMap.get(firstDnAlign.getDnScan()).getConfScores()),
                firstDnAlign.getScore());
        contigs.add(firstContig);
        //System.out.println("Debug " + firstContig.toString());

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
                //System.out.println("Debug Extended contig ");
                //System.out.println(contig.toString());
            } else {
                //If not overlapped, add the dnAligned as a new contig
                Contig newContig = new Contig(dnAligned.gettStart(), -1,
                        scanDnMap.get(dnAligned.getDnScan()).getAAs(),
                        copyConfs(scanDnMap.get(dnAligned.getDnScan()).getConfScores()),
                        dnAligned.getScore());
                //System.out.println("Debug new contig: " );
                //System.out.println(newContig.toString());
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
                if (!charEqual(contigAAs[j], dnAAs[j - startIndex])) {
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
     * @param dnAlignToLeftList The dnAlign list sorted by descending dn.tEnd, if same dn.tEnd, sort by score descending
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
                scanDnMap.get(firstDnAlign.getDnScan()).getAAs(),
                copyConfs(scanDnMap.get(firstDnAlign.getDnScan()).getConfScores()),
                firstDnAlign.getScore());
        contigs.add(firstContig);

        for (int i = 1; i < dnAlignToLeftList.size(); i++) {
            DenovoAligned dnAligned = dnAlignToLeftList.get(i);
            char[] dnAAs = scanDnMap.get(dnAligned.getDnScan()).getAAs();
            int dnTEnd = dnAligned.gettEnd();

            int overlappedContigIndex = findContigLeftOverlap(contigs, dnAAs, dnTEnd);

            //System.out.println("Debug " + dnTEnd + " " + dnAligned.toString() + " " + new String(dnAAs));

            if (overlappedContigIndex < contigs.size()) {
                Contig contig = contigs.get(overlappedContigIndex);
                char[] contigAAs = contig.getAAs();
                int startIndex = contig.gettEnd() - dnTEnd;

                //If overlapped, add up their score
                contig.setScore(contig.getScore() + dnAligned.getScore());

                //System.out.println("Debug " + contig.toString());

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
                        scanDnMap.get(dnAligned.getDnScan()).getAAs(),
                        copyConfs(scanDnMap.get(firstDnAlign.getDnScan()).getConfScores()),
                        dnAligned.getScore());
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

            if (startIndex > endIndex) {
                isOverlapped = false;
            }

            for (int j = startIndex; j < endIndex; j++) {
                if (!charEqual(contigAAs[contigAAs.length - 1 - j], dnAAs[dnAAs.length - 1 - (j - startIndex)])) {
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
     * To realize this, apply the contigs in a order of ascending scores.  Then the AA with higher score
     * will override the one with less score.
     * @param templateHooked
     * @param assembledContigs
     */
    private void old_generateCandidateTemplate(TemplateHooked templateHooked, List<Contig> assembledContigs) {
        System.out.println("Debug:" + templateHooked.getTemplateAccession());

        //TODO don't need merge.
        //List<Contig> mergedContigs = mergeContigs(assembledContigs);
        //Collections.sort(mergedContigs, Contig.cmpScore());

        Collections.sort(assembledContigs, Contig.cmpScore());
        printContigs(assembledContigs);

        char[] AAs = templateHooked.getSeq().clone();

        //for (Contig contig : mergedContigs) {
        for (Contig contig : assembledContigs) {
            int tStart = contig.gettStart() < 0 ? 0 : contig.gettStart();
            int tEnd = contig.gettEnd();

            for (int j = tStart; j <= tEnd; j++) {
                AAs[j] = contig.getAAs()[j - tStart];
            }
        }

        List<char[]> candidateTemplateSeqs = new ArrayList<>();
        candidateTemplateSeqs.add(AAs);
        templateHooked.setModifiedTemplates(candidateTemplateSeqs);
    }

    /**
     * A new candidate template generator which connect segments of template and assembled contigs.
     * Because the contigs assembled based on denovo result might contain deletion or insertion, the
     * method of generating new template is by chaining segments of old templates and the assembled contigs.
     * @param templateHooked
     * @param assembledContigs
     */
    private void generateCandidateTemplate(TemplateHooked templateHooked, List<Contig> assembledContigs) {
        System.out.println("Debug:" + templateHooked.getTemplateAccession());

        //TODO don't need merge.
        //List<Contig> mergedContigs = mergeContigs(assembledContigs);
        //Collections.sort(mergedContigs, Contig.cmpScore());

        Collections.sort(assembledContigs, Contig.cmpTStart());
        List<Character> newTemplate = new ArrayList<>();

        char[] AAs = templateHooked.getSeq().clone();
        int start = 0;

        //for (Contig contig : mergedContigs) {
        for (Contig contig : assembledContigs) {
            int tStart = contig.gettStart() < 0 ? 0 : contig.gettStart();
            int tEnd = contig.gettEnd();

            for (int i = start; i < tStart; i++) {
                newTemplate.add(AAs[i]);
            }
            for (char c : contig.getAAs()) {
                newTemplate.add(c);
            }
            start = contig.gettEnd() + 1;
        }
        for (int i = start; i < AAs.length; i++) {
            newTemplate.add(AAs[i]);
        }

        char[] newAAs = new char[newTemplate.size()];
        for (int i = 0; i < newTemplate.size(); i++) {
            newAAs[i] = newTemplate.get(i);
        }

        List<char[]> candidateTemplateSeqs = new ArrayList<>();
        candidateTemplateSeqs.add(newAAs);
        templateHooked.setModifiedTemplates(candidateTemplateSeqs);
    }


    /**
     * Merge adjacent overlapped contigs
     * @param assembledContigs The assembledContigs are sorted by descending score.
     * @param posDiffTolerance The gap size
     * @param minOverlap The minimal overlapped size required.
     * @return The extended contigs.  There might be contigs with same sequence, but
     *          different scores or same scores. The contig with higher score is added
     *          when extended with another nested contig.  The contig with lower score
     *          is not removed for the convenience of the code.  The contig with same score
     *          is the same one which is added also because of the convenience of the score.
     *          Those contigs with AA conflict but lower score are also kept.
     *          Be sure to run pickTopContigs process after this process to get correct sets of
     *          contigs.
     */
    private List<Contig> mergeContigs(List<Contig> assembledContigs,
                                      int posDiffTolerance, int minOverlap) {
        List<Contig> mergedContigList = new ArrayList<>();
        Collections.sort(assembledContigs, Contig.cmpReverseScore());
        Queue<Contig> contigQueue = new LinkedList<Contig>();
        contigQueue.addAll(assembledContigs);

        Contig contig = contigQueue.poll();
        mergedContigList.add(contig);
        while (null != (contig = contigQueue.poll())) {
            /* If there is position overlap between contig and any mergedContig,
                the contig should be tested whether it can be an extension of mergedContig.
                If it is an extension, extend it.  If not, it will be added to mergedContigList.
                Later, a process of pick topScore will filter out those overlapped with current
                mergedContigs.
             */
            Contig extendContig = null;
            int listLength = mergedContigList.size();
            for (int j = 0; j < listLength; j++) {
                Contig mergedContig = mergedContigList.get(j);
                if (contig.equals(mergedContig)) {
                    /* Set extendContig to be non null, if there is no extension,
                    it will not be added to contigQueue again.
                     */

                    extendContig = contig;
                    continue;
                }
                if ((contig.gettStart() > mergedContig.gettStart()) &&
                        (contig.gettStart() < (mergedContig.gettEnd() + posDiffTolerance))) {
                    //mergedContig is left to contig
                    extendContig = extendContig(mergedContig, contig, minOverlap);
                    if (extendContig != null) {
                        mergedContigList.set(j, extendContig);
                        /* Add the new extended Contig to contigQueue because the new
                            extended contig might overlap with other contigs.
                         */
                        contigQueue.add(extendContig);
                        break;
                    }
                }
                if ((mergedContig.gettStart() > contig.gettStart()) &&
                        (mergedContig.gettStart() < (contig.gettEnd() + posDiffTolerance))) {
                    //mergedContig is right to contig
                    extendContig = extendContig(contig, mergedContig, minOverlap);
                    if (extendContig != null) {
                        mergedContigList.set(j, extendContig);
                        /* Add the new extended Contig to contigQueue because the new
                            extended contig might overlap with other contigs.
                         */
                        contigQueue.add(extendContig);
                        break;
                    }
                }
            }

            if (extendContig == null) {
                mergedContigList.add(contig);
            }
        }
        return mergedContigList;

    }


    /**
     * Check whether contig is a right-direction extension of leftContig.
     * If yes, extend the leftContig to right and return the merged contigs.
     * The new contig will have the start of leftContig and the end of contig.
     * So there could generate deletion or insertion in the returned contig
     * @param leftContig  The contig in the mergedContig List
     * @param contig a new contig which is right to the mergedContig
     * @param minOverlapLen the minimal overlap AA numbers required to be viewed
     *                      as an extension
     * @return  null if no overlap or an extended contig
     */
    private Contig extendContig(Contig leftContig, Contig contig, int minOverlapLen) {
        Contig mergedContig = null;
        char[] seq1 = leftContig.getAAs();
        char[] seq2 = contig.getAAs();
        int overlapIndex = findOverlapIndex(seq1, seq2);
        if (overlapIndex < seq1.length && (seq1.length - overlapIndex) >= minOverlapLen) {
            char[] mergedAA = new char[overlapIndex + 1 + seq2.length];
            int[] confs = new int[overlapIndex + 1 + seq2.length];
            for (int i = 0; i < overlapIndex; i++) {
                mergedAA[i] = seq1[i];
            }
            /*
            for (int i = 0; i < seq1.length; i++) {
                confs[i] += leftContig.getConfs()[i];
            }
            */

            for (int i = 0; i < seq2.length; i++) {
                mergedAA[i + overlapIndex] = seq2[i];
            //    confs[i + overlapIndex] += contig.getConfs()[i];
            }

            /* Create the merged contigs to have the score of the summation of the two contigs.
               Although the mergedContigs might have insertion or deletion, in which case the
               length of merged contig does not equal to the summation of the two lengths, still
               use the start and end of the two merged contigs for later convenience to chain this
               assembled contigs to correct positions on the template sequence.
             */
            mergedContig = new Contig(leftContig.gettStart(), contig.gettEnd(),
                    mergedAA, confs, (leftContig.getScore() + contig.getScore()));
        }
        return mergedContig;
    }

    /* Merge overlapped contigs, choose the AAs with higher score as consensus.  there is bug in this function */
    private List<Contig> old_mergeContigs(List<Contig> assembledContigs) {
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

    private void printContigs(List<Contig> contigs) {
        for (Contig contig : contigs) {
            System.out.println(contig.toString());
        }
    }

    private List<Contig> filterContigs(List<Contig> contigs, int scoreThresh) {
        ArrayList<Contig> filteredContigs = new ArrayList<>();
        for (Contig contig : contigs) {
            if (contig.getScore() > scoreThresh * (contig.getAAs().length)) {
                filteredContigs.add(contig);
            }
        }
        return filteredContigs;
    }

    private Set<DenovoAligned> mergeDuplicateDn(Set<DenovoAligned> dnAlignSet, HashMap<String, DenovoOnly> scanDnMap) {
        Set<DenovoAligned> mergedDnAlignedSet = new HashSet<>();
        Map<String, DenovoAligned> AADnMap = new HashMap<>();
        for (DenovoAligned dnAligned : dnAlignSet) {
            String AAString = new String(scanDnMap.get(dnAligned.getDnScan()).getAAs());
            if (AADnMap.containsKey(AAString)) {
                System.out.println("merge " + AAString);
                AADnMap.get(AAString).setScore(AADnMap.get(AAString).getScore() + dnAligned.getScore());

                //tOdo add confScore
            } else {
                AADnMap.put(AAString, dnAligned);
            }
        }
        for (DenovoAligned dnAligned : AADnMap.values()) {
            mergedDnAlignedSet.add(dnAligned);
        }
        return mergedDnAlignedSet;
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
                //regionToAssemble.setDnAlignToRightSet(mergeDuplicateDn(regionToAssemble.dnAlignToRightSet, scanDnMap));
                //regionToAssemble.setDnAlignToLeftSet(mergeDuplicateDn(regionToAssemble.dnAlignToLeftSet, scanDnMap));

                assembledContigs.addAll(assembleOneRegion(regionToAssemble, scanDnMap));
            }

            /* Filter the assembled contigs according to the average score threshold
                for each AA.  A contig will be kept if its total score is greater than
                scoreThresh * length(contig).  If not, the contig is viewed as inconfident
                contig, which will be discarded.
             */
            System.out.println("Debug");
            printContigs(assembledContigs);
            int scoreThresh = 70;
            System.out.println("Filtering assembled contigs whose average conf score of AA is less than " + scoreThresh);
            List<Contig> filteredContigs = filterContigs(assembledContigs, scoreThresh);
            if (filteredContigs.size() == 0) continue;
            printContigs(filteredContigs);

            System.out.println("Extending assembled contigs according to overlap.");
            int posDiffTolerance = 10;
            int minOverlap = 3;
            List<Contig> mergedContigs = mergeContigs(filteredContigs, posDiffTolerance,
                                                        minOverlap);
            printContigs(mergedContigs);
            /* correct way currently
            System.out.println("Assembled contigs");
            generateCandidateTemplate(templateHooked, assembledContigs);
*/

            /*  There is some problem with assemble only contigs with top score because there are some other part not applied.
                Need to be considered again, how to chain confident AA or kmer other than segment.
             */

            System.out.println("Assembled contigs with top score");
            Collections.sort(mergedContigs, Contig.cmpReverseScore());
            List<Contig> topAssembledContigs = pickTopContigs(assembledContigs);
            for (Contig assembledContig : topAssembledContigs) {
                System.out.println(assembledContig.toString());
            }

            //Select contigs with max score and apply to the template to generate candidate template.
            generateCandidateTemplate(templateHooked, topAssembledContigs);



        }

    }




}

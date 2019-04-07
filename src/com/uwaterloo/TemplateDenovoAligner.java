package com.uwaterloo;
import java.util.*;

/**
 * Align denovo only peptides to templates. Build index for de novo only peptides. HashTable<KmerString, List<KmerPosition>>
 * Scan the template along.
 */
public class TemplateDenovoAligner {
    List<DenovoOnly> denovoOnlyList;
    HashMap<String, DenovoOnly> scanDnMap;

    public TemplateDenovoAligner(List<DenovoOnly> denovoOnlyList) {
        this.denovoOnlyList = denovoOnlyList;
        this.scanDnMap = buildScanDnMap();
        //this.denovoKmerIndexTable = new Hashtable<>();
    }

    /**
     * Build a map between the scan number and the DenovoOnly object.
     * @return
     */
    private HashMap<String, DenovoOnly> buildScanDnMap() {
        HashMap<String, DenovoOnly> scanDnMap = new HashMap<>();
        for (DenovoOnly dn : denovoOnlyList) {
            scanDnMap.put(dn.getScan(), dn);
        }
        return scanDnMap;
    }

    /**
     * Build a hashtable <kmer, Listof<kmerPosition>>  kmerPosition is <denovo scans, pos where kmer appears>
     * @param kmerSize the length of kmer
     */
    private Hashtable<String, List<KmerPosition>> buildKmerIndexTable(short kmerSize) {
        Hashtable<String, List<KmerPosition>> denovoKmerIndexTable = new Hashtable<>();
        int size = denovoOnlyList.size();
        for (int index = 0; index < size; index++) {
            DenovoOnly denovoOnly = denovoOnlyList.get(index);
            int dnLength = denovoOnly.getLength() - kmerSize;

            char[] AAs = denovoOnly.getAAs();

            for (short i = 0; i < dnLength; i++) {
                String kmer = "";
                for (short j = 0; j < kmerSize; j++) {
                    kmer += AAs[i + j];
                }
                KmerPosition kmerPosition = new KmerPosition(index, i);
                if (denovoKmerIndexTable.containsKey(kmer)) {
                    denovoKmerIndexTable.get(kmer).add(kmerPosition);
                } else {
                    List<KmerPosition> kmerPositionList = new ArrayList<>();
                    kmerPositionList.add(kmerPosition);
                    denovoKmerIndexTable.put(kmer, kmerPositionList);
                }
            }
        }
        return denovoKmerIndexTable;
    }

    private List<DenovoAligned> scanTemplate(TemplateHooked templateHooked, short kmerSize,
                                             Hashtable<String, List<KmerPosition>> denovoKmerIndexTable) {
        char[] seq = templateHooked.getSeq();
        /* Using set other than List is because a denovo only peptide may have
            multiple kmer overlapped with the template. After extend the kmer,
            they will generate one same alignment. Using set to keep only one.
         */
        Set<DenovoAligned> denovoAlignedSet = new HashSet<>();
        for (int i = 0; i < seq.length - kmerSize; i++) {
            String kmer = "";
            for (int j = 0; j < kmerSize; j++) {
                kmer += seq[i + j];
            }

            if (!denovoKmerIndexTable.containsKey(kmer)) {
                continue;
            }
            List<KmerPosition> kmerPositionList = denovoKmerIndexTable.get(kmer);

            denovoAlignedSet.addAll(extendDenovoMap(denovoOnlyList, kmerPositionList, templateHooked, i, kmerSize));
        }

        ArrayList<DenovoAligned> denovoAlignedList = new ArrayList<>();
        denovoAlignedList.addAll(denovoAlignedSet);
        return denovoAlignedList;
    }

    /**
     * Give a list of denovo shared kmer at position i of template, extend them to unextended.
     * @param denovoOnlyList The list of denovo only peptide
     * @param kmerPositionList A list of KmerPosition
     * @param templateHooked The template to be scanned
     * @param posOnTemplate The start position of the kmer in templateHooked
     */
    private Set<DenovoAligned> extendDenovoMap(List<DenovoOnly> denovoOnlyList, List<KmerPosition> kmerPositionList,
                                               TemplateHooked templateHooked, int posOnTemplate, short kmerSize) {
        char[] templateSeq = templateHooked.getSeq();
        Set<DenovoAligned> denovoAlignedList = new HashSet<>();

        for (KmerPosition kmerPosition : kmerPositionList) {
            DenovoOnly dn = denovoOnlyList.get(kmerPosition.getIndex());
            short dnPos = kmerPosition.getPos();
            char[] dnSeq = dn.getAAs();

            DenovoAligned denovoAligned = extendAdn(templateHooked, posOnTemplate, dn, dnPos, kmerSize);
            denovoAlignedList.add(denovoAligned);
        }
        return denovoAlignedList;
    }

    /**
     * Find the maximum sequence that a de novo peptide shared with a template
     * @param templateHooked
     * @param posOnTemplate
     * @param denovoOnly
     * @param dnPos
     */
    private DenovoAligned extendAdn(TemplateHooked templateHooked, int posOnTemplate, DenovoOnly denovoOnly, short dnPos, short kmerSize) {
        char[] templateSeq = templateHooked.getSeq();
        char[] dnSeq = denovoOnly.getAAs();
        int score = 0;

        //Extend the kmer to left until AA is different
        int leftTemp = posOnTemplate - 1;
        int leftDn = dnPos - 1;
        while ((leftDn >= 0) && (leftTemp >= 0)) {
            if ((templateSeq[leftTemp] == dnSeq[leftDn]) ||
                    ((templateSeq[leftTemp] == 'I') && (dnSeq[leftDn] == 'L')) ||
                    ((templateSeq[leftTemp] == 'L') && (dnSeq[leftDn] == 'I'))) {
                leftTemp--;
                leftDn--;
            } else {
                break;
            }
        }

        //Extend the kmer to right until AA is different.
        int rightTemp = posOnTemplate + kmerSize - 1;
        int rightDn = dnPos + kmerSize - 1;
        while ((rightTemp < templateSeq.length) && (rightDn < dnSeq.length)) {
            if ((templateSeq[rightTemp] == dnSeq[rightDn]) ||
                ((templateSeq[rightTemp] == 'I') && (dnSeq[rightDn] == 'L')) ||
                    ((templateSeq[rightTemp] == 'L') && (dnSeq[rightDn] == 'I'))) {
                rightTemp++;
                rightDn++;
            } else {
                break;
            }
        }

        for (int i = (leftDn + 1); i <= (rightDn - 1); i++) {
            score += denovoOnly.getConfScores()[i];
        }
        return new DenovoAligned(templateHooked.getTemplateId(), (leftTemp + 1), (rightTemp - 1),
                denovoOnly.getScan(), (leftDn + 1), (rightDn - 1), score);
    }

    private HashMap<String, List<DenovoAligned>> buildSortAlignList(List<DenovoAligned> denovoAlignedList) {
        HashMap<String, List<DenovoAligned>> dnAlignMap = new HashMap<>();

        //For each denovo only has overlap with template, set a list of aligned info with this scan
        for (DenovoAligned dnAligned : denovoAlignedList) {
            if (dnAlignMap.containsKey(dnAligned.getDnScan())) {
                dnAlignMap.get(dnAligned.getDnScan()).add(dnAligned);
            } else {
                List<DenovoAligned> alignListPerDn = new ArrayList<>();
                alignListPerDn.add(dnAligned);
                dnAlignMap.put(dnAligned.getDnScan(), alignListPerDn);
            }
        }

        //Sort aligned info of a dn according to the score
        for (List<DenovoAligned> dnAlignedList : dnAlignMap.values()) {
            Collections.sort(dnAlignedList, Collections.reverseOrder());
        }

        return dnAlignMap;
    }

    /**
     * According to sortedDnAlignList, build a new scanDnAlignListMap.
     * For each denovo only, only keep the alignment with highest score.
     * If there are multiple location with same score, attach the dn to
     * all of them.
     * @param sortedDnAlignListMap
     * @return
     */
    private HashMap<String, List<DenovoAligned>> extractDnAlignListWithMaxScore(HashMap<String, List<DenovoAligned>> sortedDnAlignListMap) {
        HashMap<String, List<DenovoAligned>> scanDnAlignListMap = new HashMap<>();
        for (String scan : sortedDnAlignListMap.keySet()) {
            List<DenovoAligned> dnAlignedList = sortedDnAlignListMap.get(scan);
            int maxScore = dnAlignedList.get(0).getScore();
            List<DenovoAligned> denovoAlignedListWithMaxScore = new ArrayList<>();
            for (DenovoAligned dnAligned : dnAlignedList) {
                //Attach dnAligned to template positions with maxScore
                if (dnAligned.getScore() == maxScore) {
                    denovoAlignedListWithMaxScore.add(dnAligned);
                } else {
                    break;
                }
            }
            scanDnAlignListMap.put(scan, denovoAlignedListWithMaxScore);
        }
        return scanDnAlignListMap;
    }

    /**
     * According to scanDnAlignListMap, for each denovo only, attach the de novo
     * only peptide to the template position with highest score. If there are
     * multiple location with same score, attach the dn to all of them.
     * @param templateHookedList A list of templateHooked with dnList to be added.
     * @param scanDnAlignListWithMaxScoreMap A hashmap <dnScan, dnAlignedListWithMaxScore>
     */
    private void hookDnToTemplate(List<TemplateHooked> templateHookedList,
                                  HashMap<String, List<DenovoAligned>> scanDnAlignListWithMaxScoreMap) {
        for (String scan : scanDnAlignListWithMaxScoreMap.keySet()) {
            List<DenovoAligned> dnAlignedList = scanDnAlignListWithMaxScoreMap.get(scan);
            for (DenovoAligned dnAligned : dnAlignedList) {
                int templateId = dnAligned.getTemplateId();
                /* Except hooked to same region with the template.  The denovo only
                    should also be hooked to unmatched region, so that the denovo sequence
                    can count on major vote together with the template region.
                 */
                int tStart = dnAligned.gettStart();
                int tEnd = dnAligned.gettEnd();

                int dnLength = scanDnMap.get(dnAligned.getDnScan()).getLength();

                if (dnAligned.getDnStart() == 0) {
                    //If denovo extend template to right
                    tEnd = tStart + dnLength - 1;
                    int templateLength = templateHookedList.get(templateId).getSeq().length;
                    tEnd = tEnd < templateLength ? tEnd : (templateLength - 1);
                } else {
                    //If denovo extend template to left
                    tStart = tEnd - dnLength + 1;
                    tStart = tStart > 0 ? tStart : 0;
                }
                //Attache the dnAlign to all positions it covered
                for (int i = tStart; i <= tEnd; i++) {
                    templateHookedList.get(templateId).getDnList().get(i).add(dnAligned);
                }
            }
        }
    }

    /**
     * Print denovo only peptides according to their alignment to the template.
     * e.g.
     * VTLGCLVKGYFHYSSL  2 138 147 F2:12912 0 9 656
     *    GCLVKGYFHYSSL  2 141 147 F2:9316 0 6 436
     *      LVKGYFHYSSL  2 143 147 F2:8159 0 4 374
     * VTLGCLVKGYFWNPVTL  2 138 147 F2:14482 0 9 634
     * @param templateHookedList
     */
    private void printAlignedDn(List<TemplateHooked> templateHookedList) {
        for (TemplateHooked templateHooked : templateHookedList) {
            System.out.println("Template " + templateHooked.getTemplateAccession());
            for (int i = 0; i < templateHooked.getSeq().length; i++) {
                int dnSize = templateHooked.getDnList().get(i).size();
                int psmSize = templateHooked.getMappedScanList().get(i).size();
                if (dnSize == 0) {
                    continue;
                }
                float psmDnRatio = ((float) psmSize) / dnSize;
//                if (templateHooked.getDnList().get(i).size() >= 20) {
                if (psmDnRatio < 2) {
                    String scanList = "";
                    int minStart = 500;
                    for (DenovoAligned dnA : templateHooked.getDnList().get(i)) {
                        scanList += dnA.dnScan + " ";
                        int currentStart = dnA.gettStart() - dnA.getDnStart();
                        minStart = currentStart < minStart ? currentStart : minStart;
                    }
                    System.out.println("pos " + i + "ratio: " + psmDnRatio + " db num: " + templateHooked.getMappedScanList().get(i).size() +
                            " dn num: " + templateHooked.getDnList().get(i).size()  + " " + scanList);

                    for (DenovoAligned dnA : templateHooked.getDnList().get(i)) {
                        String prefix = "";
                        for (int j = 0; j <= (dnA.gettStart() - dnA.getDnStart() - minStart); j++) {
                            prefix += " ";
                        }
                        System.out.print(prefix);
                        System.out.print(scanDnMap.get(dnA.getDnScan()).getAAs());
                        System.out.println("  " + dnA.toString());
                    }

                }
            }
        }
    }

    private HashMap<String, List<DenovoAligned>> getDnCanExtendRight(HashMap<String, List<DenovoAligned>> scanDnAlignListMap) {
        HashMap<String, List<DenovoAligned>> dnCanExtendRightMap= new HashMap<>();
        for (String scan : scanDnAlignListMap.keySet()) {
            List<DenovoAligned> dnAlignList = scanDnAlignListMap.get(scan);
            List<DenovoAligned> dnAlignCanExtend = new ArrayList<>();
            for (DenovoAligned dnAlign : dnAlignList) {
                if (dnAlign.getDnStart() == 0) {
                    dnAlignCanExtend.add(dnAlign);
                }
            }
            if (dnAlignCanExtend.size() > 0) {
                dnCanExtendRightMap.put(scan, dnAlignCanExtend);
            }
        }
        return dnCanExtendRightMap;
    }

    private HashMap<String, List<DenovoAligned>> getDnCanExtendLeft(HashMap<String, List<DenovoAligned>> scanDnAlignListMap,
                                                                    HashMap<String, DenovoOnly> scanDnMap) {
        HashMap<String, List<DenovoAligned>> dnCanExtendRightMap= new HashMap<>();
        for (String scan : scanDnAlignListMap.keySet()) {
            List<DenovoAligned> dnAlignList = scanDnAlignListMap.get(scan);
            List<DenovoAligned> dnAlignCanExtend = new ArrayList<>();
            for (DenovoAligned dnAlign : dnAlignList) {
                if (dnAlign.getDnEnd() == (scanDnMap.get(dnAlign.getDnScan()).getLength() - 1)) {
                    dnAlignCanExtend.add(dnAlign);
                }
            }
            if (dnAlignCanExtend.size() > 0) {
                dnCanExtendRightMap.put(scan, dnAlignCanExtend);
            }
        }
        return dnCanExtendRightMap;
    }

    /**
     * If a position's dn num > thresh, print both its db num and dn num on that position
     * @param templateHookedList
     */
    private void printTemplateDbvsDn(List<TemplateHooked> templateHookedList) {
        int thresh = 0;
        for (TemplateHooked templateHooked : templateHookedList) {
            int seqLength = templateHooked.getSeq().length;
            System.out.println(templateHooked.getTemplateAccession());
            System.out.println("Pos,Db,Dn,Ratio");
            for (int i = 0; i < seqLength; i++) {
                if (templateHooked.getDnList().get(i).size() > thresh) {
                    float ratio = ((float) templateHooked.getMappedScanList().get(i).size()) / templateHooked.getDnList().get(i).size();
                    System.out.println(i + "," + templateHooked.getMappedScanList().get(i).size() + "," +
                            templateHooked.getDnList().get(i).size() + "," + ratio);
                }
            }
        }
    }

    public void alignDenovoOnlyToTemplate(List<TemplateHooked> templateHookedList, short kmerSize) {
        Hashtable<String, List<KmerPosition>> denovoKmerIndexTable = buildKmerIndexTable(kmerSize);
        List<DenovoAligned> denovoAlignedList = new ArrayList<>();
        for (TemplateHooked templateHooked : templateHookedList) {
            denovoAlignedList.addAll(scanTemplate(templateHooked, kmerSize, denovoKmerIndexTable));
        }
        System.out.println("denovo align size: " + denovoAlignedList.size());
        HashMap<String, List<DenovoAligned>> sortedAlignListMap = buildSortAlignList(denovoAlignedList);
        System.out.println("Remaining denovo only: " + sortedAlignListMap.size());

        //For each denovo only, find its best matches to template
        HashMap<String, List<DenovoAligned>> scanDnAlignMap = extractDnAlignListWithMaxScore(sortedAlignListMap);
        //Hook the denovo only to the best matched template position
        //hookDnToTemplate(templateHookedList, scanDnAlignMap);

        //Find those dn whose left end are the same with the aligned template dnStart = 0

        //Find those dn whose right end are the same with the aligned template  dnEnd = denovo length - 1.
        HashMap<String, List<DenovoAligned>> dnExtendRight = getDnCanExtendRight(scanDnAlignMap);
        hookDnToTemplate(templateHookedList, dnExtendRight);

        HashMap<String, List<DenovoAligned>> dnExtendLeft = getDnCanExtendLeft(scanDnAlignMap, scanDnMap);
        hookDnToTemplate(templateHookedList, dnExtendLeft);
//        printTemplateDbvsDn(templateHookedList);
        printAlignedDn(templateHookedList);


    }
}

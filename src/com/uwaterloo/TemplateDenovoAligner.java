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
            if ((templateSeq[rightTemp] != dnSeq[rightDn]) ||
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

    private Hashtable<String, List<DenovoAligned>> buildSortAlignList(List<DenovoAligned> denovoAlignedList) {
        Hashtable<String, List<DenovoAligned>> dnAlignMap = new Hashtable<>();

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
     * According to sortedDnAlignList, for each denovo only, attach the de novo
     * only peptide to the template position with highest score. If there are
     * multiple location with same score, attach the dn to all of them.
     * @param templateHookedList A list of templateHooked with dnList to be added.
     * @param sortedDnAlignListMap A hashmap <dnScan, sortedDnAlignedList>
     */
    private void assignDnToTemplate(List<TemplateHooked> templateHookedList,
                                    Hashtable<String, List<DenovoAligned>> sortedDnAlignListMap) {
        for (String scan : sortedDnAlignListMap.keySet()) {
            List<DenovoAligned> dnAlignedList = sortedDnAlignListMap.get(scan);
            int maxScore = dnAlignedList.get(0).getScore();
            for (DenovoAligned dnAligned : dnAlignedList) {
                //Attach dnAligned to template positions with maxScore
                if (dnAligned.getScore() == maxScore) {
                   // System.out.println(dnAligned);
                    int templateId = dnAligned.getTemplateId();
                    int tStart = dnAligned.gettStart();
                    int tEnd = dnAligned.gettEnd();

                    //Attache the dnAlign to all positions it covered
                    for (int i = tStart; i <= tEnd; i++) {
                        templateHookedList.get(templateId).getDnList().get(i).add(dnAligned);
                    }
                } else {
                    break;
                }
            }
        }
        printAlignedDn(templateHookedList, sortedDnAlignListMap);

    }

    private void printAlignedDn(List<TemplateHooked> templateHookedList, Hashtable<String, List<DenovoAligned>> sortedAlignListMap) {
        for (TemplateHooked templateHooked : templateHookedList) {
            System.out.println("Template " + templateHooked.getTemplateAccession());
            for (int i = 0; i < templateHooked.getSeq().length; i++) {
                if (templateHooked.getDnList().get(i).size() >= 10) {
                    String scanList = "";
                    int minStart = 500;
                    for (DenovoAligned dnA : templateHooked.getDnList().get(i)) {
                        scanList += dnA.dnScan + " ";
                        int currentStart = dnA.gettStart() - dnA.getDnStart();
                        minStart = currentStart < minStart ? currentStart : minStart;
                    }
                    System.out.println("pos " + i + " db num: " + templateHooked.getMappedScanList().get(i).size() +
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

    public void alignDenovoOnlyToTemplate(List<TemplateHooked> templateHookedList, short kmerSize) {


        Hashtable<String, List<KmerPosition>> denovoKmerIndexTable = buildKmerIndexTable(kmerSize);
        List<DenovoAligned> denovoAlignedList = new ArrayList<>();
        for (TemplateHooked templateHooked : templateHookedList) {
            denovoAlignedList.addAll(scanTemplate(templateHooked, kmerSize, denovoKmerIndexTable));
        }
        System.out.println("denovo align size: " + denovoAlignedList.size());
        Hashtable<String, List<DenovoAligned>> sortedAlignListMap = buildSortAlignList(denovoAlignedList);
        System.out.println("Remaining denovo only: " + sortedAlignListMap.size());

        assignDnToTemplate(templateHookedList, sortedAlignListMap);

    }
}

package com.uwaterloo.ScanTemplateMapper;
import com.uwaterloo.DenovoAssembler.DenovoAligned;

import java.util.ArrayList;

public class TemplateHooked extends Template {
    ArrayList<ArrayList<String>> mappedScanList;   //The list of scans mapped to each position of the template.
    ArrayList<ArrayList<PSMAligned>> dbList;    //List of db starting at pos
    ArrayList<ArrayList<PSMAligned>> spiderList;    //List of spider starting at pos
    ArrayList<ArrayList<DenovoAligned>> dnToRightList; //List of denovo only covering pos whose left end overlapping with template
    ArrayList<ArrayList<DenovoAligned>> dnToLeftList; //List of denovo only covering pos whose right end overlapping with template

    public TemplateHooked(Template template) {
        super(template.getTemplateId(), template.getTemplateAccession(), template.getSeq());
        this.mappedScanList = new ArrayList<>();
        this.dbList = new ArrayList<>();
        this.spiderList = new ArrayList<>();
        this.dnToLeftList = new ArrayList<>();
        this.dnToRightList = new ArrayList<>();
        initializeAlignList(this.seq.length);
    }

    /* Initialize the scanList, dblist, spiderlist to have the size
        of template length
     */
    private void initializeAlignList(int size) {
        for (int i = 0; i < size; i++) {
            ArrayList<String> scanList = new ArrayList<>();
            this.mappedScanList.add(scanList);

            ArrayList<PSMAligned> dbList = new ArrayList<>();
            this.dbList.add(dbList);

            ArrayList<PSMAligned> spiderList = new ArrayList<>();
            this.spiderList.add(spiderList);

            ArrayList<DenovoAligned> dnToRightList = new ArrayList<>();
            this.dnToRightList.add(dnToRightList);

            ArrayList<DenovoAligned> dnToLeftList = new ArrayList<>();
            this.dnToLeftList.add(dnToLeftList);

        }
    }

    public ArrayList<ArrayList<PSMAligned>> getDbList() {
        return this.dbList;
    }

    public ArrayList<ArrayList<PSMAligned>> getSpiderList() {
        return this.spiderList;
    }

    public ArrayList<ArrayList<String>> getMappedScanList() {
        return this.mappedScanList;
    }

    public ArrayList<ArrayList<DenovoAligned>> getDnToLeftList() {
        return dnToLeftList;
    }

    public ArrayList<ArrayList<DenovoAligned>> getDnToRightList() {
        return dnToRightList;
    }

    /**
     * The class to store a peptide spectrum match.
     * It could be used to store db result, spider result
     * and de novo result.
     */
    public static class PSM {
        String scan;
        String peptide;
        /* The ion score of each AA. If an AA contain both ions around it, it will have 100.
            If it belongs to a seg of length n, then each AA will have int(100 / n) score.
         */
        short[] ionScores;
        int intensity;

        public PSM(String scan, String peptide) {
            this.scan = scan;
            this.peptide = peptide;

        }

        public PSM(String scan, String peptide, int intensity) {
            this.scan = scan;
            this.peptide = peptide;
            this.intensity = intensity;
        }

        public String getScan() {
            return scan;
        }

        public String getPeptide() {
            return peptide;
        }

        public short[] getIonScores() {
            return ionScores;
        }

        public void setIonScores(short[] ionScores) {
            this.ionScores = ionScores;
        }

        public int getIntensity() {
            return intensity;
        }
    }
}

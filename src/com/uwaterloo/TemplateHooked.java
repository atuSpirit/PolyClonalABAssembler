package com.uwaterloo;
import java.util.LinkedList;
import java.util.ArrayList;

/* Temorary no use */
public class TemplateHooked extends Template {
    ArrayList<ArrayList<String>> mappedScanList;   //The list of scans mapped to each position of the template.
    ArrayList<ArrayList<PSMAligned>> dbList;
    ArrayList<ArrayList<PSMAligned>> spiderList;
    //ToDo  add arrayOfDenovoAligned.

    public TemplateHooked(Template template) {
        super(template.getTemplateId(), template.getTemplateAccession(), template.getSeq());
        this.mappedScanList = new ArrayList<>();
        this.dbList = new ArrayList<>();
        this.spiderList = new ArrayList<>();
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

}

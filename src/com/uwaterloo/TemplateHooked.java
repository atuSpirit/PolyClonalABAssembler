package com.uwaterloo;
import java.util.LinkedList;
import java.util.ArrayList;

/* Temorary no use */
public class TemplateHooked extends Template {
    ArrayList<LinkedList<String>> mappedScanList;   //The list of scans mapped to each position of the template.
    ArrayList<LinkedList<PSMAligned>> dbList;
    ArrayList<LinkedList<PSMAligned>> spiderList;
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
            LinkedList<String> scanList = new LinkedList<>();
            this.mappedScanList.add(scanList);

            LinkedList<PSMAligned> dbLinkedList = new LinkedList<>();
            this.dbList.add(dbLinkedList);

            LinkedList<PSMAligned> spiderLinkedList = new LinkedList<>();
            this.spiderList.add(spiderLinkedList);
        }
    }

    public ArrayList<LinkedList<PSMAligned>> getDbList() {
        return this.dbList;
    }

    public ArrayList<LinkedList<PSMAligned>> getSpiderList() {
        return this.spiderList;
    }

    public ArrayList<LinkedList<String>> getMappedScanList() {
        return this.mappedScanList;
    }

}

package com.uwaterloo;
import java.util.LinkedList;
import java.util.ArrayList;

/* Temorary no use */
public class TemplateHooked extends Template {
    ArrayList<LinkedList<PSMAligned>> dbList;
    ArrayList<LinkedList<PSMAligned>> spiderList;
    //ToDo  add arrayOfDenovoAligned.

    public TemplateHooked(Template template) {
        super(template.getTemplateId(), template.getTemplateAccession(), template.getSeq());
        this.dbList = new ArrayList<>();
        this.spiderList = new ArrayList<>();
        initializeAlignList(this.seq.length);
    }

    private void initializeAlignList(int size) {
        for (int i = 0; i < size; i++) {
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

}

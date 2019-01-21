package com.uwaterloo;

/* Temorary no use */
public class TemplateHooked extends Template {
    PSMAligned[] arrayOfDBAligned;
    PSMAligned[] arrayOfSpiderAligned;
    //ToDo  add arrayOfDenovoAligned.

    public TemplateHooked(Template template) {
        super(template.getTemplateId(), template.getTemplateAccession(), template.getSeq());
        this.arrayOfDBAligned = new PSMAligned[seq.length];
        this.arrayOfSpiderAligned = new PSMAligned[seq.length];
    }

    public PSMAligned[] getArrayOfDBAligned() {
        return arrayOfDBAligned;
    }

    public PSMAligned[] getArrayOfSpiderAligned() {
        return arrayOfSpiderAligned;
    }
}

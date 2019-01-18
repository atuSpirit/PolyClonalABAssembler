package com.uwaterloo;

public class Template {
    int templateId;
    String templateAccession;
    char[] seq;
    char[] modifiedSeq;

    public Template(int templateId, String templateAccession, char[] seq) {
        this.templateId = templateId;
        this.templateAccession = templateAccession;
        this.seq = seq;
    }

    public int getTemplateId() {
        return templateId;
    }

    public String getTemplateAccession() {
        return templateAccession;
    }

    public char[] getSeq() {
        return seq;
    }

    public char[] getModifiedSeq() {
        return modifiedSeq;
    }
}

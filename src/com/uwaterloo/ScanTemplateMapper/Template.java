package com.uwaterloo.ScanTemplateMapper;

import java.util.List;

public class Template {
    int templateId;
    String templateAccession;
    char[] seq;
    List<char[]> modifiedTemplates;


    public Template(int templateId, String templateAccession, char[] seq) {
        this.templateId = templateId;
        this.templateAccession = templateAccession;
        this.seq = seq;
    }

    public void setModifiedTemplates(List<char[]> candidateTemplates) {
        this.modifiedTemplates = candidateTemplates;
    }

    public int getTemplateId() {
        return templateId;
    }

    public String getTemplateAccession() {
        return templateAccession;
    }

    public void setTemplateAccession(String templateAccession) {
        this.templateAccession = templateAccession;
    }

    public char[] getSeq() {
        return seq;
    }

    public List<char[]> getModifiedSeq() {
        return modifiedTemplates;
    }
}

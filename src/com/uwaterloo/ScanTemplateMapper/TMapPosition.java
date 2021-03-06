package com.uwaterloo.ScanTemplateMapper;

/* A helper class to store the start position of a peptide
    mapped to the template.
 */
public class TMapPosition {
    int templateId; //The id of the template sequence
    int start;    //The start position the peptide is mapped to the template sequence
    int end;    //The end position the peptide is mapped to the template sequence

    public TMapPosition(int templateId, int start, int end) {
        this.templateId = templateId;
        this.start = start;
        this.end = end;
    }

    public int getTemplateId() {
        return templateId;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    @Override
    public String toString() {
        return " start: " + start + " end: " + end;
    }
}

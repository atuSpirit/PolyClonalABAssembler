package com.uwaterloo;

/* A helper class to store the start position of a peptide
    mapped to the template.
 */
public class TMapPosition {
    int templateId; //The id of the template sequence
    int start;    //The start position the peptide is mapped to the template sequence

    public TMapPosition(int templateId, int start) {
        this.templateId = templateId;
        this.start = start;
    }

    public int getTemplateId() {
        return templateId;
    }

    public int getStart() {
        return start;
    }
}

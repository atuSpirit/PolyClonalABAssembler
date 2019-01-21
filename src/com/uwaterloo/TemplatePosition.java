package com.uwaterloo;

import java.util.ArrayList;
import java.util.List;

/**
 * The class storing spectra aligned to this position
 */
public class TemplatePosition {
    int pos;    //The position on the template
    char templateAA;    //The amino acid at the pos of the template
    List<PSMAligned> dbList;   //The list of spectra with db result including PTM
    List<PSMAligned> dnList;   //The List of spectra with de novo only result
    List<PSMAligned> spList;   //The list of spectra with spider result

    public TemplatePosition(int pos, char templateAA) {
        this.pos = pos;
        this.templateAA = templateAA;
        dbList = new ArrayList<>();
        dnList = new ArrayList<>();
        spList = new ArrayList<>();
    }

    public int getPos() {
        return pos;
    }

    public char getTemplateAA() {
        return templateAA;
    }

    public List<PSMAligned> getDbList() {
        return dbList;
    }

    public List<PSMAligned> getDnList() {
        return dnList;
    }

    public List<PSMAligned> getSpList() {
        return spList;
    }

}

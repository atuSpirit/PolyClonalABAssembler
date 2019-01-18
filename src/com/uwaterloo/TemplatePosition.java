package com.uwaterloo;

import java.util.ArrayList;
import java.util.List;

/**
 * The class storing spectra aligned to this position
 */
public class TemplatePosition {
    int pos;    //The position on the template
    char templateAA;    //The amino acid at the pos of the template
    List<PSM> dbList = null;   //The list of spectra with db result including PTM
    List<PSM> dnList = null;   //The List of spectra with de novo only result
    List<PSM> spList = null;   //The list of spectra with spider result

    public TemplatePosition(int pos, char templateAA) {
        this.pos = pos;
        this.templateAA = templateAA;
    }


}

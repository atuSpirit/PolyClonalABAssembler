package com.uwaterloo;

import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

/**
 * The class to store a peptide spectrum match.
 * It could be used to store db result, spider result
 * and de novo result.
 */
public class PSM {
    String scan;
    String peptide;

    /*
    public PSM(String scan, String peptide, char[] AAs, short[] confScores) {
        this.scan = scan;
        this.peptide = peptide;
        this.AAs = AAs;
        this.confScores = confScores;
        setPositionOfVariations(peptide);
        this.mapPositionList = new ArrayList<TMapPosition>();
    }
    */

    public PSM(String scan, String peptide) {
        this.scan = scan;
        this.peptide = peptide;

    }

    public String getScan() {
        return scan;
    }

    public String getPeptide() {
        return peptide;
    }


}

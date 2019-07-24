package com.uwaterloo.Tools;

import Utils.AAMass;

import java.util.List;

public class IntactMassValidater {
    public static float computeIntactMass(char[] proteinSeq) {
        float intactMass = 0.0f;

        for (char AA : proteinSeq) {
            String AAstr = "";
            AAstr += AA;
            intactMass += AAMass.AA_MASS_TABLE.get(AAstr);
        }
        return intactMass;
    }



}

package com.uwaterloo.Utils;

import com.uwaterloo.DenovoAssembler.DenovoOnly;

import java.util.ArrayList;
import java.util.List;

public class PeptideFilter {
    int confScoreThresh;
    int inConfidentAANumThresh;
    public PeptideFilter(int confScoreThresh, int inConfidentAANumThresh) {
        this.confScoreThresh = confScoreThresh;
        this.inConfidentAANumThresh = inConfidentAANumThresh;
    }

    /**
     * Filter out those denovo only peptide containing confScore less than confScoreThresh.
     * @param dnList
     * @return a new dnList with high conf score
     */
    public List<DenovoOnly> filterDnByConfScore(List<DenovoOnly> dnList) {
        List<DenovoOnly> newDnList = new ArrayList<>();
        for (DenovoOnly dn : dnList) {
            if (isConfident(dn.getConfScores())) {
                newDnList.add(dn);
            }
        }
        return newDnList;
    }

    /**
     * Filter out those psm peptide containing confScore less than confScoreThresh.
     * @param psmList
     * @return a new psmList with high conf score
     */
    public List<PSM> filterPSMByConfScore(List<PSM> psmList) {
        List<PSM> newPSMList = new ArrayList<>();
        for (PSM psm : psmList) {
            if (psm.getIonScores() == null) {
                continue;
            }
            if (isConfident(psm.getIonScores())) {
                newPSMList.add(psm);
            }
        }
        return newPSMList;
    }


    /**
     * Judge whether a peptide with confScores contains conf score less
     * than confScoreThresh
     * @param confScores
     * @return
     */
    boolean isConfident(short[] confScores) {
        boolean isConfident = true;
        int inConfidentNum = 0;
        for (short confScore : confScores) {
            if (confScore < this.confScoreThresh) {
                inConfidentNum++;
            }
            if (inConfidentNum >= this.inConfidentAANumThresh) {
                isConfident = false;
                break;
            }
        }
        return isConfident;
    }
}

package com.uwaterloo;

/**
 * The class to store a peptide spectrum match.
 * It could be used to store db result, spider result
 * and de novo result.
 */
public class PSM {
    String scan;
    String peptide;
    /* The ion score of each AA. If an AA contain both ions around it, it will have 100.
        If it belongs to a seg of length n, then each AA will have int(100 / n) score.
     */
    short[] ionScores;
    int intensity;

    public PSM(String scan, String peptide) {
        this.scan = scan;
        this.peptide = peptide;

    }

    public PSM(String scan, String peptide, int intensity) {
        this.scan = scan;
        this.peptide = peptide;
        this.intensity = intensity;
    }

    public String getScan() {
        return scan;
    }

    public String getPeptide() {
        return peptide;
    }

    public short[] getIonScores() {
        return ionScores;
    }

    public void setIonScores(short[] ionScores) {
        this.ionScores = ionScores;
    }

    public int getIntensity() {
        return intensity;
    }
}

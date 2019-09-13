package com.uwaterloo.DenovoAssembler;

public class DenovoOnly {
    String scan;
    String peptide;
    short alc;
    short length;
    char[] AAs;
    short[] confScores;
    int intensity;

    public DenovoOnly(String scan, String peptide, short alc, short length, String confString, int intensity) {
        this.scan = scan;
        this.peptide = peptide;
        this.alc = alc;
        this.length = length;
        this.AAs = convertToAA(peptide);
        this.confScores = convertToScore(confString);
        this.intensity = intensity;
    }

    public DenovoOnly(String scan, String peptide, short alc, short length, char[] AAs, short[] confScores, int intensity) {
        this.scan = scan;
        this.peptide = peptide;
        this.alc = alc;
        this.length = length;
        this.AAs = AAs;
        this.confScores = confScores;
        this.intensity = intensity;
    }

    /**
     * Convert the local confidence score into an array of scores.
     * @param confString
     * @return
     */
    private short[] convertToScore(String confString) {
        String[] scoreStrs = confString.split(" ");
        short[] scores = new short[scoreStrs.length];
        for (int i = 0; i < scoreStrs.length; i++) {
            scores[i] = Short.valueOf(scoreStrs[i]);
        }
        return scores;
    }

    /**
     * Convert the peptide to an array containing only AAs.
     * @param peptide
     * @return
     */
    private char[] convertToAA(String peptide) {
        //Remove all PTM
        peptide = peptide.replaceAll("\\(\\S(\\d)*.(\\d)+\\)", "");

        return peptide.toCharArray();
    }

    public String getScan() {
        return scan;
    }

    public String getPeptide() {
        return peptide;
    }

    public short getAlc() {
        return alc;
    }

    public short getLength() {
        return length;
    }

    public char[] getAAs() {
        return AAs;
    }

    public short[] getConfScores() {
        return confScores;
    }

    public int getIntensity() {
        return intensity;
    }

    @Override
    public String toString() {
        String dnStr = scan + " " + peptide + " " + alc + " " + length + " " + intensity;
        return dnStr;
    }
}

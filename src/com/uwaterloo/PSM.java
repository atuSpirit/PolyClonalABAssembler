package com.uwaterloo;

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
    /*The amino acid of each position. If spider peptide, score the true seq.*/
    char[] AAs;
    /* The confidence score of each AA. If peptide is a de novo only, it will be
        the ALC score. If it is a db or spider, we should find a way to measure
        whether there is a fragment peak here. If yes, it should have higher score.
     */
    short[] confScores;
    /* The list of positions where variations (substitution, insertion
        deletion) located. If not spider seq, this field is empty.
     */
    List<Integer> positionOfVariations;

    public PSM(String scan, String peptide, char[] AAs, short[] confScores) {
        this.scan = scan;
        this.peptide = peptide;
        this.AAs = AAs;
        this.confScores = confScores;

        setPositionOfVariations(peptide);
    }

    public PSM(String scan, String peptide) {
        this.scan = scan;
        this.peptide = peptide;
        this.AAs = convertToAA(peptide);
        this.confScores = convertToConfScore(AAs);
        setPositionOfVariations(peptide);
    }

    private char[] convertToAA(String peptide) {
        //Remove all PTM
        peptide = peptide.replaceAll("\\(\\S(\\d)+.(\\d)+\\)", "");

        //Remove all ins|del|sub
        peptide = peptide.replaceAll("\\(sub \\S\\)", "");
        peptide = peptide.replaceAll("\\(ins\\)", "");
        peptide = peptide.replaceAll("\\(del\\)", "");

        return peptide.toCharArray();
    }

    private short[] convertToConfScore(char[] AAs) {
        int length = AAs.length;
        short[] confScore = new short[length];
        for (int i = 0; i < length; i++) {
            confScore[i] = 99;
        }
        return confScore;
    }

    private void setPositionOfVariations(String peptide) {
        if (peptide.contains("sub") || peptide.contains("ins") || peptide.contains("del")) {
            this.positionOfVariations = new ArrayList<Integer>();

            //Remove the PTM
            peptide = peptide.replaceAll("\\(\\S(\\d)+.(\\d)+\\)", "");

            //Put the position where (sub *), (ins), (del) to positionOfVariations List.
            int i = 0;
            int length = peptide.length();
            int pos = 0;
            char c;
            while (i < length) {
                c = peptide.charAt(i);
                if (c == '(') {
                    this.positionOfVariations.add(pos - 1);
                    while (c != ')') {
                        i += 1;
                        c = peptide.charAt(i);
                    }
                }
                i += 1;
                pos += 1;
            }

        } else {
            this.positionOfVariations = null;
        }
    }

    public String getScan() {
        return scan;
    }

    public String getPeptide() {
        return peptide;
    }

    public char[] getAAs() {
        return AAs;
    }

    public short[] getConfScores() {
        return confScores;
    }

    public List<Integer> getPositionOfVariations() {
        return positionOfVariations;
    }
}

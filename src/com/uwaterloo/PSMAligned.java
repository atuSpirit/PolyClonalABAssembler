package com.uwaterloo;

import java.util.ArrayList;
import java.util.List;

public class PSMAligned extends PSM {
    int templateId; //The id of the template sequence
    int start;
    int end;

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

    /*  Based on an assumption that one psm can only map to one location one one template,
        the field is outdated: The list of positions on template the peptide of the psm is mapped to.
     * Computed by connecting PSM and protein-peptide table.
    List<TMapPosition> mapPositionList; */

    public PSMAligned(String scan, String peptide, int templateId, int start, int end) {
        super(scan, peptide);

        this.templateId = templateId;
        this.AAs = convertToAA(peptide);
        this.confScores = convertToConfScore(AAs);
        setPositionOfVariations(peptide);
        this.start = start;
        this.end = end;
    }


    private char[] convertToAA(String peptide) {
        //Remove all PTM
        peptide = peptide.replaceAll("\\(\\S(\\d)*.(\\d)+\\)", "");

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
            this.positionOfVariations = new ArrayList<>();

            //Remove the PTM
            peptide = peptide.replaceAll("\\(\\S(\\d)*.(\\d)+\\)", "");

            //Put the position where (sub *), (ins), (del) to positionOfVariations List.
            int i = 0;
            int length = peptide.length();
            int pos = 0;
            char c;
            while (i < length) {
                c = peptide.charAt(i);
                if (c == '(') {
                    pos -= 1;
                    this.positionOfVariations.add(pos);
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

    public int getTemplateId() {
        return templateId;
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

    /*
    public List<TMapPosition> getMapPositionList() {
        return mapPositionList;
    }
    */

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    @Override
    public String toString() {
        String posOfVar = "";
        for (int i : positionOfVariations) {
            posOfVar += String.valueOf(i) + " ";
        }
        String str = "Scan: " + scan + " Peptide: " + peptide +
                " char[]: " + new String(AAs) +
                " position of variations: " + posOfVar +
                " start: " + start + " end: " + end;
        return str;
    }
}

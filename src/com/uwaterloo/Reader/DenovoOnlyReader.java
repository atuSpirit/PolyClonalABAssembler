package com.uwaterloo.Reader;

import com.uwaterloo.DenovoAssembler.DenovoOnly;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/* Read in the denovo only.csv file */
public class DenovoOnlyReader extends CSVReader{
    HashMap<String, Integer> fieldIndexMap;

    public DenovoOnlyReader() {

    }

    public List<DenovoOnly> readCSVFile(String psmFile) {
        List<DenovoOnly> fieldExtracted = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(psmFile))) {
            String titleString = br.readLine();
            this.fieldIndexMap = parseTitle(titleString);

            String line;
            while ((line = br.readLine()) != null) {
                DenovoOnly dnOnly = readOneLine(line);
                if (dnOnly != null) {
                    fieldExtracted.add(dnOnly);
                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fieldExtracted;
    }

    private DenovoOnly readOneLine(String line) {
        String[] fields = line.split(",");

        String scan = fields[fieldIndexMap.get("Scan")];
        String peptide = fields[fieldIndexMap.get("Peptide")];
        short alc = Short.valueOf(fields[fieldIndexMap.get("ALC (%)")]);
        short length = Short.valueOf(fields[fieldIndexMap.get("length")]);
        String confString = fields[fieldIndexMap.get("local confidence (%)")];

        int intensity = 0;
        if (!fields[fieldIndexMap.get("Area")].equals("")) {
            intensity = Double.valueOf(fields[fieldIndexMap.get("Area")]).intValue();
        }
        return new DenovoOnly(scan, peptide, alc, length, confString, intensity);
    }

    /**
     * Fragment dn into several parts if there is AA score below threshold. Return the fragment with largest length.
     * In this way, the dn seq is truncated those ends where the conf score is below the confScoreThresh to form new denovo tag
     * @param dn
     * @param confThresh Normally use 50
     * @return
     */
    public DenovoOnly truncateDnAAByConfScore(DenovoOnly dn, int confThresh) {
        List<DenovoOnly> dnFragmentList = new ArrayList<>();
        short[] confs = dn.getConfScores();
        char[] AAs = dn.getAAs();
        List<Character> fragAAs = new ArrayList<>();
        List<Short> fragConfs = new ArrayList<>();
        int scoreSum = 0;

        List<Character> longestFragAAs = new ArrayList<>();
        List<Short> longestFragConfs = new ArrayList<>();
        int maxScoreSum = 0;

        for (int i = 0; i < AAs.length; i++) {
            if (confs[i] < confThresh) {
                if (scoreSum > maxScoreSum) {
                    longestFragAAs = fragAAs;
                    longestFragConfs = fragConfs;
                    maxScoreSum = scoreSum;
                }
                scoreSum = 0;
                fragAAs = new ArrayList<>();
                fragConfs = new ArrayList<>();
            } else {
                fragAAs.add(AAs[i]);
                fragConfs.add(confs[i]);
                scoreSum += confs[i];
            }
        }
        if (scoreSum > maxScoreSum) {
            longestFragAAs = fragAAs;
            longestFragConfs = fragConfs;
            maxScoreSum = scoreSum;
        }
        int length = longestFragAAs.size();
        AAs = new char[length];
        confs = new short[length];
        for (int i = 0; i < length; i++) {
            AAs[i] = longestFragAAs.get(i);
            confs[i] = longestFragConfs.get(i);
        }
        DenovoOnly truncatedDn = new DenovoOnly(dn.getScan(), dn.getPeptide(), dn.getAlc(), (short) length,
                AAs, confs, dn.getIntensity());
        return truncatedDn;
    }

    public List<DenovoOnly> new_filterDnByConfScore(List<DenovoOnly> dnList, int confAAThresh, short kmerSize) {
        List<DenovoOnly> newDnList = new ArrayList<>();
        for (DenovoOnly dn : dnList) {
            DenovoOnly truncatedDn = truncateDnAAByConfScore(dn, confAAThresh);
            if (truncatedDn.getLength() >= kmerSize) {
                newDnList.add(truncatedDn);
            }
        }
        return newDnList;
    }

}

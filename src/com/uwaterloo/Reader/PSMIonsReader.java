package com.uwaterloo.Reader;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import org.junit.Test;
import static org.junit.Assert.*;

public class PSMIonsReader extends CSVReader {
    HashMap<String, Integer> fieldIndexMap;

    public PSMIonsReader() {

    }

    /**
     * Read ion positions for each scan from PSM ions.csv file.     *
     * @param psmIonsFile
     * @return a hashMap of scan and a list of {0,1} denoting whether
     * there is ions at the position.
     */
    public HashMap<String, short[]> readPSMIonsFile(String psmIonsFile) {
        HashMap<String, short[]> scanIonsMap = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(psmIonsFile))) {
            String titleString = br.readLine();
            this.fieldIndexMap = parseTitle(titleString);

            String line;
            while ((line = br.readLine()) != null) {
                String[] fields = line.split(",");

                String scan = fields[fieldIndexMap.get("scan")];
                String peptide = fields[fieldIndexMap.get("Peptide")];
                int pos = Integer.valueOf(fields[fieldIndexMap.get("pos")]) - 1;    //starting from zero
                if (scanIonsMap.containsKey(scan)) {
                    scanIonsMap.get(scan)[pos] = 1;
                } else {
                    short[] poses = new short[peptide.length()];
                    poses[pos] = 1;
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return scanIonsMap;
    }

    /**
     * For each scan, transfer the ions to ion Scores.
     * @param scanIonsMap
     * @return a map of scan and its ionScores.
     */
    public HashMap<String, short[]> setIonsScore(HashMap<String, short[]> scanIonsMap) {
        HashMap<String, short[]> scanIonScoresMap = new HashMap<>();
        for (String scan : scanIonsMap.keySet()) {
            short[] ions = scanIonsMap.get(scan);
            short[] ionsScores = computeIonsScoresFromIonPos(ions);
            scanIonScoresMap.put(scan, ionsScores);
        }
        return scanIonScoresMap;
    }

    /**
     * Given the binary array of ions, compute the ion scores.
     * If there is one or more zero, the score is 100 / (length of zero + 1).
     * @param ions
     * @return
     */
    private short[] computeIonsScoresFromIonPos(short[] ions) {
        short[] ionsScores = new short[ions.length];
        int zeroStart = -1;
        int zeroEnd = 0;
        for (int i = 0; i < ions.length; i++) {
            if (ions[i] == 1) {
                if (zeroStart != -1) {
                    zeroEnd = i;
                    short score = (short) (100 / (zeroEnd - zeroStart + 1));
                    for (int j = zeroStart; j <= zeroEnd; j++) {
                        ionsScores[j] = score;
                    }
                    zeroStart = -1;
                } else {
                    ionsScores[i] = 100;
                }
            } else {
                if (zeroStart == -1) {
                    zeroStart = i;
                }
            }
        }
        if (zeroStart != -1) {
            zeroEnd = ions.length - 1;
            short score = (short) (100 / (zeroEnd - zeroStart + 1));
            for (int j = zeroStart; j <= zeroEnd; j++) {
                ionsScores[j] = score;
            }
        }
        return ionsScores;
    }

    @Test
    public void testComputeIonsScoresFromIonPos() {
        short[] ions = new short[]{0,1,1,1,0,0,0};

        short[] ionsScore = computeIonsScoresFromIonPos(ions);
        short[] trueScore = new short[]{50, 50, 100, 100, 33, 33, 33};
        assertArrayEquals(trueScore, ionsScore);

        ions = new short[]{0,0,0};
        trueScore = new short[]{33,33,33};
        ionsScore = computeIonsScoresFromIonPos(ions);
        assertArrayEquals(trueScore, ionsScore);

        ions = new short[]{1,1,1,1};
        trueScore = new short[]{100, 100, 100, 100};
        ionsScore = computeIonsScoresFromIonPos(ions);
        assertArrayEquals(trueScore, ionsScore);

        ions = new short[]{0,1,1,0,1,0,0};
        trueScore = new short[]{50, 50, 100, 50, 50, 50, 50};
        ionsScore = computeIonsScoresFromIonPos(ions);
        assertArrayEquals(trueScore, ionsScore);
    }

}

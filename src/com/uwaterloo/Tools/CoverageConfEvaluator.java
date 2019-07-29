package com.uwaterloo.Tools;

import com.uwaterloo.PSMAligned;
import com.uwaterloo.TemplateHooked;

import java.util.List;

/**
 * Evaluate the confidence of each AA on a template with DB mapped to it.
 */
public class CoverageConfEvaluator {
    public CoverageConfEvaluator() {

    }

    public int[] evaluateCoverageConf(TemplateHooked templateHooked) {
        int length = templateHooked.getSeq().length;
        int[] coverageConfs = new int[length];
        for (int i = 0; i < length; i++) {
            List<PSMAligned> dbList = templateHooked.getDbList().get(i);
            for (PSMAligned psmAligned : dbList) {
                int start = psmAligned.getStart();
                int end = psmAligned.getEnd();
                short[] ionScores = psmAligned.getIonScores();
                for (int j = start; j <= end; j++) {
                    coverageConfs[j] += ionScores[j - start];
                }
            }
        }
        return coverageConfs;
    }

    public int sumCoverageConf(int[] coverageConfs) {
        int sum = 0;
        for (int conf : coverageConfs) {
            sum += conf;
        }
        return sum;
    }


}

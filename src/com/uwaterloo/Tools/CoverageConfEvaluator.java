package com.uwaterloo.Tools;

import com.uwaterloo.ScanTemplateMapper.PSMAligned;
import com.uwaterloo.ScanTemplateMapper.TemplateHooked;

import java.util.List;

/**
 * Evaluate the confidence of each AA on a template with DB mapped to it.
 */
public class CoverageConfEvaluator {
    public CoverageConfEvaluator() {

    }

    /**
     * Coverage only consider db result, other than the spider result
     * @param templateHooked
     * @return
     */
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

    /**
     * The confidence score is the number of full fragmented ions covered in this position.
     * @param templateHooked
     * @return
     */
    public int[] evaluateFragmentationConf(TemplateHooked templateHooked) {
        int length = templateHooked.getSeq().length;
        int[] coverageConfs = new int[length];
        for (int i = 0; i < length; i++) {
            List<PSMAligned> dbList = templateHooked.getDbList().get(i);
            for (PSMAligned psmAligned : dbList) {
                int start = psmAligned.getStart();
                int end = psmAligned.getEnd();
                short[] ionScores = psmAligned.getIonScores();
                for (int j = start; j <= end; j++) {
                    if (ionScores[j - start] == 100) {
                        coverageConfs[j] += 1;
                    }

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

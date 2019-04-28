package com.uwaterloo.Tools;

import java.util.Comparator;

public class TemplateStatistics {
    int templateId;
    int psmNum;
    int coveredAANum;
    int confidentAANum;
    int moreConfidentAANum;
    int scoreSum;   //The summation of all the AA confScore

    public TemplateStatistics(int templateId, int psmNum, int coveredAANum, int confidentAANum,
                              int moreConfidentAANum, int scoreSum) {
        this.templateId = templateId;
        this.psmNum = psmNum;
        this.coveredAANum = coveredAANum;
        this.confidentAANum = confidentAANum;
        this.moreConfidentAANum = moreConfidentAANum;
        this.scoreSum = scoreSum;
    }

    public int getTemplateId() {
        return templateId;
    }

    public int getPsmNum() {
        return psmNum;
    }

    public int getCoveredAANum() {
        return coveredAANum;
    }

    public int getConfidentAANum() {
        return confidentAANum;
    }

    public int getMoreConfidentAANum() {
        return moreConfidentAANum;
    }

    public int getScoreSum() {
        return scoreSum;
    }

    @Override
    public String toString() {
        return psmNum + "\t" + coveredAANum + "\t" + confidentAANum + "\t" + moreConfidentAANum + "\t" + scoreSum;
    }

    public static Comparator<TemplateStatistics> cmpReverseConfidentAANum() {
        return new Comparator<TemplateStatistics>() {
            @Override
            public int compare(TemplateStatistics o1, TemplateStatistics o2) {
                return o2.getConfidentAANum() - o1.getConfidentAANum();
            }
        };
    }

    public static Comparator<TemplateStatistics> cmpReverseScoreSum() {
        return new Comparator<TemplateStatistics>() {
            @Override
            public int compare(TemplateStatistics o1, TemplateStatistics o2) {
                return o2.getScoreSum() - o1.getScoreSum();
            }
        };
    }
}

package com.uwaterloo;

import java.util.Comparator;

public class DenovoAligned implements Comparable<DenovoAligned> {
    int templateId;
    int tStart;
    int tEnd;
    String dnScan;
    int dnStart;
    int dnEnd;
    int score;


    public DenovoAligned(int templateId, int tStart, int tEnd,
                         String dnScan, int dnStart, int dnEnd, int score) {
        this.templateId = templateId;
        this.tStart = tStart;
        this.tEnd = tEnd;
        this.dnScan = dnScan;
        this.dnStart = dnStart;
        this.dnEnd = dnEnd;
        this.score = score;
    }


    public int getTemplateId() {
        return templateId;
    }

    public int gettStart() {
        return tStart;
    }

    public int gettEnd() {
        return tEnd;
    }

    public String getDnScan() {
        return dnScan;
    }

    public int getDnStart() {
        return dnStart;
    }

    public int getDnEnd() {
        return dnEnd;
    }

    public int getScore() {
        return score;
    }


    @Override
    public int compareTo(DenovoAligned o) {
        return (this.score - o.score);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null) {
            return false;
        }
        DenovoAligned dnAligned = (DenovoAligned) o;
        if (templateId != dnAligned.templateId) {
            return false;
        }

        if (tStart != dnAligned.tStart) {
            return false;
        }
        if (tEnd != dnAligned.tEnd) {
            return false;
        }

        if (!dnScan.equals(dnAligned.dnScan)) {
            return false;
        }
        if (dnStart != dnAligned.dnStart) {
            return false;
        }
        if (dnEnd != dnAligned.dnEnd) {
            return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + templateId;
        result = prime * result + tStart;
        result = prime * result + tEnd;
        result = prime * result + dnStart;
        result = prime * result + dnEnd;
        result = prime * result;
        char[] chars = dnScan.toCharArray();
        for (char c : chars) {
            result += c;
        }

        return result;
    }

    @Override
    public String toString() {
        String dnAlignStr = templateId + " " + tStart + " " + tEnd + " " +
                dnScan + " " + dnStart + " " + dnEnd + " " + score;
        return dnAlignStr;
    }

    public static Comparator<DenovoAligned> cmpTStart() {
        return new Comparator<DenovoAligned>() {
            @Override
            public int compare(DenovoAligned o1, DenovoAligned o2) {
                return o1.gettStart() - o2.gettStart();
            }
        };
    }

    //Sort according to tEnd descending.
    public static Comparator<DenovoAligned> cmpReverseTEnd() {
        return new Comparator<DenovoAligned>() {
            @Override
            public int compare(DenovoAligned o1, DenovoAligned o2) {
                return o2.gettEnd() - o1.gettEnd();
            }
        };
    }

}



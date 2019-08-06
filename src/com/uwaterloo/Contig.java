package com.uwaterloo;

import java.util.Comparator;

/* The structure of contig assembled from denovo only peptides */
public class Contig {
    int tStart; //The start position on the template for contig extend to right
    int tEnd;   //The start position on the template for contig extend to left
    char[] AAs;
    int[] confs;
    int score;

    public Contig(int tStart, int tEnd, char[] AAs, int[] confs, int score) {
        this.tStart = tStart;
        this.tEnd = tEnd;
        this.AAs = AAs;
        this.confs = confs;
        this.score = score;
    }

    public int gettStart() {
        return tStart;
    }

    public int gettEnd() {
        return tEnd;
    }

    public char[] getAAs() {
        return AAs;
    }

    public int[] getConfs() { return confs;}

    public int getScore() {
        return score;
    }

    public void settStart(int tStart) {
        this.tStart = tStart;
    }

    public void settEnd(int tEnd) {
        this.tEnd = tEnd;
    }

    public void setAAs(char[] AAs) {
        this.AAs = AAs;
    }

    public void setConfs(int[] confs) { this.confs = confs;}

    public void setScore(int score) {
        this.score = score;
    }

    public static Comparator<Contig> cmpScore() {
        return new Comparator<Contig>() {
            @Override
            public int compare(Contig o1, Contig o2) {
                return o1.getScore() - o2.getScore();
            }
        };
    }

    public static Comparator<Contig> cmpReverseScore() {
        return new Comparator<Contig>() {
            @Override
            public int compare(Contig o1, Contig o2) {
                return o2.getScore() - o1.getScore();
            }
        };
    }

    public static Comparator<Contig> cmpTStart() {
        return new Comparator<Contig>() {
            @Override
            public int compare(Contig o1, Contig o2) {
                return o1.gettStart() - o2.gettStart();
            }
        };
    }

    @Override
    public String toString() {
        return tStart + " " + tEnd + " " + new String(AAs) + " " + score;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || o.getClass() != this.getClass()) {
            return false;
        }
        Contig contig = (Contig) o;
        return (contig.gettStart() == this.gettStart()) &&
                (contig.gettEnd() == this.gettEnd()) &&
                (contig.getScore() == this.getScore()) &&
                (java.util.Arrays.equals(contig.getAAs(), this.getAAs()));
    }
}

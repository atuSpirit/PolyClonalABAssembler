package com.uwaterloo;

import java.util.List;

public class MutationsPattern implements Comparable {
    List<Integer> posList;
    String AAs;
    int freq;
    int score;

    public MutationsPattern(List<Integer> posList, String AAs, int freq, int score) {
        this.posList = posList;
        this.AAs = AAs;
        this.freq = freq;
        this.score = score;

    }

    public List<Integer> getPosList() {
        return posList;
    }

    public String getAAs() {
        return AAs;
    }

    public int getFreq() {
        return freq;
    }

    public int getScore() {
        return score;
    }

    public void setPosList(List<Integer> posList) {
        this.posList = posList;
    }

    public void setAAs(String AAs) {
        this.AAs = AAs;
    }

    public void setFreq(int freq) {
        this.freq = freq;
    }

    public void setScore(int score) {
        this.score = score;
    }

    @Override
    public String toString() {
        return AAs + " at " +  posList.toString() + " freq: " + freq + " score: " + score;
    }

    @Override
    public int compareTo(Object o) {
        if (!this.getAAs().equals(((MutationsPattern) o).getAAs())) {
            return -1;
        }
        if (this.getPosList().equals(((MutationsPattern) o).getPosList())) {
            return 0;
        } else {
            return -1;
        }
    }
}

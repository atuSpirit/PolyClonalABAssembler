package com.uwaterloo;

import java.util.List;

public class MutationsPattern {
    List<Integer> posList;
    String AAs;
    int freq;

    public MutationsPattern(List<Integer> posList, String AAs, int freq) {
        this.posList = posList;
        this.AAs = AAs;
        this.freq = freq;
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

    public void setPosList(List<Integer> posList) {
        this.posList = posList;
    }

    public void setAAs(String AAs) {
        this.AAs = AAs;
    }

    public void setFreq(int freq) {
        this.freq = freq;
    }

    @Override
    public String toString() {
        return AAs + " at " +  posList.toString() + " freq: " + freq;
    }
}

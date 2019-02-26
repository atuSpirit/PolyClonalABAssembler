package com.uwaterloo;

import org.junit.Test;
import org.junit.Assert.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MutationsPattern {
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
    public boolean equals(Object o) {
        if (!this.getAAs().equals(((MutationsPattern) o).getAAs())) {
            return false;
        }
        if (this.getPosList().equals(((MutationsPattern) o).getPosList())) {
            return true;
        }

        return false;
    }

    @Override
    public int hashCode() {
        int result = 7;
        result = 31 * result + posList.toString().hashCode();
        result = 31 * result + AAs.hashCode();
        return result;
    }

    public static void main(String[] args) {
        List<Integer> posList1 = new ArrayList<>();
        posList1.add(18);
        posList1.add(27);

        List<Integer> posList2 = new ArrayList<>();
        posList2.add(18);
        posList2.add(27);

        MutationsPattern pattern1 = new MutationsPattern(posList1,"DG",3, 26);
        MutationsPattern pattern2 = new MutationsPattern(posList2, "DG", 2, 28);

        HashMap<MutationsPattern, Integer> map = new HashMap<>();
        map.put(pattern1, 1);
        if (map.containsKey(pattern2)) {
            System.out.println("same");
        }

        boolean isEqual = (pattern1.equals(pattern2));
        System.out.println(isEqual);

    }

}

package com.uwaterloo.Test;

import com.uwaterloo.PSM;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class PSMTest {
    @Test
    public void testPositionOfVariations() {
        String peptide = "VL(sub T)VSS(+15.02)ASTKGPSVF";
        String scan = "F7:5719";
        PSM psm = new PSM(scan, peptide);

        List<Integer> truePosList = new ArrayList<>();
        truePosList.add(1);
        assertArrayEquals(truePosList.toArray(), psm.getPositionOfVariations().toArray());
        assertArrayEquals("VLVSSASTKGPSVF".toCharArray(), psm.getAAs());

        peptide = "M(sub T)NQVSLTCLVK";
        psm = new PSM(scan, peptide);

        List<Integer> truePosList2 = new ArrayList<>();
        truePosList2.add(0);
        assertArrayEquals(truePosList2.toArray(), psm.getPositionOfVariations().toArray());
        assertArrayEquals("MNQVSLTCLVK".toCharArray(), psm.getAAs());

        peptide = "M(ins)NQVSLTCLVK";
        psm = new PSM(scan, peptide);
        assertArrayEquals(truePosList2.toArray(), psm.getPositionOfVariations().toArray());
        assertArrayEquals("MNQVSLTCLVK".toCharArray(), psm.getAAs());
    }
}

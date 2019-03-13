package com.uwaterloo.Test;


import com.uwaterloo.PSMAligned;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class PSMAlignedTest {
    @Test
    public void testPositionOfVariations() {
        String peptide = "VL(sub T)VSS(+15.02)ASTKGPSVF";
        String scan = "F7:5719";
        int templateId = 0;
        int start = 345;
        int end = 357;
        short[] ionScores = null;
        PSMAligned psmAligned = new PSMAligned(scan, peptide, 0, templateId, start, end, ionScores);

        List<Integer> truePosList = new ArrayList<>();
        truePosList.add(1);
        assertArrayEquals(truePosList.toArray(), psmAligned.getPositionOfVariations().toArray());
        assertArrayEquals("VLVSSASTKGPSVF".toCharArray(), psmAligned.getAAs());

        peptide = "M(sub T)NQVSLTCLVK";
        psmAligned = new PSMAligned(scan, peptide, 0, templateId, start, end, ionScores);

        List<Integer> truePosList2 = new ArrayList<>();
        truePosList2.add(0);
        assertArrayEquals(truePosList2.toArray(), psmAligned.getPositionOfVariations().toArray());
        assertArrayEquals("MNQVSLTCLVK".toCharArray(), psmAligned.getAAs());

        peptide = "M(ins)NQVSLTCLVK";
        psmAligned = new PSMAligned(scan, peptide, 0, templateId, start, end, ionScores);
        assertArrayEquals(truePosList2.toArray(), psmAligned.getPositionOfVariations().toArray());
        assertArrayEquals("MNQVSLTCLVK".toCharArray(), psmAligned.getAAs());

        peptide = "EM(sub V)S(sub Q)LVESGGGV(sub L)VKPR(sub G)GSL";
        psmAligned = new PSMAligned("F3:12304", peptide, 0, templateId, 19, 36, ionScores);
        List<Integer> truePosList3 = new ArrayList<>();
        truePosList3.add(1);
        truePosList3.add(2);
        truePosList3.add(10);
        truePosList3.add(14);
        assertArrayEquals(truePosList3.toArray(), psmAligned.getPositionOfVariations().toArray());

    }
}

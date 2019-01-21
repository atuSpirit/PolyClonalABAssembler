package com.uwaterloo.Test;

import com.uwaterloo.PSM;
import com.uwaterloo.PSMAligned;
import com.uwaterloo.TMapPosition;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class PSMAlignedTest {
    @Test
    public void testPositionOfVariations() {
        String peptide = "VL(sub T)VSS(+15.02)ASTKGPSVF";
        String scan = "F7:5719";
        List<TMapPosition> tMapPositionList = new ArrayList<>();
        tMapPositionList.add(new TMapPosition(0, 345));
        PSMAligned psmAligned = new PSMAligned(scan, peptide, tMapPositionList);

        List<Integer> truePosList = new ArrayList<>();
        truePosList.add(1);
        assertArrayEquals(truePosList.toArray(), psmAligned.getPositionOfVariations().toArray());
        assertArrayEquals("VLVSSASTKGPSVF".toCharArray(), psmAligned.getAAs());

        peptide = "M(sub T)NQVSLTCLVK";
        psmAligned = new PSMAligned(scan, peptide, tMapPositionList);

        List<Integer> truePosList2 = new ArrayList<>();
        truePosList2.add(0);
        assertArrayEquals(truePosList2.toArray(), psmAligned.getPositionOfVariations().toArray());
        assertArrayEquals("MNQVSLTCLVK".toCharArray(), psmAligned.getAAs());

        peptide = "M(ins)NQVSLTCLVK";
        psmAligned = new PSMAligned(scan, peptide, tMapPositionList);
        assertArrayEquals(truePosList2.toArray(), psmAligned.getPositionOfVariations().toArray());
        assertArrayEquals("MNQVSLTCLVK".toCharArray(), psmAligned.getAAs());
    }
}

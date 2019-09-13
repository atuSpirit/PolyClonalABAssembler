package com.uwaterloo.Test;
import com.uwaterloo.DenovoAssembler.UncertainRegionAssembler;
import org.junit.Test;

import static org.junit.Assert.*;

public class UncertainRegionAssemblerTest {
    @Test
    public void testFindOverlapIndex() {
        UncertainRegionAssembler assembler = new UncertainRegionAssembler();
        char[] seq1 = "VTNMDPANTATYYCARDMLF".toCharArray();
        char[] seq2 = "DMLFDFYFDVWGQGTTVTVSSASTK".toCharArray();

        int trueOverlapIndex = 16;
        assertEquals(trueOverlapIndex, assembler.findOverlapIndex(seq1, seq2));

    }
}

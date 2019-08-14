package com.uwaterloo.Test;
import com.uwaterloo.DenovoAligned;
import com.uwaterloo.DenovoOnly;
import com.uwaterloo.UncertainRegionAssembler;
import org.junit.Before;
import org.junit.Test;

import java.util.HashMap;
import java.util.List;

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

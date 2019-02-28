package com.uwaterloo.Test;

import com.uwaterloo.TemplateCandidateBuilder;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;
import static org.junit.Assert.*;

public class TemplateCandidateBuilderTest {
    @Test
    public void testExtractAAs() {
        TemplateCandidateBuilder builder = new TemplateCandidateBuilder(null);
        char[] AAs = "ACTGACTGATG".toCharArray();
        List<Integer> posArray = new ArrayList<>();
        posArray.add(18);
        posArray.add(19);
        posArray.add(27);
        List<Integer> subPosList = new ArrayList<>();
        subPosList.add(18);
        subPosList.add(27);

        /*
        String subAAs = builder.extractAAs(subPosList, posArray, 18, AAs);
        assertEquals("AT", subAAs);

        posArray.remove(2);
        subAAs = builder.extractAAs(subPosList, posArray, 18, AAs);
        assertEquals("A", subAAs);
        */
    }
}

package com.uwaterloo.Test;

import com.uwaterloo.DenovoAssembler.DenovoOnly;
import com.uwaterloo.Reader.DenovoOnlyReader;
import org.junit.Test;

public class DenovoOnlyReaderTest {
    @Test
    public void testTruncateDnAAByConfScore() {
        DenovoOnlyReader dnReader = new DenovoOnlyReader();
        short alc = 99;
        short len = 17;
        DenovoOnly dn = new DenovoOnly("F4:14380", "QVTLRESGPALVKPTVFPTL", alc, len, "59 50 65 97 99 98 96 69 67 88 99 97 51 46 54 31 33 74 80 93", 6720000);
        DenovoOnly truncatedDn = dnReader.truncateDnAAByConfScore(dn, 50);

        System.out.println("truncatedDn: " + truncatedDn.toString());

    }
}

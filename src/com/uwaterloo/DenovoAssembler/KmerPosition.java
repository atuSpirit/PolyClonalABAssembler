package com.uwaterloo.DenovoAssembler;

/**
 * A structure to store kmers on denovo only peptide
 */
public class KmerPosition {
    int index;  //The index of the denovo peptide in denovoList
    //String scan;   //The Id of denovo only
    short pos;  //The starting position of the kmer on this denovo only peptide

    public KmerPosition(int index, short pos) {
    //    public KmerPosition(int index, String scan, short pos) {
        this.index = index;
    //    this.scan = scan;
        this.pos = pos;
    }

    public int getIndex() {
        return index;
    }

    /*
    public String getScan() {
        return scan;
    }
    */

    public short getPos() {
        return pos;
    }
}

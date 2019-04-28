package com.uwaterloo;

import java.util.Set;

/* The region has numOf(psm) / numOf(dn) less than a threshold */
public class UncertainRegion {
    int startPos;   //The start position on a template
    int endPos;
    Set<DenovoAligned> dnAlignToRightSet;
    Set<DenovoAligned> dnAlignToLeftSet;
    Set<PSMAligned> psmAlignedSet;

    public UncertainRegion(int startPos, int endPos, Set<DenovoAligned> dnAlignToRightSet,
                           Set<DenovoAligned> dnAlignToLeftSet, Set<PSMAligned> psmAlignedSet) {
        this.startPos = startPos;
        this.endPos = endPos;
        this.dnAlignToRightSet = dnAlignToRightSet;
        this.dnAlignToLeftSet = dnAlignToLeftSet;
        this.psmAlignedSet = psmAlignedSet;
    }


    public int getStartPos() {
        return startPos;
    }

    public int getEndPos() {
        return endPos;
    }

    public Set<DenovoAligned> getDnAlignToRightSet() {
        return dnAlignToRightSet;
    }

    public Set<DenovoAligned> getDnAlignToLeftSet() {
        return dnAlignToLeftSet;
    }

    public Set<PSMAligned> getPsmAlignedSet() {
        return psmAlignedSet;
    }

    public void setStartPos(int startPos) {
        this.startPos = startPos;
    }

    public void setEndPos(int endPos) {
        this.endPos = endPos;
    }

    public void setDnAlignToRightSet(Set<DenovoAligned> dnAlignToRightSet) {
        this.dnAlignToRightSet = dnAlignToRightSet;
    }

    public void setDnAlignToLeftSet(Set<DenovoAligned> dnAlignToLeftSet) {
        this.dnAlignToLeftSet = dnAlignToLeftSet;
    }
}

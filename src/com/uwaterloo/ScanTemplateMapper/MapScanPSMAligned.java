package com.uwaterloo.ScanTemplateMapper;

import java.util.HashMap;
import java.util.List;

/* Build a hash map between scan and its corresponding PSMAligned */
public class MapScanPSMAligned {
    HashMap<String, PSMAligned> scanPSMMap;

    public MapScanPSMAligned(List<PSMAligned> psmAlignedList) {
        buildScanPSMAlignedMap(psmAlignedList);
    }

    private void buildScanPSMAlignedMap(List<PSMAligned> psmAlignedList) {
        scanPSMMap = new HashMap<>();
        for (PSMAligned psmAligned : psmAlignedList) {
            String scan = psmAligned.getScan();
            scanPSMMap.put(scan, psmAligned);
        }
    }

    public HashMap<String, PSMAligned> getScanPSMMap() {
        return scanPSMMap;
    }
}

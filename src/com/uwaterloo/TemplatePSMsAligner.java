package com.uwaterloo;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

public class TemplatePSMsAligner {
    /* Each template contain a PSMAlignList storing psms mapped to it. */
    ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList;
    /**
     * Building the list of PSM aligned.
     * According to the psm's peptide sequence, find its mapping position
     * in proteinPeptideMap.
     * @param psmList   The psm list read from DB psm.csv
     * @param peptideProteinMap The HashMap<Peptide, listOfMapOnProtein> built from protein-peptide.csv
     * @return a list of aligned psms with its aligned position list added, a list of scans is stored for each position
     *          of the template.
     */
    private void buildListOfPSMAligned(int templateNum, List<PSM> psmList,
                                                   HashMap<String, List<TMapPosition>> peptideProteinMap) {

        listOfPSMAlignedList = new ArrayList<>();
        for (int i = 0; i < templateNum; i++) {
            ArrayList<PSMAligned> psmAlignedList = new ArrayList<>();
            listOfPSMAlignedList.add(psmAlignedList);
        }

        for (PSM psm : psmList) {
            String peptide = psm.getPeptide();
            List<TMapPosition> tMapPositionList = peptideProteinMap.get(peptide);
            for (TMapPosition tMapPosition : tMapPositionList) {
                int templateId = tMapPosition.getTemplateId();
                PSMAligned psmAligned = new PSMAligned(psm.getScan(), peptide, templateId,
                        tMapPosition.getStart(), tMapPosition.getEnd(), psm.ionScores);
                listOfPSMAlignedList.get(templateId).add(psmAligned);
            }
        }
    }


    /**
     * Scan each psm, hooked it to the correct position on the template.
     * @param template the template to hook
     * @param psmAlignedList the list of PSMs with the aligned position information
     * @return Each template will have a list of TemplatePosition
     */
    private TemplateHooked hookPSMsToTemplate(Template template,
                                                    List<PSMAligned> psmAlignedList) {
        TemplateHooked templateHooked = new TemplateHooked(template);

        /* Hook each psms to the template */
        for (PSMAligned psmAligned : psmAlignedList) {
            String scan = psmAligned.getScan();
            boolean isSpider = (psmAligned.getPositionOfVariations() != null);

            int start = psmAligned.getStart();

            /* Hook the psm information to corresponding dbList or spiderList */
            List<PSMAligned> psmListToBeHook = null;
            if (isSpider) {
                psmListToBeHook = templateHooked.getSpiderList().get(start);
            } else {
                psmListToBeHook = templateHooked.getDbList().get(start);
            }
            psmListToBeHook.add(psmAligned);

            /* For each position the peptide mapped, add the scan to its scanList */
            int end = psmAligned.getEnd();
            ArrayList<ArrayList<String>> mappedScanList = templateHooked.getMappedScanList();
            for (int i = start; i < end; i++) {
                mappedScanList.get(i).add(scan);
            }
        }
        return templateHooked;
    }

    public List<TemplateHooked> alignPSMstoTemplate(List<PSM> psmList, List<Template> templateList,
                                    HashMap<String, List<TMapPosition>> peptideProteinMap) {

        int templateNum = templateList.size();
        buildListOfPSMAligned(templateNum, psmList, peptideProteinMap);

        ArrayList<TemplateHooked> templateHookedList = new ArrayList<>();
        for (int i = 0; i < templateNum; i++) {
            TemplateHooked templateHooked = hookPSMsToTemplate(templateList.get(i),
                    listOfPSMAlignedList.get(i));
            templateHookedList.add(templateHooked);
        }

        return templateHookedList;
    }

    public ArrayList<ArrayList<PSMAligned>> getPsmAlignedList() {
        return listOfPSMAlignedList;
    }
}

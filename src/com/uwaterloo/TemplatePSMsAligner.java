package com.uwaterloo;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import static javax.swing.UIManager.get;

public class TemplatePSMsAligner {
    /**
     * Building the list of PSM aligned.
     * According to the psm's peptide sequence, find its mapping position
     * in proteinPeptideMap.
     * @param psmList   The psm list read from DB psm.csv
     * @param peptideProteinMap The HashMap<Peptide, listOfMapOnProtein> built from protein-peptide.csv
     * @return a list of aligned psms with its aligned position list added, a list of scans is stored for each position
     *          of the template.
     */
    private List<PSMAligned> buildListOfPSMAligned(List<PSM> psmList,
                                                   HashMap<String, List<TMapPosition>> peptideProteinMap) {
        List<PSMAligned> psmAlignedList = new ArrayList<>();
        for (PSM psm : psmList) {
            String peptide = psm.getPeptide();
            List<TMapPosition> tMapPositionList = peptideProteinMap.get(peptide);
            if (tMapPositionList == null) {
                System.out.println(psm.getScan() + " " + psm.getPeptide() + " contain no tMap position");
            }
            PSMAligned psmAligned = new PSMAligned(psm.getScan(), peptide, tMapPositionList);
            psmAlignedList.add(psmAligned);
        }
        return psmAlignedList;
    }


    /**
     * Scan each psm, hooked it to the correct templatePosition of its corresponding template.
     * @param templateList the list of templates
     * @param psmAlignedList the list of PSMs with the aligned position information
     * @return Each template will have a list of TemplatePosition
     */
    private List<TemplateHooked> hookPSMsToTemplate(List<Template> templateList,
                                                    List<PSMAligned> psmAlignedList) {
        List<TemplateHooked> listOfTemplateHooked = new ArrayList<>();
        for (Template template : templateList) {
            listOfTemplateHooked.add(new TemplateHooked(template));
        }

        /* Hook each psms to the corresponding template */
        for (PSMAligned psmAligned : psmAlignedList) {
            String scan = psmAligned.getScan();
            boolean isSpider = (psmAligned.getPositionOfVariations() != null);
            List<TMapPosition> tMapPositionList = psmAligned.getMapPositionList();
            for (TMapPosition tMapPosition : tMapPositionList) {
                int templateId = tMapPosition.getTemplateId();
                int start = tMapPosition.getStart();

                /* Hook the psm information to corresponding dbList or spiderList */
                List<PSMAligned> psmListToBeHook = null;
                if (isSpider) {
                    psmListToBeHook = listOfTemplateHooked.get(templateId).getSpiderList().get(start);
                } else {
                    psmListToBeHook = listOfTemplateHooked.get(templateId).getDbList().get(start);
                }
                psmListToBeHook.add(psmAligned);

                /* For each position the peptide mapped, add the scan to its scanList */
                int end = tMapPosition.getStart();
                ArrayList<LinkedList<String>> mappedScanList = listOfTemplateHooked.get(templateId).getMappedScanList();
                for (int i = start; i < end; i++) {
                    mappedScanList.get(i).add(scan);
                }
            }
        }
        return listOfTemplateHooked;
    }

    public List<TemplateHooked> alignPSMstoTemplate(List<PSM> psmList, List<Template> templateList,
                                    HashMap<String, List<TMapPosition>> peptideProteinMap) {
        List<PSMAligned> psmAlignedList = buildListOfPSMAligned(psmList, peptideProteinMap);
        List<TemplateHooked> templateHookedList = hookPSMsToTemplate(templateList, psmAlignedList);

        return templateHookedList;
    }

}

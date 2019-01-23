package com.uwaterloo;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class TemplatePSMsAligner {
    /**
     * Building the list of PSM aligned.
     * According to the psm's peptide sequence, find its mapping position
     * in proteinPeptideMap.
     * @param psmList   The psm list read from DB psm.csv
     * @param peptideProteinMap The HashMap<Peptide, listOfMapOnProtein> built from protein-peptide.csv
     * @return a list of aligned psms with its aligned position list added
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
            boolean isSpider = (psmAligned.getPositionOfVariations() != null);
            List<TMapPosition> tMapPositionList = psmAligned.getMapPositionList();
            for (TMapPosition tMapPosition : tMapPositionList) {
                int templateId = tMapPosition.getTemplateId();
                int start = tMapPosition.getStart();

                List<PSMAligned> psmListToBeHook = null;
                if (isSpider) {
                    psmListToBeHook = listOfTemplateHooked.get(templateId).getSpiderList().get(start);
                } else {
                    psmListToBeHook = listOfTemplateHooked.get(templateId).getDbList().get(start);
                }
                psmListToBeHook.add(psmAligned);
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

    /*
    private List<List<TemplatePosition>> hookPSMsToTemplate(int templateNum,
                                                        List<Template> templateList, List<PSMAligned> psmAlignedList) {
        List<List<TemplatePosition>> listOfTemplatePositionList = new ArrayList<List<TemplatePosition>>();
        // Initialize listOfTemplatePositionList for each template
        for (int i = 0; i < templateNum; i++) {
            listOfTemplatePositionList.add(new ArrayList<TemplatePosition>());
        }

        // Hook each psms to the corresponding template
        for (PSMAligned psmAligned : psmAlignedList) {
            boolean isSpider = (psmAligned.getPositionOfVariations() != null);
            List<TMapPosition> tMapPositionList = psmAligned.getMapPositionList();
            for (TMapPosition tMapPosition : tMapPositionList) {
                int templateId = tMapPosition.getTemplateId();
                int start = tMapPosition.getStart();
                Template template = templateList.get(templateId);
                char templateAA = template.getSeq()[start];

                TemplatePosition templatePosition = new TemplatePosition(start, templateAA);
                if (isSpider) {
                    templatePosition.getSpList().add(psmAligned);
                } else {
                    templatePosition.getDbList().add(psmAligned);
                }
                listOfTemplatePositionList.get(templateId).add(templatePosition);
            }
        }
        return listOfTemplatePositionList;
    }

*/

}

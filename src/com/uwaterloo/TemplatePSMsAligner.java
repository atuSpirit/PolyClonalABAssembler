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
     * @param proteinPeptideMap The HashMap<Peptide, listOfMapOnProtein> built from protein-peptide.csv
     * @return a list of aligned psms with its aligned position list added
     */
    private List<PSMAligned> buildListOfPSMAligned(List<PSM> psmList,
                                                   HashMap<String, List<TMapPosition>> proteinPeptideMap) {
        List<PSMAligned> psmAlignedList = new ArrayList<>();
        for (PSM psm : psmList) {
            String peptide = psm.getPeptide();
            List<TMapPosition> tMapPositionList = proteinPeptideMap.get(peptide);
            PSMAligned psmAligned = new PSMAligned(peptide, psm.getScan(), tMapPositionList);
            psmAlignedList.add(psmAligned);
        }
        return psmAlignedList;
    }


    /**
     * Scan each psm, hooked it to the correct templatePosition of its corresponding template.
     * @param templateList  the list of template with template AA information
     * @param psmAlignedList the list of PSMs with the aligned position information
     * @return Each template will have a list of TemplatePosition
     */
    private List<TemplateHooked> hookPSMsToTemplate(int templateNum, List<Template> templateList,
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

                if (isSpider) {
                    listOfTemplateHooked.get(templateId).getArrayOfSpiderAligned()[start] = psmAligned;
                } else {
                    listOfTemplateHooked.get(templateId).getArrayOfDBAligned()[start] = psmAligned;
                }
            }
        }
        return listOfTemplateHooked;
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

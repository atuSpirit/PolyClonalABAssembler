package com.uwaterloo;

import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

public class Assembler {
    public Assembler() {

    }

    public void process() {
        String psmFile = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\DB search psm.csv";
        PSMReader psmReader = new PSMReader();
        List<PSM> psmList = psmReader.readCSVFile(psmFile);

        String templateFasta = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\Nuno.2016.heavy.template.fasta";
        TemplatesLoader loader = new TemplatesLoader();
        List<Template> templateList = loader.loadTemplateFasta(templateFasta);

        ProteinPeptideReader ppReader = new ProteinPeptideReader(templateList);
        String proteinPeptideFile = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\protein-peptides.csv";
        HashMap<String, List<TMapPosition>> peptideProteinMap = ppReader.readProteinPeptideFile(proteinPeptideFile);

        TemplatePSMsAligner aligner = new TemplatePSMsAligner();
        List<TemplateHooked> templateHookedList = aligner.alignPSMstoTemplate(psmList, templateList, peptideProteinMap);


        int templateId = 0;
        /* Test one first template */
        TemplateHooked aTemplateHooked = templateHookedList.get(templateId);
        /*
        char[] seq = aTemplateHooked.getSeq();
        int size = seq.length;
        for (int i = 0; i < size; i++) {
            System.out.println(seq[i] + "\t" + aTemplateHooked.getDbList().get(i).size() + "\t" +
                    aTemplateHooked.getSpiderList().get(i).size() + "\t" +
                    aTemplateHooked.getMappedScanList().get(i).size());
        }
        */

        ArrayList<LinkedList<PSMAligned>> listOfPSMAlignedList = aligner.getPsmAlignedList();
        MapScanPSMAligned scanPSMMapper = new MapScanPSMAligned(listOfPSMAlignedList.get(templateId));
        HashMap<String, PSMAligned> scanPSMMap = scanPSMMapper.getScanPSMMap();

        MutationValidator validator = new MutationValidator();
        validator.validateMutations(aTemplateHooked, scanPSMMap);

    }

    public void testMutationValidator(TemplateHooked templateHooked) {
        MutationValidator validator = new MutationValidator();
        int maxMutationNum = validator.findMaxMutationNumPerPSM(templateHooked);
        System.out.println("Max mutation num: " + maxMutationNum);
        int mutationNum = 3;
        List<Integer> posWithMaxMutationNumList = validator.extractPositionWithMutationNum(templateHooked, mutationNum);

        for (Integer pos : posWithMaxMutationNumList) {
            List<PSMAligned> psmAlignedList = templateHooked.getSpiderList().get(pos);
            System.out.print("pos: " + pos);
            for (PSMAligned psmAligned : psmAlignedList) {
                if (psmAligned.getPositionOfVariations().size() == mutationNum) {
                    System.out.print(" " + psmAligned.getScan());
                }
            }
            System.out.println();

        }



    }
}

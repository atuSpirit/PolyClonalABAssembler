package com.uwaterloo;

import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;

import java.util.*;

public class Assembler {
    public Assembler() {

    }

    public void process() {
        String dir = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\";
        dir = "D:\\Hao\\data\\for_analysis\\PolyClonal_ab19001_SPIDER_12\\";
        String psmFile = dir + "DB search psm.csv";
        PSMReader psmReader = new PSMReader();
        List<PSM> psmList = psmReader.readCSVFile(psmFile);

        String templateFasta = dir + "Nuno.2016.heavy.template.fasta";
        templateFasta = dir + "ab19001.template_top8.fasta";
        TemplatesLoader loader = new TemplatesLoader();
        List<Template> templateList = loader.loadTemplateFasta(templateFasta);

        ProteinPeptideReader ppReader = new ProteinPeptideReader(templateList);
        String proteinPeptideFile = dir + "protein-peptides.csv";
        HashMap<String, List<TMapPosition>> peptideProteinMap = ppReader.readProteinPeptideFile(proteinPeptideFile);

        TemplatePSMsAligner aligner = new TemplatePSMsAligner();
        List<TemplateHooked> templateHookedList = aligner.alignPSMstoTemplate(psmList, templateList, peptideProteinMap);
        ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList = aligner.getPsmAlignedList();

        for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
            System.out.println("Template " + templateId);
            TemplateHooked aTemplateHooked = templateHookedList.get(templateId);
            findCandidateForOneTemplate(aTemplateHooked, templateId, listOfPSMAlignedList);
        }


    }

    private void findCandidateForOneTemplate(TemplateHooked aTemplateHooked, int templateId,
                                             ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList) {

        MapScanPSMAligned scanPSMMapper = new MapScanPSMAligned(listOfPSMAlignedList.get(templateId));
        HashMap<String, PSMAligned> scanPSMMap = scanPSMMapper.getScanPSMMap();

        MutationValidator validator = new MutationValidator();
        List<HashMap<List<Integer>, List<String>>> mutationsOnTemplateList = validator.validateMutations(aTemplateHooked, scanPSMMap);
        //  printMutationsOnTemplate(mutationsOnTemplateList);

        TemplateCandidateBuilder templateCandidateBuilder = new TemplateCandidateBuilder(mutationsOnTemplateList);
        templateCandidateBuilder.buildCandidateTemplate(aTemplateHooked);

    }

    private void printMutationsOnTemplate(List<HashMap<List<Integer>, List<String>>> mutationsOnTemplateList) {
        for (int i = 0; i < mutationsOnTemplateList.size(); i++) {
            HashMap<List<Integer>, List<String>> mutationsOnTemplate = mutationsOnTemplateList.get(i);
            for (Map.Entry entry : mutationsOnTemplate.entrySet()) {
                List<Integer> posList = (List<Integer>) entry.getKey();
                List<String> patternList = (List<String>) entry.getValue();
                for (int pos : posList) {
                    System.out.print(pos + " ");
                }
                System.out.println();
                for (String pattern : patternList) {
                    System.out.println(pattern);
                }

            }
        }

    }


    private void printMutationsSortedByPos(TreeMap<Integer, List<String>> mutationsSortedByPos) {
        for (Map.Entry entry : mutationsSortedByPos.entrySet()) {
            System.out.println(entry.getKey());
            List<String> patternList = (List<String>) entry.getValue();

            for (String pattern : patternList) {
                System.out.println(pattern);
            }

        }
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

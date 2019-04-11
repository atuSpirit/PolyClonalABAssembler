package com.uwaterloo;

import com.uwaterloo.Reader.*;

import java.io.*;
import java.util.*;

public class Assembler {
    public Assembler() {

    }

    public void process() {
        String dir = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\";
        dir = "D:\\Hao\\data\\for_analysis\\PolyClonal_ab19001_SPIDER_12\\";
        dir = "D:\\Hao\\result\\Waters_mAB_SPIDER_46\\";
        dir = "D:\\Hao\\result\\ab19001.5enzymes_SPIDER_17\\";
        dir = "D:\\Hao\\result\\ab19001.5enzymes.new_SPIDER_79\\";
        //dir = "/Users/hao/data/ab19001.5enzymes.new_SPIDER_33/";
        //dir = "D:\\Hao\\result\\Nuno2016_HC_SPIDER_70\\";
        String psmFile = dir + "DB search psm.csv";
        PSMReader psmReader = new PSMReader();
        List<PSM> psmList = psmReader.readCSVFile(psmFile);

        String psmIonsFile = dir + "PSM ions.csv";
        PSMIonsReader ionsReader = new PSMIonsReader();
        HashMap<String, short[]> scanIonPosesMap = ionsReader.readPSMIonsFile(psmIonsFile);
        HashMap<String, short[]> scanIonScoresMap = ionsReader.setIonsScore(scanIonPosesMap);
        setIonScoresForPSMList(psmList, scanIonScoresMap);

        String dnFile = dir + "de novo only peptides.csv";
        DenovoOnlyReader dnReader = new DenovoOnlyReader();
        List<DenovoOnly> dnList = dnReader.readCSVFile(dnFile);
        System.out.println(dnList.size());
        System.out.println(dnList.get(0).toString());

        String templateFasta = dir + "Nuno.2016.heavy.template.fasta";
        //templateFasta = dir + "Waters_mAB.template_top4.fasta";
        templateFasta = dir + "ab19001.5enzymes.template_top8.fasta";
        templateFasta = dir + "candidate_template2.mergedFrag.fasta";
        templateFasta = dir + "Nuno2016.HC_LC.template.top4.fasta";
        templateFasta = dir + "proteins.fasta";
        TemplatesLoader loader = new TemplatesLoader();
        List<Template> templateList = loader.loadTemplateFasta(templateFasta);

        ProteinPeptideReader ppReader = new ProteinPeptideReader(templateList);
        String proteinPeptideFile = dir + "protein-peptides.csv";
        HashMap<String, List<TMapPosition>> peptideProteinMap = ppReader.readProteinPeptideFile(proteinPeptideFile);

        //Map psms to templates
        TemplatePSMsAligner psmAligner = new TemplatePSMsAligner();
        List<TemplateHooked> templateHookedList = psmAligner.alignPSMstoTemplate(psmList, templateList, peptideProteinMap);
        ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList = psmAligner.getPsmAlignedList();

        //Map Denovo only peptides to templates
        TemplateDenovoAligner dnAligner = new TemplateDenovoAligner(dnList);
        short kmerSize = 6;
        dnAligner.alignDenovoOnlyToTemplate(templateHookedList, kmerSize);

        boolean useDenovo = false;
        if (!useDenovo) {
            //Generating candidate templates using DB and Spider PSMs
            generateCandidateTemplates(templateHookedList, listOfPSMAlignedList);
        } else {
            //Generating candidate templates using denovo only results
            float dbDnRatioThresh = 1.0f;
            dnAligner.majorityVoteInDnDominantRegion(templateHookedList, dbDnRatioThresh);;
        }

        //Export candidate templates together with contaminant sequences as a fasta file
        String contaminantFile = "D:\\Hao\\database\\contaminants.fasta";
        //contaminantFile = "/Users/hao/data/contaminants.fasta";
        String candidateTemplateWithContaminant = "D:\\Hao\\result\\database\\candidate_template_with_contaminant.fasta";
        //candidateTemplateWithContaminant = "/Users/hao/data/candidate_template_with_contaminant.fasta";
        int min_template_length = 0;  //If a template length is shorter than the min_length, don't output it.

        exportCandidateTemplates(templateHookedList, min_template_length, candidateTemplateWithContaminant, contaminantFile);
    }

    /* Export candidate templates together with contaminant sequences as a fasta file */
    private void exportCandidateTemplates(List<TemplateHooked> templateHookedList, int min_template_length,
                                          String candidateTemplateWithContaminant, String contaminantFile) {
        String seqsToExport = "";
        for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
            List<char[]> candidateTemplates = templateHookedList.get(templateId).getModifiedSeq();
            String templateAccession = templateHookedList.get(templateId).getTemplateAccession();
            for (int i = 0; i < candidateTemplates.size(); i++) {
                //Only export template longer than min_template_length to keep only heavy or light chain. Delete fragments.
                if (candidateTemplates.get(i).length < min_template_length) {
                    continue;
                }
                System.out.println(">can" + (i + 1) + "_" + templateAccession);
                System.out.println(new String(candidateTemplates.get(i)));
                seqsToExport += ">can" + (i + 1) + "_" + templateAccession + "\n";
                seqsToExport += new String(candidateTemplates.get(i)) + "\n";
            }
            //Debug
            //break;
        }

        //Attach the contaminant sequences
        String contaminantSeqs = null;
        try (BufferedReader br = new BufferedReader(new FileReader(contaminantFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                seqsToExport += line + "\n";
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        //Export the sequences to a file
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(candidateTemplateWithContaminant))) {
            bw.write(seqsToExport);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    /*
    Attach the candidate templates to the mutated templates in templateHookedList
     */
    private void generateCandidateTemplates(List<TemplateHooked> templateHookedList,
                                                    ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList) {
        for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
            System.out.println("Template " + templateId + " " + templateHookedList.get(templateId).getTemplateAccession());
            TemplateHooked aTemplateHooked = templateHookedList.get(templateId);
            List<char[]> top2CandidateTemplates = findCandidateForOneTemplate(aTemplateHooked, templateId, listOfPSMAlignedList);
            templateHookedList.get(templateId).setModifiedTemplates(top2CandidateTemplates);
            //Debug
            //break;
        }
    }

    /**
     * For each psm in psmList, try to fill in the fragment ion scores.
     * @param psmList
     * @param scanIonScoresMap
     */
    public static void setIonScoresForPSMList(List<PSM> psmList, HashMap<String,short[]> scanIonScoresMap) {
        for (int i = 0; i < psmList.size(); i++) {
            String scan = psmList.get(i).getScan();
            if (scanIonScoresMap.containsKey(scan)) {
                psmList.get(i).setIonScores(scanIonScoresMap.get(scan));
            } else {
                System.err.println("scan " + scan + " does not contain fragment ions information!");
            }
        }
    }
    /* Trim the C end of candidate template if no reads covered the C end */
    private void trimTemplateCEnd(TemplateHooked aTemplateHooked, List<char[]> topCandidateTemplates) {
        int pos = 0;
        while (pos < aTemplateHooked.getSeq().length) {
            if (aTemplateHooked.getMappedScanList().get(pos).size() > 0) {
                break;
            }
            pos++;
        }

        if (pos == 0) {
            return;
        }
        for (int i = 0; i < topCandidateTemplates.size(); i++) {
            int length = topCandidateTemplates.get(i).length - pos;
            char[] newCandidate = new char[length];
            for (int start = 0; start < length; start++) {
                newCandidate[start] = topCandidateTemplates.get(i)[start + pos];
            }
            topCandidateTemplates.set(i, newCandidate);
        }
    }

    private List<char[]> findCandidateForOneTemplate(TemplateHooked aTemplateHooked, int templateId,
                                             ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList) {

        MapScanPSMAligned scanPSMMapper = new MapScanPSMAligned(listOfPSMAlignedList.get(templateId));
        HashMap<String, PSMAligned> scanPSMMap = scanPSMMapper.getScanPSMMap();

        MutationValidator validator = new MutationValidator();
        List<HashMap<List<Integer>, List<MutationsPattern>>> mutationsOnTemplateList = validator.findSignificantMutations(aTemplateHooked, scanPSMMap);
        //  printMutationsOnTemplate(mutationsOnTemplateList);

        TemplateCandidateBuilder templateCandidateBuilder = new TemplateCandidateBuilder(mutationsOnTemplateList);
        List<char[]> topCandidateTemplates = templateCandidateBuilder.buildCandidateTemplate(aTemplateHooked, scanPSMMap);

//        trimTemplateCEnd(aTemplateHooked, topCandidateTemplates);

        return topCandidateTemplates;

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

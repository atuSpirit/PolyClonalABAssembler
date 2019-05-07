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
        dir = "D:\\Hao\\result\\ab19001.5enzymes.new_SPIDER_91\\";
        dir = "D:\\Hao\\result\\ab19001.5enzymes.4tempaltes_SPIDER_86\\";

        dir = "/Users/hao/data/ab19001.polyclonal.templateSelected_SPIDER_37/";

        dir = "C:\\Hao\\result\\ab19001.polyclonal.05.05_SPIDER_13\\";
        //dir = "C:\\Hao\\result\\ab19001.polyclonal.templateSelected_SPIDER_37\\";
        //dir = "C:\\Hao\\result\\Nuno.HC_SPIDER_19\\";
        //dir = "C:\\Hao\\result\\Nuno.LC_SPIDER_18\\";

        //dir = "/Users/hao/data/ab19001.5enzymes.new_SPIDER_33/";
        //dir = "D:\\Hao\\result\\Nuno2016_HC_SPIDER_66\\";
        //dir = "D:\\Hao\\result\\Water_mAB.clean_SPIDER_20\\";
        //dir = "D:\\Hao\\result\\Water_mAB.clean_PEAKS_19\\";
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
        System.out.println("Read in " + dnList.size() + " denovo only.");
        int confScoreThresh = 20;
        int inConfidentAANumThresh = 2;
        dnList = dnReader.filterDnByConfScore(dnList, confScoreThresh, inConfidentAANumThresh);
        System.out.println("filtered denovo size: " + dnList.size());

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

        //Transfer list of PSMAligned to list of scanPSMMap for each template
        List<HashMap<String, PSMAligned>> listOfScanPSMMap = new ArrayList<>();
        for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
            MapScanPSMAligned scanPSMMapper = new MapScanPSMAligned(listOfPSMAlignedList.get(templateId));
            HashMap<String, PSMAligned> scanPSMMap = scanPSMMapper.getScanPSMMap();
            listOfScanPSMMap.add(scanPSMMap);
        }

        HashMap<String, DenovoOnly> scanDnMap = buildScanDnMap(dnList);
        //Map Denovo only peptides to templates
        TemplateDenovoAligner dnAligner = new TemplateDenovoAligner(dnList, scanDnMap);
        short kmerSize = 6;
        dnAligner.alignDenovoOnlyToTemplate(templateHookedList, kmerSize);

        double significantThreshold = 0.2;  //Increase from 0.1 to 0.2 for Nuno data which is less accurate

        boolean useDenovo = false;
        if (!useDenovo) {
            //Generating candidate templates using DB and Spider PSMs
            generateCandidateTemplates(templateHookedList, listOfScanPSMMap, significantThreshold);
        } else {
            //Generating candidate templates using denovo only results
            float dbDnRatioThresh = 1.5f;
            UncertainRegionAssembler uncertainRegionAssembler = new UncertainRegionAssembler();
            uncertainRegionAssembler.assembleUncertainRegions(templateHookedList,
                                                listOfScanPSMMap, scanDnMap, dbDnRatioThresh);
        }

        //Export candidate templates together with contaminant sequences as a fasta file
        String contaminantFile = "D:\\Hao\\database\\contaminants.fasta";
        contaminantFile = "/Users/hao/data/contaminants.fasta";

        contaminantFile = "C:\\hao\\database\\contaminants.fasta";
        //contaminantFile = "/Users/hao/data/contaminants.fasta";
        String candidateTemplateWithContaminant = "C:\\Hao\\database\\candidate_template_with_contaminant.fasta";

        //candidateTemplateWithContaminant = "/Users/hao/data/candidate_template_with_contaminant.fasta";
        int min_template_length = 0;  //If a template length is shorter than the min_length, don't output it.

        System.out.println("Exporting candidate templates: ");
        exportCandidateTemplates(templateHookedList, min_template_length, candidateTemplateWithContaminant, contaminantFile);

    }


    /* Export candidate templates together with contaminant sequences as a fasta file */
    private void exportCandidateTemplates(List<TemplateHooked> templateHookedList, int min_template_length,
                                          String candidateTemplateWithContaminant, String contaminantFile) {
        //Export the sequences to a file
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(candidateTemplateWithContaminant))) {

            for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
                List<char[]> candidateTemplates = templateHookedList.get(templateId).getModifiedSeq();
                String templateAccession = templateHookedList.get(templateId).getTemplateAccession();
                if (candidateTemplates == null) {
                    //If no candidate, means no change need to make to the template, export the template directly
                    System.out.println(">" + templateAccession);
                    System.out.println(new String(templateHookedList.get(templateId).getSeq()));
                    bw.write(">" + templateAccession);
                    bw.write("\n");
                    bw.write(new String(templateHookedList.get(templateId).getSeq()));
                    bw.write("\n");
                    continue;
                }
                for (int i = 0; i < candidateTemplates.size(); i++) {
                    //Only export template longer than min_template_length to keep only heavy or light chain. Delete fragments.
                    if (candidateTemplates.get(i).length < min_template_length) {
                        continue;
                    }

                    if (candidateTemplates.size() < 2) {
                        System.out.println(">" + templateAccession);
                        System.out.println(new String(candidateTemplates.get(i)));
                        bw.write(">" + templateAccession);
                        bw.write("\n");
                        bw.write(new String(candidateTemplates.get(i)));
                        bw.write("\n");
                    } else {
                        // System.out.println(">can" + (i + 1) + "_" + templateAccession);
                        // System.out.println(new String(candidateTemplates.get(i)));
                        bw.write(">can" + (i + 1) + "_" + templateAccession);

                        if (candidateTemplates.get(i)[0] != '>') {
                            bw.write("\n");
                            bw.write(new String(candidateTemplates.get(i)));
                        } else {
                            String candidateTemplateWithInfo = new String(candidateTemplates.get(i));
                            String[] fields = candidateTemplateWithInfo.split(">");
                            String infoString = fields[1];
                            String templateString = fields[2];
                            bw.write("_" + infoString + "\n");
                            bw.write(templateString);
                        }

                        bw.write("\n");
                    }
                }
                //Debug
                //break;
            }

            //Attach the contaminant sequences
            String contaminantSeqs = null;
            try (BufferedReader br = new BufferedReader(new FileReader(contaminantFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    bw.write(line);
                    bw.write("\n");
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    /*
    Attach the candidate templates to the mutated templates in templateHookedList
     */
    private void generateCandidateTemplates(List<TemplateHooked> templateHookedList,
                                            List<HashMap<String, PSMAligned>> listOfScanPSMMap,
                                            double significantThreshold) {
        for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
            System.out.println("Template " + templateId + " " + templateHookedList.get(templateId).getTemplateAccession());
            TemplateHooked aTemplateHooked = templateHookedList.get(templateId);
            List<char[]> top2CandidateTemplates = findCandidateForOneTemplate(aTemplateHooked,
                    listOfScanPSMMap.get(templateId), significantThreshold);
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
                short[] ionScores = scanIonScoresMap.get(scan);
                /*
                if (psmList.get(i).getPeptide().contains("-")) {
                    ionScores = modifyIonScoreForDeletion(psmList.get(i).getPeptide(), ionScores);
                }*/
                psmList.get(i).setIonScores(ionScores);
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

    private List<char[]> findCandidateForOneTemplate(TemplateHooked aTemplateHooked,
                                                     HashMap<String, PSMAligned> scanPSMMap,
                                                     double significantThreshold) {
        MutationValidator validator = new MutationValidator();
        List<HashMap<List<Integer>, List<MutationsPattern>>> mutationsOnTemplateList = validator.findSignificantMutations(aTemplateHooked, scanPSMMap, significantThreshold);
        //  printMutationsOnTemplate(mutationsOnTemplateList);

        TemplateCandidateBuilder templateCandidateBuilder = new TemplateCandidateBuilder(mutationsOnTemplateList);
        List<char[]> topCandidateTemplates = templateCandidateBuilder.buildCandidateTemplate(aTemplateHooked, scanPSMMap, significantThreshold);

//        trimTemplateCEnd(aTemplateHooked, topCandidateTemplates);

        return topCandidateTemplates;

    }

    /**
     * Build a map between the scan number and the DenovoOnly object.
     * @return
     */
    private HashMap<String, DenovoOnly> buildScanDnMap(List<DenovoOnly> denovoOnlyList) {
        HashMap<String, DenovoOnly> scanDnMap = new HashMap<>();
        for (DenovoOnly dn : denovoOnlyList) {
            scanDnMap.put(dn.getScan(), dn);
        }
        return scanDnMap;
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

package com.uwaterloo.Statistic;

import com.uwaterloo.*;
import com.uwaterloo.Reader.PSMIonsReader;
import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static com.uwaterloo.Assembler.setIonScoresForPSMList;

public class TemplateCoverage {
    private List<TemplateHooked> hookTemplates(String dir) {
        String psmFile = dir + "DB search psm.csv";
        PSMReader psmReader = new PSMReader();
        List<PSM> psmList = psmReader.readCSVFile(psmFile);

        String psmIonsFile = dir + "PSM ions.csv";
        PSMIonsReader ionsReader = new PSMIonsReader();
        HashMap<String, short[]> scanIonPosesMap = ionsReader.readPSMIonsFile(psmIonsFile);
        HashMap<String, short[]> scanIonScoresMap = ionsReader.setIonsScore(scanIonPosesMap);
        setIonScoresForPSMList(psmList, scanIonScoresMap);

        String templateFasta = dir + "proteins.fasta";
        TemplatesLoader loader = new TemplatesLoader();
        List<Template> templateList = loader.loadTemplateFasta(templateFasta);

        ProteinPeptideReader ppReader = new ProteinPeptideReader(templateList);
        String proteinPeptideFile = dir + "protein-peptides.csv";
        HashMap<String, List<TMapPosition>> peptideProteinMap = ppReader.readProteinPeptideFile(proteinPeptideFile);

        //Map psms to templates
        TemplatePSMsAligner psmAligner = new TemplatePSMsAligner();
        List<TemplateHooked> templateHookedList = psmAligner.alignPSMstoTemplate(psmList, templateList, peptideProteinMap);
        ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList = psmAligner.getPsmAlignedList();

        return templateHookedList;
    }

    private static int getScore(PSMAligned psmAligned, int i) {
        return psmAligned.getIonScores()[i - psmAligned.getStart()];
    }

    private String statistic(List<TemplateHooked> templateHookedList, int scoreThresh) {
        String reportString = "";

        for (TemplateHooked templateHooked : templateHookedList) {
            System.out.println(templateHooked.getTemplateAccession());
            int coveredAANum = 0;
            Set<String> scanSet = new HashSet<>();
            int[] templateAAScore = new int[templateHooked.getSeq().length];
            int confidentAANum = 0;
            int moreConfidentAANum = 0;
            int scoreSum = 0;

            for (int i = 0; i < templateHooked.getSeq().length; i++) {
                List<String> scanList = templateHooked.getMappedScanList().get(i);
                if (scanList.size() == 0) {
                    continue;
                }

                //Compute the total confidence score sum of DB result at this position
                List<PSMAligned> dbAlignedList = templateHooked.getDbList().get(i);
                List<PSMAligned> spiderList = templateHooked.getSpiderList().get(i);

                if (dbAlignedList.size() > 0) {
                    for (PSMAligned psmAligned : dbAlignedList) {
                        int start = psmAligned.getStart();
                        int end = psmAligned.getEnd();
                        for (int j = start; j <= end; j++) {
                            if (psmAligned.getIonScores() == null) continue;
                            templateAAScore[j] += psmAligned.getIonScores()[j - start];
                        }
                        //Add scans with db result mapped to the template
                        scanSet.add(psmAligned.getScan());
                    }
                }

                if (spiderList.size() > 0) {
                    for (PSMAligned psmAligned : spiderList) {
                        if (psmAligned.getPeptide().contains("del")) {
                            System.err.println("Del in spider result, skip");
                            continue;
                        }

                        int start = psmAligned.getStart();
                        int end = psmAligned.getEnd();
                        for (int j = start; j <= end; j++) {
                            if (psmAligned.getIonScores() == null) continue;
                            templateAAScore[j] += psmAligned.getIonScores()[j - start];
                        }
                        //Add scans with db result mapped to the template
                        scanSet.add(psmAligned.getScan());
                    }
                }

                //Compute the coverage of db at this position
                coveredAANum += 1;

                //Add the score of this position to the total sum of scores of this template
                scoreSum += templateAAScore[i];
                //System.out.println(i + "," + templateAAScore[i]);

                if (templateAAScore[i] >= scoreThresh) {
                    confidentAANum++;
                    if (templateAAScore[i] >= scoreThresh * 2) {
                        moreConfidentAANum++;
                    }
                }

                //scanSet.addAll(scanList);

            }
            reportString += templateHooked.getTemplateAccession() + "\t" + coveredAANum + "\t" + scanSet.size() +
                    "\t" + scoreSum + "\t" + confidentAANum + "\t" + moreConfidentAANum + "\n";
            System.out.println(coveredAANum + "," + scanSet.size() + ", " + scoreSum + ", " + confidentAANum +
                    ", " + moreConfidentAANum);
        }

        return reportString;
    }

    private void exportStatitic(String dir, String reportString) {
        String filePath = dir + "statistic.txt";
        try (BufferedWriter br = new BufferedWriter(new FileWriter(filePath))) {
            br.write("Accession\tcoverage\tpsmNum\ttotalScore\tconfidentAANum\tmoreConfidentAANum\n");
            br.write(reportString);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        String dir = "D:\\Hao\\result\\ab19001.5enzymes.4tempaltes_SPIDER_65\\";
        //dir = "D:\\Hao\\result\\ab19001.5enzymes.new_SPIDER_91\\";
        //dir = "D:\\Hao\\result\\Nuno2016_HC_SPIDER_66\\";

        TemplateCoverage tc = new TemplateCoverage();
        List<TemplateHooked> templateHookedList = tc.hookTemplates(dir);

        //Draw the histogram of score of each position in a template, choose 1000 as the threshold whether this position is a confident one.
        int scoreThresh = 1000;
        String reportString = tc.statistic(templateHookedList, scoreThresh);

        tc.exportStatitic(dir, reportString);
    }


}

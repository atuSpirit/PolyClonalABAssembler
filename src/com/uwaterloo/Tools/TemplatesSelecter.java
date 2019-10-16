package com.uwaterloo.Tools;

import com.uwaterloo.Reader.PSMIonsReader;
import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;
import com.uwaterloo.ScanTemplateMapper.*;
import com.uwaterloo.Utils.PSM;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static com.uwaterloo.Assembler.setIonScoresForPSMList;

public class TemplatesSelecter {
    List<TemplateHooked> lightChainTemplates;
    List<TemplateHooked> heavyChainTemplates;
    final int MAX_LIGHT_CHAIN_LEN = 220;
    //The metrics to evaluate template
    final int METRIC_CONFIDENT_AA = 0;   //sum of score of AA passed confScoreThresh
    final int METRIC_SUM_SCORE = 1; //Summation of AA scores
    final int METRIC_FRAGCONFSUM = 2;   //Number of AA with full fragmentation in db result
    int topK;

    public TemplatesSelecter(int topK) {
        this.topK = topK;
        lightChainTemplates = new ArrayList<>();
        heavyChainTemplates = new ArrayList<>();
    }


    public static List<TemplateHooked> hookTemplates(String dir) {
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

        return templateHookedList;
    }

    private void splitHeavyLight(List<TemplateHooked> templateHookedList) {
        for (TemplateHooked templateHooked : templateHookedList) {
            if (templateHooked.getSeq().length <= MAX_LIGHT_CHAIN_LEN) {
                lightChainTemplates.add(templateHooked);
            } else {
                heavyChainTemplates.add(templateHooked);
                //System.out.println("heavy " + templateHooked.getSeq().length);
            }
        }
    }

    /**
     * Get top template according to specified criteria
     * @param templateList
     * @param scoreThresh
     * @param flag to specify the critera used to sort the template list.
     *             0: Using confidentAANum
     *             1: Using sum of AA score
     *             2: Using number of AA that have full fragmentation (fragConfsSum)
     * @return
     */
    private TemplateHooked getTopTemplate(List<TemplateHooked> templateList, int scoreThresh, int flag) {
        List<TemplateStatistics> templateStatisticsList = generateStatistics(templateList, scoreThresh);
        if (flag == METRIC_CONFIDENT_AA) {
            Collections.sort(templateStatisticsList, TemplateStatistics.cmpReverseConfidentAANum());
        } else if (flag == METRIC_FRAGCONFSUM) {
            Collections.sort(templateStatisticsList, TemplateStatistics.cmpReverseFragConfsSum());
        } else if (flag == METRIC_SUM_SCORE){
            Collections.sort(templateStatisticsList, TemplateStatistics.cmpReverseScoreSum());
        }

        int templateId = templateStatisticsList.get(0).getTemplateId();
    //    System.out.println(templateList.get(templateId).getTemplateAccession() + " " + templateStatisticsList.get(0).toString());
        return templateList.get(templateId);
    }

    private List<TemplateStatistics> generateStatistics(List<TemplateHooked> templateHookedList, int scoreThresh) {
        List<TemplateStatistics> templateStatisticsList = new ArrayList<>();

        int templateNum = templateHookedList.size();
        for (int templateId = 0; templateId < templateNum; templateId++) {
            TemplateHooked templateHooked = templateHookedList.get(templateId);
            int coveredAANum = 0;
            Set<String> scanSet = new HashSet<>();
            int[] templateAAScore = new int[templateHooked.getSeq().length];
            int confidentAANum = 0;
            int moreConfidentAANum = 0;
            int scoreSum = 0;
            int fragConfSum = 0;    //The total number of confident fragmentation AA with db result

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
                            if (psmAligned.getIonScores()[j - start] == 100) {
                                fragConfSum++;
                            }
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

            }
            TemplateStatistics templateStatistics = new TemplateStatistics(templateId, scanSet.size(),
                    coveredAANum, confidentAANum, moreConfidentAANum, scoreSum, fragConfSum);
            templateStatisticsList.add(templateStatistics);
        }
        return templateStatisticsList;
    }

    private void decreasePsmConfScores(TemplateHooked templateHooked, float decreaseRatio) {
        for (int i = 0; i < templateHooked.getSeq().length; i++) {
            List<PSMAligned> dbList = templateHooked.getDbList().get(i);
            descreaseConfScores(dbList, decreaseRatio);
            List<PSMAligned> spiderList = templateHooked.getSpiderList().get(i);
            descreaseConfScores(spiderList, decreaseRatio);
        }
    }

    private void descreaseConfScores(List<PSMAligned> psmAlignedList, float decreaseRatio) {
        for (PSMAligned psmAligned : psmAlignedList) {
            if (psmAligned.getIonScores() == null) {
                continue;
            }
            short[] ionScores = psmAligned.getIonScores();
            for (int i = 0; i < ionScores.length; i++) {
                ionScores[i] = (short) (ionScores[i] * decreaseRatio);
            }
        }
    }

    private void printStatistics(List<TemplateStatistics> templateStatisticsList) {
        for (TemplateStatistics templateStatistics : templateStatisticsList) {
            System.out.println(templateStatistics.getTemplateId() + "\t" +
                    + templateStatistics.getCoveredAANum() + "\t" + templateStatistics.getPsmNum() +
                    "\t" + templateStatistics.getScoreSum() + "\t" + templateStatistics.getConfidentAANum()
                    + "\t" + templateStatistics.getMoreConfidentAANum()
                    + "\t" + templateStatistics.getFragConfsSum());
        }
    }

    private void exportStatitic(String filePath, List<TemplateHooked> templateHookedList, int scoreThresh) {
        List<TemplateStatistics> templateStatisticsList = generateStatistics(templateHookedList, scoreThresh);
        Collections.sort(templateStatisticsList, TemplateStatistics.cmpReverseConfidentAANum());

        try (BufferedWriter br = new BufferedWriter(new FileWriter(filePath))) {
            br.write("Accession\tpsmNum\tcoveredAANum\tconfidentAANum\tmoreConfidentAANum\tfragConfSum\ttotalScore");
            br.write("\n");
            for (TemplateStatistics templateStatistics : templateStatisticsList) {
                br.write(templateHookedList.get(templateStatistics.getTemplateId()).getTemplateAccession());
                br.write("\t");
                br.write(templateStatistics.toString());
                br.write("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * Select topK templates from list of templateHooked.
     * @param templateHookedList
     * @param topK
     * @param scoreThresh
     * @param decreaseRatio
     * @param metricFlag The flag to choose metric to evaluate a template
     * @return
     */
    public List<TemplateHooked> selectTemplates(List<TemplateHooked> templateHookedList, int topK,
                                                int scoreThresh, float decreaseRatio,
                                                int metricFlag) {
        List<TemplateHooked> topKTemplates = new ArrayList<>();

        for (int i = 0; i < topK; i++) {
            TemplateHooked topHCTemplate = getTopTemplate(templateHookedList, scoreThresh, metricFlag);
            topKTemplates.add(topHCTemplate);

            decreasePsmConfScores(topHCTemplate, decreaseRatio);

        }

        return topKTemplates;
    }

    /** Outdated
    public List<TemplateHooked> selectTemplatesAccordingToFragmentConfs(List<TemplateHooked> templateHookedList,
                                                                        int topK) {
        List<TemplateHooked> topKTemplates = new ArrayList<>();
        CoverageConfEvaluator confEvaluator = new CoverageConfEvaluator();
        HashMap<Integer, List<Integer>> fragConfSumTemplateMap = new HashMap<>();
        List<Integer> fragConfSumList = new ArrayList<>();

        for (int index = 0; index < templateHookedList.size(); index++) {
            int[] fragmentConfs = confEvaluator.evaluateFragmentationConf(templateHookedList.get(index));
            int sum = 0;
            for (int conf : fragmentConfs) {
                sum += conf;
            }
            fragConfSumList.add(sum);

            if (fragConfSumTemplateMap.containsKey(sum)) {
                List<Integer> indexList = fragConfSumTemplateMap.get(sum);
                indexList.add(index);
            } else {
                List<Integer> indexList = new ArrayList<>();
                indexList.add(index);
                fragConfSumTemplateMap.put(sum, indexList);
            }
        }

        Collections.sort(fragConfSumList);

        int topIndex = 0;
        int size = fragConfSumList.size();
        while (topIndex < topK) {
            int sum = fragConfSumList.get(size - 1 - topIndex);
            System.out.println(sum);
            List<Integer> indexList = fragConfSumTemplateMap.get(sum);
            for (int index : indexList) {
                topKTemplates.add(templateHookedList.get(index));
                topIndex++;
                if (topIndex == topK) {
                    break;
                }
            }
        }
        return topKTemplates;
    }
     */

    public List<TemplateHooked> getLightChainTemplates() {
        return lightChainTemplates;
    }

    public List<TemplateHooked> getHeavyChainTemplates() {
        return heavyChainTemplates;
    }

    public void exportFasta(String filePath, List<TemplateHooked> templateHookedList) {
        try (BufferedWriter br = new BufferedWriter(new FileWriter(filePath))) {
            for (TemplateHooked templateHooked : templateHookedList) {
                br.write(">" + templateHooked.getTemplateAccession());
                br.write("\n");
                br.write(new String(templateHooked.getSeq()));
                br.write("\n");

                System.out.println(">" + templateHooked.getTemplateAccession());
                System.out.println(new String(templateHooked.getSeq()));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static void main(String[] args) {

        String dir = "C:\\Hao\\result\\ab19001.polyclonal.templateSelected_SPIDER_49\\";
        dir = "C:\\Hao\\result\\Nuno.HC_SPIDER_12\\";
        dir = "C:\\Hao\\result\\ab19001.1_SPIDER_32\\";
        dir = "C:\\Hao\\result\\Nuno.LC_SPIDER_12\\";
        dir = "C:\\hao\\result\\Hieu.mixed_data_SPIDER_22\\";
        dir = "C:\\hao\\result\\NIST_Waters.clean_SPIDER_11\\";
        dir = "C:\\hao\\result\\NIST_Waters.EThcd_SPIDER_127\\";
        dir = "C:\\hao\\result\\NIST_Waters.Peaks8_SPIDER_42\\";
        dir = "C:\\hao\\result\\NIST_Waters.Peaks8_SPIDER_63\\";
        dir = "C:\\hao\\result\\ab19001.peaks8_SPIDER_36\\";
        //dir = "C:\\hao\\result\\Nuno.HC.Peaks8_SPIDER_63\\";

        //dir = "D:\\Hao\\result\\Water_mAB.clean_SPIDER_11\\";

        int topK = 4;
        int scoreThresh = 300;  //300 to select better initial templates


        //If selecting from all antibody database, set it less than 1. If choose from template candidate, set it to 1
        float descreaseRatio = 0.0f;  //If selecting templates from antibody database, choose 0.1f.  If choose top template, use 0

        TemplatesSelecter templatesSelecter = new TemplatesSelecter(topK);
        int metricFlag = templatesSelecter.METRIC_CONFIDENT_AA;
        metricFlag = templatesSelecter.METRIC_FRAGCONFSUM;

        List<TemplateHooked> templateHookedList = hookTemplates(dir);
        templatesSelecter.splitHeavyLight(templateHookedList);

        String filePath = dir + "statistic.txt";
        templatesSelecter.exportStatitic(filePath, templateHookedList, scoreThresh);
        String hcFilePath = dir + "statistic.heavy.txt";
        templatesSelecter.exportStatitic(hcFilePath, templatesSelecter.getHeavyChainTemplates(), scoreThresh);
        String lcFilePath = dir + "statistic.light.txt";
        templatesSelecter.exportStatitic(lcFilePath, templatesSelecter.getLightChainTemplates(), scoreThresh);

        List<TemplateHooked> topKHeavyTemplates = templatesSelecter.selectTemplates(templatesSelecter.getHeavyChainTemplates(),
                                        topK, scoreThresh, descreaseRatio, metricFlag);
        String hcFastaFile = dir + "heavy.top" + topK + ".fasta";
        templatesSelecter.exportFasta(hcFastaFile, topKHeavyTemplates);

        List<TemplateHooked> topKLightTemplates = templatesSelecter.selectTemplates(templatesSelecter.getLightChainTemplates(),
                                        topK, scoreThresh, descreaseRatio, metricFlag);
        String lcFastaFile = dir + "light.top" + topK + ".fasta";
        templatesSelecter.exportFasta(lcFastaFile, topKLightTemplates);


    }

}

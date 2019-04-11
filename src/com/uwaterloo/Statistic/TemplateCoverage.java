package com.uwaterloo.Statistic;

import com.uwaterloo.*;
import com.uwaterloo.Reader.PSMIonsReader;
import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;

import java.util.*;

import static com.uwaterloo.Assembler.setIonScoresForPSMList;

public class TemplateCoverage {
    public static void main() {
        String dir = "D:\\Hao\\result\\ab19001.5enzymes.new_SPIDER_79\\";
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
        //ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList = psmAligner.getPsmAlignedList();

        for (TemplateHooked templateHooked : templateHookedList) {
            System.out.println(templateHooked.getTemplateAccession());
            int coveredAANum = 0;
            Set<String> scanSet = new HashSet<>();
            for (int i = 0; i < templateHooked.getSeq().length; i++) {
                List<String> scanList = templateHooked.getMappedScanList().get(i);
                if (scanList.size() > 0) {
                    coveredAANum += 1;
                    scanSet.addAll(scanList);
                }
            }
            System.out.println(coveredAANum + "," + scanSet.size());
        }

    }
}

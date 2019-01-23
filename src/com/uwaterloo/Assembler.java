package com.uwaterloo;

import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;

import java.util.HashMap;
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

        TemplateHooked templateHooked1 = templateHookedList.get(0);
        char[] seq = templateHooked1.getSeq();
        int size = seq.length;
        for (int i = 0; i < size; i++) {
            System.out.println(seq[i] + "\t" + templateHooked1.getDbList().get(i).size() + "\t" +
                    templateHooked1.getSpiderList().get(i).size());
        }
    }
}

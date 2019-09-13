package com.uwaterloo.Tools;

import com.uwaterloo.Utils.AAMass;
import com.uwaterloo.Reader.TemplatesLoader;
import com.uwaterloo.ScanTemplateMapper.Template;

import java.util.List;

public class IntactMassValidater {
    public static float computeIntactMass(char[] proteinSeq) {
        float intactMass = 0.0f;

        for (char AA : proteinSeq) {
            String AAstr = "";
            AAstr += AA;
            if (!AAMass.AA_MASS_TABLE.containsKey(AAstr)) {
                //System.out.println(AAstr);
                continue;
            }
            intactMass += AAMass.AA_MASS_TABLE.get(AAstr);

        }
        return intactMass;
    }

    public static void main(String[] args) {
        String templateFasta = "C:\\hao\\database\\abysis.fasta";
        templateFasta = "C:\\hao\\result\\NIST_Waters.1_SPIDER_147\\proteins.fasta";
        TemplatesLoader loader = new TemplatesLoader();
        List<Template> templateList = loader.loadTemplateFasta(templateFasta);
        for (Template template : templateList) {
            float intactMass = IntactMassValidater.computeIntactMass(template.getSeq());
            System.out.println(template.getTemplateAccession() + "," + intactMass);
        }
    }



}

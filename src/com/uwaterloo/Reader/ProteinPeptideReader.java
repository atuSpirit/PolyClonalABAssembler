package com.uwaterloo.Reader;

import com.sun.xml.internal.ws.api.message.ExceptionHasMessage;
import com.uwaterloo.TMapPosition;
import com.uwaterloo.Template;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/* Read in the protein-peptide.csv result of Peaks. */
public class ProteinPeptideReader extends CSVReader {
    HashMap<String, Integer> fieldIndexMap;
    //A peptide could map to multiple templates
    HashMap<String, List<TMapPosition>> peptideProteinMap;

    HashMap<String, Integer> proteinAccessionsIdMap;

    public ProteinPeptideReader(List<Template> templateList) {
        peptideProteinMap = new HashMap<>();
        this.proteinAccessionsIdMap = mapProteinAccessionId(templateList);
    }

    public HashMap<String, List<TMapPosition>> readProteinPeptideFile(String proteinPeptideFile) {
        try (BufferedReader br = new BufferedReader(new FileReader(proteinPeptideFile))) {
            String titleString = br.readLine();
            this.fieldIndexMap = parseTitle(titleString);

            String line;
            while ((line = br.readLine()) != null) {
                readOneLine(line);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return this.peptideProteinMap;
    }

    private void readOneLine(String line) {
        String[] fields = line.trim().split(",");

        String peptide = removeProteinEnd(fields[fieldIndexMap.get("Peptide")]);
        int start = Integer.valueOf(fields[fieldIndexMap.get("Start")]) - 1;
        int end = Integer.valueOf(fields[fieldIndexMap.get("End")]) - 1;
        String templateAccession = fields[fieldIndexMap.get("Protein Accession")];
        int templateId = this.proteinAccessionsIdMap.get(templateAccession);

        TMapPosition tMapPosition = new TMapPosition(templateId, start, end);
        if (peptideProteinMap.get(peptide) == null) {
            List<TMapPosition> tMapPositionList = new ArrayList<>();
            tMapPositionList.add(tMapPosition);
            this.peptideProteinMap.put(peptide, tMapPositionList);
        } else {
            this.peptideProteinMap.get(peptide).add(tMapPosition);
        }

/*
        if (peptideProteinMap.get("KGFYPSDIAVEWESN(+.98)GQPEN(+.98)NY") != null) {
            System.out.println();
        }
        */
    }

    /**
     * Helper function to remove the first AA and the last AA in the peptide
     * sequence in protein-peptide.csv.  eg. Q.FE(sub N)WYM(sub V)DGVEVHNAK.T
     * remove Q.  and .T
     * @param peptide The peptide containing the AA next to the peptide sequence
     *                in the protein sequence
     * @return a peptide in the middle
     */
    private String removeProteinEnd(String peptide) {
        Pattern p = Pattern.compile("^(\\w\\.)(.+)(\\.\\w)$");
        Matcher m = p.matcher(peptide);
        if (m.find()) {
            peptide = m.group(2);
        } else {
            p = Pattern.compile("^(\\w\\.)(.+)");
            m = p.matcher(peptide);
            if (m.find()) {
                peptide = m.group(2);
            } else {
                p = Pattern.compile("(.+)(\\.\\w)$");
                m = p.matcher(peptide);
                if (m.find()) {
                    peptide = m.group(1);
                } else {
                    System.err.println(peptide + " does not match.");
                }
            }
        }
        return peptide;
    }


    /* Map the protein accession to protein id */
    public HashMap<String, Integer> mapProteinAccessionId(List<Template> templateList) {
        HashMap<String, Integer> proteinAccessionsIdMap = new HashMap<>();
        for (Template t : templateList) {
            proteinAccessionsIdMap.put(t.getTemplateAccession(), t.getTemplateId());
        }
        return proteinAccessionsIdMap;
    }

    /* Test the class */
    public static void main(String[] args) {
        String templateFasta = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\Nuno.2016.heavy.template.fasta";
        TemplatesLoader templatesLoader = new TemplatesLoader();
        List<Template> templateList = templatesLoader.loadTemplateFasta(templateFasta);
        ProteinPeptideReader ppReader = new ProteinPeptideReader(templateList);
        String proteinPeptideFile = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\protein-peptides.csv";
        HashMap<String, List<TMapPosition>> peptideProteinMap = ppReader.readProteinPeptideFile(proteinPeptideFile);
        System.out.println("Total loaded " + peptideProteinMap.size());
        List<TMapPosition> tMapPositionList = peptideProteinMap.get("D.S(+27.99)VKGRFTISR.D");
        System.out.println("D.S(+27.99)VKGRFTISR.D" + " is "
                + "mapped to protein " + tMapPositionList.toString());

    }


}

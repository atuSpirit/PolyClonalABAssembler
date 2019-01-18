package com.uwaterloo.Reader;

import com.uwaterloo.TMapPosition;
import com.uwaterloo.Template;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

/* Read in the protein-peptide.csv result of Peaks. */
public class ProteinPeptideReader extends CSVReader {
    HashMap<String, Integer> fieldIndexMap;
    HashMap<String, TMapPosition> peptideProteinMap;

    HashMap<String, Integer> proteinAccessionsIdMap;

    public ProteinPeptideReader(List<Template> templateList) {
        peptideProteinMap = new HashMap<String, TMapPosition>();
        setProteinAccessionIdMap(templateList);
    }

    public void readProteinPeptideFile(String proteinPeptideFile, List<Template> templateList) {
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
    }

    private void readOneLine(String line) {
        String[] fields = line.trim().split(",");

        String peptide = fields[fieldIndexMap.get("Peptide")];
        int start = Integer.valueOf(fields[fieldIndexMap.get("Start")]);
        String templateAccession = fields[fieldIndexMap.get("Protein Accession")];
        int templateId = this.proteinAccessionsIdMap.get(templateAccession);

        TMapPosition tMapPosition = new TMapPosition(templateId, start);
        peptideProteinMap.put(peptide, tMapPosition);
    }

    public void setProteinAccessionIdMap(List<Template> templateList) {
        for (Template t : templateList) {
            proteinAccessionsIdMap.put(t.getTemplateAccession(), t.getTemplateId());
        }
    }


}

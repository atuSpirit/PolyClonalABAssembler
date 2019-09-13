package com.uwaterloo.Reader;

import com.uwaterloo.ScanTemplateMapper.TemplateHooked;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/* Read the DB Search psm.csv file */
public class PSMReader extends CSVReader{
    HashMap<String, Integer> fieldIndexMap;

    public PSMReader() {
    }

    public List<TemplateHooked.PSM> readCSVFile(String psmFile) {
        List<TemplateHooked.PSM> fieldExtracted = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(psmFile))) {
            String titleString = br.readLine();
            this.fieldIndexMap = parseTitle(titleString);

            String line;
            while ((line = br.readLine()) != null) {
                TemplateHooked.PSM psm = readOneLine(line);
                if (psm != null) {
                    fieldExtracted.add(psm);
                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fieldExtracted;
    }

    private TemplateHooked.PSM readOneLine(String line) {
        String[] fields = line.split(",");

        String scan = fields[fieldIndexMap.get("Scan")];
        String peptide = fields[fieldIndexMap.get("Peptide")];
        int intensity = 0;

        if (!fields[fieldIndexMap.get("Area")].equals("")) {
            intensity = Double.valueOf(fields[fieldIndexMap.get("Area")]).intValue();
        }


        if (fieldIndexMap.get("Accession") >= fields.length) {
            return null;
        }
        String proteinAccession = fields[fieldIndexMap.get("Accession")];

        if (proteinAccession.equals("")) {
            return null;
        }

        return new TemplateHooked.PSM(scan, peptide, intensity);
    }


}

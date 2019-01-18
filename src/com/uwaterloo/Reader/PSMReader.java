package com.uwaterloo.Reader;

import com.uwaterloo.PSM;

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

    public List<PSM> readCSVFile(String psmFile) {
        List<PSM> fieldExtracted = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(psmFile))) {
            String titleString = br.readLine();
            this.fieldIndexMap = parseTitle(titleString);

            String line;
            while ((line = br.readLine()) != null) {
                fieldExtracted.add(readOneLine(line));
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fieldExtracted;
    }

    private PSM readOneLine(String line) {
        line = line.trim();
        String[] fields = line.split(",");

        String scan = fields[fieldIndexMap.get("Scan")];
        String peptide = fields[fieldIndexMap.get("Peptide")];

        return new PSM(scan, peptide);
    }


}

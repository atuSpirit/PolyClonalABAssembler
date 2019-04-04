package com.uwaterloo.Reader;

import com.uwaterloo.DenovoOnly;
import com.uwaterloo.PSM;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/* Read in the denovo only.csv file */
public class DenovoOnlyReader extends CSVReader{
    HashMap<String, Integer> fieldIndexMap;

    public DenovoOnlyReader() {

    }

    public List<DenovoOnly> readCSVFile(String psmFile) {
        List<DenovoOnly> fieldExtracted = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(psmFile))) {
            String titleString = br.readLine();
            this.fieldIndexMap = parseTitle(titleString);

            String line;
            while ((line = br.readLine()) != null) {
                DenovoOnly dnOnly = readOneLine(line);
                if (dnOnly != null) {
                    fieldExtracted.add(dnOnly);
                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fieldExtracted;
    }

    private DenovoOnly readOneLine(String line) {
        String[] fields = line.split(",");

        String scan = fields[fieldIndexMap.get("Scan")];
        String peptide = fields[fieldIndexMap.get("Peptide")];
        short alc = Short.valueOf(fields[fieldIndexMap.get("ALC (%)")]);
        short length = Short.valueOf(fields[fieldIndexMap.get("length")]);
        String confString = fields[fieldIndexMap.get("local confidence (%)")];

        int intensity = 0;
        if (!fields[fieldIndexMap.get("Area")].equals("")) {
            intensity = Double.valueOf(fields[fieldIndexMap.get("Area")]).intValue();
        }
        return new DenovoOnly(scan, peptide, alc, length, confString, intensity);
    }

}

package com.uwaterloo.Reader;

import java.util.HashMap;

public class CSVReader {

    public CSVReader() {
    }


    public HashMap<String, Integer> parseTitle(String titleString) {
        String[] fields = titleString.trim().split(",");
        HashMap<String, Integer> fieldIndexMap = new HashMap<>();

        int length = fields.length;
        for (int i = 0; i < length; i++) {
            fieldIndexMap.put(fields[i].trim(), i);
        }

        return fieldIndexMap;
    }


}

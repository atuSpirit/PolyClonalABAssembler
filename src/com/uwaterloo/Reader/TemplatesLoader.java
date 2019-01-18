package com.uwaterloo.Reader;

import com.uwaterloo.Template;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;


public class TemplatesLoader {
    public List<Template> loadTemplateFasta(String templateFastaFile) {
        List<Template> templateList = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(templateFastaFile))) {
            String line;
            String accession = "";
            int templateId = -1;
            String protein_seq = "";
            boolean isContaminant = false;

            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (!protein_seq.equals("") && !isContaminant) {
                        templateId += 1;
                        Template template = new Template(templateId, accession,
                                protein_seq.toCharArray());
                        templateList.add(template);
                    }

                    if (line.contains("#CONTAM#")) {
                        isContaminant = true;
                    } else {
                        isContaminant = false;
                    }

                    Pattern p = Pattern.compile(">(\\S)+");
                    Matcher m = p.matcher(line);
                    if (m.find()) {
                        accession = m.group(0);
                    }
                    protein_seq = "";
                } else {
                    protein_seq += line;
                }
            }

            if (!isContaminant) {
                templateId += 1;
                Template template = new Template(templateId, accession,
                        protein_seq.toCharArray());
                templateList.add(template);
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return templateList;
    }

    public static void main(String[] args) {
        String templateFasta = "D:\\Hao\\data\\for_analysis\\polyclonalAssemblerData\\Nuno.2016.heavy.template.fasta";
        TemplatesLoader loader = new TemplatesLoader();
        List<Template> templateList = loader.loadTemplateFasta(templateFasta);
        System.out.println(templateList.size());
        System.out.println(templateList.get(2).getSeq());
    }
}

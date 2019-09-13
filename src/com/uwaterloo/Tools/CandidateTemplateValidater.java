package com.uwaterloo.Tools;

import com.uwaterloo.Utils.Vertex;
import com.uwaterloo.*;
import com.uwaterloo.Reader.PSMIonsReader;
import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;
import com.uwaterloo.ScanTemplateMapper.*;
import com.uwaterloo.SignificantMutationsFinder.MutationsPattern;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.uwaterloo.Assembler.setIonScoresForPSMList;

/**
 * Given a group of candidateTemplates with the score when construction them and the new PEAKS
 * result mapping to them. For those the score increased, keep them. Otherwise, they are not
 * correct path combination.
 */
public class CandidateTemplateValidater {

    public static void main(String[] args) {
        CandidateTemplateValidater candidateTemplateValidater = new CandidateTemplateValidater();

        String dir = "C:\\hao\\result\\ab19001.polyclonal.05.05_SPIDER_16\\";
        dir = "C:\\Hao\\result\\NIST_Waters.EThcd_SPIDER_84\\";
        String contaminantFile = "C:\\hao\\database\\contaminants.fasta";
        String candidateTemplateWithContaminant = "C:\\Hao\\database\\candidate_template_with_contaminant.fasta";

        TemplatePSMsAligner psmsAligner = candidateTemplateValidater.initialize(dir);
        List<TemplateHooked> templateHookedList = psmsAligner.getTemplateHookedList();
        ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList = psmsAligner.getPsmAlignedList();

        String templateFile = dir + "proteins.fasta";
        List<MutationsPattern> templatePatternList = candidateTemplateValidater.parseTemplateAccessionInfo(templateFile);
        System.out.println("Candidate Template number: " + templatePatternList.size());

        List<TemplateHooked> validatedTemplateHookedList = new ArrayList<>();

        //Transfer list of PSMAligned to list of scanPSMMap for each template
        List<HashMap<String, PSMAligned>> listOfScanPSMMap = new ArrayList<>();
        for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
            MapScanPSMAligned scanPSMMapper = new MapScanPSMAligned(listOfPSMAlignedList.get(templateId));
            HashMap<String, PSMAligned> scanPSMMap = scanPSMMapper.getScanPSMMap();
            listOfScanPSMMap.add(scanPSMMap);
        }


        for (int i = 0; i < templateHookedList.size(); i++) {
            TemplateHooked templateHooked = templateHookedList.get(i);
            MutationsPattern templatePattern = templatePatternList.get(i);
            HashMap<String, PSMAligned> scanPSMMap = listOfScanPSMMap.get(i);
            int scoreDiff = candidateTemplateValidater.validateTemplate(templateHooked, scanPSMMap, templatePattern);
            if (scoreDiff >= 0) {
                templateHooked.setTemplateAccession(templateHooked.getTemplateAccession() + "_scoreIncrease_" +  scoreDiff);
                validatedTemplateHookedList.add(templateHooked);
            }
        }
        System.out.println("Validated candidate template number: " + validatedTemplateHookedList.size());

        for (TemplateHooked templateHooked : validatedTemplateHookedList) {
            System.out.println(templateHooked.getTemplateAccession());
        }
        candidateTemplateValidater.exportValidTemplates(validatedTemplateHookedList, candidateTemplateWithContaminant, contaminantFile);
    }

    private TemplatePSMsAligner initialize(String dir) {
        String psmFile = dir + "DB search psm.csv";
        PSMReader psmReader = new PSMReader();
        List<TemplateHooked.PSM> psmList = psmReader.readCSVFile(psmFile);

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
        return psmAligner;
    }

    private List<MutationsPattern> parseTemplateAccessionInfo(String templateFastaFile) {
        List<MutationsPattern> templatePatterns = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(templateFastaFile))) {
            String line;

            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    String[] fields = line.split("_");
                    String info = fields[fields.length - 1];
                    MutationsPattern templatePattern = parseInfo(info);
                    templatePatterns.add(templatePattern);
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return templatePatterns;

    }


    private MutationsPattern parseInfo(String info) {
        Pattern p = Pattern.compile("(\\S+)\\sat\\s\\[(.+)\\]\\sfreq:\\s(\\d+)\\sscore:\\s(\\d+)\\sintensity:\\s(\\d+)");
        Matcher m = p.matcher(info);
        if (m.find()) {
            String AAs = m.group(1);
            String[] posStrings = m.group(2).split(", ");
            List<Integer> posList = new ArrayList<>();
            for (String posString : posStrings) {
                posList.add(Integer.valueOf(posString));
            }

            int freq = Integer.valueOf(m.group(3));
            int score = Integer.valueOf(m.group(4));
            long intensity = Long.valueOf(m.group(5));

            MutationsPattern pattern = new MutationsPattern(posList, AAs, freq, score);
            return pattern;

        }
        return null;

    }

    /**
     * Return the new path score - old path score.  If decrease, it for sure is not a good path combination, discard.
     * @param templateHooked
     * @param scanPSMMap
     * @param templatePattern
     * @return
     */
    private int validateTemplate(TemplateHooked templateHooked, HashMap<String, PSMAligned> scanPSMMap,
                                     MutationsPattern templatePattern) {
        /* In the case the template is not diverged, there will be no templatePattern. We still need to
           put the sequence in, because this is the correct sequence, which don't need update any more.
         */
        if (templatePattern == null) {
            return 1;
        }
        List<Integer> posArray = templatePattern.getPosList();

        TemplateCandidateBuilder templateCandidateBuilder = new TemplateCandidateBuilder();
        TreeMap<Integer, List<MutationsPattern>> extendedMutations = templateCandidateBuilder.extendPatterns(posArray, templateHooked, scanPSMMap);

        TreeMap<Integer, Map<Character, Vertex>> verticesMap = buildTemplateVertice(templatePattern);
        List<List<Vertex>> verticesList = templateCandidateBuilder.buildEdges(extendedMutations, verticesMap, templateHooked);
        templateCandidateBuilder.addConnection(verticesList);
        //templateCandidateBuilder.printEdges(verticesMap);
        List<MutationsPattern> pathCombination = templateCandidateBuilder.generatePathCombination(verticesList);

        //The pathCombination should contain only one mutationPatter, the templateMutationPattern
        return (pathCombination.get(0).getScore() - templatePattern.getScore());
    }

    /**
     * Build vertice only for the templatePattern.
     * @param templateMutationPattern
     * @return
     */
    private TreeMap<Integer, Map<Character, Vertex>> buildTemplateVertice(MutationsPattern templateMutationPattern) {
        TreeMap<Integer, Map<Character, Vertex>> vertexesMap = new TreeMap<>();
        for (int i = 0; i < templateMutationPattern.getPosList().size(); i++) {
            int pos = templateMutationPattern.getPosList().get(i);
            Map<Character, Vertex> vertices = new HashMap<>();
            char AA = templateMutationPattern.getAAs().charAt(i);
            Vertex v = new Vertex(pos, AA);
            vertices.put(AA, v);

            vertexesMap.put(pos, vertices);
        }
        return vertexesMap;
    }


    private void exportValidTemplates(List<TemplateHooked> templateHookedList,
                                      String candidateTemplateWithContaminant, String contaminantFile) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(candidateTemplateWithContaminant))) {

            for (int templateId = 0; templateId < templateHookedList.size(); templateId++) {
                bw.write(">" + templateHookedList.get(templateId).getTemplateAccession());
                bw.write("\n");
                bw.write(new String(templateHookedList.get(templateId).getSeq()));
                bw.write("\n");
            }

            //Attach the contaminant sequences
            String contaminantSeqs = null;
            try (BufferedReader br = new BufferedReader(new FileReader(contaminantFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    bw.write(line);
                    bw.write("\n");
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        } catch (IOException e) {
        e.printStackTrace();
        }

    }
}

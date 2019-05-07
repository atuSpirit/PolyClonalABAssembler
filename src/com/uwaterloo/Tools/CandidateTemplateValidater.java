package com.uwaterloo.Tools;

import com.uwaterloo.*;
import com.uwaterloo.Reader.PSMIonsReader;
import com.uwaterloo.Reader.PSMReader;
import com.uwaterloo.Reader.ProteinPeptideReader;
import com.uwaterloo.Reader.TemplatesLoader;

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

        String dir = "C:\\Hao\\result\\ab19001.polyclonal.templateSelected_SPIDER_49\\";
        TemplatePSMsAligner psmsAligner = candidateTemplateValidater.initialize(dir);
        List<TemplateHooked> templateHookedList = psmsAligner.getTemplateHookedList();
        ArrayList<ArrayList<PSMAligned>> listOfPSMAlignedList = psmsAligner.getPsmAlignedList();

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
            HashMap<String, PSMAligned> scanPSMMap = listOfScanPSMMap.get(i);
            boolean isValidateTemplate = candidateTemplateValidater.validateTemplate(templateHooked, scanPSMMap);
            if (isValidateTemplate) {
                validatedTemplateHookedList.add(templateHooked);
            }
        }

    }

    private TemplatePSMsAligner initialize(String dir) {
        String psmFile = dir + "DB search psm.csv";
        PSMReader psmReader = new PSMReader();
        List<PSM> psmList = psmReader.readCSVFile(psmFile);

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

    private MutationsPattern parseInfo(String info) {
        Pattern p = Pattern.compile("(\\S+)\\sat\\s[(\\.+)]\\sfreq:\\s(\\d+)\\sscore:\\s(\\d+)\\sintensity:\\s(\\d+)");
        Matcher m = p.matcher(info);
        if (m.find()) {
            String AAs = m.group(1);
            String[] posStrings = m.group(2).split(",");
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

    private boolean validateTemplate(TemplateHooked templateHooked, HashMap<String, PSMAligned> scanPSMMap) {
        String templateAccession = templateHooked.getTemplateAccession();
        String[] fields = templateAccession.split("_");
        String info = fields[fields.length - 1];

        MutationsPattern templatePattern = parseInfo(info);
        List<Integer> posArray = templatePattern.getPosList();

        TemplateCandidateBuilder templateCandidateBuilder = new TemplateCandidateBuilder();
        TreeMap<Integer, List<MutationsPattern>> extendedMutations = templateCandidateBuilder.extendPatterns(posArray, templateHooked, scanPSMMap);

        TreeMap<Integer, Map<Character, Vertex>> verticesMap = buildTemplateVertice(templatePattern);
        List<List<Vertex>> verticesList = templateCandidateBuilder.buildEdges(extendedMutations, verticesMap, templateHooked);
        templateCandidateBuilder.printEdges(verticesMap);
        List<MutationsPattern> pathCombination = templateCandidateBuilder.generatePathCombination(verticesList);

        for (MutationsPattern mutationsPattern : pathCombination) {
            System.out.println(mutationsPattern.toString());
        }

        return pathCombination.get(0).getScore() > templatePattern.getScore();
    }

    /**
     * Build vertice only for the templatePattern.
     * @param templateMutationPattern
     * @return
     */
    private TreeMap<Integer, Map<Character, Vertex>> buildTemplateVertice(MutationsPattern templateMutationPattern) {
        TreeMap<Integer, Map<Character, Vertex>> vertexesMap = new TreeMap<>();
        for (int pos : templateMutationPattern.getPosList()) {
            Map<Character, Vertex> vertices = new HashMap<>();
            char AA = templateMutationPattern.getAAs().charAt(pos);
            Vertex v = new Vertex(pos, AA);
            vertices.put(AA, v);

            vertexesMap.put(pos, vertices);
        }
        return vertexesMap;
    }


}

package com.uwaterloo;

import java.util.*;

public class TemplateCandidateBuilder {
    List<HashMap<List<Integer>, List<String>>> mutationsList;

    public TemplateCandidateBuilder(List<HashMap<List<Integer>, List<String>>> mutationsList) {
        this.mutationsList = mutationsList;
    }

    /* Get the mutation positions set from mutationList*/
    private List<Integer> getPosSet() {
        Set<Integer> posSet = new HashSet<>();
        for (int i = 0; i < mutationsList.size(); i++) {
            HashMap<List<Integer>, List<String>> mutations = mutationsList.get(i);
            for (List<Integer> posList : mutations.keySet()) {
                for (int pos : posList) {
                    posSet.add(pos);
                }
            }
        }
        List<Integer> posList = new ArrayList<>(posSet);
        Collections.sort(posList);

        return posList;
    }

    private TreeMap<Integer, List<MutationsPattern>> sortMutationsAccordingPos() {
        TreeMap<Integer, List<MutationsPattern>> mutationsSorted = new TreeMap<>();
        for (int i = 0; i < mutationsList.size(); i++) {
            HashMap<List<Integer>, List<String>> mutations = mutationsList.get(i);
            for (List<Integer> posList : mutations.keySet()) {
                List<String> patterns = mutations.get(posList);
                int pos = posList.get(0);
                for (String pattern : patterns) {
                    String[] items = pattern.split("_");
                    MutationsPattern mutationsPattern = new MutationsPattern(posList, items[0],
                            Integer.valueOf(items[1]));
                    if (!mutationsSorted.containsKey(pos)) {
                        List<MutationsPattern> newPatternList = new ArrayList<>();
                        mutationsSorted.put(pos, newPatternList);
                    }
                    mutationsSorted.get(pos).add(mutationsPattern);
                }
            }
        }
        return mutationsSorted;
    }
    /* For some pattern, they does not contain the consecutive position from the posList,
        add the AA on template to them to make them consecutive for later procession.
        For exmaple, pos List contain {18, 19, 27}, but one pattern only has {18, 27}.
        Then the AA at 19 on template will be added to the middle.
     */
    private void makePatternsConsecutive(TemplateHooked templateHooked, List<Integer> posArray) {

    }

    private void printMutationsAlongPos(TreeMap<Integer, List<MutationsPattern>> posMutationMap) {
        for (int pos : posMutationMap.keySet()) {
            System.out.println(pos);
            List<MutationsPattern> patterns = posMutationMap.get(pos);
            for (MutationsPattern pattern : patterns) {
                System.out.println(pattern.toString());
            }
        }
    }

    private TreeMap<Integer, List<MutationsPattern>> assembleMutationPatterns(TreeMap<Integer, List<MutationsPattern>> mutationSorted) {
        
    }

    public void buildCandidateTemplate(TemplateHooked templateHooked) {
        List<Integer> posArray = getPosSet();
        System.out.println(posArray.toString());

        TreeMap<Integer, List<MutationsPattern>> mutationSorted = sortMutationsAccordingPos();
        printMutationsAlongPos(mutationSorted);

        TreeMap<Integer, List<MutationsPattern>> assembledMutations = assembleMutationPatterns(mutationSorted);


        //makePatternsConsecutive(templateHooked, posArray);
    }




}

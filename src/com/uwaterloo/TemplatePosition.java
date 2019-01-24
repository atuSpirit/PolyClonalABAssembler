package com.uwaterloo;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * The class storing the amino acide matched to this position and
 * their corresponding scan number and score.
 */
public class TemplatePosition {
    char templateAA;    //The amino acid at the pos of the template
    //The amino acid and
    HashMap<Character, LinkedList<String>> mappedSpectrum;

    public TemplatePosition(char templateAA) {
        this.templateAA = templateAA;
        this.mappedSpectrum = new HashMap<>();
    }


    public char getTemplateAA() {
        return templateAA;
    }

    public HashMap<Character, LinkedList<String>> getMappedSpectrum() {
        return this.mappedSpectrum;
    }
}

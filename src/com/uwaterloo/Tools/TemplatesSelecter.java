package com.uwaterloo.Tools;

import com.uwaterloo.Statistic.TemplateCoverage;
import com.uwaterloo.TemplateHooked;
import java.util.List;

public class TemplatesSelecter {

    public static void main(String[] args) {
        String dir = "D:\\Hao\\result\\Waters_mAB_SPIDER_9\\";
        List<TemplateHooked> templateHookedList = TemplateCoverage.hookTemplates(dir);

    }

}

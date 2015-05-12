package org.bioinfo.babelomics.utils;

import org.bioinfo.commons.io.utils.IOUtils;

import java.io.IOException;
import java.util.*;

/**
 * Created by ralonso on 12/27/14.
 */
public class GOManager {


    private String babelomicsHome = System.getenv("BABELOMICS_HOME");
    private Map<String, List<String>> terms;

    public GOManager() {
        terms = new HashMap<String, List<String>>();
        String annotationFile = this.babelomicsHome + "/conf/annotations/go.txt";
        try {
            List<String> raw = IOUtils.readLines(annotationFile);
            for (String goTerm: raw){
                String fields[] = goTerm.split("\t");
                List<String> values = new ArrayList<String>();
                values.add(fields[1]);
                values.add(fields[2]);
                values.add(fields[3]);
                terms.put(fields[0],values);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    public Map<String,List<String>> getTerms(){
        return terms;
    }



}

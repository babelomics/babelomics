package org.bioinfo.babelomics.utils;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.sun.jersey.api.client.Client;
import com.sun.jersey.api.client.ClientResponse;
import com.sun.jersey.api.client.WebResource;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

/**
 * Created by ralonso on 12/27/14.
 */
public class AnnotationManager {

    private String species;
    /**
     * file system databases *
     */
    private Set<String> fsDb;

    /****/
    private String babelomicsHome = System.getenv("BABELOMICS_HOME");

    public AnnotationManager(String species) {
        this.init(species);
    }


    private void init(String species) {

        this.species = species;

        /** This should be red from a file **/
        this.fsDb = new HashSet<String>();
        this.fsDb.add("recon");
        this.fsDb.add("biological_process");
        this.fsDb.add("cellular_component");
        this.fsDb.add("molecular_function");
        this.fsDb.add("go_slim");
        this.fsDb.add("go");
        this.fsDb.add("interpro");
    }

    public Map<String, List<String>> getIdAnnotationsGenome(String db){
        String annotationFile = this.babelomicsHome + "/conf/annotations/" + this.species + "/" + db + ".txt";
        List<String> rawAnnots = new ArrayList<String>();
        try {
            rawAnnots = IOUtils.readLines(annotationFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<String, List<String>> featAnnot = parseFeatureAnnotations(rawAnnots);
        return featAnnot;
    }
    public Map<String, List<String>> getAnnotationIdsGenome(String db){
        String annotationFile = this.babelomicsHome + "/conf/annotations/" + this.species + "/" + db + ".txt";
        List<String> rawAnnots = new ArrayList<String>();
        try {
            rawAnnots = IOUtils.readLines(annotationFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Map<String, List<String>> featAnnot = parseFeatureAnnotations(rawAnnots);
        Map<String, List<String>> annotFeat = parseAnnotationsOriginalId(featAnnot);
        return annotFeat;
    }

    private Map<String, List<String>> parseAnnotationsOriginalId(Map<String, List<String>> oriIdAnnotation) {
        Map<String, List<String>> annotFeat = new HashMap<String, List<String>>();

        for (String id : oriIdAnnotation.keySet()) {
            List<String> annotations = oriIdAnnotation.get(id);
            for (String annot : annotations) {
                List<String> orisId = new ArrayList<String>();
                if (annotFeat.containsKey(annot)) {
                    orisId = annotFeat.get(annot);
                }
                orisId.add(id);
                annotFeat.put(annot, orisId);
            }
        }
        return annotFeat;
    }



    private Map<String, List<String>> parseFeatureAnnotations(List<String> rawAnnots) {
        Map<String, List<String>> featureAnnot = new HashMap<String, List<String>>();
        for (String raw : rawAnnots) {
            String fields[] = raw.split("\t");
            if (raw.startsWith("#") || !raw.contains("\t") || fields.length < 2) {
                continue;
            }
            String feature = fields[0];
            String annot = fields[1];
            List<String> annotations = new ArrayList<String>();
            if (!featureAnnot.containsKey(feature)) {
                featureAnnot.put(feature, annotations);
            } else {
                annotations = featureAnnot.get(feature);

            }
            annotations.add(annot);
            featureAnnot.put(feature, annotations);

        }
        return featureAnnot;
    }

}

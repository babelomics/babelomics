package org.bioinfo.babelomics.utils;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.sun.jersey.api.client.Client;
import com.sun.jersey.api.client.ClientResponse;
import com.sun.jersey.api.client.WebResource;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.babelomics.utils.GOManager;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

/**
 * Created by ralonso on 12/27/14.
 */
public class XrefManager {

    private List<String> list;
    private String species;
    private Map<String, List<String>> listXref;
    private String db;
    /**
     * file system databases *
     */
    private Set<String> fsDb;

    /****/
    private String uri;
    private int requests;
    private String babelomicsHome = System.getenv("BABELOMICS_HOME");

    public XrefManager(String id, String species) {
        this.list = new ArrayList<String>();
        this.list.add(id);
        this.init(species);
    }

    public XrefManager(List<String> list, String species) {
        this.list = list;
        this.init(species);
    }

    private void init(String species) {

        this.species = species;
        this.listXref = new HashMap<String, List<String>>();

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

    public void fillXref(StringBuilder batch, String db) {
        Client client = Client.create();

        int lastIndexOfComma = batch.lastIndexOf(",");
        batch.deleteCharAt(lastIndexOfComma);

        String uri = this.uri + batch.toString() + "/xref?dbname=" + db;
        System.out.println("uri = " + uri);
//      WebResource webResource = client.resource("https://www.ebi.ac.uk/cellbase/webservices/rest/v3/hsapiens/feature/id/" + batch.toString() + "/xref?dbname=ensembl_transcript");
        WebResource webResource = client.resource(uri);
        ClientResponse response = webResource.get(ClientResponse.class);
        String resp = response.getEntity(String.class);
        ObjectMapper mapper = new ObjectMapper();
        JsonNode actualObj = null;
        try {
            actualObj = mapper.readTree(resp).get("response");
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (JsonNode jnode : actualObj) {
            if (listXref.containsKey(jnode.get("id")))
                continue;
            List<String> ids = new ArrayList<String>();
            listXref.put(jnode.get("id").asText(), ids);
            Iterator<JsonNode> it = jnode.get("result").iterator();
            while (it.hasNext()) {
                JsonNode node = it.next();
                ids.add(node.get("id").asText());
            }
        }
    }

    public Map<String, List<String>> getXrefs(String db) {
        this.db = db;
        this.resolveUri();
        String oriDb = db;
        if (this.fsDb.contains(db)) {
            db = "ensembl_transcript";
        }
        int numberIds = 0;
        StringBuilder batch = new StringBuilder();
        for (String id : list) {
            if (id.contains("/") || id.contains("\"") || id.contains("\'") || id.contains(" ")) {

                continue;
            }
            if (numberIds == requests) {
                this.fillXref(batch, db);
                numberIds = 0;
                batch = new StringBuilder();
            }
            batch.append(id);
            batch.append(",");
            numberIds++;
        }
        this.fillXref(batch, db);

        /** if we want go, kegg... keep going**/
        if (this.fsDb.contains(oriDb)) {
            listXref = this.getXrefsFileSystem(oriDb);
        }
        return listXref;
    }

    public Map<String, List<String>> getXrefsFileSystem(String db) {
        String annotationFile = this.babelomicsHome + "/conf/annotations/" + this.species + "/" + db + ".txt";
        Map<String, List<String>> oriIdAnnotation = null;
        try {
            List<String> rawAnnots = IOUtils.readLines(annotationFile);

            /** Get map with this structure: [featureId in ENSTRANSCRIPT: anotationsList] **/
            Map<String, List<String>> featAnnot = parseFeatureAnnotations(rawAnnots);

            /** Get map with the annotations related to the original id with this structure: [annotationId: originalIdList ]**/
//            annotFeat = this.parseAnnotationsOriginalId(listXref, featAnnot);
            oriIdAnnotation = this.parseOriginalIdAnnotation(listXref, featAnnot);

            /** Filter annotations **/
//            FeatureList<AnnotationItem> annotations = this.filter(annotFeat);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return oriIdAnnotation;
    }

    public Map<String, List<String>> parseOriginalIdAnnotation(Map<String, List<String>> list1Xref, Map<String, List<String>> featAnnot) {
        Map<String, List<String>> annotOriId = new HashMap<String, List<String>>();

        for (String id : list1Xref.keySet()) {
            List<String> xrefs = list1Xref.get(id);
            for (String xref : xrefs) {
                if (featAnnot.containsKey(xref)) {
                    Set<String> annotationsSet = new HashSet<String>();
                    if (annotOriId.containsKey(id))
                        annotationsSet.addAll(annotOriId.get(id));
                    annotationsSet.addAll(featAnnot.get(xref));
                    List<String> annotations = new ArrayList<String>();
                    annotations.addAll(annotationsSet);
                    annotOriId.put(id, annotations);
                }

            }
        }
        return annotOriId;
    }

    public Map<String, List<String>> parseAnnotationsOriginalId(Map<String, List<String>> oriIdAnnotation) {
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

    public FeatureList<AnnotationItem> filter(Map<String, List<String>> oriIdAnnotation, FunctionalFilter filter) {
        Map<String, List<String>> annotFeature = parseAnnotationsOriginalId(oriIdAnnotation);
//
//        for (String id: annotFeature.keySet()){
//            System.out.println("id = " + id);
//            System.out.println("annotFeature = " + annotFeature.get(id));
//        }


        FeatureList<AnnotationItem> annotations = new FeatureList<AnnotationItem>();
        for (String annot : annotFeature.keySet()) {
            List<String> features = annotFeature.get(annot);
            /** Size filter **/
            if (filter.getMinNumberGenes() <= features.size() && features.size() <= filter.getMaxNumberGenes()) {
                for (String feature : features) {
                    annotations.add(new AnnotationItem(feature, annot));

                }
            }
        }

        return annotations;
    }

    public Map<String, List<String>> parseFeatureAnnotations(List<String> rawAnnots) {
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


    private void resolveUri() {
        Properties prop = new Properties();
        InputStream input = null;

        try {
            String file = this.babelomicsHome + "/conf/webservices.conf";
            input = new FileInputStream(file);
            prop.load(input);
            this.uri = prop.getProperty("HOST");
            String auxSpecies = species;

            /** This should be in a file **/
            if (species.equalsIgnoreCase("hsa"))
                auxSpecies = "hsapiens";
            if (species.equalsIgnoreCase("dre"))
                auxSpecies = "drerio";
            if (species.equalsIgnoreCase("mmu"))
                auxSpecies = "mmusculus";
            if (species.equalsIgnoreCase("rno"))
                auxSpecies = "rnorvegicus";
            if (species.equalsIgnoreCase("dme"))
                                auxSpecies = "dmelanogaster";
            if (species.equalsIgnoreCase("sce"))
                auxSpecies = "scerevisiae";
            if (species.equalsIgnoreCase("cel"))
               auxSpecies = "celegans";
            if (species.equalsIgnoreCase("ath"))
              auxSpecies = "athaliana";

            this.uri += auxSpecies + "/feature/id/";

            this.requests = Integer.parseInt(prop.getProperty("REQUESTS"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

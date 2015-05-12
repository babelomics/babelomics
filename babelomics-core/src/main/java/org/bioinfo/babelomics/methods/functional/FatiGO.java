package org.bioinfo.babelomics.methods.functional;

import java.io.IOException;

import org.bioinfo.babelomics.utils.XrefManager;
import org.bioinfo.babelomics.utils.AnnotationManager;
import org.bioinfo.babelomics.utils.GOManager;

import java.lang.annotation.Annotation;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.*;

import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.log.Logger;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.GOSlimFilter;
import org.bioinfo.infrared.funcannot.filter.InterproFilter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.stats.inference.FisherExactTest;
import org.bioinfo.babelomics.utils.filters.ReconFilter;


public class FatiGO {

    public static final int REMOVE_NEVER = 0;
    public static final int REMOVE_EACH = 1;
    public static final int REMOVE_REF = 2;
    public static final int REMOVE_GENOME = 3;
    public static final int REMOVE_ALL = 4;

    private enum modeTypes {
        list2list, list2genome, list2rest
    }

    private String mode;

    // input params
    private List<String> list1;
    private List<String> list2;
    private FunctionalFilter filter;
    //    private DBConnector dbConnector;
    private int testMode;
    private int duplicatesMode;
    private boolean isYourAnnotations;
    private Logger logger;

    private String species;
    private String db;

    // test
    private TwoListFisherTest fisher;

    // results
    private List<TwoListFisherTestResult> results;
    private FeatureList<AnnotationItem> annotations;

    // Term sizes
    private Map<String, Integer> termSizes;

    // summary
    // list 1
    private int list1AnnotatedCounter;
    private double list1MeanAnnotationsPerId;
    private int list1SizeBeforeDuplicates;
    private int list1SizeAfterDuplicates;
    // list 2
    private int list2AnnotatedCounter;
    private double list2MeanAnnotationsPerId;
    private int list2SizeBeforeDuplicates;
    private int list2SizeAfterDuplicates;

    private String babelomicsHome = System.getenv("BABELOMICS_HOME");


    // Two list constructor
    public FatiGO(List<String> list1, List<String> list2, FunctionalFilter filter, String species, int testMode, int duplicatesMode) {
        this.list1 = list1;
        this.list2 = list2;
        this.filter = filter;
        this.species = species;
        this.testMode = testMode;
        this.duplicatesMode = duplicatesMode;
        this.isYourAnnotations = false;
        this.mode = modeTypes.list2list.name();
    }

    // One list against Genome constructor
    public FatiGO(List<String> list1, FunctionalFilter filter, String species) {
        this.list1 = list1;
//		this.list2 = InfraredUtils.getGenome(dbConnector);
        /** Leer el fichero del genoma **/
        try {

            this.list2 = IOUtils.readLines(this.babelomicsHome + "/conf/annotations/" + species + "/genome.txt");
            this.mode = modeTypes.list2genome.name();
//            this.list2 = Files.readAllLines(Paths.get(this.babelomicsHome + "/conf/data/genome.txt"), Charset.defaultCharset());

        } catch (IOException e) {
            e.printStackTrace();
        }
        this.filter = filter;
        this.species = species;
        this.testMode = FisherExactTest.GREATER;
        this.duplicatesMode = REMOVE_GENOME;
        this.isYourAnnotations = false;
    }

    // Your annotations two list constructor
    public FatiGO(List<String> list1, List<String> list2, FeatureList<AnnotationItem> annotations, int testMode, int duplicatesMode, String species) {
        this.list1 = list1;
        this.list2 = list2;
        this.annotations = annotations;
        this.testMode = testMode;
        this.duplicatesMode = duplicatesMode;
        this.isYourAnnotations = true;
        this.mode = modeTypes.list2list.name();
        this.species = species;
        System.out.println("");
    }

    // Your annotations one list constructor
    public FatiGO(List<String> list1, FeatureList<AnnotationItem> annotations) {
        this.list1 = list1;
        this.list2 = getAnnotationIds(annotations);
        this.annotations = annotations;
        this.testMode = FisherExactTest.GREATER;
        this.duplicatesMode = REMOVE_REF;
        this.isYourAnnotations = true;
        this.mode = modeTypes.list2rest.name();
    }


    private List<String> getAnnotationIds(FeatureList<AnnotationItem> annotations) {
        List<String> idList = new ArrayList<String>();
        HashMap<String, Boolean> idHash = new HashMap<String, Boolean>();
        for (AnnotationItem item : annotations) {
            if (!idHash.containsKey(item.getId())) {
                idList.add(item.getId());
                idHash.put(item.getId(), true);
            }
        }
        return idList;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    public void setDb(String db) {
        this.db = db;
    }

    public void run() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, IOException {

        if (logger == null) logger = new Logger("FatiGO");

        logger.println("Starting FatiGO test...");

        // duplicates managing
        logger.print("removing duplicates...");
        removeDuplicates();
        logger.println("OK");

        logger.print("preparing list union...");
        List<String> all = new ArrayList<String>();

        all.addAll(list1);
        if (!this.mode.equalsIgnoreCase(modeTypes.list2genome.name())) {
            all.addAll(list2);
        }

        //logger.println("list2_annotation.size: " + InfraredUtils.getAnnotations(dbConnector, list2, filter).size());

        // annotation
        logger.println("getting annotations from file system");
        String db = "";
        if (!isYourAnnotations) {


            if (filter instanceof GOFilter) {
                db = ((GOFilter) filter).getNamespace();
            }
            if (filter instanceof ReconFilter) {
                db = "recon";
            }
            if (filter instanceof GOSlimFilter) {
                db = "go_slim";
            }
            if (filter instanceof InterproFilter) {
                db = "interpro";
            }

            System.out.println("db = " + db);


            XrefManager xrefManager = new XrefManager(all, this.species);
            Map<String, List<String>> xrefs = xrefManager.getXrefs(db);
//            System.out.println("xrefs = " + xrefs);
            if (this.mode.equalsIgnoreCase(modeTypes.list2genome.name())) {
                Set<String> list2set = new HashSet<String>();
                list2set.addAll(list2);
                AnnotationManager annoManager = new AnnotationManager(this.species);
                Map<String, List<String>> annots = annoManager.getIdAnnotationsGenome(db);
                for (String s : annots.keySet()) {
                    if (list2set.contains(s))
                        xrefs.put(s, annots.get(s));
                }
            }

            annotations = xrefManager.filter(xrefs, filter);
//            System.out.println("annotations = " + annotations);


            /** Set term sizes **/
//            annotations = InfraredUtils.getAnnotations(dbConnector, all, filter);
        }
        if (!db.equals("") && !this.species.equals("none")) {
//            this.termSizes = new HashMap<String, Integer>();
            this.termSizes = getAnnotationsTermSizesInGenome(db, this.species);
        }

        /** Add go name to id **/
//        Map<String, List<String>> goTerms = new HashMap<String, List<String>>();
//        if (db.equalsIgnoreCase("biological_process") || db.equalsIgnoreCase("cellular_component") || db.equalsIgnoreCase("molecular_function")) {
//            GOManager goManager = new GOManager();
//            goTerms = goManager.getTerms();
//        }
//        for (AnnotationItem annot : annotations) {
//            if (goTerms.containsKey(annot.getFunctionalTermId())) {
//                String descTerm = goTerms.get(annot.getFunctionalTermId()).get(0) + "(" + annot + ")";
//                annot.setFunctionalTermId(descTerm);
//
//            }
//        }

        logger.println("OK");

        computeAnnotateds();

        // run test
        fisher = new TwoListFisherTest();

        logger.print("executing fisher test...");
        fisher.test(list1, list2, annotations, testMode, termSizes);
        logger.println("OK");

        // get result
        logger.print("getting results...");
        results = fisher.getResults();
        logger.println("OK");


        logger.println("end of FatiGO test...");
    }

    protected Map<String, Integer> getAnnotationsTermSizesInGenome(String db, String species) {
        Map<String, Integer> termSizes = new HashMap<String, Integer>();
        AnnotationManager annotationManager = new AnnotationManager(species);
        Map<String, List<String>> annots = annotationManager.getAnnotationIdsGenome(db);
        String termId;
        for (AnnotationItem annotation : annotations) {
            termId = annotation.getFunctionalTermId();
            if (annots.containsKey(termId)) {
                termSizes.put(termId, annots.get(termId).size());
            }
        }

        logger.println("OK");
        return termSizes;
    }

    public void removeDuplicates() {
        // before
        list1SizeBeforeDuplicates = list1.size();
        list2SizeBeforeDuplicates = list2.size();

        // each list
        if (duplicatesMode != REMOVE_NEVER) {
            list1 = ListUtils.unique(list1);
            list2 = ListUtils.unique(list2);
        }
        //complementary
        if (duplicatesMode == REMOVE_REF) {
            for (String id : list1) {
                if (list2.contains(id)) {
                    list2.remove(id);
                    //list1.remove(id);
                }
            }
        }
        //genome
        if (duplicatesMode == REMOVE_GENOME) {
            list1 = ListUtils.unique(list1);
            XrefManager xrefManager = new XrefManager(list1, this.species);
            Map<String, List<String>> xrefs = xrefManager.getXrefs("ensembl_transcript");
            List<String> ensemblList1 = new ArrayList<String>();
            for (String xref : xrefs.keySet()) {
                ensemblList1.addAll(xrefs.get(xref));
            }
            ensemblList1 = ListUtils.unique(ensemblList1);
            for (String id : ensemblList1) {
                if (list2.contains(id)) {
                    list2.remove(id);
                }
            }
        }
        // all
        if (duplicatesMode == REMOVE_ALL) {
            list1 = ListUtils.unique(list1);
            list2 = ListUtils.unique(list2);
            ArrayList<String> toRemove = new ArrayList<String>();
            for (String id : list1) {
                if (list2.contains(id)) {
                    //list1.remove(id);
                    toRemove.add(id);
                    list2.remove(id);
                }
            }
            for (String item : toRemove) {
                list1.remove(item);
            }
        }

        // after
        list1SizeAfterDuplicates = list1.size();
        list2SizeAfterDuplicates = list2.size();
    }

    private void computeAnnotateds() {
        // list 1
        HashMap<String, Integer> list1Annotations = new HashMap<String, Integer>();
        for (String id : list1) {
            list1Annotations.put(id, 0);
        }
        // list 1
        HashMap<String, Integer> list2Annotations = new HashMap<String, Integer>();
        for (String id : list2) {
            list2Annotations.put(id, 0);
        }
        // run annotations
        int count;
        for (AnnotationItem annot : annotations) {
            String id = annot.getId();
            if (list1Annotations.containsKey(id)) {
                count = list1Annotations.get(id);
                //System.err.print("vale " + count);
                count++;
                list1Annotations.put(id, count);
                //System.err.println(" y lo paso a " + count + " " + list1Annotations.get(id));
            }
            if (list2Annotations.containsKey(id)) {
                count = list2Annotations.get(id);
                count++;
                list2Annotations.put(id, count);
            }
        }
        // counts
        // list 1
        Iterator<String> it1 = list1Annotations.keySet().iterator();
        list1AnnotatedCounter = 0;
        int list1Total = 0;
        String id;
        while (it1.hasNext()) {
            id = it1.next();
            count = list1Annotations.get(id);
            if (count > 0) {
                list1AnnotatedCounter++;
            }
            list1Total += count;
        }
        list1MeanAnnotationsPerId = (double) list1Total / (double) list1SizeAfterDuplicates;
        // list 2
        Iterator<String> it2 = list2Annotations.keySet().iterator();
        list2AnnotatedCounter = 0;
        int list2Total = 0;
        while (it2.hasNext()) {
            id = it2.next();
            count = list2Annotations.get(id);
            if (count > 0) {
                list2AnnotatedCounter++;
            }
            list2Total += count;
        }
        list2MeanAnnotationsPerId = (double) list2Total / (double) list2SizeAfterDuplicates;
    }

    public List<TwoListFisherTestResult> getSignificant(double threshold) {
        if (fisher != null) return fisher.getSignificantResults(threshold);
        return null;
    }

    /**
     * @return the annotations
     */
    public FeatureList<AnnotationItem> getAnnotations() {
        return annotations;
    }


    /**
     * @param annotations the annotations to set
     */
    public void setAnnotations(FeatureList<AnnotationItem> annotations) {
        this.annotations = annotations;
    }


    /**
     * @return the results
     */
    public List<TwoListFisherTestResult> getResults() {
        return results;
    }


    /**
     * @param results the results to set
     */
    public void setResults(List<TwoListFisherTestResult> results) {
        this.results = results;
    }

    /**
     * @return the list1
     */
    public List<String> getList1() {
        return list1;
    }

    /**
     * @param list1 the list1 to set
     */
    public void setList1(List<String> list1) {
        this.list1 = list1;
    }

    /**
     * @return the list2
     */
    public List<String> getList2() {
        return list2;
    }

    /**
     * @param list2 the list2 to set
     */
    public void setList2(List<String> list2) {
        this.list2 = list2;
    }

    /**
     * @return the logger
     */
    public Logger getLogger() {
        return logger;
    }

    /**
     * @param logger the logger to set
     */
    public void setLogger(Logger logger) {
        this.logger = logger;
    }


    /**
     * @return the list1SizeBeforeDuplicates
     */
    public int getList1SizeBeforeDuplicates() {
        return list1SizeBeforeDuplicates;
    }

    /**
     * @param list1SizeBeforeDuplicates the list1SizeBeforeDuplicates to set
     */
    public void setList1SizeBeforeDuplicates(int list1SizeBeforeDuplicates) {
        this.list1SizeBeforeDuplicates = list1SizeBeforeDuplicates;
    }

    /**
     * @return the list1SizeAfterDuplicates
     */
    public int getList1SizeAfterDuplicates() {
        return list1SizeAfterDuplicates;
    }

    /**
     * @param list1SizeAfterDuplicates the list1SizeAfterDuplicates to set
     */
    public void setList1SizeAfterDuplicates(int list1SizeAfterDuplicates) {
        this.list1SizeAfterDuplicates = list1SizeAfterDuplicates;
    }

    /**
     * @return the list2SizeBeforeDuplicates
     */
    public int getList2SizeBeforeDuplicates() {
        return list2SizeBeforeDuplicates;
    }

    /**
     * @param list2SizeBeforeDuplicates the list2SizeBeforeDuplicates to set
     */
    public void setList2SizeBeforeDuplicates(int list2SizeBeforeDuplicates) {
        this.list2SizeBeforeDuplicates = list2SizeBeforeDuplicates;
    }

    /**
     * @return the list2SizeAfterDuplicates
     */
    public int getList2SizeAfterDuplicates() {
        return list2SizeAfterDuplicates;
    }

    /**
     * @param list2SizeAfterDuplicates the list2SizeAfterDuplicates to set
     */
    public void setList2SizeAfterDuplicates(int list2SizeAfterDuplicates) {
        this.list2SizeAfterDuplicates = list2SizeAfterDuplicates;
    }

    /**
     * @return the list1AnnotatedCounter
     */
    public int getList1AnnotatedCounter() {
        return list1AnnotatedCounter;
    }

    /**
     * @param list1AnnotatedCounter the list1AnnotatedCounter to set
     */
    public void setList1AnnotatedCounter(int list1AnnotatedCounter) {
        this.list1AnnotatedCounter = list1AnnotatedCounter;
    }


    /**
     * @return the list2AnnotatedCounter
     */
    public int getList2AnnotatedCounter() {
        return list2AnnotatedCounter;
    }

    /**
     * @param list2AnnotatedCounter the list2AnnotatedCounter to set
     */
    public void setList2AnnotatedCounter(int list2AnnotatedCounter) {
        this.list2AnnotatedCounter = list2AnnotatedCounter;
    }

    /**
     * @return the list1MeanAnnotationsPerId
     */
    public double getList1MeanAnnotationsPerId() {
        return list1MeanAnnotationsPerId;
    }

    /**
     * @param list1MeanAnnotationsPerId the list1MeanAnnotationsPerId to set
     */
    public void setList1MeanAnnotationsPerId(double list1MeanAnnotationsPerId) {
        this.list1MeanAnnotationsPerId = list1MeanAnnotationsPerId;
    }

    /**
     * @return the list2MeanAnnotationsPerId
     */
    public double getList2MeanAnnotationsPerId() {
        return list2MeanAnnotationsPerId;
    }

    /**
     * @param list2MeanAnnotationsPerId the list2MeanAnnotationsPerId to set
     */
    public void setList2MeanAnnotationsPerId(double list2MeanAnnotationsPerId) {
        this.list2MeanAnnotationsPerId = list2MeanAnnotationsPerId;
    }

    /**
     * @return the termSizes
     */
    public Map<String, Integer> getTermSizes() {
        return termSizes;
    }

    /**
     * @param termSizes the termSizes to set
     */
    public void setTermSizes(Map<String, Integer> termSizes) {
        this.termSizes = termSizes;
    }


}

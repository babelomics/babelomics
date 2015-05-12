package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bioinfo.babelomics.exception.EmptyAnnotationException;
import org.bioinfo.commons.log.Logger;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.babelomics.utils.XrefManager;
import org.bioinfo.infrared.funcannot.filter.GOFilter;

public abstract class GeneSetAnalysis {

	public static final int ASCENDING_SORT = 1;
	public static final int DESCENDING_SORT = 2;
	public static final double DEFAULT_PVALUE_THRESHOLD = 0.05;
	public static final String FATISCAN = "FATISCAN";
	public static final String LOGISTIC = "LOGISTIC";
	
	protected String method;
	
	// input params
	protected List<String> idList;	
	protected List<Double> statistic;
	protected FeatureData rankedList;
	protected FunctionalFilter filter;
//	protected DBConnector dbConnector;
	protected int order;
	protected boolean isYourAnnotations;
	
	// results
	protected FeatureList<AnnotationItem> annotations;
	protected List<GeneSetAnalysisTestResult> results;
	
	// term sizes
	protected Map<String,Integer> termSizes;
	protected Map<String,List<String>> geneMap;
	
	// summary
	private int annotatedCounter;
	private double meanAnnotationsPerId;
	protected Logger logger;
	
	public void prepare(String species) throws InvalidIndexException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, EmptyAnnotationException{
				
		// init logger
		if (logger==null) logger = new Logger();
		logger.setStdOutputActive(true);
		
		logger.print("preparing gene lists...");
		prepareLists();
		logger.println("OK");
		
		// annotation		
		logger.print("getting annotations...");
		if(!isYourAnnotations) {
			String db ="";
			if (filter instanceof GOFilter) {
				db = ((GOFilter) filter).getNamespace();
			}
			XrefManager xrefManager = new XrefManager(idList, species);
			Map<String, List<String>> xrefs = xrefManager.getXrefs(db);
			annotations = xrefManager.filter(xrefs, filter);
//			annotations = InfraredUtils.getAnnotations(dbConnector, idList, filter);
		}
		if(annotations==null || annotations.size()==0) throw new EmptyAnnotationException();
		logger.println(annotations.size() + " annotations found...OK");
		
		logger.print("loading gene map...");
		loadGeneMap();
		logger.println("OK");
		
		results = new ArrayList<GeneSetAnalysisTestResult>();
	}
		
	public void prepareLists() throws InvalidIndexException{
		// id list
		idList = rankedList.getDataFrame().getRowNames();//.getColumn(0);
		
		// statistic		
		statistic = ArrayUtils.toList(rankedList.getDataFrame().getColumnAsDoubleArray(0));
		
		// order ranked list		
		int[] sortIndex = ListUtils.order(statistic);
		statistic = ListUtils.ordered(statistic,sortIndex);
		idList = ListUtils.ordered(idList,sortIndex);
		//if(order==ASCENDING_SORT){
		if(order==DESCENDING_SORT){
			Collections.reverse(idList);
			Collections.reverse(statistic);
		}		
	}
	
	private void loadGeneMap(){
		
		// init annotation counters
		meanAnnotationsPerId = 0;
		annotatedCounter = 0;
		HashMap<String,Integer> listAnnotations = new HashMap<String, Integer>();
		HashMap<String,Boolean> idHash = new HashMap<String, Boolean>();
		for(String id: idList){
			listAnnotations.put(id, 0);
			idHash.put(id, true);
		}
		
		// process annotations
		System.err.println("loading geneMap (found " + annotations.size() + " annotations)");
		geneMap = new HashMap<String, List<String>>();
		int count;
		
		for(AnnotationItem annotation:annotations){
//			System.err.print("annotation: " + annotation.getFunctionalTermId());
//			System.err.println(annotation.getId());
//			System.err.println("listAnnotations.get(annotation.getId()): " + annotation.getId() + " : " + listAnnotations.get(annotation.getId()));
			if(!geneMap.containsKey(annotation.getFunctionalTermId())){
				geneMap.put(annotation.getFunctionalTermId(), new ArrayList<String>());
			}
			if(idHash.containsKey(annotation.getId())){
				geneMap.get(annotation.getFunctionalTermId()).add(annotation.getId());
			}
			if(listAnnotations.containsKey(annotation.getId())){
				count = listAnnotations.get(annotation.getId()) + 1;				
			} else {
				count = 0;
			}
			listAnnotations.put(annotation.getId(), count);		
		}
		System.err.println("loading geneMap (found " + annotations.size() + " annotations)");
		
		// compute counters
		Set<String> geneSet = listAnnotations.keySet();
		Iterator<String> geneIterator = geneSet.iterator();
		String gene;
		int totalAnnotations = 0;
		while(geneIterator.hasNext()){
			gene = geneIterator.next();
			if(listAnnotations.get(gene)>0){
				annotatedCounter++;
				totalAnnotations+= listAnnotations.get(gene);
			}			
		}
		meanAnnotationsPerId = totalAnnotations/(double)geneSet.size();
		System.err.println("totalAnnotations " + totalAnnotations);
		System.err.println("meanAnnotationsPerId " + meanAnnotationsPerId);
		System.err.println("annotatedCounter " + annotatedCounter);
	}
	
	// run (must be override)
	public abstract void run() throws Exception;
	
	
	// get significants
	public List<GeneSetAnalysisTestResult> getSignificant(){
		return getSignificant(DEFAULT_PVALUE_THRESHOLD);
	}
	public List<GeneSetAnalysisTestResult> getSignificant(double threshold){
		List<GeneSetAnalysisTestResult> significant = new ArrayList<GeneSetAnalysisTestResult>();
		for(GeneSetAnalysisTestResult result: this.results){			
			if(result.getAdjPValue()<threshold) significant.add(result);
		}
		return significant;
	}

	// get result list as a String list
	protected abstract List<String> resultListToStringList(List<GeneSetAnalysisTestResult> results, boolean header);
	
	// get test results as string list
	public List<String> resultsToStringList(){
		return resultListToStringList(results,true);
	}
	public List<String> resultsToStringList(boolean header){
		return resultListToStringList(results,header);
	}		
	public List<String> significantResultsToStringList(List<GeneSetAnalysisTestResult> significant, double threshold){
		return resultListToStringList(getSignificant(threshold),true);
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
	public List<GeneSetAnalysisTestResult> getResults() {
		return results;
	}

	/**
	 * @param results the results to set
	 */
	public void setResults(List<GeneSetAnalysisTestResult> results) {
		this.results = results;
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

	/**
	 * @return the method
	 */
	public String getMethod() {
		return method;
	}

	/**
	 * @param method the method to set
	 */
	public void setMethod(String method) {
		this.method = method;
	}

	/**
	 * @return the idList
	 */
	public List<String> getIdList() {
		return idList;
	}

	/**
	 * @param idList the idList to set
	 */
	public void setIdList(List<String> idList) {
		this.idList = idList;
	}

	/**
	 * @return the rankedList
	 */
	public FeatureData getRankedList() {
		return rankedList;
	}

	/**
	 * @param rankedList the rankedList to set
	 */
	public void setRankedList(FeatureData rankedList) {
		this.rankedList = rankedList;
	}

	/**
	 * @return the annotatedCounter
	 */
	public int getAnnotatedCounter() {
		return annotatedCounter;
	}

	/**
	 * @param annotatedCounter the annotatedCounter to set
	 */
	public void setAnnotatedCounter(int annotatedCounter) {
		this.annotatedCounter = annotatedCounter;
	}

	/**
	 * @return the meanAnnotationsPerId
	 */
	public double getMeanAnnotationsPerId() {
		return meanAnnotationsPerId;
	}

	/**
	 * @param meanAnnotationsPerId the meanAnnotationsPerId to set
	 */
	public void setMeanAnnotationsPerId(double meanAnnotationsPerId) {
		this.meanAnnotationsPerId = meanAnnotationsPerId;
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
	
}

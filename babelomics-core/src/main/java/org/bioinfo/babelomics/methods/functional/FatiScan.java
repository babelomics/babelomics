package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.exception.EmptyAnnotationException;
import org.bioinfo.commons.log.Logger;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.math.stats.MultipleTestCorrection;


public class FatiScan extends GeneSetAnalysis {


	public static final int SHORT_FORMAT = 1;
	public static final int LONG_FORMAT = 2;
	public final static int DEFAULT_NUMBER_OF_PARTITIONS = 30;
	public static final int REMOVE_NEVER = 0;
	public static final int REMOVE_DUPLICATES = 1;
	
	
	// input params	
	private int testMode;	
	private int numberOfPartitions;
	private int outputFormat;
	private int duplicatesMode;

	protected DBConnector dbConnector;

	// test
	private TwoListFisherTest fisher;
		
	
	// Two list constructor
	public FatiScan(FeatureData rankedList, FunctionalFilter filter, DBConnector dbConnector, int numberOfPartitions, int testMode, int outputFormat, int order, int duplicatesMode) {
		this.rankedList = rankedList;
		this.filter = filter;
		this.dbConnector = dbConnector;		
		this.order = order;
		this.isYourAnnotations = false;
		// fatiscan specific
		this.numberOfPartitions = numberOfPartitions;
		this.testMode = testMode;
		this.outputFormat = outputFormat;
		
		// set analysis type id
		this.method = FATISCAN;
		this.duplicatesMode = duplicatesMode; 
	}
	
	public FatiScan(FeatureData rankedList, FeatureList<AnnotationItem> annotations, int numberOfPartitions, int testMode, int outputFormat, int order, int duplicatesMode) {
		
		this.rankedList = rankedList;
		this.annotations = annotations;		
		this.order = order;
		this.isYourAnnotations = true;
		
		// fatiscan specific
		this.numberOfPartitions = numberOfPartitions;
		this.testMode = testMode;
		this.outputFormat = outputFormat;
		
		// set analysis type id
		this.method = FATISCAN;
		this.duplicatesMode = duplicatesMode;
	}
	
	@Override
	public void run() throws InvalidIndexException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, EmptyAnnotationException {
		
		if(logger==null) logger = new Logger("Fatiscan");
				
		// duplicates managing
		logger.print("removing duplicates...");		
		removeDuplicates();		
		logger.println("OK");
		
		logger.print("removing duplicates...");		
		removeDuplicates();		
		logger.println("OK");
		// prepare list
		prepare(this.dbConnector.getSpecies());
		
		int thresholdPosition;
		List<String> list1,list2;
		
		double inc = -(double)(statistic.get(0)-statistic.get(statistic.size()-1))/(numberOfPartitions+1);
		double acum = statistic.get(0) + inc;
		
		System.err.println("test mode: " + testMode);
		// test each partition	
		for(int i=0; i<numberOfPartitions; i++){
			
			thresholdPosition = getThresholdPosition(acum);
			
			System.err.println(i + ": threshold = " + acum + " (" + thresholdPosition + ") ");
			
			// top
			list1 = new ArrayList<String>(idList.subList(0, thresholdPosition));
						
			if(thresholdPosition<(idList.size()-1)){
				
				// bottom
				//list2 = new ArrayList<String>(idList.subList(thresholdPosition + 1, idList.size()-1));
				list2 = new ArrayList<String>(idList.subList(thresholdPosition, idList.size()));
				
				// list1 and list2 have length > 0
				if(list2.size()>0){
					// run test
					fisher = new TwoListFisherTest();
					fisher.test(list1,list2,annotations,testMode,termSizes);							
					// get result
					results.addAll(toGeneSetAnalysisTestResult(fisher.getResults()));					
				}
						
			}
						
			acum+=inc;
			
		}

		//MultipleTestCorrection.BHCorrection(result);
		
		if(outputFormat == SHORT_FORMAT) {			
			HashMap<String,GeneSetAnalysisTestResult> resultsMap = new HashMap<String,GeneSetAnalysisTestResult>();
			// unique term
			for(GeneSetAnalysisTestResult testResult: results){
				if(resultsMap.containsKey(testResult.getTerm())){					
					if(resultsMap.get(testResult.getTerm()).getAdjPValue()>testResult.getAdjPValue()){
						resultsMap.remove(testResult.getTerm());
						resultsMap.put(testResult.getTerm(), testResult);
					}
				} else {
					resultsMap.put(testResult.getTerm(), testResult);
				}
			}			
			// update results
			results.clear();
			results.addAll(resultsMap.values());
			
		}
				
		System.err.println("final results.size: " + results.size());
		
	}
		
	// from TwoListFisherTest result to GeneSetAnalysis result
	private List<GeneSetAnalysisTestResult> toGeneSetAnalysisTestResult(List<TwoListFisherTestResult> twoListFisherTest){
		List<GeneSetAnalysisTestResult> result = new ArrayList<GeneSetAnalysisTestResult>(twoListFisherTest.size());
		for(TwoListFisherTestResult test: twoListFisherTest){			
			GeneSetAnalysisTestResult gseaTest = new GeneSetAnalysisTestResult(test);
			result.add(gseaTest);
		}
		return result;
	}
	
	
	public void removeDuplicates(){
		
		// each list
		if(duplicatesMode==REMOVE_DUPLICATES){			
			this.rankedList = (FeatureData) ListUtils.unique(this.rankedList.getDataFrame().getRowNames());
		}
		
	}
	
	
	
	private int getThresholdPosition(double acum){
		int position = 0;
		for(int i=0; i<statistic.size(); i++){
			if( (order==ASCENDING_SORT && statistic.get(i)>=acum) || (order==DESCENDING_SORT && statistic.get(i)<acum) ) {				
				position = i;
				break;
			}
		}
		return position;
	}


	@Override	
	public List<String> resultListToStringList(List<GeneSetAnalysisTestResult> resultList, boolean header){
		List<String> results = new ArrayList<String>();
		if(header) results.add(GeneSetAnalysisTestResult.fatiScanHeader());
		for(int i=0; i<resultList.size(); i++){			
			results.add(resultList.get(i).toFatiScanString());
		}
		return results;
	}
	
	/**
	 * @return the statistic
	 */
	public List<Double> getStatistic() {
		return statistic;
	}

	/**
	 * @param statistic the statistic to set
	 */
	public void setStatistic(List<Double> statistic) {
		this.statistic = statistic;
	}

}

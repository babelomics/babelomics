package org.bioinfo.babelomics.methods.functional;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.math.data.IntegerMatrix;
import org.bioinfo.math.result.FisherTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.MultipleTestCorrection;
import org.bioinfo.math.stats.inference.FisherExactTest;

public class TwoListFisherTest extends FunctionalTest {
	
	public final static double DEFAULT_PVALUE_THRESHOLD = 0.05;
	
	private List<TwoListFisherTestResult> results;
	
	@Override
	public void test(List<String> list1, List<String> list2, FeatureList<AnnotationItem> annotations, int testMode) {
		test(list1, list2, annotations, testMode, null);
	}
		
	public void test(List<String> list1, List<String> list2, FeatureList<AnnotationItem> annotations, int testMode, Map<String,Integer> termSizes) {
		
		// init term hashes
		HashMap<String, Boolean> mList1 = new HashMap<String,Boolean>();
		for(String id0: list1){
			mList1.put(id0, true);
		}
		HashMap<String, Boolean> mList2 = new HashMap<String,Boolean>();
		for(String id0: list2){
			mList2.put(id0,true);
		}
		//System.err.println("list1: " + list1);
		
		// init counters
		Map<String, Integer> map1 = new HashMap<String, Integer>();
		Map<String, Integer> map2 = new HashMap<String, Integer>();
		Map<String, List<String>> list1Positives = new HashMap<String, List<String>>();
		Map<String, List<String>> list2Positives = new HashMap<String, List<String>>();
		List<String> terms = new ArrayList<String>();
		HashMap<String,Boolean> mTerms = new HashMap<String,Boolean>();
		String term,id;
				
		// count annotations for every gene
		for(AnnotationItem annotation: annotations) {
			term = annotation.getFunctionalTermId();
			id = annotation.getId();
			
			// backup term name
			if(!mTerms.containsKey(term)) {
				terms.add(term);
				mTerms.put(term, true);
			}
			
			// list 1 count
			  // init entry
			if(!map1.containsKey(term)) {
				map1.put(term, 0);
				list1Positives.put(term, new ArrayList<String>());
			}
			  // add term
			//if(list1.contains(id)) {
			if(mList1.containsKey(id)){
				map1.put(term, map1.get(term) + 1);
				list1Positives.get(term).add(id);
			}
			
			// list 2 count
			  // init entry
			if(!map2.containsKey(term)) {
				map2.put(term, 0);
				list2Positives.put(term, new ArrayList<String>());
			}
			  // add term
			//if(list2.contains(id)) {
			if(mList2.containsKey(id)){
				map2.put(term, map2.get(term) + 1);
				list2Positives.get(term).add(id);
			}
		}

		// count summary		
		if(terms!= null && terms.size()>0){			
			IntegerMatrix fisherCounts = new IntegerMatrix(terms.size(), 4);
			for(int i=0 ; i<terms.size() ; i++) {	
				fisherCounts.set(i, 0, map1.get(terms.get(i)));
				fisherCounts.set(i, 1, list1.size()-map1.get(terms.get(i)));
				fisherCounts.set(i, 2, map2.get(terms.get(i)));
				fisherCounts.set(i, 3, list2.size()-map2.get(terms.get(i)));
			}			
			TestResultList<FisherTestResult> testResult = new FisherExactTest().fisherTest(fisherCounts, testMode);
			
			
			// p-value adjustment
			MultipleTestCorrection.BHCorrection(testResult);
			
			// fill results
			results = new ArrayList<TwoListFisherTestResult>(testResult.size());
			int termSizeInGenome;
			for(int i=0; i<testResult.size(); i++){
				termSizeInGenome = 0;
				if(termSizes!=null && termSizes.containsKey(terms.get(i))) {
					termSizeInGenome = termSizes.get(terms.get(i));
				}
				results.add(new TwoListFisherTestResult(terms.get(i),fisherCounts.get(i,0)+fisherCounts.get(i,2),termSizeInGenome,fisherCounts.get(i,0),fisherCounts.get(i,1),fisherCounts.get(i,2),fisherCounts.get(i,3),list1Positives.get(terms.get(i)),list2Positives.get(terms.get(i)),testResult.get(i).getOddRatio(),testResult.get(i).getPValue(),testResult.get(i).getAdjPValue()));
				TwoListFisherTestResult mytest =  new TwoListFisherTestResult(terms.get(i),fisherCounts.get(i,0)+fisherCounts.get(i,2),termSizeInGenome,fisherCounts.get(i,0),fisherCounts.get(i,1),fisherCounts.get(i,2),fisherCounts.get(i,3),list1Positives.get(terms.get(i)),list2Positives.get(terms.get(i)),testResult.get(i).getOddRatio(),testResult.get(i).getPValue(),testResult.get(i).getAdjPValue());
				//if(mytest.getPValue()<0.1) System.err.println(mytest.getList1Positives() + " " + mytest.getList1Negatives() + " " + mytest.getList2Positives() + " " + mytest.getList2Negatives() + " : " + mytest.getPValue() + "/" + mytest.getAdjPValue());
			}
			
		} else {
			//FIXME thrown an exception
			System.out.println("\nNo annotations was found\n");
			results = null;
		}
			
	}

	public List<TwoListFisherTestResult> getSignificantResults(){
		return getSignificantResults(DEFAULT_PVALUE_THRESHOLD);
	}
	
	public List<TwoListFisherTestResult> getSignificantResults(double threshold){
		List<TwoListFisherTestResult> significant = new ArrayList<TwoListFisherTestResult>();
		for(TwoListFisherTestResult result: this.results){
			if(result.getAdjPValue()<threshold) significant.add(result);
		}
		return significant;
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



}

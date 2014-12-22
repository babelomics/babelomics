package org.bioinfo.babelomics.tools.functional;

import java.sql.SQLException;
import java.util.List;

import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.FisherTestResult;
import org.bioinfo.math.result.TestResultList;

@Deprecated
public class FatigoTest {

	public static final int REMOVE_NEVER = 0;
	public static final int REMOVE_EACH = 1;	
	public static final int REMOVE_REF = 2;
	public static final int REMOVE_ALL = 3;
	
	// input params
	private List<String> list1;
	private List<String> list2;
	private Filter filter;
	DBConnector dbConnector;
	private int testMode;
	private int duplicatesMode;
	

	// results
	TestResultList<FisherTestResult> result;
	FeatureList<AnnotationItem> annotations;

	
	public FatigoTest(List<String> list1, List<String> list2, Filter filter, DBConnector dbConnector, int testMode, int duplicatesMode ) {
		this.list1 = list1;
		this.list2 = list2;
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.testMode = testMode;
		this.duplicatesMode = duplicatesMode;
	}
	
		
	public void run() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException {

		
//		List<String> alls = new ArrayList<String>(list1);
//		alls.addAll(list2);	
//		
//		// duplicates managing
//		if(duplicatesMode!=REMOVE_NEVER){
//			list1 = ListUtils.unique(list1);
//			list2 = ListUtils.unique(list2);
//		}
//		if(duplicatesMode==REMOVE_REF){
//			for (String id:list1) {
//				if(list2.contains(id)){
//					list2.remove(id);
//				}
//			}
//		}
//		if(duplicatesMode==REMOVE_ALL){
//			list1 = ListUtils.unique(list1);
//			list2 = ListUtils.unique(list2);
//			for (String id:list1) {
//				if(list2.contains(id)) {
//					list1.remove(id);
//					list2.remove(id);
//				}	
//			}
//		}
//
//		
//		// annotation
//		AnnotationDBManager annotationMng = new AnnotationDBManager(this.dbConnector);
//		if(filter instanceof GOFilter){
//			annotations = annotationMng.getGOAnnotation(alls, (GOFilter) filter);
//		} else if(filter instanceof KeggFilter){
//			annotations = annotationMng.getKeggAnnotation(alls, (KeggFilter) filter);
//		}
//
//				
//		Map<String, Integer> map1 = new HashMap<String, Integer>();
//		Map<String, Integer> map2 = new HashMap<String, Integer>();
//
//						
//	
//			
//
//		names = new ArrayList<String>();
//
//		String term;
//		
//		for(AnnotationItem annotation: annotations) {
//			term = annotation.getFunctionalTermId();
//
//			if ( !names.contains(term) ) names.add(term);
//			if ( !map1.containsKey(term) ) map1.put(term, 0);
//			if ( !map2.containsKey(term) ) map2.put(term, 0);
//
//			if ( list1.contains(item.getId()) ) map1.put(term, map1.get(term) + 1);
//			if ( list2.contains(item.getId()) ) map2.put(term, map2.get(term) + 1);
//		}
//
//		if ( names != null && names.size() > 0 ) {
//			
//			System.out.println("GO terms (" + StringUtils.arrayToString(names, "\t") + ")");
//			System.out.println("fisher matrix (" + names.size() + ", 4):");
//			
//			matrix = new IntegerMatrix(names.size(), 4);
//			for(int i=0 ; i<names.size() ; i++) {
//				//System.out.println(i + ", " + names.get(i) + ": " + map1.get(names.get(i)) + ", " + map2.get(names.get(i)) + ", " + (list1.size() - map1.get(names.get(i))) + ", " + (list1.size() - map1.get(names.get(i))) + ")");
//				matrix.set(i, 0, map1.get(names.get(i)));
//				matrix.set(i, 1, map2.get(names.get(i)));
//				matrix.set(i, 2, list1.size() - map1.get(names.get(i)));
//				matrix.set(i, 3, list2.size() - map2.get(names.get(i)));
//			}
//
//			System.out.println(getMatrix().toString());
//
//			System.out.println("\nbegin fisher test\n");
//
//			result = new FisherExactTest().fisherTest(matrix, mode);
//		} else {
//			System.out.println("\nannotation is nullllll\n");
//			result = null;
//		}
//		System.out.println("\nend fisher test\n");
//
//		// p-value adjustment
//		MultipleTestCorrection.BHCorrection(result);
//		
//		return result;
	}


	/**
	 * @return the result
	 */
	public TestResultList<FisherTestResult> getResult() {
		return result;
	}


	/**
	 * @param result the result to set
	 */
	public void setResult(TestResultList<FisherTestResult> result) {
		this.result = result;
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
	
}

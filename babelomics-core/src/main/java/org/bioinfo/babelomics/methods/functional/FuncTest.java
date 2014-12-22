package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
import org.bioinfo.math.data.IntegerMatrix;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.FisherTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.FisherExactTest;

public class FuncTest {

	private List<String> list1;
	private List<String> list2;
	private Filter filter;
	private String db;

	private List<String> names;
	private IntegerMatrix matrix;
	TestResultList<FisherTestResult> result;

	DBConnector dbConnector;
	AnnotationDBManager annotationMng;

	public static final String GO = "go";
	public static final String KEGG = "kegg";

	public FuncTest(List<String> list1, List<String> list2) {
		this(list1, list2, "go", "hsa");
	}

	public FuncTest(List<String> list1, List<String> list2, String db) {
		this(list1, list2, db, "hsa");
	}

	public FuncTest(List<String> list1, List<String> list2, String db, String specie) {
		this(list1, list2, db, new DBConnector(specie));
	}

	public FuncTest(List<String> list1, List<String> list2, String db, DBConnector dbConnector) {
		this.list1 = list1;
		this.list2 = list2;
		this.db = db;
		this.dbConnector = dbConnector;
		this.annotationMng = new AnnotationDBManager(this.dbConnector);
	}

	public TestResultList<FisherTestResult> run() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException {
		return run(FisherExactTest.TWO_SIDED);
	}

	private TestResultList<FisherTestResult> run(int mode) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException {
		String term;
		List<String> alls = new ArrayList<String>(list1);
		alls.addAll(list2);

		//System.out.println("list1 " + list1.size() + ":\n" + StringUtils.arrayToString(list1, "\n"));				
		//System.out.println("list2 " + list2.size() + ":\n" + StringUtils.arrayToString(list2, "\n"));				
		//System.out.println("alls " + alls.size() + ":\n" + StringUtils.arrayToString(alls, "\n"));				

		Map<String, Integer> map1 = new HashMap<String, Integer>();
		Map<String, Integer> map2 = new HashMap<String, Integer>();

		//System.out.println("input :\n" + StringUtils.arrayToString(alls, "\n"));		

		//System.out.println("alls size = " + alls.size());				

		FeatureList<AnnotationItem> items;
		if ( KEGG.equalsIgnoreCase(db) ) {
			System.out.println("begin of getKeggAnnotation\n");		
			items = annotationMng.getKeggAnnotation(alls, (KeggFilter) filter);
			System.out.println("end of getKeggAnnotation\n");		
		} else {
			System.out.println("begin of getGOAnnotation\n");		
			items = annotationMng.getGOAnnotation(alls, (GOFilter) filter);			
			System.out.println("end of getGOAnnotation\n");		
		}

		System.out.println("size (input, annotation) = (" + alls.size() + ", " + items.size() + ")");

		names = new ArrayList<String>();
		for(AnnotationItem item: items) {
			//System.out.println("item (id, functional term id) -> " + (item.getId() + ", " + item.getFunctionalTermId() + ")"));

			term = item.getFunctionalTermId();

			if ( !names.contains(term) ) names.add(term);
			if ( !map1.containsKey(term) ) map1.put(term, 0);
			if ( !map2.containsKey(term) ) map2.put(term, 0);

			if ( list1.contains(item.getId()) ) map1.put(term, map1.get(term) + 1);
			if ( list2.contains(item.getId()) ) map2.put(term, map2.get(term) + 1);
		}

		if ( names != null && names.size() > 0 ) {
			
			System.out.println("GO terms (" + ListUtils.toString(names, "\t") + ")");
			System.out.println("fisher matrix (" + names.size() + ", 4):");
			
			matrix = new IntegerMatrix(names.size(), 4);
			for(int i=0 ; i<names.size() ; i++) {
				//System.out.println(i + ", " + names.get(i) + ": " + map1.get(names.get(i)) + ", " + map2.get(names.get(i)) + ", " + (list1.size() - map1.get(names.get(i))) + ", " + (list1.size() - map1.get(names.get(i))) + ")");

				matrix.set(i, 0, map1.get(names.get(i)));
				matrix.set(i, 1, map2.get(names.get(i)));
				matrix.set(i, 2, list1.size() - map1.get(names.get(i)));
				matrix.set(i, 3, list2.size() - map2.get(names.get(i)));
			}

			System.out.println(getMatrix().toString());

			System.out.println("\nbegin fisher test\n");

			result = new FisherExactTest().fisherTest(matrix, mode);
		} else {
			System.out.println("\nannotation is nullllll\n");
			result = null;
		}
		System.out.println("\nend fisher test\n");

		return result;
	}

	public List<String> getNames() {
		return names;
	}

	public void setNames(List<String> names) {
		this.names = names;
	}

	public IntegerMatrix getMatrix() {
		return matrix;
	}

	public void setMatrix(IntegerMatrix matrix) {
		this.matrix = matrix;
	}

	public TestResultList<FisherTestResult> getResult() {
		return result;
	}

	public void setResult(TestResultList<FisherTestResult> result) {
		this.result = result;
	}

	public void setFilter(Filter filter) {
		this.filter = filter;
	}

	public Filter getFilter() {
		return filter;
	}

	public void setDb(String db) {
		this.db = db;
	}

	public String getDb() {
		return db;
	}	
}

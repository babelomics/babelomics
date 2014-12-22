package org.bioinfo.babelomics.methods.functional;


import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.log.Logger;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.math.exception.InvalidParameterException;


public class GeneCodis {

	
	public static final int REMOVE_NEVER = 0;
	public static final int REMOVE_EACH = 1;	
	public static final int REMOVE_REF = 2;
	public static final int REMOVE_GENOME = 3;
	public static final int REMOVE_ALL = 4;
	
	
	public enum correctionFactor{
		fdr, 
		permutation, 
		none
		};
		
	public enum testFactor{
		hypergeometric,
		chiSquare,
		both
		};
		
	public enum analysisFactor{
		concurrence,
		singular
		};
	
	// input params
	private List<String> list1;
	private List<String> list2;
	private FunctionalFilter filter;
	private DBConnector dbConnector;
	private int testMode;
	private int duplicatesMode;
	private boolean isYourAnnotations;
	private int RFactor; //my list
	private int rFactor;//referenced list
	private String binPath;
	private Logger logger;
	// test
	
	// results	
	List<String> results;	
	private FeatureList<AnnotationItem> annotations;
	private String outdir;
	private String name;
	private int support = 3;
	private int supportRandom = 3;
	private correctionFactor correction;
	private testFactor test;
	
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
	
	//sig terms
	
	private int significantTerms = -1;

//two list
	
	public GeneCodis(String binPath,String outdir,String name, List<String> list1, List<String> list2, FunctionalFilter filter, DBConnector dbConnector,int duplicatesMode, int support,int supportRandom, correctionFactor correction, testFactor test) {
	
		this.filter = filter;
		this.list1 = list1;
		this.list2 = list2;
		this.dbConnector = dbConnector;
		this.duplicatesMode = duplicatesMode;
		this.isYourAnnotations = false;
		this.binPath = binPath;
		this.outdir= outdir;
		this.name=name;
		this.support= support;
		this.supportRandom = supportRandom;
		this.correction = correction;
		this.test=test;
		this.setRFactor(list1.size());
		this.setrFactor(list2.size());
	}

	
	
	// Your annotations two list constructor
	public GeneCodis(String binPath, String outdir, String name,List<String> list1, List<String> list2,	FeatureList<AnnotationItem> yourAnnotations,DBConnector dbConnector, int duplicatesMode, int support,int supportRandom, correctionFactor correction, testFactor test) {
		this.list1 = list1;
		this.list2 = list2;
		//this.filter = filter;
		this.dbConnector = dbConnector;
		this.duplicatesMode = duplicatesMode;
		this.isYourAnnotations = true;
		this.binPath = binPath;
		this.outdir= outdir;
		this.name=name;
		this.support= support;
		this.supportRandom = supportRandom;
		this.correction = correction;
		this.test=test;
		this.setRFactor(list1.size());
		this.setrFactor(list2.size());
		System.err.println("list2---------------"+list2.size());
		
	}
	
	//one list agains genome
	public GeneCodis(String binPath,String outdir,String name, List<String> list1, FunctionalFilter filter, DBConnector dbConnector, int duplicatesMode, int support, int supportRandom, correctionFactor correction, testFactor test) {
		this.list1 = list1;
		this.list2 = InfraredUtils.getGenome(dbConnector);
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.duplicatesMode = REMOVE_GENOME;
		this.isYourAnnotations = false;
		this.binPath = binPath;
		this.outdir= outdir;
		this.name=name;
		this.support= support;
		this.supportRandom = supportRandom;
		this.correction = correction;
		this.test=test;
		this.setRFactor(list1.size());
		this.setrFactor(list2.size());
		System.err.println("list2---------------"+list2.size());
	}
	
	

	public GeneCodis(String binPath, String outdir, String name,List<String> list1, FeatureList<AnnotationItem> yourAnnotations, DBConnector dbConnector, int duplicatesMode, int support,int supportRandom, correctionFactor correction, testFactor test) {
		this.list1 = list1;
		this.list2 = getAnnotationIds(yourAnnotations);
		//this.filter = filter;
		this.dbConnector = dbConnector;
		this.duplicatesMode = duplicatesMode;
		this.isYourAnnotations = true;
		this.binPath = binPath;
		this.outdir= outdir;
		this.name=name;
		this.support= support;
		this.supportRandom = supportRandom;
		this.correction = correction;
		this.test=test;
		
		this.setRFactor(list1.size());
		
		System.err.println("list2---------------"+list2.size());
	}
	
	private List<String> getAnnotationIds(FeatureList<AnnotationItem> annotations){
		List<String> idList = new ArrayList<String>();
		HashMap<String,Boolean> idHash = new HashMap<String,Boolean>();
		for(AnnotationItem item: annotations){
			if(!idHash.containsKey(item.getId())){
				idList.add(item.getId());
				idHash.put(item.getId(), true);
			}			
		}
		return idList;
	}

	public void run() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException {
		if(logger==null) logger = new Logger("geneCodis");
		logger.print("removing duplicates...");		
		removeDuplicates();
		logger.println("OK");
		

		logger.print("preparing list union...");
		List<String> all = new ArrayList<String>(list1.size()+list2.size());
		all.addAll(this.list1);
		all.addAll(this.list2);
		logger.println("OK");
		
		
		Hashtable<String, Boolean> myHash  = new Hashtable<String, Boolean>();		
		for (int i=0; i< list1.size(); i++){
			myHash.put(list1.get(i), true);
		}
		
		
		// annotation
		logger.print("getting annotations from infrared");		
		if(!isYourAnnotations) annotations = InfraredUtils.getAnnotations(dbConnector, all, filter);	
		logger.println("OK");
		
		
		logger.println("computeAnnotateds()");
		computeAnnotateds();
		logger.println("OK");
		
		//create annot hash
		Map<String, List<String>> myHashAnnotation  = new LinkedHashMap<String, List<String>>();
		
		//fill annotation Hash : <id of gene>, <annotation compact>
		for (int i=0; i< this.annotations.size(); i++){
			if(!myHashAnnotation.containsKey(annotations.get(i).getId())) {
				myHashAnnotation.put(annotations.get(i).getId(), new ArrayList<String>());
			}
			myHashAnnotation.get(annotations.get(i).getId()).add(annotations.get(i).getFunctionalTermId());
		}
		
		List<String> itemsContent = new ArrayList<String>(list1.size()+list2.size());

		logger.print("preparing file for genecodis call...");
		String isRef;
		Iterator<String> e = myHashAnnotation.keySet().iterator();
			String key;
			int autoNum = 0;
			  while (e.hasNext()) {
			     key = (String)e.next();
			     //if id exist in my first hash list, fill 3er column with id of gene
			     isRef = "0";
			     if (myHash.get(key) != null){
			    	 isRef = (String) key;
			    	 }
			     itemsContent.add(autoNum + "\t" + ListUtils.toString(myHashAnnotation.get(key), ",") + "\t\t" + isRef);
			     autoNum ++;
			  }
		setResults(itemsContent);
		
		logger.println("------------doSingleProcess for sing analisis");
		doSingleProcess(this.binPath ,this.outdir+ "/"+name+"_WellFormedInput",analysisFactor.singular, outdir + "/"+ name+ "_singular.txt");
		
		logger.println("------------doSingleProcess for cong analisis");
		doSingleProcess(this.binPath ,this.outdir+ "/"+name+"_WellFormedInput",analysisFactor.concurrence, outdir + "/"+ name+ "_concurrence.txt");
		
		logger.println("-------------OK");
	}
	
	public void doSingleProcess(String binPath, String inputAbsolutPath, analysisFactor analysis, String outputAbsolutPath){		
		
		try {
			IOUtils.write(this.outdir+ "/"+name+"_WellFormedInput", getResults());
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		int ANALISIS = 1;
		int CORRECTION = 0;
		int TEST= 1;
		
		//ANALISIS
		switch(analysis){
		case concurrence:  ANALISIS= 1;break;
		case singular : ANALISIS = 2;break;
		default:
			System.out.println("ANALISIS format not valid");
		}
		
		//CORRECTION
		switch(this.correction){
		case fdr:CORRECTION= -1;break;
		case permutation: CORRECTION= 1;break;
		case none: CORRECTION =0;break;
		default:
			System.out.println("CORRECTION format not valid");
		}
		
		//TEST
		switch(this.test){
		case chiSquare:TEST = 1;break;
		case hypergeometric: TEST = 0;break;
		case both: TEST = 2;break;
		default:
			System.out.println("TEST format not valid");
		}
	
		String cmdStr = binPath +" "+" "+inputAbsolutPath +" "+ this.support+" "+ " -a" + ANALISIS + " -i" + this.supportRandom +" -r"+this.getrFactor()+ " -R"+this.getRFactor() + " -s"+ CORRECTION + " -t" + TEST+ " -o "+outputAbsolutPath;
		System.err.println("genecodis command: " + cmdStr);
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runAsync();
		sp.waitFor();
		
		//adding tag in header
		try {
			String content = IOUtils.toString(outputAbsolutPath);
			IOUtils.write(outputAbsolutPath, "#" + content);
			setSignificantTerms(IOUtils.countLines(outputAbsolutPath)-1);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		///end_execution
		System.err.println();
		System.err.println("Output result:");
		System.err.println(sp.getRunnableProcess().getOutput());		
		System.err.println("Error:");
		System.err.println(sp.getRunnableProcess().getError());
		System.err.println("Exception:");
		System.err.println(sp.getRunnableProcess().getException());
		System.err.println("=====================> END OF LOCAL genecodis\n");
	}
	
	public void removeDuplicates(){
		// before
		list1SizeBeforeDuplicates = list1.size();
		list2SizeBeforeDuplicates = list2.size();
		
		// each list
		
		System.err.println("--------------duplicatesMode---------------"+duplicatesMode);
		
		if(duplicatesMode!=REMOVE_NEVER){
			list1 = ListUtils.unique(list1);
			list2 = ListUtils.unique(list2);
		}
		//complementary
		if(duplicatesMode==REMOVE_REF){
			for (String id:list1) {
				if(list2.contains(id)){
					list2.remove(id);
				}
			}
		}
		//genome
		if(duplicatesMode==REMOVE_GENOME){
			list1 = ListUtils.unique(list1);
			List<String> ensemblList1 = InfraredUtils.toEnsemblId(dbConnector, list1);
			for (String id:ensemblList1) {
				if(list2.contains(id)){
					list2.remove(id);
				}
			}
		}
		// all
		if(duplicatesMode==REMOVE_ALL){
			list1 = ListUtils.unique(list1);
			list2 = ListUtils.unique(list2);
			for (String id:list1) {
				if(list2.contains(id)) {
					list1.remove(id);
					list2.remove(id);
				}	
			}
		}
		
		// after
		list1SizeAfterDuplicates = list1.size();
		list2SizeAfterDuplicates = list2.size();
	}
	
	
	private void computeAnnotateds(){
		// list 1
		HashMap<String,Integer> list1Annotations = new HashMap<String, Integer>();		
		for(String id: list1){
			list1Annotations.put(id, 0);			
		}
		// list 1
		HashMap<String,Integer> list2Annotations = new HashMap<String, Integer>();		
		for(String id: list2){
			list2Annotations.put(id, 0);
		}
		// run annotations
		int count;
		for(AnnotationItem annot: annotations){
			String id = annot.getId();			
			if(list1Annotations.containsKey(id)){				
				count = list1Annotations.get(id);
				//System.err.print("vale " + count);
				count++;
				list1Annotations.put(id,count);
				//System.err.println(" y lo paso a " + count + " " + list1Annotations.get(id));
			}
			if(list2Annotations.containsKey(id)){
				count = list2Annotations.get(id);
				count++;
				list2Annotations.put(id,count);
			}
		}
		// counts
		  // list 1
		Iterator<String> it1 = list1Annotations.keySet().iterator();
		list1AnnotatedCounter = 0;
		int list1Total = 0;
		String id;
		while(it1.hasNext()){
			id = it1.next();
			count = list1Annotations.get(id);			
			if(count>0) {
				list1AnnotatedCounter++;
			}
			list1Total+=count;
		}
		list1MeanAnnotationsPerId = (double)list1Total/(double)list1SizeAfterDuplicates;
		  // list 2
		Iterator<String> it2 = list2Annotations.keySet().iterator();
		list2AnnotatedCounter = 0;
		int list2Total = 0;
		while(it2.hasNext()){
			id = it2.next();
			count = list2Annotations.get(id);			
			if(count>0) {
				list2AnnotatedCounter++;
			}
			list2Total+=count;
		}
		list2MeanAnnotationsPerId = (double)list2Total/(double)list2SizeAfterDuplicates;
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
	public List<String> getResults() {
		return results;
	}


	/**
	 * @param results the results to set
	 */
	public void setResults(List<String> results) {
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

	public void setRFactor(int rFactor) {
		RFactor = rFactor;
	}

	public int getRFactor() {
		return RFactor;
	}

	public void setrFactor(int rFactor) {
		this.rFactor = rFactor;
	}

	public int getrFactor() {
		return rFactor;
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



	public void setSignificantTerms(int significantTerms) {
		this.significantTerms = significantTerms;
	}



	public int getSignificantTerms() {
		return significantTerms;
	}	
}

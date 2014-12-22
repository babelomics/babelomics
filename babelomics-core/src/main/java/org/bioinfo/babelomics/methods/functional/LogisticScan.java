package org.bioinfo.babelomics.methods.functional;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.exception.EmptyAnnotationException;
import org.bioinfo.babelomics.utils.RCommand;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;

public class LogisticScan extends GeneSetAnalysis{
	
	// working paths
	private String scriptPath;
	private String workingDirectory;	
	
	// files
	private String annotationFileName;
	private String rankedListFileName;
	private String outputFileName;
	

	// Infrared annotation constructor
	public LogisticScan(FeatureData rankedList, FunctionalFilter filter, DBConnector dbConnector, int order) {
			
		// params
		this.rankedList = rankedList;
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.order = order;
		this.isYourAnnotations = false;
	
		// set analysis type id
		this.method = LOGISTIC;
		
		// configure working paths
		scriptPath = System.getenv("BABELOMICS_HOME") + "/bin/logistic_univariate.r";
		workingDirectory = "/tmp";
		initTemporalFiles();
		
		
	}
	
	// Your annotations constructor
	public LogisticScan(FeatureData rankedList, FeatureList<AnnotationItem> annotations, int order) {
			
		// params
		this.rankedList = rankedList;
		this.annotations = annotations;
		this.order = order;
		this.isYourAnnotations = true;
		
		// set analysis type id
		this.method = LOGISTIC;
		
		// configure working paths
		scriptPath = System.getenv("BABELOMICS_HOME") + "/bin/logistic_univariate.r";
		workingDirectory = "/tmp";
		initTemporalFiles();
	
	}
	
	// init temporal files
	private void initTemporalFiles(){
		String random = StringUtils.randomString(30);
		this.annotationFileName = random + "_annotations.txt";
		this.rankedListFileName = random + "_ranked_list.txt";
		this.outputFileName = random + "_result.txt";		
	}
	
	@Override
	public void run() throws InvalidIndexException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, IOException, EmptyAnnotationException {
				
		// prepare list
		prepare();
		
		// annotation		
		if(!isYourAnnotations) annotations = InfraredUtils.getAnnotations(dbConnector, idList, filter);
	
		// prepare files		
		saveRankedList();
		saveAnnotations();
		
		// RComand
		RCommand rCommand = new RCommand(scriptPath,workingDirectory);
		// input
		rCommand.addParam("rankingfile",workingDirectory + "/" + rankedListFileName);
		rCommand.addParam("annotationfile",workingDirectory + "/" + annotationFileName);
		// output
		rCommand.addParam("outfile", outputFileName);

		// RUN
		rCommand.exec();
		
		// Load results
		loadResults();

	}
		
	
	private void saveRankedList() throws IOException{
		IOUtils.write(workingDirectory + "/" + rankedListFileName, rankedList.toString());
	}
	
	private void saveAnnotations() throws IOException{
		IOUtils.write(workingDirectory + "/" + annotationFileName, annotations.toString());
	}
	
	private void loadResults() throws IOException{
		
		// init results
		results = new ArrayList<GeneSetAnalysisTestResult>();
		
		// load result file
		FileUtils.checkFile(workingDirectory + "/" + outputFileName);
		String resultsContent = IOUtils.toString(workingDirectory + "/" + outputFileName);
		List<String> terms = StringUtils.toList(resultsContent,"\n");
		
		// init params
		String[] fields;
		String term;
		double adjPValue,logRatio;
		boolean converged;
		int size;
		List<String>list1Ids;
		List<String>list2Ids;
		
		// run rows
		int termSize, termSizeInGenome;
		List<String> genes;
		for(String row: terms){		
			if(!row.startsWith("#") && row.contains("\t")){
				// split row
				fields = row.split("\t");		
				// parse columns
				term = fields[0];
				size = Integer.parseInt(fields[1]);
				if(fields[2].trim().equals("1")) converged = true;
				else converged = false;
				logRatio = Double.parseDouble(fields[3]);
				adjPValue = Double.parseDouble(fields[4]);
				// compute term sizes and gene list 
				if(geneMap.containsKey(term)){
					if(term.equalsIgnoreCase("GO:0004004")){
						System.err.println("genes: " + geneMap.get(term) + " : " + geneMap.get(term).size());
					}
					termSize = geneMap.get(term).size();
					genes = geneMap.get(term);
				} else {
					termSize = 0;
					genes = new ArrayList<String>();
				}
				if(termSizes!=null && termSizes.containsKey(term)){					
					termSizeInGenome = termSizes.get(term);					
				} else {
					termSizeInGenome = 0;
				}				
				GeneSetAnalysisTestResult gseaTest = new GeneSetAnalysisTestResult(term,termSize,termSizeInGenome,genes,converged,logRatio,adjPValue);
				results.add(gseaTest);
			}
		}
	}
	
	@Override
	public List<String> resultListToStringList(List<GeneSetAnalysisTestResult> resultList, boolean header){
		List<String> results = new ArrayList<String>();
		if(header) results.add(GeneSetAnalysisTestResult.LogisticHeader());
		for(int i=0; i<resultList.size(); i++){			
			results.add(resultList.get(i).toLogisticString());
		}
		return results;
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
	 * @return the scriptPath
	 */
	public String getScriptPath() {
		return scriptPath;
	}

	/**
	 * @param scriptPath the scriptPath to set
	 */
	public void setScriptPath(String scriptPath) {
		this.scriptPath = scriptPath;
	}

	/**
	 * @return the workingDirectory
	 */
	public String getWorkingDirectory() {
		return workingDirectory;
	}

	/**
	 * @param workingDirectory the workingDirectory to set
	 */
	public void setWorkingDirectory(String workingDirectory) {
		this.workingDirectory = workingDirectory;
	}

}

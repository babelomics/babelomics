package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.exception.EmptyAnnotationException;
import org.bioinfo.babelomics.methods.functional.FatiScan;
import org.bioinfo.babelomics.methods.functional.GeneSetAnalysis;
import org.bioinfo.babelomics.methods.functional.GeneSetAnalysisTestResult;
import org.bioinfo.babelomics.methods.functional.LogisticScan;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

import es.blast2go.prog.graph.GoGraphException;

public class FatiScanTool  extends FunctionalProfilingTool{
		
	// list1
	private enum Method {FatiScan, Logistic};	
	private Method method;
	private FeatureData rankedList;
	private int numberOfPartitions;
	private int outputFormat;
	private int order;
	protected int duplicatesMode;
	
	
	
	public FatiScanTool(){
		initOptions();
	}

	@Override
	public void initOptions() {
		// parent options
		super.initOptions();
		// mode
		options.addOption(OptionFactory.createOption("method", "Gene set analysis method [fatiscan, logistic]", false, true));
		// list 1
		options.addOption(OptionFactory.createOption("ranked-list", "the feature data containig the ranked list"));
		// test
		options.addOption(OptionFactory.createOption("partitions", "the number of partitions",false,true));
		// extras
		options.addOption(OptionFactory.createOption("output-format", "short (just most significant partition) or long (term results for all partitions) [short]",false,true));
		options.addOption(OptionFactory.createOption("order", "ascend or descend [ascend]",false,true));
		options.addOption(OptionFactory.createOption("higher-label", "label for condition with higher statistical values",false,true));
		options.addOption(OptionFactory.createOption("lower-label", "label for condition with lower statistical values",false,true));
		options.addOption(OptionFactory.createOption("duplicates", "to remove duplicated IDs", false));
	}
		
	@Override
	public void prepare() throws IOException, ParseException, InvalidIndexException {
		super.prepare();
		
		
		String duplicates = commandLine.getOptionValue("duplicates", "never");		
		if(duplicates.equalsIgnoreCase("never")) duplicatesMode = FatiScan.REMOVE_NEVER;
		if(duplicates.equalsIgnoreCase("remove_duplicates")) duplicatesMode = FatiScan.REMOVE_DUPLICATES;
		
		// your annotations

		// method 
		method = Method.FatiScan;
		
		if(commandLine.hasOption("method") && commandLine.getOptionValue("method").equalsIgnoreCase("logistic")) {
			method = Method.Logistic;
		}
		
		// ranked list		
		
		rankedList = new FeatureData(new File(commandLine.getOptionValue("ranked-list")), true);
		
		// test
		numberOfPartitions = Integer.parseInt(commandLine.getOptionValue("partitions","" + FatiScan.DEFAULT_NUMBER_OF_PARTITIONS));
		
		// output format
		if(commandLine.hasOption("output-format") && commandLine.getOptionValue("output-format").equalsIgnoreCase("long")) {
			outputFormat = FatiScan.LONG_FORMAT;
		} else {
			outputFormat = FatiScan.SHORT_FORMAT;
		}
		// sort order
		if(commandLine.hasOption("order") && commandLine.getOptionValue("order").equalsIgnoreCase("ascend")) {
			order = FatiScan.ASCENDING_SORT;
		} else {
			order = FatiScan.DESCENDING_SORT;
		}

	}
	
	@Override
	protected void execute() {
		try {
			
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");			
			
			// infrared connector			
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));		
			
			// prepare params
			prepare();	
	
			// run fatigo's
			if(filterList.size()==0  && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				
				// save id lists
				IOUtils.write(outdir + "/ranked_list.txt", rankedList.toString());
				result.addOutputItem(new Item("ranked_list","ranked_list.txt","Ranked list",Item.TYPE.FILE,Arrays.asList("RANKED_LIST","CLEAN"),new HashMap<String,String>(),"Input data"));
								
				// order ranked list and save genes and rank separately
				FatiScan fatiscan = new FatiScan(rankedList,null,dbConnector,FatiScan.DEFAULT_NUMBER_OF_PARTITIONS,testMode,outputFormat,order,duplicatesMode);
				fatiscan.prepareLists();
				IOUtils.write(outdir + "/id_list.txt", fatiscan.getIdList());
				result.addOutputItem(new Item("id_list","id_list.txt","Id list (sorted)",Item.TYPE.FILE,Arrays.asList("IDLIST","SORTED"),new HashMap<String,String>(),"Input data"));
				IOUtils.write(outdir + "/statistic.txt", ListUtils.toStringArray(fatiscan.getStatistic()));
				result.addOutputItem(new Item("statistic","statistic.txt","Statistic (sorted)",Item.TYPE.FILE,Arrays.asList("STATISTIC","SORTED"),new HashMap<String,String>(),"Input data"));
				
							
				// significant count
				significantCount = new ArrayList<StringBuilder>();
				for(int i=0; i<DEFAULT_PVALUES.length; i++){
					significantCount.add(new StringBuilder("#DB\tNº of significant terms\n"));
				}
				
				// annotation report
				annotationReport = new StringBuilder();
				annotationReport.append("#DB").append("\t").append("Number of annotations").append("\n");
				
				
				// do fatiscans
				double progress = 20;
				double inc = 60.0/filterList.size();
				for(FunctionalFilter filter: filterList) {
					doTest(rankedList,filter,dbConnector,method);					
					jobStatus.addStatusMessage("" + progress, "Executing test");
					progress+=inc;
				}
				
				// update status					
				jobStatus.addStatusMessage("90", "Executing test");
				
				if(isYourAnnotations){
					doYourAnnotationsTest(rankedList,yourAnnotations,method);
				}
				
				// update status					
				jobStatus.addStatusMessage("95", "Saving results");
				
				saveAnnotationsReport();
				
				// save significant count table
				for(int i=0; i<DEFAULT_PVALUES.length; i++){
					IOUtils.write(outdir + "/significant_count_" + pvalueFormatter.format(DEFAULT_PVALUES[i]) + ".txt", significantCount.get(i).toString());
				}
				result.getOutputItems().add(3, new Item("significant","significant_count_${pvalue}.txt","Number of significant terms per DB",Item.TYPE.FILE,Arrays.asList("TABLE,SUMMARY_TABLE,SIGNIFICANT_COUNT_TABLE"),new HashMap<String,String>(),"Significant Results"));				
				result.addMetaItem(new Item("flags","SHOW_PVALUES","",TYPE.MESSAGE));
				
				 
			}
					
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

	
	
	/*
	 * 
	 * FATISCAN
	 * 
	 * 
	 */
	
	
	private void doTest(FeatureData rankedList,FunctionalFilter filter,DBConnector dbConnector, Method method) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		
		// db attributes
		FunctionalDbDescriptor filterInfo = new FunctionalDbDescriptor(filter);

		// get term sizes
		AnnotationDBManager adbm = new AnnotationDBManager(dbConnector);
		Map<String, Integer> termSizes = adbm.getAnnotationTermsSize(filterInfo.getPrefix());
				
		logger.info(filterInfo.getTitle() + "...\n");

		GeneSetAnalysis gsea;
		
		if(method==Method.Logistic){			
			gsea = new LogisticScan(rankedList,filter,dbConnector,order);
			((LogisticScan)gsea).setWorkingDirectory(outdir);
		}else {
			gsea = new FatiScan(rankedList,filter,dbConnector,numberOfPartitions,testMode,outputFormat,order, duplicatesMode);	
		}
			
		// run
		try{

			gsea.setLogger(logger);
			
			// set term sizes
			gsea.setTermSizes(termSizes);
			
			gsea.run();
							
			logger.print("saving results...");
			// save results
			saveGeneSetAnalysisResults(gsea,filterInfo);
			logger.println("OK");
			
		} catch (EmptyAnnotationException ene){
			result.addOutputItem(new Item("annot_" + filterInfo.getName(),"No annotation was found for " + filterInfo.getTitle() + " ids","Annotations for " + filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Annotation files"));
		} catch (Exception ene){
			ene.printStackTrace();
			result.addOutputItem(new Item("annot_" + filterInfo.getName(),"Unexpected error was found while running analysis for " + filterInfo.getTitle() + " ids",filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("ERROR"),new HashMap<String,String>(),"All results"));
		}
		
		logger.info("...end of " + filterInfo.getTitle());
		
	}
	
	private void doYourAnnotationsTest(FeatureData rankedList, FeatureList<AnnotationItem> annotations, Method method) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, InvalidIndexException{
		
		// db attributes
		FunctionalDbDescriptor filterInfo = new FunctionalDbDescriptor("your_annotation","Your annotations", "your_annotations","Your annotations");
		
		// get term sizes
		Map<String, Integer> termSizes = getYourAnnotationsTermSizes();
		
		logger.info(filterInfo.getTitle() + "...\n");
		
		GeneSetAnalysis gsea;
		
		if(method==Method.Logistic){			
			gsea = new LogisticScan(rankedList, annotations,order);
		}else {
			gsea = new FatiScan(rankedList, annotations,numberOfPartitions,testMode,outputFormat,order, duplicatesMode);
		}
		
		// run
		try{
			
			gsea.setLogger(logger);
			
			// set term sizes
			gsea.setTermSizes(termSizes);
			
			gsea.run();
				
			// save results
			saveGeneSetAnalysisResults(gsea,filterInfo);
			
		} catch (EmptyAnnotationException ene){
			result.addOutputItem(new Item("annot_" + filterInfo.getName(),"No annotation was found for " + filterInfo.getTitle() + " ids","Annotations for " + filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Annotation files"));
		} catch (Exception ene){
			result.addOutputItem(new Item("annot_" + filterInfo.getName(),"Unexpected error was found while running analysis for " + filterInfo.getTitle() + " ids",filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("ERROR"),new HashMap<String,String>(),"All results"));
			ene.printStackTrace();
		}
		
		logger.info("...end of " + filterInfo.getTitle());
		
	}
	
	private void saveGeneSetAnalysisResults(GeneSetAnalysis gsea, FunctionalDbDescriptor filterInfo) throws IOException{
		
		String fileName = filterInfo.getName() + ".txt";
		String annotFileName = filterInfo.getName() + ".annot";
		
		// save statistic results				
		IOUtils.write(outdir + "/" + fileName, ListUtils.toString(gsea.resultsToStringList(),"\n"));
		result.addOutputItem(new Item(filterInfo.getName(),fileName,filterInfo.getTitle(),Item.TYPE.FILE,Arrays.asList(filterInfo.getPrefix().toUpperCase() + "_TERM"),new HashMap<String,String>(),"All results"));
						
		// save significant results
		int numberOfSignificantTerms;
		List<GeneSetAnalysisTestResult> significant;
		String formattedPValue;
		for(int i=0; i<DEFAULT_PVALUES.length; i++){
			formattedPValue = pvalueFormatter.format(DEFAULT_PVALUES[i]);
			significant = gsea.getSignificant(DEFAULT_PVALUES[i]);
			if(significant!=null && significant.size()>0){
				IOUtils.write(outdir + "/significant_" + filterInfo.getName() + "_" + formattedPValue + ".txt", gsea.significantResultsToStringList(significant, DEFAULT_PVALUES[i]));					
				numberOfSignificantTerms = significant.size();					
			} else {
				numberOfSignificantTerms = 0;
			}
			significantCount.get(i).append(filterInfo.getTitle()).append("\t").append(numberOfSignificantTerms).append("\n");
			if(numberOfSignificantTerms>0){
				// go graph
				if(filterInfo.getPrefix().equalsIgnoreCase("go")) {
					try {
						String limitMessage = "";
						boolean outOfbounds = createGseaGoGraph(significant,DEFAULT_PVALUES[i],filterInfo);
						if(outOfbounds) limitMessage = " just most 100 significant terms";										
						Item item = new Item("go_graph_significant_" + filterInfo.getName() + "_" + formattedPValue,"go_graph_" + filterInfo.getName() + "_" + formattedPValue + "_graphimage.png",filterInfo.getTitle() + " DAG (significant terms, pvalue<" + formattedPValue + ") " + limitMessage,Item.TYPE.IMAGE,Arrays.asList("SIGNIFICANT,THUMBNAIL,GO_GRAPH_VIZ_JNLP"),new HashMap<String,String>(),"Significant Results." + filterInfo.getTitle());
						item.setContext("pvalue==" + formattedPValue);
						result.getOutputItems().add(4, item);
					} catch(GoGraphException gge){
						Item item = new Item("go_graph_significant_" + filterInfo.getName() + "_" + formattedPValue,"Graph not found",filterInfo.getTitle() + " DAG (significant terms, pvalue<" + formattedPValue + ")",Item.TYPE.MESSAGE,Arrays.asList("ERROR"),new HashMap<String,String>(),"Significant Results." + filterInfo.getTitle());
						item.setContext("pvalue==" + formattedPValue);
						result.getOutputItems().add(4, item);
					}
				}
				// table
				Item item = new Item("significant_" + filterInfo.getName(),"significant_" + filterInfo.getName() + "_" + formattedPValue + ".txt",filterInfo.getTitle() + " significant terms (pvalue<" + formattedPValue + ")",Item.TYPE.FILE,Arrays.asList("SIGNIFICANT","TABLE",gsea.getMethod().toUpperCase() + "_TABLE",filterInfo.getPrefix().toUpperCase() + "_TERM"),new HashMap<String,String>(),"Significant Results." + filterInfo.getTitle());
				item.setContext("pvalue==" + formattedPValue);
				result.getOutputItems().add(4, item);
				
			} else {
				Item item = new Item("significant_" + filterInfo.getName(),"No significant terms for current pvalue " + formattedPValue,filterInfo.getTitle() + " significant terms (pvalue<" + formattedPValue + ")",Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Significant Results." + filterInfo.getTitle()); 
				item.setContext("pvalue==" + formattedPValue);
				result.getOutputItems().add(4, item);
			}
			
			
		}			
		
		//annotations report
		addAnnotationReport(gsea,filterInfo.getTitle());
		
		// save annotations
		if(gsea.getAnnotations()!=null && gsea.getAnnotations().size()>0){
			IOUtils.write(outdir + "/" + annotFileName, gsea.getAnnotations().toString());
			result.addOutputItem(new Item("annot_" + filterInfo.getName(),annotFileName,"Annotations for " + filterInfo.getTitle(),Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
		} else {
			result.addOutputItem(new Item("annot_" + filterInfo.getName(),"no annotations found for input ids","Annotations for " + filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Annotation files"));
		}	
	
	
	}
	
		

	

	
	
	/*
	 * Annotations count
	 */
	private void addAnnotationReport(GeneSetAnalysis fatiscan, String dbName){
		DecimalFormat formatter = new DecimalFormat("#######.##");
		double list1Percentage = ((double)(fatiscan.getAnnotatedCounter())/(double)fatiscan.getIdList().size())*100.0;		
		String listMessage = fatiscan.getAnnotatedCounter() + " of " + fatiscan.getIdList().size() + " (" + formatter.format(list1Percentage) + "%) " + formatter.format(fatiscan.getMeanAnnotationsPerId()) + " annotations/id";
		annotationReport.append(dbName).append("\t").append(listMessage).append("\n");
	}
	private void saveAnnotationsReport() throws IOException{
		IOUtils.write(outdir + "/annotations_per_db.txt", annotationReport.toString());
		result.getOutputItems().add(1,new Item("annotations_per_db","annotations_per_db.txt","Id annotations per DB",Item.TYPE.FILE,Arrays.asList("TABLE","SUMMARY_TABLE"),new HashMap<String,String>(),"Summary"));
	}
	
	
}

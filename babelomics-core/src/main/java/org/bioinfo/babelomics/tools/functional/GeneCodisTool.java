package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.FatiGO;
import org.bioinfo.babelomics.methods.functional.GeneCodis;
import org.bioinfo.babelomics.methods.functional.GeneCodis.correctionFactor;
import org.bioinfo.babelomics.methods.functional.GeneCodis.testFactor;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class GeneCodisTool extends FunctionalProfilingTool{
	
	// list1
	private FeatureData list1;
	
	// list2
	private FeatureData list2;		
	
	
	// other capabilities
	protected int duplicatesMode;
	private StringBuilder duplicatesReport;
	private StringBuilder annotationReport;
	private String list2label;
	
	//private List<Chromosome> chromosomes;
	
	private int support;   //Minimum number of genes required that have to be implicated in a rule to take it into account (default: 3).
	private int supportRandom;//i . i support for random:  Min support taked into account when the algorithm is correcting p-values (usually the same than Support).
	
	private testFactor test;//-t- [1,2] statistical test, 0: Hypergeometric test (default), 1: Chi Square test
	private correctionFactor correction; // -s 0: none	any number < 0: FDR method	any number > 0: Number of permutations to correct p-values with permutations method
	List<String> significantDbs;
	
	@Override
	public void initOptions() {
		// parent options
		super.initOptions();
		
		options.addOption(OptionFactory.createOption("datalist", "the ranked list"));
		options.addOption(OptionFactory.createOption("support", "Minimum number of genes",true));
		options.addOption(OptionFactory.createOption("support-for-random", "Minimum number of genes for correcting p-values",true));
		options.addOption(OptionFactory.createOption("hypergeometric", "",false));
		options.addOption(OptionFactory.createOption("chi-square", "",false));
		options.addOption(OptionFactory.createOption("correction", "default correction",true));
		
		OptionGroup input2 = new OptionGroup();
		input2.setRequired(false);		
			input2.addOption(OptionFactory.createOption("list2", "the file containig the list #2 of genes, or the feature data file",false));
			input2.addOption(OptionFactory.createOption("genome", "compares list #1 with the rest of genome",false,false));
			input2.addOption(OptionFactory.createOption("chromosomes", "the chromosome region to compare with list #1. Use chromosome-start and chromosome-end options tp delimite the region within the chromosome",false));
		options.addOptionGroup(input2);
		
		// extras
		options.addOption(OptionFactory.createOption("duplicates", "to remove duplicated IDs, valid values: never, each, ref. By default, never", false));
		
	}

	
	@Override
	public void prepare() throws IOException, ParseException, InvalidIndexException {
		super.prepare();		

		
		// list 1		
		list1 = new FeatureData(new File(commandLine.getOptionValue("datalist")), true);
		
		// list 2: list to compare
		if(commandLine.hasOption("list2")){
			list2 = new FeatureData(new File(commandLine.getOptionValue("list2")), true);
		} else if(commandLine.hasOption("genome")) {
			this.setRestOfGenome(true);
		} else if(commandLine.hasOption("chromosomes")) {
			throw new ParseException("chromosome reading not yet implemented");
		} else {
			throw new ParseException("No comparison provided, use list2, genome or chromosome options to set your comparison");
		}
		
		//command parameters
		support = Integer.parseInt(commandLine.getOptionValue("support"));  
		supportRandom = Integer.parseInt(commandLine.getOptionValue("support-for-random")); //i
		
		test = testFactor.hypergeometric;		
		if(commandLine.hasOption("chi-square")){
			if (commandLine.hasOption("hypergeometric")){
				test = testFactor.both;
			} else {
				test = testFactor.chiSquare;
			}
		} 
				
		correction = correctionFactor.none; 
		if (commandLine.getOptionValue("correction").equalsIgnoreCase("correction_fdr")){
			correction = correctionFactor.fdr;
			
		} else if (commandLine.getOptionValue("correction").equalsIgnoreCase("correction_permutation")){
			correction = correctionFactor.permutation;
		
		}
		
		// extras
		String duplicates = commandLine.getOptionValue("duplicates", "never");
		if(duplicates.equalsIgnoreCase("each")) duplicatesMode = GeneCodis.REMOVE_EACH;
		if(duplicates.equalsIgnoreCase("ref")) duplicatesMode = GeneCodis.REMOVE_REF;
		if(duplicates.equalsIgnoreCase("all")) duplicatesMode = GeneCodis.REMOVE_ALL;		
		if(duplicates.equalsIgnoreCase("never")) duplicatesMode = GeneCodis.REMOVE_NEVER;

	}
	
	@Override
	protected void execute() {
		try {
			logger.println("Starting GeneCodis...");
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");
			
			// infrared connector			
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));		
			
			// prepare params		
			prepare();
			
			GeneCodis genecodis = null;
			list2label = "List 2";
			
			jobStatus.addStatusMessage("" + ("40"), "preparing Genecodis execution");
			
			
			// list 1			
			List<String> idList1 = list1.getDataFrame().getRowNames();//list1.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list1.getDataFrame().getColumn(0));
			
			// list 2
			List<String> idList2 = null;
			addInputParams();
			if(filterList.size()==0 && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				significantDbs = new ArrayList<String>();
				if(list2!=null)idList2 = list2.getDataFrame().getRowNames();
				annotationReport = new StringBuilder();
				annotationReport.append("#DB").append("\t").append("List1").append("\t").append(list2label).append("\n");
				// run genecodis				
				for(FunctionalFilter filter: filterList) {
					genecodis = doTest(idList1,idList2,filter,dbConnector);
				}				
				if(isYourAnnotations){		
					genecodis = doTest(idList1,idList2,yourAnnotations,dbConnector);
				}
				
				logger.println("Starting saveDuplicatesReport(genecodis);...");
				saveDuplicatesReport(genecodis);
				logger.println("OK...");
				saveAnnotationsReport();
				
				}
			
			}
		catch (Exception e) {
		}
	}
	
	
	private GeneCodis doTest(List<String> idList1, List<String> idList2,	FeatureList<AnnotationItem> yourAnnotations, DBConnector dbConnector) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, IOException {
		
		// db attributes
		FunctionalDbDescriptor filterInfo = new FunctionalDbDescriptor("your_annotation","Your annotations", "your_annotations","Your annotations");
			
		logger.println("Starting doing test your_annotations...");
		
		GeneCodis genecodis = null;
		if(list2!=null) genecodis = new GeneCodis(babelomicsHomePath + "/bin/genecodis/genecodis_bin ",outdir,filterInfo.getName(),idList1, idList2, yourAnnotations, dbConnector, duplicatesMode, support, supportRandom,correction, test);
		else if(isRestOfGenome()) genecodis = new GeneCodis(babelomicsHomePath + "/bin/genecodis/genecodis_bin ",outdir,filterInfo.getName(),idList1, yourAnnotations, dbConnector, duplicatesMode, support, supportRandom, correction, test);
                                                 
		runAndSave(genecodis,filterInfo);
		
		logger.println("addAnnotationReport...");
		significantDbs.add(filterInfo.getPrefix().toUpperCase() + "_TERM");
		addAnnotationReport(genecodis,"Your annotations");
		return genecodis;
	}


	private GeneCodis doTest(List<String> idList1, List<String> idList2, FunctionalFilter filter, DBConnector dbConnector) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, IOException{

		// db attributes				
		FunctionalDbDescriptor filterInfo= new FunctionalDbDescriptor(filter);
		
		//File inputFile = new File(outdir + "/"+name+"_WellFormedInput");
		logger.println("Starting doing test 2 list and filter...");
		
		GeneCodis genecodis = null;
		if(list2!=null) genecodis = new GeneCodis(babelomicsHomePath + "/bin/genecodis/genecodis_bin ",outdir,filterInfo.getName(),idList1, idList2, filter, dbConnector, duplicatesMode, support, supportRandom,correction, test);
		else if(isRestOfGenome()) genecodis = new GeneCodis(babelomicsHomePath + "/bin/genecodis/genecodis_bin ",outdir,filterInfo.getName(),idList1, filter, dbConnector, duplicatesMode, support, supportRandom, correction, test);
		
		runAndSave(genecodis,filterInfo);
		
		significantDbs.add(filterInfo.getPrefix().toUpperCase() + "_TERM");
		logger.println("addAnnotationReport...");
		addAnnotationReport(genecodis,filterInfo.getTitle());
		
		return genecodis;
	}
	
	private void runAndSave(GeneCodis genecodis,FunctionalDbDescriptor filterInfo) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, IOException{
		
		//run run...
		genecodis.run();
		
		//save results		
		saveGenecodisResults(genecodis,filterInfo);
		
	}
	
	private void addAnnotationReport(GeneCodis genecodis, String dbName){
		logger.println("doing addAnnotationReport...");
		DecimalFormat formatter = new DecimalFormat("#######.##");
		double list1Percentage = ((double)(genecodis.getList1AnnotatedCounter())/(double)genecodis.getList1SizeBeforeDuplicates())*100.0;
		double list2Percentage = ((double)(genecodis.getList2AnnotatedCounter())/(double)genecodis.getList2SizeBeforeDuplicates())*100.0;
		String list1Message = genecodis.getList1AnnotatedCounter() + " of " + genecodis.getList1SizeAfterDuplicates() + " (" + formatter.format(list1Percentage) + "%) " + formatter.format(genecodis.getList1MeanAnnotationsPerId()) + " annotations/id";
		String list2Message = genecodis.getList2AnnotatedCounter() + " of " + genecodis.getList2SizeAfterDuplicates() + " (" + formatter.format(list2Percentage) +"%) " + formatter.format(genecodis.getList2MeanAnnotationsPerId()) + " annotations/id";
		annotationReport.append(dbName).append("\t").append(list1Message).append("\t").append(list2Message).append("\n");
		logger.println("ok...");
	}
	
	private void saveAnnotationsReport() throws IOException{
		IOUtils.write(outdir + "/annotations_per_db.txt", annotationReport.toString());
		result.getOutputItems().add(1,new Item("annotations_per_db","annotations_per_db.txt","Id annotations per DB",Item.TYPE.FILE,Arrays.asList("TABLE","SUMMARY_TABLE"),new HashMap<String,String>(),"Summary"));
	}
	
	private void addInputParams(){
		
		// species		
		addInputParam("species", "Species", species);		
		
		// duplicates		
		HashMap<Integer,String> duplicates = new HashMap<Integer, String>();		
		duplicates.put(FatiGO.REMOVE_NEVER, "Never remove");
		duplicates.put(FatiGO.REMOVE_EACH, "Remove separately on each list");
		duplicates.put(FatiGO.REMOVE_REF, "Remove list 1 ids from list 2 (complementary list)");
		duplicates.put(FatiGO.REMOVE_GENOME, "Remove list 1 ids from genome");
		duplicates.put(FatiGO.REMOVE_ALL, "Remove all duplicates (owned and shared duplicates)");		
		addInputParam("duplicates", "Duplicates management", duplicates.get(duplicatesMode));
		
		// others params laterality		
		addInputParam("support", "Minimum number of genes", String.valueOf(support));		
		addInputParam("supportRandom", "Minimum number of genes for correcting p-values", String.valueOf(supportRandom));		
		addInputParam("test", "Statistical test", test.toString());		
		addInputParam("correction", "Correction", correction.toString());
		
	}	
	private void addInputParam(String id, String label, String value){
		result.addOutputItem(new Item(id,value,label,Item.TYPE.MESSAGE,Arrays.asList("INPUT_PARAM"),new HashMap<String,String>(),"Input params"));
	}
	
	private void saveDuplicatesReport(GeneCodis genecodis) throws IOException{		
		DecimalFormat formatter = new DecimalFormat("##.##");		
		
		// list 1	
		int list1Duplicates = genecodis.getList1SizeBeforeDuplicates()-genecodis.getList1SizeAfterDuplicates();		
		double list1DuplicatesPercentage = (double)(list1Duplicates)/(double)genecodis.getList1SizeBeforeDuplicates();		
		String list1DuplicatesMessage = list1Duplicates + " of " + genecodis.getList1SizeBeforeDuplicates() + " (" + formatter.format(list1DuplicatesPercentage) + "%)";		
		
		// list 2		
		int list2Duplicates = genecodis.getList2SizeBeforeDuplicates()-genecodis.getList2SizeAfterDuplicates();		
		double list2DuplicatesPercentage = (double)(list2Duplicates)/(double)genecodis.getList2SizeBeforeDuplicates();		
		String list2DuplicatesMessage = list2Duplicates + " of " + genecodis.getList2SizeBeforeDuplicates() + " (" + formatter.format(list2DuplicatesPercentage) + "%)";
		
		// report		
		duplicatesReport = new StringBuilder();		
		duplicatesReport.append("#Detail").append("\t").append("List 1").append("\t").append(list2label).append("\n");		
		duplicatesReport.append("Number of duplicates").append("\t").append(list1DuplicatesMessage).append("\t").append(list2DuplicatesMessage).append("\n");		
		duplicatesReport.append("Number of finally used ids").append("\t").append(genecodis.getList1SizeAfterDuplicates()).append("\t").append(genecodis.getList2SizeAfterDuplicates()).append("\n");
		IOUtils.write(outdir + "/duplicates.txt", duplicatesReport.toString());
		result.addOutputItem(new Item("duplicates","duplicates.txt","Duplicates management",Item.TYPE.FILE,Arrays.asList("TABLE","SUMMARY_TABLE"),new HashMap<String,String>(),"Summary"));
		
	}
	
	private void saveGenecodisResults(GeneCodis genecodis,FunctionalDbDescriptor filterInfo) throws IOException{
		
		String fileName = filterInfo.getName() + ".txt";
		String annotFileName = filterInfo.getName() + ".annot";
		updateJobStatus("80", "saving results");
				String legend = "Id: the id of the result returned<br />";
				legend += "Items: contains the enriched items<br />";
				legend += "S: contains the number of genes in the input list that have been matched in the enriched set and the total number of genes in the list.<br />";
				legend += "TS: contains the number of genes in the reference list that match with the enriched set and total number of genes in the reference list.<br />";
				legend += "Hyp: pvalue calculated by the hypergeometric function.<br />";
				legend += "Chi: pvalue calculated by the chi-square method.<br />";				
				legend += "Genes: the genes of the input list that match with the enriched items.<br/>";
		if (genecodis.getSignificantTerms()>0){
			String correctionParam = (correction== correctionFactor.none)?"CORRECTEDNONE":"CORRECTED";
			result.addOutputItem(new Item(filterInfo.getName() + "_concurrence.txt", filterInfo.getName() + "_concurrence.txt", "Co-ocurrence annotations results", TYPE.FILE, Arrays.asList("TABLE","GENECODIS_TABLE_"+test.toString().toUpperCase()+"_"+correctionParam,filterInfo.getPrefix().toUpperCase() + "_TERM"), new HashMap<String, String>(), "Co-ocurrence annotations results"));
			result.addOutputItem(new Item(filterInfo.getName() + "_singular.txt", filterInfo.getName() + "_singular.txt", "Singular annotations results", TYPE.FILE, Arrays.asList("TABLE","GENECODIS_TABLE_"+test.toString().toUpperCase()+"_"+correctionParam,filterInfo.getPrefix().toUpperCase() + "_TERM"), new HashMap<String, String>(), "Singular annotations results"));				
			result.addOutputItem(new Item(filterInfo.getName() + "_description", filterInfo.getDescription(),filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("MINI_COMMENT"),new HashMap<String,String>(),"Co-ocurrence annotations results"));
			result.addOutputItem(new Item(filterInfo.getName() + "_description", filterInfo.getDescription(),filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("MINI_COMMENT"),new HashMap<String,String>(),"Singular annotations results"));
		}else if (genecodis.getSignificantTerms()==0){
			result.addOutputItem(new Item(filterInfo.getName() + "_genecodis_file","no significant terms found","Annotations results ",Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Annotations results"));
		} else if (genecodis.getSignificantTerms()==-1){
			result.addOutputItem(new Item(filterInfo.getName() +"_genecodis_file","Error running tool","Annotations results " + filterInfo.getName(),Item.TYPE.MESSAGE,Arrays.asList("ERROR"),new HashMap<String,String>(),"Significant terms"));
		}
		result.addOutputItem(new Item(filterInfo.getName() +"_genecodis_file",legend,"",Item.TYPE.MESSAGE,"Legend"));
		// save annotation
		IOUtils.write(outdir + "/" + annotFileName, genecodis.getAnnotations().toString());
		result.addOutputItem(new Item("annot_" + filterInfo.getName(),annotFileName,"Annotations for " + filterInfo.getTitle(),Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
					
	}
	
}

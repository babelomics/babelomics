package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.GeneSetAnalysisTestResult;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.BiocartaFilter;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.GOSlimFilter;
import org.bioinfo.infrared.funcannot.filter.InterproFilter;
import org.bioinfo.infrared.funcannot.filter.JasparFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
import org.bioinfo.infrared.funcannot.filter.MiRnaTargetFilter;
import org.bioinfo.infrared.funcannot.filter.OregannoFilter;
import org.bioinfo.infrared.funcannot.filter.ReactomeFilter;
import org.bioinfo.math.stats.inference.FisherExactTest;
import org.bioinfo.tool.OptionFactory;

import es.blast2go.prog.graph.GetGraphApi;
import es.blast2go.prog.graph.GoGraphException;
import es.blast2go.prog.graph.MakeGraphDot;


public abstract class FunctionalProfilingTool extends BabelomicsTool {

	public final static double DEFAULT_PVALUE_THRESHOLD = 0.05;
	public final static double[] DEFAULT_PVALUES = {0.1,0.05,0.01,0.005};
	public DecimalFormat pvalueFormatter = new DecimalFormat("#.####");
	public final static int MAX_GOS = 100;
	
//	// JT.2009.09.07
//	private String removeDuplicates;
	
	// input data		
	protected boolean restOfGenome;

	// reports
	protected StringBuilder duplicatesReport;
	protected StringBuilder annotationReport;
	
	// test
	protected int testMode;
		
	// dbs
	protected List<FunctionalFilter> filterList;
	
	// your annotations
	protected boolean isYourAnnotations;
	protected FeatureList<AnnotationItem> yourAnnotations;

	protected List<StringBuilder> significantCount;
	
	// working
	protected int numberOfResults;
	
	public FunctionalProfilingTool (){
		this.numberOfResults = 0;
	}
	
	@Override
	public void initOptions() {

		// commons options
		getOptions().addOption(OptionFactory.createOption("fisher", "the Fisher test mode, valid values: less, greater, two-tailed. By default, two-tailed", false));
		getOptions().addOption(OptionFactory.createOption("go-dag", "Compute DAG for significant GO terms", false));
		
		// GO biological process options
		addGOOptions("bp");
		addGOOptions("cc");
		addGOOptions("mf");

		addGenericOptions("go-slim");
		addGenericOptions("interpro");
		addGenericOptions("kegg");
		addGenericOptions("reactome");
		addGenericOptions("biocarta");
		addGenericOptions("mirna");
		addGenericOptions("jaspar");
		addGenericOptions("oreganno");

		// your annotations
		getOptions().addOption(OptionFactory.createOption("annotations", "Your own annotations",false,true));
		
		
		
	}


	private void addGOOptions(String namespace){
		String namespaceTitle = "biological process";
		if(namespace.equalsIgnoreCase("cc")) namespaceTitle = "cellular component";
		if(namespace.equalsIgnoreCase("mf")) namespaceTitle = "molecular function";
		getOptions().addOption(OptionFactory.createOption("go-" + namespace, "GO " + namespaceTitle + " database",false,false));
//		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-inclusive", "GO " + namespaceTitle + ", inclusive analysis (one per GO level), otherwise joins GO levels",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-min-level", "GO " + namespaceTitle + ", min go level to take into account, default 3",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-max-level", "GO " + namespaceTitle + ", max GO level to take into account, default 15",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-min-num-genes", "GO " + namespaceTitle + ", min number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-max-num-genes", "GO " + namespaceTitle + ", max number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-nannot-domain", "GO " + namespaceTitle + ", computes the number of annotated genes from all genome, otherwise from your input list",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-keyword-filter", "GO " + namespaceTitle + ", enables keywords filter",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-keyword-list", "GO " + namespaceTitle + ", keywords filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-keyword-operator", "GO " + namespaceTitle + ", keywords filter logic: all or any",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-propagation", "GO " + namespaceTitle + ", direct vs propagated annotation",false));
	}
	
	protected void addGenericOptions(String db){
		getOptions().addOption(OptionFactory.createOption(db, db, false, false));
		getOptions().addOption(OptionFactory.createOption(db + "-min-num-genes", db + " min number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption(db + "-max-num-genes", db + " max number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption(db + "-nannot-domain", db + " computes the number of annotated genes from all genome, otherwise from you input list",false,false));
	}
	
	public void prepare() throws IOException, ParseException, InvalidIndexException{

		// species
		setSpecies(commandLine.getOptionValue("species","hsa"));
		
		// fisher test		
		String testMode = commandLine.getOptionValue("fisher", "two-tailed");
		if(testMode.equals("less") || testMode.equals("lower")) setTestMode(FisherExactTest.LESS);
		if(testMode.equals("greater")) setTestMode(FisherExactTest.GREATER);
		if(testMode.equals("two-tailed") || testMode.equals("two_sided") ) setTestMode(FisherExactTest.TWO_SIDED);
		if(getTestMode()==0){
			logger.warn("Unknown fisher test mode. Set Two tailed instead.");
			setTestMode(FisherExactTest.TWO_SIDED);			
		}
		
		filterList = new ArrayList<FunctionalFilter>();
		
		// go bp
		if(commandLine.hasOption("go-bp")) {
			parseGODb(commandLine,"bp");
		}
		// go cc
		if(commandLine.hasOption("go-cc")) {
			parseGODb(commandLine,"cc");
		}
		// go mf
		if(commandLine.hasOption("go-mf")) {
			parseGODb(commandLine,"mf");
		}

		parseGenericDb(commandLine,"go-slim");
		parseGenericDb(commandLine,"interpro");
		parseGenericDb(commandLine,"kegg");
		parseGenericDb(commandLine,"reactome");
		parseGenericDb(commandLine,"biocarta");
		parseGenericDb(commandLine,"mirna");
		parseGenericDb(commandLine,"jaspar");
		parseGenericDb(commandLine,"oreganno");
				
		// species must be provided if any db is selected
		if(commandLine.hasOption("go-bp") || commandLine.hasOption("go-mf") || commandLine.hasOption("go-cc") || commandLine.hasOption("go-slim") || commandLine.hasOption("interpro") || commandLine.hasOption("kegg") || commandLine.hasOption("reactome") || commandLine.hasOption("biocarta") || commandLine.hasOption("mirna") || commandLine.hasOption("jaspar") || commandLine.hasOption("oreganno")){
			if(!commandLine.hasOption("species")) throw new ParseException("species is mandatory if any database specified");
		}
		
		// your annotations
		if(commandLine.hasOption("annotations") && !commandLine.getOptionValue("annotations").equalsIgnoreCase("") && !commandLine.getOptionValue("annotations").equalsIgnoreCase("none")) {
			isYourAnnotations = true;
			FileUtils.checkFile(new File(commandLine.getOptionValue("annotations")));
			
			// read annotation file
			List<String> annots = IOUtils.readLines(commandLine.getOptionValue("annotations"));
			
			// init annotation object
			yourAnnotations = new FeatureList<AnnotationItem>(annots.size());
			
			// init map to detect duplicated annotations
			Map<String,Boolean> uniqAnnots = new HashMap<String, Boolean>();
			
			// read annotations
			String[] fields;
			for(String annot: annots){
				if(!uniqAnnots.containsKey(annot)){
					if(!annot.startsWith("#") && annot.contains("\t")){
						fields = annot.split("\t");
						yourAnnotations.add(new AnnotationItem(fields[0],fields[1]));
						uniqAnnots.put(annot, true);
					}
				}
			}
		
		}
		
		config.load(new FileInputStream(new File(babelomicsHomePath + "/conf/blast2go.properties")));
		
	}

	
	
	public void parseGODb(CommandLine cmdLine, String namespace) throws ParseException{
		if(cmdLine.hasOption("go-" + namespace + "")) {
//			if(cmdLine.hasOption("go-" + namespace + "-inclusive")) {				
//				int min = Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","5"));
//				int max = Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","12"));
//				for(int i=min; i<=max; i++){
//					GOFilter goBpFilter = parseGOFilter(cmdLine, namespace);
//					goBpFilter.setMinLevel(i);
//					goBpFilter.setMaxLevel(i);					
//					filterList.add(goBpFilter);					
//				}
//			} else {
				GOFilter goBpFilter = parseGOFilter(cmdLine, namespace);				
				filterList.add(goBpFilter);
//			}			
		}		
	}
	
	public GOFilter parseGOFilter(CommandLine cmdLine, String namespace) throws ParseException{
		
		// ontology
		String infraredNamespace = "biological_process";		
		if(namespace.equalsIgnoreCase("cc")) infraredNamespace = "cellular_component";
		if(namespace.equalsIgnoreCase("mf")) infraredNamespace = "molecular_function";
		GOFilter goFilter = new GOFilter(infraredNamespace);
		
		// direct vs propagated annotation
		goFilter.setPropagated(true);
		if(cmdLine.hasOption("go-" + namespace + "-propagation") && !cmdLine.getOptionValue("go-" + namespace + "-propagation").equalsIgnoreCase("propagate")){				
			goFilter.setPropagated(false);
		}
		
		// levels
		goFilter.setMinLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","3")));
		goFilter.setMaxLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-level","9")));
		
		// number of annotations
		goFilter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-num-genes","5")));	
		goFilter.setMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-num-genes","500")));
		goFilter.setGenomicNumberOfGenes(true);
		System.err.println("estas? " + cmdLine.getOptionValue("go-" + namespace + "-nannot-domain"));
		if(cmdLine.hasOption("go-" + namespace + "-nannot-domain") && !cmdLine.getOptionValue("go-" + namespace + "-nannot-domain").equals("genome")){
			goFilter.setGenomicNumberOfGenes(false);
		}		
		
		// filter by keywords
		if(cmdLine.hasOption("go-" + namespace + "-keyword-filter")) {
			if(!cmdLine.hasOption("go-" + namespace + "-keyword-list")){
				throw new ParseException("go-" + namespace + "-keyword-list is not defined");
			} else {
				goFilter.addKeywords(StringUtils.toList(cmdLine.getOptionValue("go-" + namespace + "-keyword-list", "")));
				goFilter.setLogicalOperator(cmdLine.getOptionValue("go-" + namespace + "-keyword-operator","all"));
				// 	add childre???
			}
		}
		

		
		return goFilter;
	}
	
	public void parseGenericDb(CommandLine cmdLine, String db){
		if(commandLine.hasOption(db)) {
			FunctionalFilter filter = null;
			
			if(db.equalsIgnoreCase("go-slim")){
				filter = new GOSlimFilter();			
			}			
			if(db.equalsIgnoreCase("interpro")){
				filter = new InterproFilter();		
			}
			if(db.equalsIgnoreCase("kegg")){
				filter = new KeggFilter();						
			}
			if(db.equalsIgnoreCase("reactome")){
				filter = new ReactomeFilter();						
			}
			if(db.equalsIgnoreCase("biocarta")){
				filter = new BiocartaFilter();						
			}
			if(db.equalsIgnoreCase("mirna")){
				filter = new MiRnaTargetFilter();
			}
			if(db.equalsIgnoreCase("jaspar")){
				filter = new JasparFilter();						
			}
			if(db.equalsIgnoreCase("oreganno")){
				filter = new OregannoFilter();						
			}
			if(filter!=null){
				filter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue(db + "-min-num-genes","2")));	
				filter.setMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue(db + "-max-num-genes","500")));		
				if(cmdLine.getOptionValue(db + "-nannot-domain","genome").equals("genome")) filter.setGenomicNumberOfGenes(true); 
				else filter.setGenomicNumberOfGenes(false);
				filterList.add(filter);
			}
			
		}
	}

//	protected String getDBName(Filter filter){
//		String name = StringUtils.randomString();
//		if(filter instanceof GOFilter) {						
//			GOFilter goFilter = (GOFilter) filter.clone();
//			name = "go_" + goFilter.getNamespace() + "_" + goFilter.getMinLevel() + "_" + goFilter.getMaxLevel();
//		} else if(filter instanceof GOSlimFilter) {
//			name = "go-slim";
//		} else if(filter instanceof KeggFilter) {
//			name = "kegg";
//		} else if(filter instanceof InterproFilter) {
//			name = "interpro";
//		} else if(filter instanceof ReactomeFilter) {
//			name = "reactome";
//		} else if(filter instanceof BiocartaFilter) {
//			name = "biocarta";
//		} else if(filter instanceof JasparFilter) {
//			name = "jaspar";
//		} else if(filter instanceof OregannoFilter) {
//			name = "oreganno";
//		}
//		return name;
//	}
//	
//	protected String getDBPrefix(Filter filter){
//		String name = StringUtils.randomString();
//		if(filter instanceof GOFilter) {		
//			name = "go";
//		} else if(filter instanceof GOSlimFilter) {
//			name = "go-slim";
//		} else if(filter instanceof KeggFilter) {
//			name = "kegg";
//		} else if(filter instanceof InterproFilter) {
//			name = "interpro";
//		} else if(filter instanceof ReactomeFilter) {
//			name = "reactome";
//		} else if(filter instanceof BiocartaFilter) {
//			name = "biocarta";
//		} else if(filter instanceof JasparFilter) {
//			name = "jaspar";
//		} else if(filter instanceof OregannoFilter) {
//			name = "oreganno";
//		}
//		return name;
//	}
//	
//	protected String getDBTitle(Filter filter){
//		String title = "Untitled",levels;
//		if(filter instanceof GOFilter) {						
//			GOFilter goFilter = (GOFilter) filter.clone();
//			if(goFilter.getMinLevel()==goFilter.getMaxLevel()) {
//				levels = "(level " + goFilter.getMinLevel() + ")";
//			} else{
//				levels = "(levels from " + goFilter.getMinLevel() + " to " + goFilter.getMaxLevel() + ")"; 
//			}						
//			title = "GO " + goFilter.getNamespace().replace("_", " ") + " " + levels;
//		} else if(filter instanceof KeggFilter) {
//			title = "GOSlim GOA";
//		} else if(filter instanceof KeggFilter) {
//			title = "Kegg";
//		} else if(filter instanceof InterproFilter) {
//			title = "InterPro";
//		}  else if(filter instanceof ReactomeFilter) {
//			title = "Reactome";
//		}	else if(filter instanceof BiocartaFilter) {
//			title = "Biocarta";
//		}	else if(filter instanceof JasparFilter) {
//			title = "Jaspar";
//		}	else if(filter instanceof OregannoFilter) {
//			title = "Oreganno";
//		}	
//		return title;
//	}

	
	protected List<String> testResultToStringList(List<TwoListFisherTestResult> testResult){
		return testResultToStringList(testResult,true);
	}
	
	protected List<String> testResultToStringList(List<TwoListFisherTestResult> testResult, boolean header){
		List<String> result = new ArrayList<String>();
		if(header) result.add(TwoListFisherTestResult.header());
		for(int i=0; i<testResult.size(); i++){			
			result.add(testResult.get(i).toString());
		}
		return result;
	}
		
	protected boolean createGoGraph(List<TwoListFisherTestResult> raw, double pvalue, FunctionalDbDescriptor filterInfo) throws GoGraphException{
		List<GeneSetAnalysisTestResult> result = new ArrayList<GeneSetAnalysisTestResult>(raw.size());
		for(TwoListFisherTestResult test: raw){			
			GeneSetAnalysisTestResult gseaTest = new GeneSetAnalysisTestResult(test);
			result.add(gseaTest);
		}		
		return createGseaGoGraph(result,pvalue,filterInfo);
	}
	
	protected boolean createGseaGoGraph(List<GeneSetAnalysisTestResult> significants, double pvalue, FunctionalDbDescriptor filterInfo) throws GoGraphException{
		boolean outOfBounds = false;
		
		DecimalFormat pvalueLabelFormatter = new DecimalFormat("#.####E0");
		String prefix = "go_graph_" + filterInfo.getName() + "_" + pvalueFormatter.format(pvalue);
				
		// preparing association file
		StringBuilder association = new StringBuilder();		
		double ratio,intensity;
		
		// get pvalue rank
		List<Double> pvalues = new ArrayList<Double>();
		for(GeneSetAnalysisTestResult result: significants){
			pvalues.add(result.getAdjPValue());
		}		
		int[] pvalueOrder = ListUtils.order(pvalues);
		
		// sort and filter		
		List<GeneSetAnalysisTestResult> sortedSignificants = new ArrayList<GeneSetAnalysisTestResult>();
		int pos;
		for(int i=0; i<pvalueOrder.length; i++){
			pos = pvalueOrder[i];
			sortedSignificants.add(significants.get(pos));
//			System.err.println(significants.get(pos).getAdjPValue());
			if(i>MAX_GOS){
				outOfBounds=true;
				break;
			}
		}
				
		if(outOfBounds) System.err.println("max gos exceeded");
		
		for(GeneSetAnalysisTestResult result: sortedSignificants){			
			
			//color=""+ (result.getList1Percentage()/(result.getList1Percentage()+result.getList2Percentage()));
			ratio = result.getAdjPValue()/pvalue;
			if(result.getList1Percentage()> result.getList2Percentage()){
				intensity = (0.5 + (1-ratio)/2.0);	
			} else {
				intensity = ratio/2.0;
			}
			//System.err.println(result.getList1Percentage() + ":" + result.getList2Percentage() + " " + result.getAdjPValue() + " " + intensity);
			
			association.append(result.getTerm()).append("\t").append(result.getTerm()).append("\t").append(intensity).append("\t").append("adj.pvalue=").append(pvalueLabelFormatter.format(result.getAdjPValue())).append("\n");
		}
		
		try {
			
			// create graph directory
			FileUtils.createDirectory(outdir + "/graphs");
			
			// save association file
			IOUtils.write(outdir + "/graphs/" + prefix  + "_association.txt", association.toString());
			
			// namespace
			String namespace = "b";
			if(filterInfo.getName().contains("molecular_function")) namespace = "m";
			if(filterInfo.getName().contains("cellular_component")) namespace = "c";
			System.err.println("namespace: " + namespace);
			
			// init graph api
			GetGraphApi graph = new GetGraphApi(outdir + "/graphs/",prefix,prefix  + "_association.txt",namespace,0,MakeGraphDot.BYDESCWITHLABEL,0.6,0,"orange",12,filterInfo.getTitle());

			// setting server params
			graph.setDownloader(config.getProperty("JNLP_DOWNLOADER_HOST_NAME"));
			graph.setDataBase(config.getProperty("BLAST2GO_HOST_NAME"),config.getProperty("BLAST2GO_DB_NAME"),config.getProperty("BLAST2GO_DB_USER"), config.getProperty("BLAST2GO_DB_PASSWORD"));
			
			System.err.println("graph : " + config.getProperty("JNLP_DOWNLOADER_HOST_NAME") +" " + config.getProperty("BLAST2GO_HOST_NAME") +" " + config.getProperty("BLAST2GO_DB_NAME") +" " + config.getProperty("BLAST2GO_DB_USER") +" " + config.getProperty("BLAST2GO_DB_PASSWORD"));
			// run
			graph.run();
			
			// copy files
			String imagePrefix = "go_graph_" + filterInfo.getName() + "_" + pvalueFormatter.format(pvalue) + "_graphimage";
				// png
			FileUtils.touch(new File(outdir + "/" + imagePrefix + ".png"));
			FileUtils.copy(outdir + "/graphs/" + imagePrefix + ".png", outdir + "/" + imagePrefix + ".png");
				// svg
			FileUtils.touch(new File(outdir + "/" + imagePrefix + ".png.svg"));
			FileUtils.copy(outdir + "/graphs/" + imagePrefix + ".svg", outdir + "/" + imagePrefix + ".png.svg");

		} catch (Exception e) {
			e.printStackTrace();
			throw new GoGraphException(e.getMessage());
		}	
		
		return outOfBounds;
	}

	
	protected Map<String, Integer> getYourAnnotationsTermSizes(){
		logger.print("Computing your annotations term sizes...");
		Map<String, Integer> termSizes = new HashMap<String, Integer>();
		String termId;
		int cont;
		for(AnnotationItem annotation: yourAnnotations){
			termId = annotation.getFunctionalTermId();
			if(termSizes.containsKey(termId)){
				cont = termSizes.get(termId);
				termSizes.put(termId,cont+1);
			} else {
				termSizes.put(termId,1);
			}
		}
		logger.println("OK");
		return termSizes;
	}
	
	
	public void setRestOfGenome(boolean restOfGenome) {
		this.restOfGenome = restOfGenome;
	}

	public boolean isRestOfGenome() {
		return restOfGenome;
	}

	public void setTestMode(int testMode) {
		this.testMode = testMode;
	}

	public int getTestMode() {
		return testMode;
	}


}

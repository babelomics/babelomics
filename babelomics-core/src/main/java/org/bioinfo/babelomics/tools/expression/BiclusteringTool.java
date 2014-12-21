package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.bioinfo.babelomics.methods.expression.clustering.Biclustering;
import org.bioinfo.babelomics.methods.expression.clustering.Biclustering.MISSING_VALUES_FILLING;
import org.bioinfo.babelomics.methods.expression.clustering.Biclustering.MISSING_VALUES_SCHEME;
import org.bioinfo.babelomics.methods.expression.clustering.Biclustering.SORTING_CRITERIA;
import org.bioinfo.babelomics.methods.expression.clustering.Biclustering.SORTING_ORDER;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

import smadeira.biclustering.Bicluster;
import smadeira.biclustering.CCC_Bicluster;
import smadeira.biclustering.CCC_Bicluster_SignChanges;

public class BiclusteringTool extends BabelomicsTool {

	private PrintWriter out;
	private PrintWriter stats;
	
	private MISSING_VALUES_SCHEME missingValuesScheme;
	private MISSING_VALUES_FILLING missingValuesFillingScheme;
	private int missingValuesNeighbors;
	
	private boolean signChanges;
	private SORTING_CRITERIA sortingCriteria;
	private SORTING_ORDER sortingOrder;
	
	private boolean filterByconstantPattern;	
	private boolean filterByMinNumberOfGenes;
	private boolean filterByMinNumberOfConditions;
	private boolean filterByMaxMSR;
	private boolean filterByMaxPValue;
	private boolean filterByOverlapping;
	private int minNumberOfGenes;
	private int minNumberOfConditions;
	private float maxMSR;
	private float maxPValue;
	private float maxOverlappingPercentage;
	
	private Biclustering biclustering;
	
	private boolean generateCharts;
	
	@Override
	public void initOptions() {
		
		// data
		getOptions().addOption(OptionFactory.createOption("data", "the data matrix"));
		
		// missing values
		getOptions().addOption(OptionFactory.createOption("missing-values", "[remove,fill,jump] missing values treatment (default is jump)",false,true));
		getOptions().addOption(OptionFactory.createOption("missing-values-filling", "[neighbors,gene] fill missing values with the neighbors mean or the gene mean (default gene)",false,true));
		getOptions().addOption(OptionFactory.createOption("missing-values-filling-neighbors", "[int] number of neighbors used to fill (default 1)",false,true));
		
		// sign changes
		getOptions().addOption(OptionFactory.createOption("sign-changes", "compute signchanges",false,false));
		
		// filtering results
		getOptions().addOption(OptionFactory.createOption("filter-constant-pattern", "filter constant patterns (genes that do not change over time)",false,false));		
		getOptions().addOption(OptionFactory.createOption("min-number-of-genes", "filter biclusters with less than min number of genes",false,true));
		getOptions().addOption(OptionFactory.createOption("min-number-of-conditions", "filter biclusters with less than min number of conditions",false,true));
		getOptions().addOption(OptionFactory.createOption("max-msr", "filter biclusters with msr higher than max MSR",false,true));
		getOptions().addOption(OptionFactory.createOption("max-pvalue", "filter biclusters with pvalue higher than max pvalue",false,true));
		getOptions().addOption(OptionFactory.createOption("filter-overlapping", "filter biclusters overlapping more than a percentage",false,true));
		
		// sorting results
		getOptions().addOption(OptionFactory.createOption("sort", "[genes,conditions,size,pvalue,msr] sorting criteria (default pvalue)",false,true));
		getOptions().addOption(OptionFactory.createOption("order", "[ascending,descencing] sorting order (default ascending, just for genes,conditions and size values of sorting param)",false,true));
		
		// charts
		getOptions().addOption(OptionFactory.createOption("charts", "generate charts",false,false));
		
	}

	
	@SuppressWarnings("rawtypes")
	@Override
	protected void execute() {
						
		try {

			// init
			out = new PrintWriter(new File(outdir + "/biclusters.txt"));
			stats = new PrintWriter(new File(outdir + "/stats.txt"));
			
			// input data
			String geneExpressionMatrixFile = commandLine.getOptionValue("data");
			if(geneExpressionMatrixFile==null || !(new File(geneExpressionMatrixFile).exists())) throw new Exception("data file not found");
			
			// missing values			
			missingValuesScheme = MISSING_VALUES_SCHEME.JUMP;
			if(commandLine.hasOption("missing-values")){
				// missing values scheme
				if(commandLine.getOptionValue("missing-values").equalsIgnoreCase("remove")) missingValuesScheme = MISSING_VALUES_SCHEME.REMOVE;
				else if(commandLine.getOptionValue("missing-values").equalsIgnoreCase("jump")) missingValuesScheme = MISSING_VALUES_SCHEME.JUMP;
				else if(commandLine.getOptionValue("missing-values").equalsIgnoreCase("fill")) 	{
					// filling scheme
					missingValuesScheme = MISSING_VALUES_SCHEME.FILL;
					missingValuesFillingScheme = MISSING_VALUES_FILLING.GENE;
					if(commandLine.hasOption("missing-value-filling")) {
						if(commandLine.getOptionValue("missing-value-filling").equalsIgnoreCase("neighbors")) {
							//neighbor filling
							missingValuesFillingScheme = MISSING_VALUES_FILLING.NEIGHBORS;
							missingValuesNeighbors = 1;
							if(commandLine.hasOption("missing-value-filling-neighbors")) {
								try{
									missingValuesNeighbors = Integer.parseInt(commandLine.getOptionValue("missing-value-filling-neighbors"));
								} catch (NumberFormatException e){
									throw new Exception("invalid number of neighbors");
								}
							}
						} else if(commandLine.getOptionValue("missing-value-filling").equalsIgnoreCase("gene")) {
							missingValuesFillingScheme = MISSING_VALUES_FILLING.GENE;
						} else throw new Exception("invalid filling value for missing values [neighbors,gene]");
					}	
					
				} else throw new Exception("invalid missing values scheme [remove,fill,jump]");
			}
			
			
			// missing values
			if(commandLine.hasOption("sign-changes")) signChanges = true;
			
			//Post-Processing Results: Sorting (default or user defined) -> Filtering (optional and user defined)
			
			// sorting
			  // sorting criteria
			sortingCriteria = SORTING_CRITERIA.BY_PVALUE;
			if(commandLine.hasOption("sort")){
				String sort = commandLine.getOptionValue("sort");
				if (sort.equalsIgnoreCase("genes")) sortingCriteria = SORTING_CRITERIA.BY_NUMBER_OF_GENES;
				else if (sort.equalsIgnoreCase("conditions")) sortingCriteria = SORTING_CRITERIA.BY_NUMBER_OF_CONDITIONS;
				else if (sort.equalsIgnoreCase("size")) sortingCriteria = SORTING_CRITERIA.BY_SIZE;
				else if(sort.equalsIgnoreCase("msr")) sortingCriteria = SORTING_CRITERIA.BY_MSR;
				else if(sort.equalsIgnoreCase("pvalue")) sortingCriteria = SORTING_CRITERIA.BY_PVALUE;				
				else throw new Exception("invalid sorting criteria [genes,conditions,size,msr,pvalue]");
			}			
			  // sorting order
			sortingOrder = SORTING_ORDER.ASCENDING;
			if(commandLine.hasOption("order")){
				String order = commandLine.getOptionValue("order");
				if (order.equalsIgnoreCase("ascending")) sortingOrder = SORTING_ORDER.ASCENDING;
				else if (order.equalsIgnoreCase("descending")) sortingOrder = SORTING_ORDER.DESCENDING;
				else throw new Exception("invalid sorting order [ascending,descending]");
				if((sortingOrder==SORTING_ORDER.DESCENDING) && (sortingCriteria==sortingCriteria.BY_MSR ||sortingCriteria==sortingCriteria.BY_PVALUE)){
					logger.println("you don't kwown anything (go home): order must be ascending");
				}
			}
			
			// filtering
			  // constant pattern
			if(commandLine.hasOption("filter-constant-pattern")){
				filterByconstantPattern = true;
			}			
			  // min number of genes
			if(commandLine.hasOption("min-number-of-genes")){
				try {
					int value = Integer.parseInt(commandLine.getOptionValue("min-number-of-genes"));
					minNumberOfGenes = value;
					filterByMinNumberOfGenes = true;
				} catch(NumberFormatException e){
					throw new Exception("invalid min-number-of-genes value");
				}
			}
			  // min number of conditions
			if(commandLine.hasOption("min-number-of-conditions")){
				try {
					int value = Integer.parseInt(commandLine.getOptionValue("min-number-of-conditions"));
					minNumberOfConditions = value;
					filterByMinNumberOfConditions = true;
				} catch(NumberFormatException e){
					throw new Exception("invalid min-number-of-conditions value");
				}
			}
			  // min number of conditions
			if(commandLine.hasOption("max-msr")){
				try {
					float value = Float.parseFloat(commandLine.getOptionValue("max-msr"));
					maxMSR = value;
					filterByMaxMSR = true;
				} catch(NumberFormatException e){
					throw new Exception("invalid max-msr value");
				}
			}
			  // pvalue
			if(commandLine.hasOption("max-pvalue")){
				try {
					float value = Float.parseFloat(commandLine.getOptionValue("max-pvalue"));
					if(value<0 || value>1) throw new Exception("pvalue must be between 0 and 1");
					maxPValue = value;
					filterByMaxPValue = true;
				} catch(NumberFormatException e){
					throw new Exception("invalid max-pvalue value");
				}  
			}
			  // overlapping
			if(commandLine.hasOption("filter-overlapping")){
				try {
					float value = Float.parseFloat(commandLine.getOptionValue("filter-overlapping"));
					maxOverlappingPercentage = value;
					filterByOverlapping = true;
				} catch(NumberFormatException e){
					throw new Exception("invalid overlapping percentage value");
				}	
			}
			// charts
			if(commandLine.hasOption("charts")) generateCharts = true;
			
			
			addInputItems();
			
			// init biclustering
			biclustering = new Biclustering(geneExpressionMatrixFile,missingValuesScheme,signChanges);
			
			// missing values
			if(missingValuesScheme==MISSING_VALUES_SCHEME.FILL){
				biclustering.setMissingValuesFillingScheme(missingValuesFillingScheme);
				if(missingValuesFillingScheme==MISSING_VALUES_FILLING.NEIGHBORS) biclustering.setMissingValuesNeighbors(missingValuesNeighbors);
			}
						
			  // sorting 
			biclustering.setSortingCriteria(sortingCriteria);
			biclustering.setSortingOrder(sortingOrder);
			
			
			ArrayList<ArrayList> filters = new ArrayList<ArrayList>();
			
			// filtering
			if(filterByconstantPattern) 
			{
				biclustering.setFilterByconstantPattern(true);
				ArrayList filter = new ArrayList();
				filter.add("Constant pattern");
				filter.add("applied");
				filters.add(filter);
			
			}
			if(filterByMinNumberOfGenes) {
				biclustering.setFilterByMinNumberOfGenes(true);
				biclustering.setMinNumberOfGenes(minNumberOfGenes);
				
				ArrayList filter = new ArrayList();
				filter.add("Minimum number of genes");
				filter.add(minNumberOfGenes);
				filters.add(filter);
			}
			if(filterByMinNumberOfConditions) {
				biclustering.setFilterByMinNumberOfConditions(true);
				biclustering.setMinNumberOfConditions(minNumberOfConditions);
				
				
				ArrayList filter = new ArrayList();
				filter.add("Minimum number of condition");
				filter.add(minNumberOfConditions);
				filters.add(filter);
			}
			if(filterByMaxMSR) {
				biclustering.setFilterByMaxMSR(true);
				biclustering.setMaxMSR(maxMSR);
			}			
			if(filterByMaxPValue) {
				biclustering.setFilterByMaxPValue(true);
				biclustering.setMaxPValue(maxPValue);
				
				ArrayList filter = new ArrayList();
				filter.add("Maximum p-value");
				filter.add(maxPValue);
				filters.add(filter);
			}
			if(filterByOverlapping) {
				biclustering.setFilterByOverlapping(true);
				biclustering.setMaxOverlappingPercentage(maxOverlappingPercentage);
				ArrayList filter = new ArrayList();
				filter.add("Maximum overlapping percentage");
				filter.add(maxOverlappingPercentage);
				filters.add(filter);
				
			}
			
			// set logger
			biclustering.setLogger(logger);
			
			// run biclustering
			try{
				biclustering.run();
			} catch (Exception e){
				if(biclustering.getBiclusters()!=null && biclustering.getBiclusters().size()==0 && biclustering.getNumberOfOriginalBiclusters()>0){
					result.addOutputItem(new Item("message","you run out biclusters: filtering options are too restrictive","Error",Item.TYPE.MESSAGE,Arrays.asList("ERROR"),new HashMap(),"Output"));
					throw new Exception("you run out biclusters: filtering options are too restrictive");
				} else {
					e.printStackTrace();
					result.addOutputItem(new Item("message","Sorry, a critical error has ocurred. Please check data format and filtering params","Error",Item.TYPE.MESSAGE,Arrays.asList("ERROR"),new HashMap(),"Output"));
					throw new Exception("you have a big problem");
				}
			}
						
			
			
			//input params
			
			biclustering.setSortingCriteria(sortingCriteria);
			biclustering.setSortingOrder(sortingOrder);
			
			
			//input params : filter
			for (ArrayList arrayList : filters) {
			
				result.addOutputItem(new Item(arrayList.get(0).toString() ,arrayList.get(1).toString(),arrayList.get(0).toString(),Item.TYPE.MESSAGE,Arrays.asList(""),new HashMap(),"Input parameters"));
			}
			
			result.addOutputItem(new Item("sorting_input" ,String.valueOf(sortingCriteria), "Sorting Criteria",Item.TYPE.MESSAGE,Arrays.asList(""),new HashMap(),"Input parameters"));
			
			
			
			
			// output results
			stats.println("Number of biclusters: " + biclustering.getNumberOfBiclusters());			
			stats.println("Number of original biclusters: " + biclustering.getNumberOfOriginalBiclusters());			
			stats.println("Number of filtered biclusters: " + biclustering.getNumberOfFilteredBiclusters());
		
			logger.print("saving biclusters...");
			
			
			  // number of biclusters
			String numberOfBiclustersMessage = "" + biclustering.getNumberOfBiclusters();
			//if(filterByconstantPattern || filterByMinNumberOfGenes || filterByMinNumberOfConditions || filterByMaxMSR || filterByMaxPValue || filterByOverlapping){
			//	numberOfBiclustersMessage += " of "  + biclustering.getNumberOfOriginalBiclusters() + " unfiltered biclusters (filtered " + StringUtils.decimalFormat(100.0*((double)biclustering.getNumberOfFilteredBiclusters()/(double)biclustering.getNumberOfOriginalBiclusters()), "#0.00") + "%)";         
			//} 
			
			
			
			String original = String.valueOf(biclustering.getNumberOfOriginalBiclusters());
			String filtered = String.valueOf(biclustering.getNumberOfFilteredBiclusters()) + " ( " + StringUtils.decimalFormat(100.0*((double)biclustering.getNumberOfFilteredBiclusters()/(double)biclustering.getNumberOfOriginalBiclusters()), "#0.00") + "%)";
			
			
			result.addOutputItem(new Item("number_of_biclusters",numberOfBiclustersMessage,"Number of biclusters",Item.TYPE.MESSAGE,Arrays.asList(""),new HashMap(),"Summary"));
			result.addOutputItem(new Item("number_of_original_biclusters", original,"Number of original biclusters",Item.TYPE.MESSAGE,Arrays.asList(""),new HashMap(),"Summary"));
			result.addOutputItem(new Item("number_of_filtered_biclusters", filtered,"Number of filtered biclusters",Item.TYPE.MESSAGE,Arrays.asList(""),new HashMap(),"Summary"));
			
			
			
			
			
			  // biclusters
			result.addOutputItem(new Item("biclusters","biclusters.txt","Biclusters file",Item.TYPE.FILE,Arrays.asList("BICLUSTERING"),new HashMap(),"Biclusters"));
			
			
			
			
			
			for(Bicluster bicluster:biclustering.getBiclusters()){
				printBicluster((CCC_Bicluster)bicluster);
			}
			logger.println("OK");

			if(generateCharts){
				logger.print("saving charts...");
				for(Bicluster bicluster:biclustering.getBiclusters()){					
					saveBiclusterCharts((CCC_Bicluster)bicluster);
					result.addOutputItem(new Item("bicluster_expression_" + bicluster.getID(),"bicluster_expression_" + bicluster.getID() + ".png","Bicluster " + bicluster.getID(),Item.TYPE.IMAGE,Arrays.asList(""),new HashMap(),"Expression Charts"));
					result.addOutputItem(new Item("bicluster_pattern_" + bicluster.getID(),"bicluster_pattern_" + bicluster.getID() + ".png","Bicluster " + bicluster.getID(),Item.TYPE.IMAGE,Arrays.asList(""),new HashMap(),"Pattern Charts"));
				}
				logger.println("OK");
			}
						
			out.close();
			stats.close();
			
			logger.println("");
			logger.println("");

			
			
			
			
		} catch (Exception e) {			
			e.printStackTrace();
		}  

	}
	
	
	private void addInputItems(){
		
		//result.addOutputItem(new Item())
		
	}
	
	private void saveBiclusterCharts(CCC_Bicluster bicluster){		
		bicluster.printColorChartGeneExpressionBiclusterConditions_ToFile(outdir,"bicluster_expression_" + bicluster.getID() + ".png",400,400,true,false);		
		bicluster.printBiclusterPatternChart_ToFile(outdir,"bicluster_pattern_" + bicluster.getID() + ".png",400,400,true,false);
	}
	
	private void printBicluster(CCC_Bicluster bicluster){
				
		// get biclustering data
		String[] geneNames = bicluster.getGenesNames();
		String[] conditionNames = bicluster.getConditionsNames();
		float[][] expressionMatrix = bicluster.getBiclusterExpressionMatrix();
		
		// print bicluster
		out.println(">>Bicluster " + bicluster.getID());
		//out.println("+ " + bicluster.getNumberOfGenes() + " genes X " + bicluster.getNumberOfConditions() + " conditions (size = " + (bicluster.getNumberOfGenes()*bicluster.getNumberOfConditions()) + ")");
		  // basic info
		out.println("+ genes = " + bicluster.getNumberOfGenes());
		out.println("+ conditions = " + bicluster.getNumberOfConditions());
		out.println("+ size = " + (bicluster.getNumberOfGenes()*bicluster.getNumberOfConditions()));
		  // pattern
		if(signChanges){
			out.println("+ pattern = " + getPrettyPattern(bicluster.getBiclusterExpressionPattern()) + "/" + getPrettyPattern(((CCC_Bicluster_SignChanges)bicluster).getBiclusterExpressionPattern_SignChanges()));
		} else {
			out.println("+ pattern = " + getPrettyPattern(bicluster.getBiclusterExpressionPattern()));
		}
		  // pvalue
		if(sortingCriteria==SORTING_CRITERIA.BY_PVALUE || filterByMaxPValue) out.println("+ pvalue = " + bicluster.getBiclusterPatternWithColumns_pValue(1));		
		if(sortingCriteria==SORTING_CRITERIA.BY_MSR || filterByMaxMSR) out.println("+ MSR = " + bicluster.computeMSR());
		
		out.println();		
		out.println("GENES\t" + ArrayUtils.toString(conditionNames, "\t"));
		for(int i=0; i<bicluster.getNumberOfGenes(); i++){
			out.print(geneNames[i] + "\t");
			for(int j=0; j<bicluster.getNumberOfConditions(); j++){
				 out.print(expressionMatrix[i][j]);
				 if(j<(bicluster.getNumberOfConditions()-1)) out.print("\t");				 
			}
			out.println();
		}		
		out.println();
		out.println();
			
		
	}
	
	private String getPrettyPattern(char[] pattern){
		return Arrays.toString(pattern).replaceAll(",", "");
	}
	
}

package org.bioinfo.babelomics.methods.expression.clustering;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.log.Logger;

import smadeira.biclustering.Bicluster;
import smadeira.biclustering.CCC_Biclustering;
import smadeira.biclustering.CCC_Biclustering_Jumping_MissingValues;
import smadeira.biclustering.CCC_Biclustering_SignChanges;
import smadeira.biclustering.CCC_Biclustering_SignChanges_Jumping_MissingValues;
import smadeira.biclustering.DiscretizedExpressionMatrix;
import smadeira.biclustering.ExpressionMatrix;
import smadeira.biclustering.PostProcessedBiclusteringInDiscretizedMatrix;
import smadeira.biclustering.PreProcessedExpressionMatrix;

public class Biclustering {

	public static final float DEFAULT_MISSING_VALUE = Integer.MIN_VALUE;
	public static enum MISSING_VALUES_FILLING {NEIGHBORS,GENE};
	public static enum MISSING_VALUES_SCHEME {REMOVE,FILL,JUMP};
		
	// sort
	public static enum SORTING_CRITERIA {BY_NUMBER_OF_GENES,BY_NUMBER_OF_CONDITIONS,BY_SIZE,BY_PVALUE,BY_MSR};
	public static enum SORTING_ORDER {ASCENDING,DESCENDING};
		
	// input as file
	private String matrixFile;
	
	// input as data
	private String[] geneNames;
	private String[] conditionNames;
	private float[][] expressionMatrix;
	
	// missing values
	private MISSING_VALUES_SCHEME missingValuesScheme;
	private MISSING_VALUES_FILLING missingValuesFillingScheme;
	private int missingValuesNeighbors;
	
	// sign changes
	private boolean signChanges;
	
	// result
	private ArrayList<Bicluster> biclusters;
	
	// post processing
	private SORTING_CRITERIA sortingCriteria;
	private SORTING_ORDER sortingOrder;
	
	// filter
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
	
	// others	
	private int numberOfOriginalBiclusters;
	private int numberOfBiclusters;
	private int numberOfFilteredBiclusters;
	private Logger logger;
	
	public Biclustering(String matrixFile, MISSING_VALUES_SCHEME missingValuesScheme, boolean signChanges){
		this.matrixFile = matrixFile;
		this.missingValuesScheme = missingValuesScheme;
		this.signChanges = signChanges;
		initDefaultParams();		
	}
	
	public Biclustering(String[] geneNames, String[] conditionNames, float[][] expressionMatrix, MISSING_VALUES_SCHEME missingValuesScheme, boolean signChanges){
		this.geneNames = geneNames;
		this.conditionNames = conditionNames;
		this.expressionMatrix = expressionMatrix;
		this.missingValuesScheme = missingValuesScheme;
		this.signChanges = signChanges;
		initDefaultParams();
	}
	
	private void initDefaultParams(){
		this.sortingCriteria = SORTING_CRITERIA.BY_PVALUE;
		this.sortingOrder = SORTING_ORDER.ASCENDING;	
		this.minNumberOfGenes = 2;
		this.minNumberOfConditions = 1;
		this.maxMSR = 0.5f;
		this.maxPValue = 0.05f;	
		this.numberOfOriginalBiclusters = 0;
		this.missingValuesFillingScheme = MISSING_VALUES_FILLING.GENE;
	}
	
	public void run() throws Exception{
	
		// init logger
		if(logger==null) logger = new Logger();
		
		// load data
		if(matrixFile!=null){
			loadData();
		}
							
//        float[][] expressionMatrix2 = {
//                {
//                Integer.MIN_VALUE, (float) 0.73, (float) - 0.54, (float) 0.45,
//                (float) 0.25}
//                , {
//                (float) - 0.34, (float) 0.46, (float) 0.5, (float) 0.76,
//                (float) - 0.44}
//                , {
//                Integer.MIN_VALUE, (float) 0.17, Integer.MIN_VALUE, (float) 0.44,
//                (float) - 0.11}
//                , {
//                (float) 0.70, 2.5f, (float) - 0.41, (float) 0.33,
//                (float) 0.35}
//                , {
//                (float) 0.70, Integer.MIN_VALUE, (float) 0.70, (float) - 0.33,
//                (float) 0.75}
//            };
//        
//        String[] genesNames2 = {
//	            "G1", "G2", "G3", "G4", "G5"};
//	    String[] conditionsNames2 = {
//	            "C1", "C2", "C3", "C4", "C5"};
//		
//	    expressionMatrix=expressionMatrix2;
//	    geneNames = genesNames2;
//	    conditionNames = conditionsNames2;
//		
		// expression matrix
		logger.print("Loading expression matrix...");
		ExpressionMatrix m = new ExpressionMatrix(geneNames, conditionNames, expressionMatrix,DEFAULT_MISSING_VALUE,null,null,null);
		logger.println("OK");
	
		// preprocessed matrix
		logger.print("Preprocessing matrix...");
		PreProcessedExpressionMatrix pm = new PreProcessedExpressionMatrix(m);			
		if(missingValuesScheme==MISSING_VALUES_SCHEME.REMOVE) {
			pm.removeGenesWithMissingValues();	
			if(pm.getNumberOfGenes()==0){
				throw new Exception("run out of genes: try fill or jump missing values");
			}
		} else if(missingValuesScheme==MISSING_VALUES_SCHEME.FILL){
			if(missingValuesFillingScheme==MISSING_VALUES_FILLING.GENE){
				pm.fillMissingsAverageGeneExpression();	
			} else {
				pm.fillMissingsNeighborsAverageGeneExpression(missingValuesNeighbors);
			}
		}
		pm.normalizeGeneExpressionMatrixByGene();
		logger.println("OK");
		
		// discretized matrix
		logger.print("Discretizing matrix...");
		DiscretizedExpressionMatrix dm = new DiscretizedExpressionMatrix(pm);			
		dm.computeDiscretizedMatrix_VariationsBetweenTimePoints((float) 1.0, 'D','N','U');
		logger.println("OK");
		
		// biclustering
		logger.print("Performing biclustering...");
		CCC_Biclustering biclustering = null;
		if(missingValuesScheme==MISSING_VALUES_SCHEME.JUMP) {
		 	if(signChanges){
				biclustering = new CCC_Biclustering_SignChanges_Jumping_MissingValues(dm);
			} else {
				biclustering = new CCC_Biclustering_Jumping_MissingValues(dm);
			}
		} else {			
			if(signChanges){
				biclustering = new CCC_Biclustering_SignChanges(dm);
			} else {
				biclustering = new CCC_Biclustering(dm);
			}		
		}
		
		biclustering.computeBiclusters();
		numberOfOriginalBiclusters = biclustering.getBiclusters().size();		
		logger.println("OK");

		// list biclusters
		logger.println("found " + numberOfOriginalBiclusters + " biclusters");
		
		
		if(biclustering.getBiclusters().size()>0) {
									
			// sorting results
			switch(sortingCriteria){
			case BY_NUMBER_OF_GENES:
				logger.print("Sorting biclusters by number of genes...");
				if(sortingOrder==SORTING_ORDER.ASCENDING){
					PostProcessedBiclusteringInDiscretizedMatrix.sort_GENES_ASC(biclustering);
				} else {
					PostProcessedBiclusteringInDiscretizedMatrix.sort_GENES_DESC(biclustering);
				}				
				break;
			case BY_NUMBER_OF_CONDITIONS:
				logger.print("Sorting biclusters by number of conditions...");
				if(sortingOrder==SORTING_ORDER.ASCENDING){
					PostProcessedBiclusteringInDiscretizedMatrix.sort_CONDITIONS_ASC(biclustering);
				} else {
					PostProcessedBiclusteringInDiscretizedMatrix.sort_CONDITIONS_DESC(biclustering);
				}
				break;
			case BY_SIZE:
				logger.print("Sorting biclusters by bicluster size (genes X conditions)...");
				if(sortingOrder==SORTING_ORDER.ASCENDING){
					PostProcessedBiclusteringInDiscretizedMatrix.sort_SIZE_ASC(biclustering);
				} else {
					PostProcessedBiclusteringInDiscretizedMatrix.sort_SIZE_DESC(biclustering);
				}
				break;
			case BY_MSR:
				logger.print("Sorting biclusters by mean squared residue score (MSR)...");
				PostProcessedBiclusteringInDiscretizedMatrix.sort_MSR_ASC(biclustering);	
				break;
			case BY_PVALUE:
				logger.print("Sorting biclusters by pattern pvalue...");
				PostProcessedBiclusteringInDiscretizedMatrix.sort_BiclusterPatternWithColumns_pValue_ASC(biclustering,1);			
				break;
			}

			logger.println("OK");
		
			
			// filtering results		
			if(filterByMinNumberOfGenes) {
				logger.println("Filtering by min number of genes (" + minNumberOfGenes + ")...");				
				if(m.getNumberOfGenes()<minNumberOfGenes) throw new Exception("the number of input genes is lower than the min number of genes");
				PostProcessedBiclusteringInDiscretizedMatrix.filter_GENES(biclustering, minNumberOfGenes);			
			}			
			if(filterByMinNumberOfConditions) {
				logger.println("Filtering by min number of conditions (" + minNumberOfConditions + ")...");
				if(m.getNumberOfConditions()<minNumberOfConditions) throw new Exception("the number of input conditions is lower than the min number of conditions");
				PostProcessedBiclusteringInDiscretizedMatrix.filter_CONDITIONS(biclustering, minNumberOfConditions);
			}
			if(filterByconstantPattern) {
				logger.println("Filtering constant patterns...");
				PostProcessedBiclusteringInDiscretizedMatrix.filter_CONSTANT_EXPRESSION_PATTERN(biclustering);		
			}
			if(filterByMaxMSR) {
				logger.println("Filtering by max MSR (" + maxMSR + ")...");
				PostProcessedBiclusteringInDiscretizedMatrix.filter_MSR(biclustering, maxMSR);
			}
			if(filterByMaxPValue) {
				logger.println("Filtering by max pvalue (" + maxPValue + ")...");
				PostProcessedBiclusteringInDiscretizedMatrix.filter_BiclusterPatternWithColumns_pValue(biclustering, maxPValue,1);
			}
			if(filterByOverlapping) {
				logger.println("Filtering by overlapping (" + maxOverlappingPercentage + "%)...");
				PostProcessedBiclusteringInDiscretizedMatrix.filter_OVERLAPPING_ROWS_COLUMNS(biclustering, (float)(maxOverlappingPercentage/100.0));
			}

			
			logger.println("OK");

		} 
			
		biclusters = biclustering.getBiclusters();
		
		numberOfBiclusters = biclusters.size();
		numberOfFilteredBiclusters = numberOfOriginalBiclusters - numberOfBiclusters;
		
		logger.println("found " + biclusters.size() + " after filtering (" + numberOfFilteredBiclusters + " deleted)");
				
		
		
	}

	
	
	private void loadData() throws IOException{
		
		// read file
		ArrayList<String> lines = (ArrayList<String>) IOUtils.readLines(matrixFile);
		
		// condition names
		String[] first  = lines.get(0).split("\t");
		conditionNames = Arrays.copyOfRange(first,1,first.length);
		
		// init gene names
		geneNames = new String[lines.size()-1];
		
		// init expression matrix
		expressionMatrix = new float[geneNames.length][conditionNames.length];
		
		// init working variables
		String[] row;
		String line;
		String separator = "\t";
		
		// fill data
		for(int i=0 ; i<geneNames.length; i++){
			
			// read line and avoid missing values		
			line = lines.get(i+1) + "<";			
			line = line.replaceAll(separator + "<",separator + "@");
//			line = line.replaceAll(separator + separator + separator, separator + "@" + separator + "@" + separator);
//			line = line.replaceAll(separator + separator, separator + "@" + separator);	
//			line = line.replaceAll("@" + separator + separator + "@", "@" + separator + "@" + separator + "@");
			String last;
			do {
				last = line;
				line = line.replaceAll(separator + separator, separator + "@" + separator);
			} while(!last.equals(line));
			
			line = line.replaceAll("<", "");
			row = line.split("\t");
			
			// gene name
			geneNames[i] = row[0];
					
//			System.out.println("el line: " + line.replaceAll("\t", "_tab_"));
			
			// row values
			for(int j=0; j<conditionNames.length; j++){
				
//				System.out.println(row[j+1]);
				
				if(row[j+1].equals("@")){
					expressionMatrix[i][j] = DEFAULT_MISSING_VALUE;
				} else {					
//					if (!row[j+1].isEmpty())
//					{
						expressionMatrix[i][j] = Float.parseFloat(row[j+1]);
//					}
				}
			}			
			
		}		
		
	}
	
//	private void loadData() throws IOException{		
//		Dataset dataset = new Dataset(new File(matrixFile));
//		dataset.load();
//		expressionMatrix = dataset.getFloatMatrix();
//		geneNames = (String[])dataset.getFeatureNames().toArray();
//		conditionNames = (String[])dataset.getSampleNames().toArray();
//	}
	
	/**
	 * @return the geneNames
	 */
	public String[] getGeneNames() {
		return geneNames;
	}

	/**
	 * @param geneNames the geneNames to set
	 */
	public void setGeneNames(String[] geneNames) {
		this.geneNames = geneNames;
	}

	/**
	 * @return the conditionNames
	 */
	public String[] getConditionNames() {
		return conditionNames;
	}

	/**
	 * @param conditionNames the conditionNames to set
	 */
	public void setConditionNames(String[] conditionNames) {
		this.conditionNames = conditionNames;
	}

	/**
	 * @return the expressionMatrix
	 */
	public float[][] getExpressionMatrix() {
		return expressionMatrix;
	}

	/**
	 * @param expressionMatrix the expressionMatrix to set
	 */
	public void setExpressionMatrix(float[][] expressionMatrix) {
		this.expressionMatrix = expressionMatrix;
	}

	/**
	 * @return the biclusters
	 */
	public ArrayList<Bicluster> getBiclusters() {
		return biclusters;
	}

	/**
	 * @param biclusters the biclusters to set
	 */
	public void setBiclusters(ArrayList<Bicluster> biclusters) {
		this.biclusters = biclusters;
	}

	/**
	 * @return the logger
	 */
	public Logger getLogger() {
		return logger;
	}

	/**
	 * @param logger the logger to set
	 */
	public void setLogger(Logger logger) {
		this.logger = logger;
	}

	/**
	 * @return the sortingCriteria
	 */
	public SORTING_CRITERIA getSortingCriteria() {
		return sortingCriteria;
	}

	/**
	 * @param sortingCriteria the sortingCriteria to set
	 */
	public void setSortingCriteria(SORTING_CRITERIA sortingCriteria) {
		this.sortingCriteria = sortingCriteria;
	}

	/**
	 * @return the sortingOrder
	 */
	public SORTING_ORDER getSortingOrder() {
		return sortingOrder;
	}

	/**
	 * @param sortingOrder the sortingOrder to set
	 */
	public void setSortingOrder(SORTING_ORDER sortingOrder) {
		this.sortingOrder = sortingOrder;
	}

	/**
	 * @return the filterByconstantPattern
	 */
	public boolean isFilterByconstantPattern() {
		return filterByconstantPattern;
	}

	/**
	 * @param filterByconstantPattern the filterByconstantPattern to set
	 */
	public void setFilterByconstantPattern(boolean filterByconstantPattern) {
		this.filterByconstantPattern = filterByconstantPattern;
	}

	/**
	 * @return the filterByMinNumberOfGenes
	 */
	public boolean isFilterByMinNumberOfGenes() {
		return filterByMinNumberOfGenes;
	}

	/**
	 * @param filterByMinNumberOfGenes the filterByMinNumberOfGenes to set
	 */
	public void setFilterByMinNumberOfGenes(boolean filterByMinNumberOfGenes) {
		this.filterByMinNumberOfGenes = filterByMinNumberOfGenes;
	}

	/**
	 * @return the filterByMinNumberOfConditions
	 */
	public boolean isFilterByMinNumberOfConditions() {
		return filterByMinNumberOfConditions;
	}

	/**
	 * @param filterByMinNumberOfConditions the filterByMinNumberOfConditions to set
	 */
	public void setFilterByMinNumberOfConditions(
			boolean filterByMinNumberOfConditions) {
		this.filterByMinNumberOfConditions = filterByMinNumberOfConditions;
	}

	/**
	 * @return the filterByMaxMSR
	 */
	public boolean isFilterByMaxMSR() {
		return filterByMaxMSR;
	}

	/**
	 * @param filterByMaxMSR the filterByMaxMSR to set
	 */
	public void setFilterByMaxMSR(boolean filterByMaxMSR) {
		this.filterByMaxMSR = filterByMaxMSR;
	}

	/**
	 * @return the minNumberOfGenes
	 */
	public int getMinNumberOfGenes() {
		return minNumberOfGenes;
	}

	/**
	 * @param minNumberOfGenes the minNumberOfGenes to set
	 */
	public void setMinNumberOfGenes(int minNumberOfGenes) {
		this.minNumberOfGenes = minNumberOfGenes;
	}

	/**
	 * @return the minNumberOfConditions
	 */
	public int getMinNumberOfConditions() {
		return minNumberOfConditions;
	}

	/**
	 * @param minNumberOfConditions the minNumberOfConditions to set
	 */
	public void setMinNumberOfConditions(int minNumberOfConditions) {
		this.minNumberOfConditions = minNumberOfConditions;
	}

	/**
	 * @return the maxMSR
	 */
	public float getMaxMSR() {
		return maxMSR;
	}

	/**
	 * @param maxMSR the maxMSR to set
	 */
	public void setMaxMSR(float maxMSR) {
		this.maxMSR = maxMSR;
	}

	/**
	 * @return the matrixFile
	 */
	public String getMatrixFile() {
		return matrixFile;
	}

	/**
	 * @param matrixFile the matrixFile to set
	 */
	public void setMatrixFile(String matrixFile) {
		this.matrixFile = matrixFile;
	}

	/**
	 * @return the signChanges
	 */
	public boolean isSignChanges() {
		return signChanges;
	}

	/**
	 * @param signChanges the signChanges to set
	 */
	public void setSignChanges(boolean signChanges) {
		this.signChanges = signChanges;
	}

	/**
	 * @return the maxPValue
	 */
	public float getMaxPValue() {
		return maxPValue;
	}

	/**
	 * @param maxPValue the maxPValue to set
	 */
	public void setMaxPValue(float maxPValue) {
		this.maxPValue = maxPValue;
	}

	/**
	 * @return the filterByMaxPValue
	 */
	public boolean isFilterByMaxPValue() {
		return filterByMaxPValue;
	}

	/**
	 * @param filterByMaxPValue the filterByMaxPValue to set
	 */
	public void setFilterByMaxPValue(boolean filterByMaxPValue) {
		this.filterByMaxPValue = filterByMaxPValue;
	}

	/**
	 * @return the filterByOverlapping
	 */
	public boolean isFilterByOverlapping() {
		return filterByOverlapping;
	}

	/**
	 * @param filterByOverlapping the filterByOverlapping to set
	 */
	public void setFilterByOverlapping(boolean filterByOverlapping) {
		this.filterByOverlapping = filterByOverlapping;
	}

	/**
	 * @return the maxOverlappingPercentage
	 */
	public float getMaxOverlappingPercentage() {
		return maxOverlappingPercentage;
	}

	/**
	 * @param maxOverlappingPercentage the maxOverlappingPercentage to set
	 */
	public void setMaxOverlappingPercentage(float maxOverlappingPercentage) {
		this.maxOverlappingPercentage = maxOverlappingPercentage;
	}

	/**
	 * @return the numberOfOriginalBiclusters
	 */
	public int getNumberOfOriginalBiclusters() {
		return numberOfOriginalBiclusters;
	}

	/**
	 * @param numberOfOriginalBiclusters the numberOfOriginalBiclusters to set
	 */
	public void setNumberOfOriginalBiclusters(int numberOfOriginalBiclusters) {
		this.numberOfOriginalBiclusters = numberOfOriginalBiclusters;
	}

	/**
	 * @return the missingValuesScheme
	 */
	public MISSING_VALUES_SCHEME getMissingValuesScheme() {
		return missingValuesScheme;
	}

	/**
	 * @param missingValuesScheme the missingValuesScheme to set
	 */
	public void setMissingValuesScheme(MISSING_VALUES_SCHEME missingValuesScheme) {
		this.missingValuesScheme = missingValuesScheme;
	}

	/**
	 * @return the missingValuesFillingScheme
	 */
	public MISSING_VALUES_FILLING getMissingValuesFillingScheme() {
		return missingValuesFillingScheme;
	}

	/**
	 * @param missingValuesFillingScheme the missingValuesFillingScheme to set
	 */
	public void setMissingValuesFillingScheme(
			MISSING_VALUES_FILLING missingValuesFillingScheme) {
		this.missingValuesFillingScheme = missingValuesFillingScheme;
	}

	/**
	 * @return the missingValuesNeighbors
	 */
	public int getMissingValuesNeighbors() {
		return missingValuesNeighbors;
	}

	/**
	 * @param missingValuesNeighbors the missingValuesNeighbors to set
	 */
	public void setMissingValuesNeighbors(int missingValuesNeighbors) {
		this.missingValuesNeighbors = missingValuesNeighbors;
	}

	/**
	 * @return the numberOfBiclusters
	 */
	public int getNumberOfBiclusters() {
		return numberOfBiclusters;
	}

	/**
	 * @param numberOfBiclusters the numberOfBiclusters to set
	 */
	public void setNumberOfBiclusters(int numberOfBiclusters) {
		this.numberOfBiclusters = numberOfBiclusters;
	}

	/**
	 * @return the numberOfFilteredBiclusters
	 */
	public int getNumberOfFilteredBiclusters() {
		return numberOfFilteredBiclusters;
	}

	/**
	 * @param numberOfFilteredBiclusters the numberOfFilteredBiclusters to set
	 */
	public void setNumberOfFilteredBiclusters(int numberOfFilteredBiclusters) {
		this.numberOfFilteredBiclusters = numberOfFilteredBiclusters;
	}


}


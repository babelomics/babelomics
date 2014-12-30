package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.XYLineChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.mlpr.classifier.GenericClassifier;
import org.bioinfo.mlpr.classifier.Knn;
import org.bioinfo.mlpr.classifier.RForest;
import org.bioinfo.mlpr.classifier.Svm;
import org.bioinfo.mlpr.classifier.result.ClassificationResult;
import org.bioinfo.mlpr.evaluation.KFoldCrossValidation;
import org.bioinfo.mlpr.selection.AbstractFeatureSelector;
import org.bioinfo.mlpr.selection.CfsFeatureSelector;
import org.bioinfo.mlpr.selection.ListFeatureSelector;
import org.bioinfo.mlpr.selection.PcaFeatureSelector;
import org.bioinfo.mlpr.utils.InstancesBuilder;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

import weka.core.Attribute;
import weka.core.Instances;

public class Predictor extends BabelomicsTool {

	private double progressStep;
	private double progressCurrent = 0;
	private StringBuilder combinedTable;
	private int numberOfBestSelectedClassifications = 5;

	private List<GenericClassifier> selectedClassifiers;
	private List<GenericClassifier> trainedClassifiers;

	private List<String> sampleNames;
	private List<List<Double>> correctSampleRatio;
	private List<String> bestClassiferParamList;
	private Instances instances;
	private Instances testInstances;

	public Predictor() {		
		initOptions();
	}


	@Override
	public void initOptions() {
		// data
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("dataset-arff", "dataset is in arff format",false,false));

		// data test
		options.addOption(OptionFactory.createOption("test", "the data test",false));
		options.addOption(OptionFactory.createOption("test-arff", "dataset test is in arff format",false,false));

		// class attribute
		options.addOption(OptionFactory.createOption("class", "corresponding class attribute in the dataset",false,true));
		options.addOption(OptionFactory.createOption("class-file", "class variable file",false,true));

		// sample and feature filters
		options.addOption(OptionFactory.createOption("sample-filter", "Sample filter", false));
		options.addOption(OptionFactory.createOption("feature-filter", "Feature filter", false));

		// classifiers algorithms
		OptionGroup classifiers = new OptionGroup();
		classifiers.setRequired(false);

		// Validation
		options.addOption(OptionFactory.createOption("loo", "Perform leaving-one-out cross validation analysis",false,false));
		options.addOption(OptionFactory.createOption("kfold", "Perform kfold cross validation analysis",false,false));
		options.addOption(OptionFactory.createOption("folds", "Number of folds in cross validation evaluation",false,true));
		options.addOption(OptionFactory.createOption("repeats", "Number of repeat each randomization",false,true));

		options.addOption(OptionFactory.createOption("feature-selection", "Feature selection",false,true));			

		// KNN
		classifiers.addOption(OptionFactory.createOption("knn", "Classify dataset with a KNN classifier",false,false));
		options.addOption(OptionFactory.createOption("knn-tune", "Perform automated number of neighbors tunning",false,false));
		options.addOption(OptionFactory.createOption("knn-neighbors", "Knn number of neighbors (5 by default)", false, true));

		// SVM
		classifiers.addOption(OptionFactory.createOption("svm", "Classify dataset with a SVM classifier", false,false));
		options.addOption(OptionFactory.createOption("svm-tune", "Perform automated parameter tunning",false,false));
		options.addOption(OptionFactory.createOption("svm-cost", "----", false));

		// Random forest
		classifiers.addOption(OptionFactory.createOption("random-forest", "Classify dataset with a Random Forest classifier", false,false));
		options.addOption(OptionFactory.createOption("random-forest-tune", "Perform automated number of tree tunning",false,false));
		options.addOption(OptionFactory.createOption("random-forest-trees", "Number of trees", false));

		//			// DLDA
		//			classifiers.addOption(OptionFactory.createOption("dlda", "Classify dataset with a DLDA classifier", false,false));
		//			options.addOption(OptionFactory.createOption("dlda-tune", "Perform automated parameter tunning",false,false));
		//			
		//			// SOM
		//			classifiers.addOption(OptionFactory.createOption("som", "Classify dataset with a SOM classifier", false,false));
		//			options.addOption(OptionFactory.createOption("som-tune", "Perform automated parameter tunning",false,false));
		//			
		//			// PAM
		//			classifiers.addOption(OptionFactory.createOption("pam", "Classify dataset with a PAM classifier", false,false));
		//			options.addOption(OptionFactory.createOption("pam-tune", "Perform automated parameter tunning",false,false));

		options.addOptionGroup(classifiers);

		selectedClassifiers = new ArrayList<GenericClassifier>();
		trainedClassifiers = new ArrayList<GenericClassifier>();

		// feature selection (gene selection), and other options
		//		options.addOption(OptionFactory.createOption("gene-selection", "the gene selection, valid values: f-ratio, wilcoxon", false));
		//options.addOption(OptionFactory.createOption("trainning-size", "number of genes to use in trainning separated by commas, default:2,5,10,20,35,50", false));		

	}

	@Override
	public void execute() {

		try {

			logger.info("Welcome to prophet...");

			// init status
			combinedTable = new StringBuilder();
			initStatus();

			// load instances
			File datasetFile = new File(commandLine.getOptionValue("dataset"));
			FileUtils.checkFile(datasetFile);
			loadInstances(datasetFile);

			// init sample names
			sampleNames = new ArrayList<String>(instances.numInstances());
			for(int i=0; i<instances.numInstances(); i++){
				sampleNames.add(instances.attribute("sample_name").value(i));
			}
			correctSampleRatio = new ArrayList<List<Double>>(30);
			bestClassiferParamList = new ArrayList<String>(30);

			// classifiers
			if(commandLine.hasOption("knn")) executeKnn(instances,testInstances);
			if(commandLine.hasOption("svm")) executeSvm(instances,testInstances);			
			if(commandLine.hasOption("random-forest")) executeRandomForest(instances,testInstances);		
			//if(commandLine.hasOption("dlda")) executeDlda(instances);
			//if(commandLine.hasOption("naive-bayes")) executeNaiveBayes(instances);s

			// save best classifiers table
			String header = GenericClassifier.getResultsTableHeader();
			String featureSelectionTag = "";
			if(commandLine.hasOption("feature-selection") && !commandLine.getOptionValue("feature-selection").equalsIgnoreCase("none")) {
				header = GenericClassifier.getFeatureSelectionResultsTableHeader();
				featureSelectionTag = "_FEATURE_SELECTION";
			}
			IOUtils.write(new File(outdir + "/best_classifiers_table.txt"), header + "\n" +  combinedTable.toString());	
			result.getOutputItems().add(0, new Item("combined_table", "best_classifiers_table.txt", " Combined results (best " + numberOfBestSelectedClassifications + " per classifier)", TYPE.FILE,Arrays.asList("TABLE","PREDICTOR" + featureSelectionTag + "_TABLE"),new HashMap<String,String>(),"Train.Summary",""));


			// test			
			if(testInstances!=null && selectedClassifiers.size()>0){

				// init selected names
				List<String> selectedNames = new ArrayList<String>(selectedClassifiers.size());
				for(int i=0; i<selectedClassifiers.size(); i++){
					selectedNames.add(selectedClassifiers.get(i).getClassifierName() + " " + selectedClassifiers.get(i).getParams());
				}

				// init test result table				
				List<List<String>> testResultTable = new ArrayList<List<String>>(testInstances.numInstances());
				for(int i=0; i<testInstances.numInstances(); i++){
					List<String> classifications = new ArrayList<String>(selectedClassifiers.size());
					for(int j=0; j<selectedClassifiers.size(); j++) {
						classifications.add("-");
					}					
					testResultTable.add(classifications);
				}

				// fill test result table 
				List<String> sampleNames = new ArrayList<String>(testInstances.numInstances());
				testInstances.setClassIndex(testInstances.numAttributes()-1);
				
				String predictedClass;
				for(int i=0; i<selectedClassifiers.size(); i++){
					
					GenericClassifier testClassifier = selectedClassifiers.get(i);
										
					//System.err.println("getting  feature selector: " + testClassifier.getFeatureSelector().getSelectedAttributeNames());
					
					// execute (train and) test
					List<ClassificationResult> testResult = testClassifier.test(instances, testInstances);
					
					// sample name column
					if(i==0){
						for(int j=0; j<testInstances.numInstances(); j++){
							sampleNames.add(testResult.get(j).getInstanceName());
						}
					}
					
					// values
					for(int j=0; j<testInstances.numInstances(); j++){
						predictedClass = testResult.get(j).getClassName();
						testResultTable.get(j).set(i,predictedClass);
					}
				}

				// save to disk
				StringBuilder testResultString = new StringBuilder();
				testResultString.append("#Sample_names").append("\t").append(ListUtils.toString(selectedNames, "\t")).append("\n");
				for(int i=0; i<testInstances.numInstances(); i++){
					testResultString.append(sampleNames.get(i)).append("\t").append(ListUtils.toString(testResultTable.get(i))).append("\n");
				}
				IOUtils.write(outdir + "/test_result.txt", testResultString.toString());
				result.addOutputItem(new Item("test_result", "test_result.txt", "Test result table", TYPE.FILE, Arrays.asList("TABLE"), new HashMap<String, String>(),"Test.Test result"));
			}

			//
			savecorrectSampleRatioTable();

		} catch(Exception e){
			e.printStackTrace();
		}

	}


	private void executeKnn(Instances instances, Instances test) throws Exception {

				
		// update status
		updateStatus(progressCurrent,"executing KNN classifier");

		// init classifier
		Knn knn = new Knn();

		// params
		if(commandLine.hasOption("knn-tune")) {
			knn.setTuneParameters(true);
		} else {
			if(commandLine.hasOption("knn-neighbors")) knn.setKnn(Integer.parseInt(commandLine.getOptionValue("knn-neighbors")));	
		}

		// validation
		setFeatureSelector(knn);
		setCrossValidation(knn);

		// train
		knn.train(instances);

		// select the best classifiers
		int[] best = knn.getBestClassifiers(numberOfBestSelectedClassifications);
		bestClassiferParamList.addAll(ListUtils.subList(knn.getParamList(),best));

		// select best classifiers for testing
		if(test!=null){			
			trainedClassifiers.add(knn);
			for(int i=0; i<best.length; i++){				
				Knn testKnn;
				if(knn.getFeatureSelector()!=null){
					System.err.println("preparing knn:" + knn.getKnnValues().get(best[i]) + " features: " + knn.getSelectedFeatureList().get(best[i]));
					ListFeatureSelector featureSelector = new ListFeatureSelector(knn.getSelectedFeatureList().get(best[i]));
					testKnn = new Knn(knn.getKnnValues().get(best[i]),featureSelector);
				} else {
					testKnn = new Knn(knn.getKnnValues().get(best[i]));
				}
				selectedClassifiers.add(testKnn);
			}			
		}

		// acum correct classification sample ratios
		double[][] ratios = knn.getEvaluationResultList().getCorrectClassificationRatio(sampleNames,best);
		addCorrectClassificationRatios(ratios);

//		Canvas canvas = knn.getEvaluationResultList().generateHeatmap(ratios,sampleNames,bestClassiferParamList);
//		canvas.save(outdir + "/ratios.png");
//		result.addOutputItem(new Item("ratios", "ratios.png", "Sample correct classification ratio (100% blue, 0% red)", TYPE.IMAGE, Arrays.asList(""), new HashMap<String, String>(),"Summary"));
		
		// save results
		saveClassifierResults(knn);


	}

	private void addCorrectClassificationRatios(double[][] ratios){
		for(int j=0; j<ratios[0].length; j++){
			List<Double> singleRatio = new ArrayList<Double>(ratios.length);
			for(int i=0; i<ratios.length; i++){			
				singleRatio.add(ratios[i][j]);
			}
			correctSampleRatio.add(singleRatio);
		}		
	}

	private void executeSvm(Instances instances, Instances test) throws Exception {

		// update status
		updateStatus(progressCurrent,"executing SVM classifier");

		// init classifier
		Svm svm = new Svm();

		// init params
		if(commandLine.hasOption("svm-tune")) svm.setTuneParameters(true);
		else {
			if(commandLine.hasOption("svm-cost")) svm.setCost(Integer.parseInt(commandLine.getOptionValue("svm-cost")));						
		}

		// validation
		setFeatureSelector(svm);
		setCrossValidation(svm);

		// train
		svm.train(instances);

		// select the best classifiers
		int[] best = svm.getBestClassifiers(numberOfBestSelectedClassifications);
		bestClassiferParamList.addAll(ListUtils.subList(svm.getParamList(),best));
		
		// select best classifiers for testing
		if(test!=null){			
			trainedClassifiers.add(svm);
			for(int i=0; i<best.length; i++){				
				Svm testSvm;
				if(svm.getFeatureSelector()!=null){
					System.err.println("preparing svm:" + svm.getCostValues().get(best[i]) + " features: " + svm.getSelectedFeatureList().get(best[i]));
					ListFeatureSelector featureSelector = new ListFeatureSelector(svm.getSelectedFeatureList().get(best[i]));
					testSvm = new Svm(svm.getCostValues().get(best[i]),featureSelector);
				} else {
					testSvm = new Svm(svm.getCostValues().get(best[i]));
				}
				selectedClassifiers.add(testSvm);
			}			
		}
		
		// acum correct classification sample ratios
		double[][] ratios = svm.getEvaluationResultList().getCorrectClassificationRatio(sampleNames,best);
		addCorrectClassificationRatios(ratios);
		
		// save results
		saveClassifierResults(svm);
		
	}

	private void executeRandomForest(Instances instances, Instances test) throws Exception {

		// update status
		updateStatus(progressCurrent,"executing Random Forest classifier");

		// init classifier
		RForest randomForest = new RForest();

		// init params
		if(commandLine.hasOption("random-forest-tune")) randomForest.setTuneParameters(true);
		else {
			if(commandLine.hasOption("random-forest-trees")) randomForest.setNumTrees(Integer.parseInt(commandLine.getOptionValue("random-forest-trees")));						
		}

		// validation
		setFeatureSelector(randomForest);
		setCrossValidation(randomForest);

		// train
		randomForest.train(instances);

		// select the best classifiers
		int[] best = randomForest.getBestClassifiers(numberOfBestSelectedClassifications);
		bestClassiferParamList.addAll(ListUtils.subList(randomForest.getParamList(),best));
		
		// select best classifiers for testing
		if(test!=null){			
			trainedClassifiers.add(randomForest);
			for(int i=0; i<best.length; i++){				
				RForest testRandomForest;
				if(randomForest.getFeatureSelector()!=null){
					System.err.println("preparing random forest:" + randomForest.getNumTreesValues().get(best[i]) + " features: " + randomForest.getSelectedFeatureList().get(best[i]));
					ListFeatureSelector featureSelector = new ListFeatureSelector(randomForest.getSelectedFeatureList().get(best[i]));
					testRandomForest = new RForest(randomForest.getNumTreesValues().get(best[i]),featureSelector);
				} else {
					testRandomForest = new RForest(randomForest.getNumTreesValues().get(best[i]));
				}
				selectedClassifiers.add(testRandomForest);
				System.err.println("getting  feature selector: " + testRandomForest.getFeatureSelector().getSelectedAttributeNames());
			}			
		}
		
		// acum correct classification sample ratios
		double[][] ratios = randomForest.getEvaluationResultList().getCorrectClassificationRatio(sampleNames,best);
		addCorrectClassificationRatios(ratios);
		
		// save results
		saveClassifierResults(randomForest);

	}

	private String testResultToString(List<ClassificationResult> classificationResult){
		StringBuilder testResult = new StringBuilder();
		testResult.append("#Sample_name\tPredicted_class\n");
		for(ClassificationResult classification: classificationResult){
			testResult.append(classification.getInstanceName()).append("\t").append(classification.getClassName()).append("\n");
		}
		return testResult.toString();
	}

	private void saveClassifierResults(GenericClassifier classifier){

//		int best = classifier.getEvaluationResultList().getBestAreaUnderRocIndex();
		String name = classifier.getClassifierName();

//		// best case
//		try {
//			IOUtils.write(new File(outdir + "/" + name + ".txt"), "Best AUC classification\n\nParameters: " + classifier.getParamList().get(best) + "\n\n" + classifier.getEvaluationResultList().getBestAreaUnderRoc().toString());			
//			result.addOutputItem(new Item(name + "_result_file", name + ".txt", name + " result file", TYPE.FILE,name + " results",""));
//		} catch (IOException e) {
//			printError("ioexception_execute_" + name + "_predictor", "Error saving " + name + " results", "Error saving " + name + " results");
//		}

		// classification table
		String featureSelectionTag = "";
		if(commandLine.hasOption("feature-selection") && !commandLine.getOptionValue("feature-selection").equalsIgnoreCase("none")) featureSelectionTag = "_FEATURE_SELECTION";
		try {
			IOUtils.write(new File(outdir + "/" + name + "_table.txt"), classifier.getResultsTable(true));
			result.addOutputItem(new Item(name + "_table", name + "_table.txt", name + " classifications", TYPE.FILE,Arrays.asList("TABLE","PREDICTOR" + featureSelectionTag + "_TABLE"),new HashMap<String,String>(),"Train." + name + " results",""));
		} catch (IOException e) {
			printError("ioexception_" + name + "_predictor", "Error saving " + name + " classifications table", "Error saving " + name + " classifications table");
		}

		// Comparative plot
		try {
			XYLineChart comparativeplot = classifier.getComparativeMetricsPlot();
			comparativeplot.save(outdir + "/" + name + "_comparative_plot.png",500,350,"png");
			result.addOutputItem(new Item(name + "_comparative_plot", name + "_comparative_plot.png", name + " comparative plot", TYPE.IMAGE, Arrays.asList(""),new HashMap<String,String>(),"Train." + name + " results",""));
		} catch (IOException e) {
			printError("ioexception_" + name + "_predictor", "Error saving " + name + " comparative plot", "Error saving " + name + " comparative plot");
			e.printStackTrace();
		}

		combinedTable.append(classifier.getSortedResultsTable(numberOfBestSelectedClassifications));
	}


	private void setFeatureSelector(GenericClassifier classifier) throws Exception{
		
				
		// feature selection		
		if(commandLine.hasOption("feature-selection") && !commandLine.getOptionValue("feature-selection").equalsIgnoreCase("none")){
			
			AbstractFeatureSelector featureSelector = null;			
			
			String featureSelectionMethod = commandLine.getOptionValue("feature-selection");
			
			// CFS subset
			if("cfs".equalsIgnoreCase(featureSelectionMethod)){
				featureSelector = new CfsFeatureSelector();
				logger.println("Setting CFS feature selector");
			}
			
			// PCA 
			else if("pca".equalsIgnoreCase(featureSelectionMethod)){
				featureSelector = ((AbstractFeatureSelector)new PcaFeatureSelector());
				classifier.initNumberOfFeaturesRange(instances.numAttributes()-2);				
				logger.println("Setting PCA feature selector");
			}
			
			// Genetic algorithm
			else if("ga".equalsIgnoreCase(featureSelectionMethod)){
				featureSelector = new CfsFeatureSelector();
				logger.println("Setting Genetic Algorithm feature selector");
			} 
			
			// error
			else {				
				throw new Exception("ERROR: feature selection method " + featureSelectionMethod + " is undefined");
			}
			
			classifier.setFeatureSelector(featureSelector);
			
		}
		
	}
	
	private void setCrossValidation(GenericClassifier classifier) throws Exception{
		// validation
		if(commandLine.hasOption("loo")){
			classifier.setClassifierEvaluation(new KFoldCrossValidation(1, instances.numInstances()-1));
			logger.println("Setting Leave-one-out cross-validation");	
		} else {
			int repeats = Integer.parseInt(commandLine.getOptionValue("repeats", "10"));
			int folds = Integer.parseInt(commandLine.getOptionValue("folds", "5"));
			classifier.setClassifierEvaluation(new KFoldCrossValidation(repeats, folds));
			logger.println("Setting Kfold cross-validation (repeats=" + repeats + ",folds=" + folds + ")");
		}
	}

	private void executeDlda(Instances instances) {		
		printError("executeDlda_predictor", "DLDA is not implemented yet", "DLDA is not implemented yet");
	}

	private void executeNaiveBayes(Instances instances) {
		printError("executeDlda_predictor", "Naive Bayes is not implemented yet", "Naive Bayes is not implemented yet");
	}
	

	private void savecorrectSampleRatioTable() throws IOException{
		DecimalFormat percentageFormatter = new DecimalFormat("##.##");
		
		// correct sample classification ratio
		StringBuilder ratiosTable = new StringBuilder();
		StringBuilder ratiosHtmlTable = new StringBuilder();
		
		// print header
		
		ratiosTable.append("#Sample").append("\t");
		ratiosHtmlTable.append("#Sample").append("\t");
		
		for(int j=0; j<bestClassiferParamList.size(); j++){
			
			ratiosTable.append(bestClassiferParamList.get(j));
			ratiosHtmlTable.append(bestClassiferParamList.get(j));
			
			if(j<(bestClassiferParamList.size()-1)){
				
				ratiosTable.append("\t");
				ratiosHtmlTable.append("\t");
				
			}
		}
		
		ratiosTable.append("\n");
		ratiosHtmlTable.append("\n");
		
		//// print values
		for(int i=0; i<sampleNames.size(); i++){
			
			ratiosTable.append(sampleNames.get(i)).append("\t");
			ratiosHtmlTable.append(sampleNames.get(i)).append("\t");
			
			for(int j=0; j<correctSampleRatio.size(); j++){
				
				ratiosTable.append(percentageFormatter.format(correctSampleRatio.get(j).get(i)*100.0));
				ratiosHtmlTable.append(renderCell(correctSampleRatio.get(j).get(i)));
				
				if(j<(correctSampleRatio.size()-1)){
					
					ratiosTable.append("\t");
					ratiosHtmlTable.append("\t");
					
				}		
			}
			
			ratiosTable.append("\n");
			ratiosHtmlTable.append("\n");
			
		}
		IOUtils.write(outdir + "/ratios.txt", ratiosTable.toString());
		IOUtils.write(outdir + "/ratios.html", ratiosHtmlTable.toString());
		result.addOutputItem(new Item("ratios", "ratios.html", "Percentage of correct classification per sample/classifier", TYPE.FILE, Arrays.asList("TABLE", "NO_LINK"), new HashMap<String, String>(),"Train.Summary"));
	}

	private String renderCell(double percentage){
		DecimalFormat percentageFormatter = new DecimalFormat("##.##");
		double red = (1-percentage)*200.0;
		double green = (percentage)*200.0;
		String color = "color: rgb(" + (int)Math.floor(red) + "," + (int)Math.floor(green) + ",0)";		
		return "<font style='font-weight: bold; " + color + "'>" + percentageFormatter.format(percentage*100.0) + "%</font>";
	}
	
	private void loadInstances(File datasetFile) throws Exception{

		// update status
		updateStatus(progressCurrent,"Reading dataset");

		// data set loading 
		if(commandLine.hasOption("dataset-arff")){
			instances = InstancesBuilder.getInstancesFromArrfFile(datasetFile,"sample_name");
			instances.setClassIndex(instances.numAttributes()-1);
		} else {
			// convert Dataset to Instances format (data is trasposed!!)
			Dataset dataset = new Dataset(datasetFile, true);
			List<List<String>> data = new ArrayList<List<String>>(dataset.getColumnDimension());
			for(int i = 0 ; i<dataset.getColumnDimension() ; i++) {
				data.add(ArrayUtils.toStringList(dataset.getDoubleMatrix().getColumn(i)));
			}
			// class values
			List<String> classValues = null;
			if(commandLine.hasOption("class-file")){
				String classFileName = commandLine.getOptionValue("class-file");
				FileUtils.checkFile(classFileName);
				classValues = StringUtils.toList(IOUtils.toString(classFileName));
			} else {
				if(commandLine.hasOption("class")){
					String className = commandLine.getOptionValue("class");
					if(dataset.getVariables()!=null && dataset.getVariables().getByName(className)!=null && dataset.getVariables().getByName(className).getValues()!=null){
						classValues = dataset.getVariables().getByName(className).getValues();
					} else {
						throw new Exception("class not found in dataset");
					}
				} else {
					throw new Exception("class or class-file is not defined"); 
				}
			}
			System.err.println("class values: " + classValues);
			// create instances
			instances = InstancesBuilder.getInstancesFromDataList(data, dataset.getFeatureNames(), dataset.getSampleNames(), classValues);
			// define class attribute
			Attribute classAttr = instances.attribute("class");
			instances.setClassIndex(classAttr.index());

			if(commandLine.hasOption("test") && !commandLine.getOptionValue("test").equalsIgnoreCase("none")){
				File testFile = new File(commandLine.getOptionValue("test"));	
				// data set loading 
				if(commandLine.hasOption("test-arff")){
					instances = InstancesBuilder.getInstancesFromArrfFile(datasetFile,"sample_name");				
				} else {
					Dataset datasetTest = new Dataset(testFile, true);
					List<List<String>> dataTest = new ArrayList<List<String>>(datasetTest.getColumnDimension());
					for(int i = 0 ; i<datasetTest.getColumnDimension() ; i++) {
						dataTest.add(ArrayUtils.toStringList(datasetTest.getDoubleMatrix().getColumn(i)));
					}
					// create instances
					testInstances = InstancesBuilder.getTestInstancesFromDataList(dataTest, datasetTest.getFeatureNames(), datasetTest.getSampleNames(),classValues);
					//				Attribute classAttr = test.attribute("class");
					//				test.setClassIndex(classAttr.index());
				}
			}
		}
	}


	/*
	 * 
	 * STATUS MANAGEMENT (TO REMOVE!!!!!!!!)
	 * 
	 */

	private void initStatus() {
		int steps = 2; // it includes reading dataset and done

		if(commandLine.hasOption("svm")) steps++;
		if(commandLine.hasOption("knn")) steps++;
		if(commandLine.hasOption("random-forest")) steps++;
		if(commandLine.hasOption("dlda")) steps++;
		if(commandLine.hasOption("som")) steps++;
		if(commandLine.hasOption("pam")) steps++;

		progressStep = 100.0 / steps;
		progressCurrent = progressStep;
	}

	private void updateStatus(double progressCurrent, String message){
		try {
			jobStatus.addStatusMessage(StringUtils.decimalFormat(progressCurrent),message);
			progressCurrent += progressStep;
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_executeknn_predictor", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
	}

}

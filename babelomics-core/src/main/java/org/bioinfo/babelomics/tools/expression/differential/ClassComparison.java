package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.methods.expression.differential.Limma;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.dataset.Variables;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.result.AnovaTestResult;
import org.bioinfo.math.result.LimmaTestResult;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.AnovaTest;
import org.bioinfo.math.stats.inference.FoldChangeTest;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class ClassComparison extends BabelomicsTool {

	private Dataset dataset;
	private String test;
	private String className = null;
	private List<String> classValues = null;
	private String correction;
	private double pValue = 0.05;
	private double foldChangeValue = 2;

	private int minDisplay = 10;
	private int maxDisplay = Integer.MAX_VALUE;
	private String msg = "";

	public ClassComparison() {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("class-name", "class"));
		options.addOption(OptionFactory.createOption("class-values", "value"));
		options.addOption(OptionFactory.createOption("test", "test: t, limma, anova, fold_change"));
		options.addOption(OptionFactory.createOption("correction", "Multiple-test correction: fdr, bh, by, bonferroni, hochberg, hold", false));
		options.addOption(OptionFactory.createOption("p-value", "p-value for significative genes", false));
		options.addOption(OptionFactory.createOption("fold-change-value", "fold-change for significative genes", false));
		//options.addOption(OptionFactory.createOption("batch", "class class variable"));
	}

	@Override
	public void execute() {

		List<String> values = null;

		// init
		//
		test = commandLine.getOptionValue("test", null);
		className = commandLine.getOptionValue("class-name", null);
		values = (commandLine.hasOption("class-values") ? StringUtils.toList(commandLine.getOptionValue("class-values", null), ",") : null);
		correction = commandLine.getOptionValue("correction", "fdr");
		String pValueParam = commandLine.getOptionValue("p-value", "0.05");
		String foldChangeValueParam = commandLine.getOptionValue("fold-change-value", "2");
		String datasetParam = commandLine.getOptionValue("dataset");

		if (className != null) {
			if (values == null) {
				values = ListUtils.unique(dataset.getVariables().getByName(className).getValues());
			} else {
				values = ListUtils.unique(values);
			}

			classValues = new ArrayList<String>();
			for(String val: values) {
				if ( val != null && val.trim().length() > 0 ) {
					classValues.add(val.trim());
				}
			}
		}

		if ("fold-change".equalsIgnoreCase(test) || "fold_change".equalsIgnoreCase(test)) {
			try {
				foldChangeValue = Double.parseDouble(foldChangeValueParam);
			} catch (NumberFormatException e) {
				foldChangeValue = 0.05;
			}
		} else {
			try {
				pValue = Double.parseDouble(pValueParam);
				if (pValue > 1 || pValue < 0) {
					pValue = 0.05;
				}
			} catch (NumberFormatException e) {
				pValue = 0.05;
			}
		}

		// input parameters
		//
		result.addOutputItem(new Item("dataset_input_param", (datasetParam == null ? "" : new File(datasetParam).getName()), "Dataset file name", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("test_input_param", test, "Test", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("class_input_param", (className == null ? "" : className) + " [" + ListUtils.toString(classValues, ", ") + "]", "Class", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		if ("fold_change".equalsIgnoreCase(test) || "fold-change".equalsIgnoreCase(test)) {
			result.addOutputItem(new Item("fold_change_value_input_param", foldChangeValueParam, "Fold-change value", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		} else {
			result.addOutputItem(new Item("correction_input_param", correction, "Multiple-test correction", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("pvalue_input_param", pValueParam, "Adjusted p-value", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		}

		// check input parameters
		//
		if ( datasetParam == null ) {
			abort("missingdataset_execute_classcomparison", "Missing dataset", "Missing dataset", "Missing dataset");
		}

		if ( className == null ) {
			abort("classnamemissing_execute_classcomparison", "class name missing", "class name missing", "class name missing");
		}

		if ( classValues == null ) {
			abort("classvaluesmissing_execute_classcomparison", "class values missing", "class values missing", "class values missing");
		}

		if ( test == null ) {
			abort("testmissing_execute_classcomparison", "class comparison test missing", "class comparison test missing", "class comparison test missing");
		}

		// reading dataset
		//
		File datasetFile = new File(datasetParam);
		try {
			dataset = new Dataset(datasetFile);
			if (!dataset.load() && !dataset.validate()) {
				abort("exception_execute_classcomparison", "Error", "Error loading dataset " + datasetFile.getName() + ": " + dataset.getMessages().getErrorString(""), "");
			}
		} catch (Exception e) {
			abort("exception_execute_clustering", "Error", "Error reading dataset " + datasetFile.getName(), "");
		}

		// executing test
		//
		updateJobStatus("40", "computing " + test);
		if ( "limma".equalsIgnoreCase(test) ) {
			executeLimma();
		} else if ( "t".equalsIgnoreCase(test) ) {
			if ( classValues.size() == 2 ) {
				executeT();
			} else {
				abort("testmismatched_execute_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");
			}
		} else if ("fold_change".equalsIgnoreCase(test) || "fold-change".equalsIgnoreCase(test)) {
			if ( classValues.size() == 2 ) {
				executeFoldChange();
			} else {
				abort("testmismatched_execute_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");
			}
		} else if ( "anova".equalsIgnoreCase(test) ) {
			if ( classValues.size() > 2 ) {
				executeAnova();
			} else {
				abort("testmismatched_execute_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");
			}
		} else {
			abort("testunknown_execute_classcomparison", "unknown test " + test, "unknown test " + test, "unknown test " + test);
		}
	}

	/**
	 *
	 */
	public void executeT() {

		int[] cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(0));
		DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);

		cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(1));
		DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);

		try {
			Dataset subDataset = dataset.getSubDataset(className, classValues);

			// apply test and multiple test correction according
			//
			TTest tTest = new TTest();
			TestResultList<TTestResult> res = tTest.tTest(sample1, sample2);
			DiffExpressionUtils.multipleTestCorrection(res, correction);

			// create output file
			//
			int[] columnOrder = ListUtils.order(subDataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

			DataFrame dataFrame = new DataFrame(subDataset.getFeatureNames().size(), 0);
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
			dataFrame.setRowNames(ListUtils.ordered(subDataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);
			File file = new File(outdir + "/t.txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item("tfile", file.getName(), "T-test output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "T-test output files"));
			}

			// getting significative genes
			//
			DiffExpressionUtils.addSignificativeResults(subDataset, test, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues(), "p-value", res.getPValues(), null, null, null, null, className, columnOrder, pValue, maxDisplay, this);
			DiffExpressionUtils.createFatiScanRedirection(dataFrame, test, "statistic", result, outdir);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("e.getMessage() = " + e.getMessage());
			abort("exception_executet_classcomparison", "ERROR", e.getMessage(), "");
		}
	}

	/**
	 *
	 */
	public void executeFoldChange() {

		int[] cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(0));
		DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);

		cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(1));
		DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);

		try {

			Dataset subDataset = dataset.getSubDataset(className, classValues);

			FoldChangeTest foldChange = new FoldChangeTest();
			double[] logRes  = foldChange.logFoldChange(sample1, sample2);
			double[] diffRes = foldChange.diffFoldChange(sample1, sample2);

			updateJobStatus("80", "saving results");

			setFoldChangeResults(subDataset, "log", "Log", logRes, className);
			setFoldChangeResults(subDataset, "diff", "Diff", diffRes, className);

//			String test, testLabel;
//			DataFrame dataFrame;
//			FeatureData featureData;
//			File file;
//
//			int[] columnOrder;
//			int[] rowOrder;
//			Canvas heatmap;
//			String heatmapFilename;
//
//			// log fold change
//			//
//			test = "log";
//			testLabel = "Log";
//
//			columnOrder = ListUtils.order(subDataset.getVariables().getByName(className).getValues());
//			rowOrder = ListUtils.order(ArrayUtils.toList(logRes), true);
//
//			dataFrame = new DataFrame(subDataset.getFeatureNames().size(), 0);
//			dataFrame.addColumn(test, ListUtils.ordered(ArrayUtils.toStringList(logRes), rowOrder));
//			dataFrame.setRowNames(ListUtils.ordered(subDataset.getFeatureNames(), rowOrder));	
//			featureData = new FeatureData(dataFrame);
//			file = new File(outdir + "/" + test + "_foldchange.txt");
//			IOUtils.write(file, dataFrame.toString(true, true));
//			//featureData.save(file);
//			if ( file.exists() ) {
//				result.addOutputItem(new Item(test + "_foldchange", file.getName(), testLabel + " fold-change output file", TYPE.FILE, StringUtils.toList("TABLE," + test.toUpperCase() + "_FOLD_CHANGE_TABLE", ","), new HashMap<String, String>(2), testLabel + " fold-change.Output files"));											
//			}		
//
//			// generating heatmap
//			//
//			heatmap = DiffExpressionUtils.generateHeatmap(subDataset, className, columnOrder, rowOrder, test + " fold-change", logRes, null, null);			
//			heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
//			try {
//				heatmap.save(heatmapFilename);
//				if ( new File(heatmapFilename).exists() ) {
//					result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", testLabel + ". fold-change heatmap with all terms", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), testLabel + " fold-change.Heatmap image"));
//				}
//			} catch (IOException e) {
//				printError("ioexception_executet_classcomparison", "error generating heatmap", e.toString(), e);
//			}
//
//
//			test = "diff";
//			testLabel = "Diff";
//
//			columnOrder = ListUtils.order(subDataset.getVariables().getByName(className).getValues());
//			rowOrder = ListUtils.order(ArrayUtils.toList(logRes), true);
//
//			dataFrame = new DataFrame(subDataset.getFeatureNames().size(), 0);
//			dataFrame.addColumn(test, ListUtils.ordered(ArrayUtils.toStringList(diffRes), rowOrder));
//			dataFrame.setRowNames(ListUtils.ordered(subDataset.getFeatureNames(), rowOrder));	
//			featureData = new FeatureData(dataFrame);
//			file = new File(outdir + "/" + test + "_foldchange.txt");
//			IOUtils.write(file, dataFrame.toString(true, true));
//			//featureData.save(file);
//			if ( file.exists() ) {
//				result.addOutputItem(new Item(test + "_foldchange", file.getName(), testLabel + " fold-change output file", TYPE.FILE, StringUtils.toList("TABLE," + test.toUpperCase() + "_FOLD_CHANGE_TABLE", ","), new HashMap<String, String>(2), testLabel + " fold-change.Output files"));											
//			}		
//
//			// generating heatmap
//			//
//			heatmap = DiffExpressionUtils.generateHeatmap(subDataset, className, columnOrder, rowOrder, test + " fold-change", diffRes, null, null);			
//			heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
//			try {
//				heatmap.save(heatmapFilename);
//				if ( new File(heatmapFilename).exists() ) {
//					result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", testLabel + ". fold-change heatmap with all terms", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), testLabel + " fold-change.Heatmap image"));
//				}
//			} catch (IOException e) {
//				printError("ioexception_executet_classcomparison", "error generating heatmap", e.toString(), e);
//			}

		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_executefoldchange_classcomparison", "error running fold-change", "error running fold-change: " + e.getMessage(), "error running fold-change: " + e.getMessage());
		}
	}


	private void setFoldChangeResults(Dataset subDataset, String test, String testLabel, double[] res, String className) throws InvalidIndexException, IOException {


//		int[] columnOrder = ListUtils.order(subDataset.getVariables().getByName(className).getValues());
		int[] rowOrder = ListUtils.order(ArrayUtils.toList(res), true);

		DataFrame dataFrame = new DataFrame(subDataset.getFeatureNames().size(), 0);
		dataFrame.addColumn(test, ListUtils.ordered(ArrayUtils.toStringList(res), rowOrder));
		dataFrame.setRowNames(ListUtils.ordered(subDataset.getFeatureNames(), rowOrder));

        FeatureData featureData = new FeatureData(dataFrame);
		File file = new File(outdir + "/" + test + "_foldchange.txt");
		IOUtils.write(file, dataFrame.toString(true, true));

		/** Get significative terms, babelomics 5 **/
		List<String> featureNames = new ArrayList<String>();
		List<String> featureValues = new ArrayList<String>();

		List<String> featureNamesUp = new ArrayList<String>();
		List<String> featureNamesDown = new ArrayList<String>();


        for (String rowName : dataFrame.getRowNames()) {
            List<String> row = dataFrame.getRow(rowName);
            double stats = Double.parseDouble(row.get(0));
            if (stats >= foldChangeValue) {
                featureNames.add(rowName);
                featureValues.add(row.get(0));
                if (stats >= 0)
                    featureNamesUp.add(rowName);
                else
                    featureNamesDown.add(rowName);
            }
        }

		dataFrame = new DataFrame(featureNames.size(), 0);
		dataFrame.addColumn(test, featureValues);
		dataFrame.setRowNames(featureNames);


        List<Integer> sigRowIndexes = new ArrayList<Integer>();
        for(String feat: featureNames){
            int idx = 0;
            for(String featSub: subDataset.getFeatureNames()){
                if(feat.equalsIgnoreCase(featSub))
                    sigRowIndexes.add(idx);
                idx++;
            }
        }
        DoubleMatrix doubleMatrix = new DoubleMatrix(dataFrame.getRowDimension(), subDataset.getColumnDimension());
        for (int i = 0; i < sigRowIndexes.size(); i++) {
            doubleMatrix.setRow(i, subDataset.getDoubleMatrix().getRow(sigRowIndexes.get(i)));
        }
        File fileXX = new File(outdir + "/" + test + "_foldchange_significative_dataset.txt");
        Dataset sigDataset = new Dataset(subDataset.getSampleNames(), featureNames, doubleMatrix);
        sigDataset.save(fileXX);


		File fileAux = new File(outdir + "/" + test + "_foldchange_significative_table.txt");
		IOUtils.write(fileAux, dataFrame.toString(true, true));

		fileAux = new File(outdir + "/" + test + "_foldchange_significative_table_up.txt");
		IOUtils.write(fileAux, featureNamesUp);

		fileAux = new File(outdir + "/" + test + "_foldchange_significative_table_down.txt");
		IOUtils.write(fileAux, featureNamesDown);

//        fileAux = new File(outdir + "/" + test + "_foldchange_significative_dataset.txt");
//        IOUtils.write(fileAux, featureNamesDown);


//
//		List<String> featuresUp = new ArrayList<String>();
//		List<String> featuresDown = new ArrayList<String>();
//		for (String rowName : dataFrame.getRowNames()) {
//                List<String> row = dataFrame.getRow(rowName);
//                double stats = Double.parseDouble(row.get(0));
//                if(stats>=foldChangeValue){
//
//
////                    featuresUp.add(rowName);
////                else
////                    featuresDown.add(rowName);
//				}
//            }
//			file = new File(tool.getOutdir() + "/" + test + "_significative_table_up.txt");
//			IOUtils.write(file, featuresUp);
//			file = new File(tool.getOutdir() + "/" + test + "_significative_table_down.txt");
//			IOUtils.write(file, featuresDown);


		//featureData.save(file);
		if ( file.exists() ) {
			result.addOutputItem(new Item(test + "_foldchange", file.getName(), testLabel + " fold-change output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), testLabel + " fold-change.Output files"));

			String json = "{\\\"paramfilename\\\": \\\"input_params.txt\\\", \\\"testfilename\\\": \\\"" + file.getName() + "\\\"}";
			result.addOutputItem(new Item("diff_expr_" + StringUtils.randomString(8), json, "Significative results", TYPE.FILE, StringUtils.toList("DIFF_EXPRESSION_VIEWER"), new HashMap<String, String>(), testLabel + " fold-change. Significative results"));

			DiffExpressionUtils.createFatiScanRedirection(dataFrame, test, test, result, outdir, testLabel + " fold-change.");
		}

		
/*		
		List<Double> orderedRes = ListUtils.ordered(ArrayUtils.toList(res), rowOrder);
		int posValues = 0;
		int negValues = 0;
		for(int i=0 ; i<orderedRes.size() ; i++) {
			if (Math.abs(orderedRes.get(i))>foldChangeValue) {
				if (orderedRes.get(i)>0) {
					posValues++;
				} else {
					negValues++;
				}
			}
		}
		
		if (posValues + negValues == 0) {
			result.addOutputItem(new Item("no_sig_results", "No significative results (fold-change value = " + foldChangeValue + ")", "Significative results", TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(2), testLabel + " fold-change.Significative results"));															
			return;
		}
		
		int halfDisplay = maxDisplay/2;
		int posValuesToDisplay = posValues;
		int negValuesToDisplay = negValues;
		if (posValues + negValues > maxDisplay) {
			if (Math.min(posValues, negValues)>halfDisplay) {
				posValuesToDisplay = halfDisplay;
				negValuesToDisplay = halfDisplay;
			} else {
				posValuesToDisplay = posValues>negValues ? (maxDisplay-negValues) : posValues;
				negValuesToDisplay = negValues>posValues ? (maxDisplay-posValues) : negValues;
			}
		}
		
		int nbToDisplay = posValuesToDisplay + negValuesToDisplay;
		DoubleMatrix doubleMatrix = new DoubleMatrix(nbToDisplay, subDataset.getColumnDimension());
		int rowIndex = 0;
		int negLimit = rowOrder.length-negValuesToDisplay;
		List<Integer> sigRowIndexes = new ArrayList<Integer>();
		for(int i=0 ; i<rowOrder.length ; i++) {
			if (i<posValuesToDisplay || i>=negLimit) {
				doubleMatrix.setRow(rowIndex++, subDataset.getDoubleMatrix().getRow(rowOrder[i]));
				//System.out.println(subDataset.getFeatureNames().get(sigOrder[i]));
				sigRowIndexes.add(rowOrder[i]);
			}
		}
		
		file = new File(outdir + "/" + test +"_fold_change_significative_dataset.txt");
		Dataset sigDataset = new Dataset(subDataset.getSampleNames(), ListUtils.subList(subDataset.getFeatureNames(), ListUtils.toIntArray(sigRowIndexes)), doubleMatrix);
		sigDataset.setVariables(subDataset.getVariables());
		sigDataset.validate();
		sigDataset.save(file);
		if (file.exists()) {
			String tags = "datamatrix,expression";
			result.addOutputItem(new Item(test + "_sig_dataset", file.getName(), "Significative values dataset (fold-change value = " + foldChangeValue + ")", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(2),  testLabel + " fold-change.Significative results"));											

			File redirectionFile = new File(outdir + "/clustering.redirection");
			DiffExpressionUtils.createClusteringRedirectionFile(redirectionFile, file);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Clustering tool...)";
				result.addOutputItem(new Item(test + "_sig_dataset", file.getName(), "Significative values dataset (fold-change value = " + foldChangeValue + ")", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2),  testLabel + " fold-change.Significative results"));											
			}
		}
		
		rowOrder = ListUtils.order(ListUtils.subList(ArrayUtils.toList(res), ListUtils.toIntArray(sigRowIndexes)), true);	
				
		DataFrame sigDataFrame = new DataFrame(sigDataset.getFeatureNames().size(), 0);
		
		sigDataFrame.addColumn(test, ListUtils.toStringList(ListUtils.ordered(ListUtils.subList(ArrayUtils.toList(res), ListUtils.toIntArray(sigRowIndexes)), rowOrder)));
		sigDataFrame.setRowNames(ListUtils.ordered(ListUtils.subList(subDataset.getFeatureNames(), ListUtils.toIntArray(sigRowIndexes)), rowOrder));

		// adding table to results
		//
		file = new File(outdir + "/" + test + "fold_change_significative_table.txt");
		IOUtils.write(file, sigDataFrame.toString(true, true));
		if ( file.exists() ) {
			result.addOutputItem(new Item(test + "fold_change_table", file.getName(), "Significative values table (fold-change value = " + foldChangeValue + ")", TYPE.FILE, StringUtils.toList("TABLE," + test.toUpperCase() + "_FOLD_CHANGE_TABLE", ","), new HashMap<String, String>(2), testLabel + " fold-change.Significative results"));											
		}

		// adding heatmap to results
		//
		Canvas sigHeatmap = DiffExpressionUtils.generateHeatmap(sigDataset, className, columnOrder, rowOrder, testLabel, ListUtils.toDoubleArray(ListUtils.subList(ArrayUtils.toList(res), ListUtils.toIntArray(sigRowIndexes))), null, null);
		if (sigHeatmap == null) {
			printError("ioexception_execute_fold_change_classcomparison", "ERROR", "Error generating " + test + " fold-change heatmap image");
		} else {
			try {
				File sigHeatmapFile = new File(outdir + "/" + test + "fold_change_heatmap_significative.png");
				sigHeatmap.save(sigHeatmapFile.getAbsolutePath());
				if (sigHeatmapFile.exists()) {
					result.addOutputItem(new Item(test + "_fold_change_heatmap_significative", sigHeatmapFile.getName(), testLabel + " fold-change heatmap with significative values (fold-change value = " + foldChangeValue + ")", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), testLabel + " fold-change.Significative results"));
				}
			} catch (IOException e) {
				printError("ioexception_execute_fold_change_classcomparison", "ERROR", "Error saving " + test + " fold-change heatmap image");
			}
		}
		DiffExpressionUtils.createFatiGoRedirection(dataFrame.getRowNames(), dataFrame.getColumnAsDoubleArray(test), test, result, outdir, testLabel + " fold-change.");
		DiffExpressionUtils.createFatiScanRedirection(sigDataFrame, test, test, result, outdir, testLabel + " fold-change.");
*/
	}

	/**
	 *
	 */
	public void executeAnova() {

		DoubleMatrix matrix = null;
		List<String> vars = new ArrayList<String>();
		List<Integer> indices = new ArrayList<Integer>();
		List<String> values = dataset.getVariables().getByName(className).getValues();

		if ( values.size() == classValues.size() ) {
			matrix = dataset.getDoubleMatrix();
			vars = values;
		} else {
			for(int i=0 ; i<values.size() ; i++) {
				if ( classValues.contains(values.get(i)) ) {
					indices.add(i);
					vars.add(values.get(i));
				}
			}
			matrix = dataset.getSubMatrixByColumns(ListUtils.toIntArray(indices));
		}

		try {
			Dataset subDataset = dataset.getSubDataset(className, classValues);

			// apply test and multiple test correction according
			//
			AnovaTest anova = new AnovaTest(matrix, vars);
			TestResultList<AnovaTestResult> res = anova.compute();
			DiffExpressionUtils.multipleTestCorrection(res, correction);

			// create output file
			//
			int[] columnOrder = ListUtils.order(subDataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

			DataFrame dataFrame = new DataFrame(subDataset.getFeatureNames().size(), 0);
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
			dataFrame.setRowNames(ListUtils.ordered(subDataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);
			File file = new File(outdir + "/" + test + ".txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item(test + "file", file.getName(), "Anova output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Anova output files"));
			}

			// getting significative genes
			//
			DiffExpressionUtils.addSignificativeResults(subDataset, test, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues(), "p-value", res.getPValues(), null, null, null, null, className, columnOrder, pValue, maxDisplay, this);
			DiffExpressionUtils.createFatiScanRedirection(dataFrame, test, "statistic", result, outdir);
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_executeanova_classcomparison", "error running anova", "error running anova: " + e.getMessage(), "error running anova: " + e.getMessage());
		}
	}

	/**
	 * @throws InvalidColumnIndexException
	 * @throws IOException
	 *
	 */
	public void executeLimma() {

		Limma limma = null;
		if ( classValues.size() > 2 ) {
			limma = new Limma(babelomicsHomePath + "/bin/diffexp/limma_multiclasses.r");
		} else if ( classValues.size() == 2 ) {
			limma = new Limma(babelomicsHomePath + "/bin/diffexp/limma_twoclasses.r");
		} else if ( classValues.size() == 1 ) {
			limma = new Limma(babelomicsHomePath + "/bin/diffexp/limma_oneclass.r");
		} else {
			abort("testmismatched_executelimma_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");
		}

		//System.out.println("dataset = " + dataset.toString());
		System.out.println("class name = " + className);

		limma.setInputFilename(dataset.getDatasetFile().getAbsolutePath());
		limma.setClasses(dataset.getVariables().getByName(className).getValues());
		limma.setContrast(classValues);

		try {
			Dataset subDataset = dataset.getSubDataset(className, classValues);

			// apply test and multiple test correction according
			//
			TestResultList<LimmaTestResult> res = limma.compute();
			DiffExpressionUtils.multipleTestCorrection(res, correction);

			// create output file
			//
			int[] columnOrder = ListUtils.order(subDataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

			DataFrame dataFrame = new DataFrame(subDataset.getFeatureNames().size(), 0);
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
			dataFrame.setRowNames(ListUtils.ordered(subDataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);
			File file = new File(outdir + "/" + test + ".txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item(test + "file", file.getName(), "Limma output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Limma output files"));
			}

			// getting significative genes
			//
			DiffExpressionUtils.addSignificativeResults(subDataset, test, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues(), "p-value", res.getPValues(), null, null, null, null, className, columnOrder, pValue, maxDisplay, this);
			DiffExpressionUtils.createFatiScanRedirection(dataFrame, test, "statistic", result, outdir);
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_executelimma_classcomparison", "error running limma", "error running limma: " + e.toString(), "error running limma: " + e.toString());
		}
	}
}

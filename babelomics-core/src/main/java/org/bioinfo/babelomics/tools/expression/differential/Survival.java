package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.math.result.CoxTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.survival.CoxTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Survival extends BabelomicsTool {

	private double pValue = 0.05;
	private int maxDisplay = Integer.MAX_VALUE;
	
	public Survival() {
	}

	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("test", "the test, possible values: cox"));
		options.addOption(OptionFactory.createOption("time-class", "class variable", false));
		options.addOption(OptionFactory.createOption("censored-class", "class class variable", false));
		options.addOption(OptionFactory.createOption("correction", "Multiple-test correction: fdr, bh, by, bonferroni, hochberg, hold"));
		options.addOption(OptionFactory.createOption("p-value", "p-value for significative genes", false));
	}

	@Override
	public void execute() {
		
		// reading dataset
		//
		Dataset dataset = initDataset(new File(commandLine.getOptionValue("dataset")));
				
		String test = commandLine.getOptionValue("test", null);
		String timeClass = commandLine.getOptionValue("time-class", null);
		String censoredClass = commandLine.getOptionValue("censored-class", null);
		String correction = commandLine.getOptionValue("correction", "fdr");
		String pValueParam = commandLine.getOptionValue("p-value", "0.05");
		
		List<String> timeValues = null;
		try {
			timeValues = dataset.getVariables().getByName(timeClass).getValues();
		} catch (Exception e) {
			timeValues = new ArrayList<String>();
		}

		List<String> censoredValues = null;
		try {
			censoredValues = dataset.getVariables().getByName(censoredClass).getValues();
		} catch (Exception e) {
			censoredValues = new ArrayList<String>();
		}

		try {
			pValue = Double.parseDouble(pValueParam);
			if (pValue > 1 || pValue < 0) {
				pValue = 0.05;
				printWarning("pvalueinvalid_execute_cox_survival", "Warning", "The adjusted p-value " + pValueParam + " has been set to " + pValue);
			}
		} catch (NumberFormatException e) {
			pValue = 0.05;
			printWarning("pvalueinvalid_execute_cox_survival", "Warning", "The adjusted p-value " + pValueParam + " has been set to " + pValue);
		}

		// input parameters
		//
		//result.addOutputItem(new Item("dataset_input_param", (datasetParam == null ? "" : new File(datasetParam).getName()), "Dataset file name", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));								
		result.addOutputItem(new Item("test_input_param", test, "Test", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("correction_input_param", correction, "Multiple-test correction", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("timeclass_input_param", timeClass + " [" + (timeValues.size() > 0 ? ListUtils.toString(timeValues, ",") : "") + "]", "Time class", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("censoredclass_input_param", censoredClass + " [" + (censoredValues.size() > 0 ? ListUtils.toString(censoredValues, ",") : "") + "]", "Censored class", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("pvalue_input_param", pValueParam, "Adjusted p-value", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		
		
		if ( ! "cox".equalsIgnoreCase(test) ) {
			abort("unknowntest_execute_survival", "unknown test (" + test + ")", "unknown test (" + test + ")", "unknown test (" + test + ")");
		}

		String aux = "";
		List<Double> timeVars = new ArrayList<Double>();
		if (timeValues.size()>0) {
			try {
				for(int i=0 ; i<timeValues.size() ; i++) {
					aux = timeValues.get(i);
					timeVars.add(Double.parseDouble(aux));
				}
			} catch (NumberFormatException e) {
				abort("exception_execute_survival", "Error", "Invalid time-class value: " + aux + ". Time-class values must be numeric for " + test + " test", "");
			}
		} else {
			abort("exception_execute_survival", "Error", "Missing time-class values for class " + timeClass, "");			
		}
		
		//vars = dataset.getVariables().getByName(censoredClass).getValues();
		List<Double> censoredVars = new ArrayList<Double>();
		if (censoredValues.size()>0) {
			try {
				for(int i=0 ; i<censoredValues.size() ; i++) {
					aux = censoredValues.get(i);
					if ("1".equalsIgnoreCase(aux) || "0".equalsIgnoreCase(aux)) {
						censoredVars.add(Double.parseDouble(aux));
					} else {
						abort("exception_execute_survival", "Error", "Invalid censored-class value: " + aux + ". Censored-class values must be 0 or 1 for " + test + " test", "");
					}
				}
			} catch (NumberFormatException e) {
				abort("exception_execute_survival", "Error", "Invalid censored-class value: " + aux + ". Censored-class values must be 0 or 1 for " + test + " test", "");		}
		} else {
			abort("exception_execute_survival", "Error", "Missing censored-class values for class " + censoredClass, "");			
		}

		int[] columnOrder = ListUtils.order(timeVars);
		
		
		
		// cox test
		//
		updateJobStatus("40", "computing cox test");
		
		try {
			//Dataset subDataset = dataset.getSubDataset(className, classValues);

			// apply test and multiple test correction according
			//
			CoxTest coxtest = new CoxTest();
			TestResultList<CoxTestResult> res = coxtest.compute(dataset.getDoubleMatrix(), timeVars, censoredVars);
			DiffExpressionUtils.multipleTestCorrection(res, correction);

			// create output file
			//
			//int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getCoefs()), true);
			
			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("coeff.", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getCoefs()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);
			File file = new File(outdir + "/" + test + ".txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item(test + "_file", file.getName(), "Cox output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Cox output files"));							
			}

			// getting significative genes
			//
			DiffExpressionUtils.addSignificativeResults(dataset, test, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues(), "p-value", res.getPValues(), "coefs.", res.getCoefs(), null, null, timeClass, columnOrder, pValue, maxDisplay, this);
			DiffExpressionUtils.createFatiScanRedirection(dataFrame, test, "statistic", result, outdir);			
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_execute_survival", "ERROR", "Error running Cox survival test", "");
		}		
		
		
		
//		CoxTest coxtest = new CoxTest();
//		TestResultList<CoxTestResult> res = null;	
//		try {
//			res = coxtest.compute(dataset.getDoubleMatrix(), timeVars, censoredVars);
//			
//			// apply multiple test correction according to input correction
//			DiffExpressionUtils.multipleTestCorrection(res, correction);						
//		} catch (Exception e) {
//			abort("exception_run_cox", "error running cox test", e.toString(), StringUtils.getStackTrace(e));
//		}
//
//		//int[] columnOrder = ListUtils.order(vars);
//		int[] columnOrder = ListUtils.order(timeVars);
//		
////		System.out.println("time vars = " + ListUtils.toString(timeVars, ","));
////		System.out.println("column order = " + ArrayUtils.toString(columnOrder, ","));
//		
//		int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);
//
//		// generating heatmap
//		//
//		updateJobStatus("60", "generating heatmap");
//		Canvas heatmap = null, sigHeatmap = null;
//		try {
//			//heatmap = DiffExpressionUtils.generateHeatmap(dataset, timeClass, columnOrder, rowOrder, "coeff.", res.getCoefs(), "adj. p-value", res.getAdjPValues());
//			sigHeatmap = DiffExpressionUtils.generateSigHeatmap(dataset, timeClass, columnOrder, "coeff.", res.getCoefs(), "adj. p-value", res.getAdjPValues(), pValue);
//		} catch (Exception e) {
//			printError("ioexception_executesurvival_survival", "error generating heatmaps for " + test + " test", e.toString(), e);
//		}
//
//		// saving data
//		//
//		updateJobStatus("80", "saving results");
//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//
//		try {
//			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
//			dataFrame.addColumn("coeff.", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getCoefs()), rowOrder)));
//			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
//			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
//
//			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
//			FeatureData featureData = new FeatureData(dataFrame);
//			featureData.save(new File(getOutdir() + "/cox.txt"));
//			result.addOutputItem(new Item("cox_file", "cox.txt", "Cox output file", TYPE.FILE, StringUtils.toList("TABLE,COX_TABLE", ","), new HashMap<String, String>(2), "Output files"));
////			result.addOutputItem(new Item("cox_file", "cox.txt", "Cox output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));
////			
////			IOUtils.write(new File(getOutdir() + "/cox_table.txt"), dataFrame.toString(true, true));			
////			result.addOutputItem(new Item("cox_table", "cox_table.txt", "Cox table", TYPE.FILE, StringUtils.toList("TABLE,COX_TABLE", ","), new HashMap<String, String>(2), "Output files"));
//			
//			
//		} catch (Exception e) {
//			printError("ioexception_execute_survival", "error saving results", "error saving results");
//		}
//
//		String sigHeatmapFilename = getOutdir() + "/" + test + "_heatmap_significative.png";
//		if ( sigHeatmap == null ) {
//			result.addOutputItem(new Item(test + "_heatmap_significative", "None significative terms", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.TEXT, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
//		} else {
//			try {
//				sigHeatmap.save(sigHeatmapFilename);
//				if ( new File(sigHeatmapFilename).exists() ) {
//					result.addOutputItem(new Item(test + "_heatmap_significative", test + "_heatmap_significative.png", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
//				}
//			} catch (IOException e) {
//				printError("ioexception_execute_survival", "error generating heatmap", "error generating heatmap");
//			}
//		}
//
//		try {
//			DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);
//		} catch (Exception e) {
//			printError("ioexception_execute_survival", "error", "error creating redirection files");
//		}
//
////		String heatmapFilename = getOutdir() + "/cox_heatmap.png";
////		try {
////			heatmap.save(heatmapFilename);
////			if ( new File(heatmapFilename).exists() ) {
////				result.addOutputItem(new Item("cox_heatmap", "cox_heatmap.png", "Cox heatmap", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
////			}
////		} catch (IOException e) {
////			printError("ioexception_cox_cox", "error generating heatmap", e.toString(), e);
////		}
	}
}

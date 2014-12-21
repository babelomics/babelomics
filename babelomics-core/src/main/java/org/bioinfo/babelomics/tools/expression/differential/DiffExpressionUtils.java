package org.bioinfo.babelomics.tools.expression.differential;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.data.list.NamedArrayList;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.result.TestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.MultipleTestCorrection;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.bioinfo.tool.result.Result;

public class DiffExpressionUtils {

	public static Canvas generateHeatmap(Dataset dataset, String className, int[] columnOrder, int[] rowOrder, String infoName1, double[] infoList1, String infoName2, double[] infoList2) {

		int xHeatMap = 2;				
		int yHeatMap = 2;
		int rowDimension = dataset.getRowDimension();
		int columnDimension = dataset.getColumnDimension();
		int cellSide = 16;
		int rowLabelsWidth = 70;
		int colLabelsWidth = 140;
		int infoWidth = 140;
		//double min = dataset.getDoubleMatrix().getMinValue();
		//double max = dataset.getDoubleMatrix().getMaxValue();
		double offset, min, max, mean, deviation, standard;
		//System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + ")");

		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);

		GridPanel gridPanel = new GridPanel("", xHeatMap, yHeatMap, (columnDimension * cellSide) + rowLabelsWidth + infoWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), ((rowDimension+1) * cellSide) + colLabelsWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator());
		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
		gridTrack.setRowLabels(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
		List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
		for(int i=0 ; i <columnOrder.length ; i++) {
			columnLabels.add(dataset.getVariables().getByName(className).getValues().get(columnOrder[i]) + ": " +  dataset.getSampleNames().get(columnOrder[i]));
		}
		gridTrack.setColumnLabels(columnLabels);
		//gridTrack.setName();
		gridTrack.setTopRegion(colLabelsWidth);
		gridTrack.setLeftRegion(rowLabelsWidth);
		gridTrack.setRightRegion(infoWidth);

		double[] values = new double[columnOrder.length]; 
		int row, column;

		ScoreFeature feature;
		for(int i=0 ; i<rowOrder.length ; i++) {
			row = rowOrder[i];
			//			System.out.println("row = " + Arrays.toString(dataset.getDoubleMatrix().getRow(row)));
			//			System.out.println("row mean = " + MathUtils.mean(dataset.getDoubleMatrix().getRow(row)));
			//			System.out.println("row deviation = " + MathUtils.standardDeviation(dataset.getDoubleMatrix().getRow(row)));
			//			System.exit(-1);

			mean = dataset.getDoubleMatrix().getRowMean(row);
			deviation = dataset.getDoubleMatrix().getRowStdDeviation(row);
			min = Double.MAX_VALUE;
			max = Double.MIN_NORMAL;
			for(int j=0 ; j<columnOrder.length ; j++) {
				column = columnOrder[j];
				values[column] = (deviation == 0) ? Double.NaN : (dataset.getDoubleMatrix().get(row, column)-mean)/(deviation);
				if ( min > values[column] ) min = values[column];
				if ( max < values[column] ) max = values[column];
			}

			offset = ( min <= 0 ) ? Math.abs(min) : (-1 * min);
			//System.out.println("mean = " + mean + ", deviation = " + deviation + ", min = " + min + ", max = " + max + ", offset = " + offset);
			for(int j=0 ; j<columnOrder.length ; j++) {
				column = columnOrder[j];
				//System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
				//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
				//standard = (deviation == 0) ? Double.NaN : (dataset.getDoubleMatrix().get(row, column)-mean)/(deviation);
				standard = (values[column] + offset) / ( max + offset);
				//System.out.println("(value, standard) = (" + dataset.getDoubleMatrix().get(row, column) + ", " + standard + ")");
				feature = new ScoreFeature("name, " + dataset.getDoubleMatrix().get(row, column), "", 0, 0, standard);
				//feature = new ScoreFeature("name (" + row + ", " + column + ")", "", 0, 0, dataset.getDoubleMatrix().get(row, column));
				//feature.setJsFunction("http://www.cipf.es");
				//				gridTrack.setFeature(row, column, feature);
				gridTrack.setFeature(i, j, feature);
			}
		}

		List<String> aux = new ArrayList<String>();
		DecimalFormat df = new DecimalFormat("##0.0000");
		if (infoList1 != null) {
			aux.clear();
			for(int i=0 ; i<infoList1.length ; i++) {
				aux.add(df.format(infoList1[i]));
			}
			gridTrack.addInfo(new NamedArrayList(infoName1, ListUtils.ordered(aux, rowOrder)));
		}
		if (infoList2 != null) {
			aux.clear();
			for(int i=0 ; i<infoList2.length ; i++) {
				aux.add(df.format(infoList2[i]));
			}
			gridTrack.addInfo(new NamedArrayList(infoName2, ListUtils.ordered(aux, rowOrder)));
		}

		gridPanel.add(gridTrack);
		canvas.addPanel(gridPanel);		
		canvas.render();

		return canvas;
	}

	public static Canvas generateSigHeatmap(Dataset dataset, String className, int[] columnOrder, String statisticsLabel, double[] statistics, String adjPValuesLabel, double[] adjPValues, double pValue) throws IOException, InvalidIndexException {
		List<Integer> filterRowIndexes = new ArrayList<Integer>();
		List<Integer> sigRowIndexes = new ArrayList<Integer>();
		for(int i=0 ; i<adjPValues.length ; i++) {
			if ( adjPValues[i] <= pValue ) { 
				sigRowIndexes.add(i);
			} else {
				filterRowIndexes.add(i);
			}
		}

		if ( sigRowIndexes.size() > 0  ) {				
			List<Double> sigStatistics = ListUtils.subList(ArrayUtils.toList(statistics), ListUtils.toIntArray(sigRowIndexes));
			List<Double> sigAdjPValues = ListUtils.subList(ArrayUtils.toList(adjPValues), ListUtils.toIntArray(sigRowIndexes));

			int[] sigRowOrder = ListUtils.order(sigStatistics, true);

			if ( filterRowIndexes.size() == 0 ) {
				return DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, sigRowOrder, statisticsLabel, ListUtils.toDoubleArray(sigStatistics), adjPValuesLabel, ListUtils.toDoubleArray(sigAdjPValues));				
			} else {
				Dataset sigDataset = dataset.filterRows(filterRowIndexes);
				sigDataset.validate();
				//			sigDataset.save("/tmp/sig_dataset.txt");
				//			System.out.println("column dimension = " + sigDataset.getColumnDimension());
				//			System.out.println("row dimension = " + sigDataset.getRowDimension());
				//			System.out.println("sig row indexes size = " + sigRowIndexes.size());
				//			System.out.println("sig statistics size = " + sigStatistics.size());
				//			System.out.println("sig adj p-values size = " + sigAdjPValues.size());
				//			System.out.println("sig row indexes = " + ListUtils.toString(sigRowIndexes));


				return DiffExpressionUtils.generateHeatmap(sigDataset, className, columnOrder, sigRowOrder, statisticsLabel, ListUtils.toDoubleArray(sigStatistics), adjPValuesLabel, ListUtils.toDoubleArray(sigAdjPValues));
			}
		}
		return null;
	}

	/**
	 * 
	 * @param res
	 */
	public static void multipleTestCorrection(TestResultList<? extends TestResult> res, String correction) {
		if ( "bonferroni".equalsIgnoreCase(correction) ) {
			MultipleTestCorrection.BonferroniCorrection(res);	
		} else if ( "bh".equalsIgnoreCase(correction) ) {
			MultipleTestCorrection.BHCorrection(res);	
		} else if ( "by".equalsIgnoreCase(correction) ) {
			MultipleTestCorrection.BYCorrection(res);	
		} else if ( "hochberg".equalsIgnoreCase(correction) ) {
			MultipleTestCorrection.HochbergCorrection(res);	
		} else if ( "holm".equalsIgnoreCase(correction) ) {
			MultipleTestCorrection.HolmCorrection(res);	
		} else {
			MultipleTestCorrection.FDRCorrection(res);	
		}
	}

	public static void generateRankedList(List<String> names, List<String> statistics, String label, File outFile) throws IOException, InvalidIndexException {
		// preparing list to fatiscan
		//
		DataFrame df = new DataFrame(names.size(), StringUtils.toList(label));
		if ( names.size() != statistics.size() ) {
			throw new InvalidIndexException("mismatch size, names size = " + names.size() + " and statistics size = " + statistics.size());
		}

		for(int i=0 ; i<statistics.size() ; i++) {
			df.addRow(names.get(i), StringUtils.toList(statistics.get(i)));
		}
		FeatureData outFeatureData = new FeatureData(df);
		outFeatureData.save(outFile);
	}

	public static void generateSignificantLists(List<String> names, List<String> adjPvalues, double pvalue, File topFile, File bottomFile) throws IOException, InvalidIndexException {
		// preparing list to fatigo
		//
		if ( names.size() != adjPvalues.size() ) {
			throw new InvalidIndexException("mismatch size, names size = " + names.size() + " and adj p-values size = " + adjPvalues.size());
		}
		int size = names.size();
		List<String> topNames = new ArrayList<String> ();
		List<String> bottomNames = new ArrayList<String> ();
		boolean endTop = false;
		boolean endBottom = false;
		for(int i=0 ; i<size ; i++) {

			// top list
			//
			if ( !endTop && Double.parseDouble(adjPvalues.get(i)) <= pvalue ) {
				topNames.add(names.get(i));
			} else {
				endTop = true;
			}

			// bottom list
			//
			if ( !endBottom && Double.parseDouble(adjPvalues.get(size - (i + 1))) <= pvalue ) {
				bottomNames.add(names.get(size - (i + 1)));
			} else {
				endBottom = true;
			}
		}
		if ( topNames.size() > 0 ) IOUtils.write(topFile, topNames);
		if ( bottomNames.size() > 0 ) IOUtils.write(bottomFile, bottomNames);		
	}

	public static void addOutputLists(DataFrame dataFrame, String test, String colName, Result result, String outDir) throws IOException, InvalidIndexException {	
		// preparing ranked list
		//
		String tags;
		File redirectionFile;
		File rankedListFile = new File(outDir + "/" + test + "_ranked_list.txt");
		DiffExpressionUtils.generateRankedList(dataFrame.getRowNames(), dataFrame.getColumn(colName), colName, rankedListFile);
		if ( rankedListFile.exists() ) {
			redirectionFile = new File(outDir + "/fatiscan.redirection");
			createFatiScanRedirectionFile(redirectionFile, rankedListFile);
			if ( redirectionFile.exists() ) {
				//				tags = "DATA,RANKED,REDIRECTION(" + redirectionFile.getName() + ":Send to FatiScan tool...)";
				//				result.addOutputItem(new Item(test + "_ranked_list_file", rankedListFile.getName(), "Send ranked list to FatiScan tool", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Continue processing"));
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiScan tool...)";
				result.addOutputItem(new Item(test + "_fatiscan", "", "Send ranked list to FatiScan tool", TYPE.TEXT, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Continue processing"));
			}
		}				

		// preparing significant list (top and bottom)
		//
		File topListFile = new File(outDir + "/" + test + "_top_list.txt");
		File bottomListFile = new File(outDir + "/" + test + "_bottom_list.txt");			
		DiffExpressionUtils.generateSignificantLists(dataFrame.getRowNames(), dataFrame.getColumn("adj. p-value"), 0.005, topListFile, bottomListFile);
		if ( topListFile.exists() || bottomListFile.exists() ) {
			redirectionFile = new File(outDir + "/fatigo.redirection");
			createFatiGORedirectionFile(redirectionFile, topListFile, bottomListFile);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiGO tool...)";
				result.addOutputItem(new Item(test + "_fatigo", "", "Send significative results to FatiGO tool", TYPE.TEXT, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Continue processing"));
			}
		}
	}

	public static void createFatiGoRedirection(List<String> names, double[] statistic, String test, Result result, String outDir) throws IOException, InvalidIndexException {
		createFatiGoRedirection(names, statistic, test, result, outDir, "");		
	}

	public static void createFatiGoRedirection(List<String> names, double[] statistic, String test, Result result, String outDir, String groupPrefix) throws IOException, InvalidIndexException {
		// preparing top and bottom list
		//
		String tags;

		int size = names.size();
		if (size != statistic.length) {
			throw new InvalidIndexException("mismatch size, names size = " + names.size() + " and statistic size = " + statistic.length);
		}

		//		for(int i=0 ; i<size ; i++) {
		//			System.out.println(names.get(i) + "\t" + statistic[i]);
		//		}

		List<String> topNames = new ArrayList<String> ();
		List<String> bottomNames = new ArrayList<String> ();
		for(int i=0 ; i<size ; i++) {
			if (statistic[i] > 0) {
				topNames.add(names.get(i));
			} else {
				bottomNames.add(names.get(i));
			}
		}

		File topListFile = 	new File(outDir + "/" + test + "_top_list.txt");
		if ( topNames.size() > 0 ) {
			IOUtils.write(topListFile, topNames);
		}
		File bottomListFile = new File(outDir + "/" + test + "_bottom_list.txt");			
		if ( bottomNames.size() > 0 ) {
			IOUtils.write(bottomListFile, bottomNames);		
		}

		if (topListFile.exists() && bottomListFile.exists()) {
			File redirectionFile = new File(outDir + "/" + test + "_top_bottom.fatigo.redirection");		
			createFatiGoRedirectionFile(redirectionFile, topListFile, bottomListFile);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiGO tool...)";
				result.addOutputItem(new Item(test + "_top_genome_fatigo", "", "Send top list vs bottom list to FatiGO tool", TYPE.TEXT, StringUtils.toList(tags, ","), new HashMap<String, String>(2), groupPrefix + "Continue processing"));
			}						
		}
		if (topListFile.exists()) {
			File redirectionFile = new File(outDir + "/" + test + "_top_genome.fatigo.redirection");		
			createFatiGoRedirectionFile(redirectionFile, topListFile);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiGO tool...)";
				result.addOutputItem(new Item(test + "_top_genome_fatigo", "", "Send top list vs genome to FatiGO tool", TYPE.TEXT, StringUtils.toList(tags, ","), new HashMap<String, String>(2), groupPrefix + "Continue processing"));
			}			
		}
		if (bottomListFile.exists()) {
			File redirectionFile = new File(outDir + "/" + test + "_bottom_genome.fatigo.redirection");		
			createFatiGoRedirectionFile(redirectionFile, bottomListFile);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiGO tool...)";
				result.addOutputItem(new Item(test + "_bottom_genome_fatigo", "", "Send bottom list vs genome to FatiGO tool", TYPE.TEXT, StringUtils.toList(tags, ","), new HashMap<String, String>(2), groupPrefix + "Continue processing"));
			}						
		}
	}

	public static void createFatiGoRedirectionFile(File redirectionFile, File listFile) {
		List<String> redirectionInputs = new ArrayList<String>();
		if ( listFile.exists() ) {
			redirectionInputs.add("comparison=list2genome");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + listFile.getName() + " (list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/" + listFile.getName());
			redirectionInputs.add("tool=fatigo");
			redirectionInputs.add("jobname=fatigo");
			redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
			try {
				IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public static void createFatiGoRedirectionFile(File redirectionFile, File topListFile, File bottomListFile) {

		List<String> redirectionInputs = new ArrayList<String>();

		if ( topListFile.exists() && bottomListFile.exists() ) {
			redirectionInputs.add("comparison=list2list");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + topListFile.getName() + " (top list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/" + topListFile.getName());

			redirectionInputs.add("list2_wum_data=true");
			redirectionInputs.add("list2_databox=" + bottomListFile.getName() + " (bottom list from job $JOB_NAME)");
			redirectionInputs.add("list2=$JOB_FOLDER/" + bottomListFile.getName());

			redirectionInputs.add("tool=fatigo");
			redirectionInputs.add("jobname=fatigo");
			redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
			try {
				IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	public static void createFatiScanRedirection(DataFrame dataFrame, String test, String colName, Result result, String outDir) throws IOException, InvalidIndexException {
		createFatiScanRedirection(dataFrame, test, colName, result, outDir, "");
	}

	public static void createFatiScanRedirection(DataFrame dataFrame, String test, String colName, Result result, String outDir, String groupPrefix) throws IOException, InvalidIndexException {
		// preparing ranked list
		//
		String tags;
		File redirectionFile;
		File rankedListFile = new File(outDir + "/" + test + "_ranked_list.txt");
		DiffExpressionUtils.generateRankedList(dataFrame.getRowNames(), dataFrame.getColumn(colName), colName, rankedListFile);
		if ( rankedListFile.exists() ) {
			redirectionFile = new File(outDir + "/fatiscan.redirection");
			createFatiScanRedirectionFile(redirectionFile, rankedListFile);
			if ( redirectionFile.exists() ) {
				//				tags = "DATA,RANKED,REDIRECTION(" + redirectionFile.getName() + ":Send to FatiScan tool...)";
				//				result.addOutputItem(new Item(test + "_ranked_list_file", rankedListFile.getName(), "Send ranked list to FatiScan tool", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Continue processing"));
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiScan tool...)";
				result.addOutputItem(new Item(test + "_fatiscan", "", "Send ranked list to FatiScan tool", TYPE.TEXT, StringUtils.toList(tags, ","), new HashMap<String, String>(2), groupPrefix + "Continue processing"));
			}
		}				
	}

	public static void createFatiScanRedirectionFile(File redirectionFile, File rankedListFile) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=fatiscan");
		redirectionInputs.add("jobname=fatiscan");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("ranked_list_databox=" + rankedListFile.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("ranked_list=$JOB_FOLDER/" + rankedListFile.getName());
		redirectionInputs.add("ranked_list_wum_data=true");
		redirectionInputs.add("method=fatiscan");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}

	public static void createFatiGORedirectionFile(File redirectionFile, File topListFile, File bottomListFile) {

		List<String> redirectionInputs = new ArrayList<String>();

		if ( topListFile.exists() && bottomListFile.exists() ) {
			redirectionInputs.add("comparison=list2list");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + topListFile.getName() + " (top list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/" + topListFile.getName());

			redirectionInputs.add("list2_wum_data=true");
			redirectionInputs.add("list2_databox=" + bottomListFile.getName() + " (bottom list from job $JOB_NAME)");
			redirectionInputs.add("list2=$JOB_FOLDER/" + bottomListFile.getName());
		} else if ( topListFile.exists() ) {
			redirectionInputs.add("comparison=list2genome");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + topListFile.getName() + " (top list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/" + topListFile.getName());
		} else if ( bottomListFile.exists() ) {
			redirectionInputs.add("comparison=list2genome");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + bottomListFile.getName() + " (top list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/" + bottomListFile.getName());
		}		

		redirectionInputs.add("tool=fatigo");
		redirectionInputs.add("jobname=fatigo");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}

	public static DataFrame getSigDataFrame(DataFrame dataFrame, double pValue, int maxSize) {
		DataFrame df = new DataFrame();
		//		for(int i=0 ; i<adjPValues.length ; i++) {
		//			if ( adjPValues[i] <= pValue ) { 
		//				sigRowIndexes.add(i);
		//			}
		//		}
		return df;
	}

	public static int getNumberOfSigValues(double[] adjPValues, double pValue) {
		int count = 0;
		for(int i=0 ; i<adjPValues.length ; i++) {
			if ( adjPValues[i] <= pValue ) { 
				count++;
			}
		}
		return count;
	}

	public static List<Integer> getSigIndexes(double[] adjPValues, double pValue) {
		List<Integer> indexes = new ArrayList<Integer>();
		for(int i=0 ; i<adjPValues.length ; i++) {
			if ( adjPValues[i] <= pValue ) { 
				indexes.add(i);
			}
		}
		return indexes;
	}

	public static List<Integer> getNoSigIndexes(double[] adjPValues, double pValue) {
		List<Integer> indexes = new ArrayList<Integer>();
		for(int i=0 ; i<adjPValues.length ; i++) {
			if ( adjPValues[i] > pValue ) { 
				indexes.add(i);
			}
		}
		return indexes;
	}


	public static void addSignificativeResults(Dataset subDataset, String test,
			String statLabel, double[] statistics,
			String adjPValLabel, double[] adjPValues, 
			String pValLabel, double[] pValues,
			String label1, double[] values1,
			String label2, double[] values2,
			String className, int[] columnOrder,
			double pValue, int maxDisplay, BabelomicsTool tool) throws IOException, InvalidIndexException {

		boolean isCox = false, isRegression = false, isCorrelation = false;
		if ("cox".equalsIgnoreCase(test) && label1 != null && values1 != null) {
			isCox = true;
		} else if (("pearson".equalsIgnoreCase(test) || "spearman".equalsIgnoreCase(test)) && label1 != null && values1 != null) {
			isCorrelation = true;
		} else if ("regression".equalsIgnoreCase(test) && label1 != null && values1 != null && label2 != null && values2 != null) {
			isRegression = true;
		}
		System.out.println("is cox ? " + isCox);
		System.out.println("is regression ? " + isRegression);
		System.out.println("is correlation ? " + isCorrelation);

		List<Integer> sigIndexes = getSigIndexes(adjPValues, pValue);
		int numberSigValues = (sigIndexes == null ? 0 : sigIndexes.size());
		if (numberSigValues > 0) {

			//System.out.println("number of sig. values = " + numberSigValues);

			int nbToDisplay = numberSigValues;
			if (numberSigValues > maxDisplay) {
				nbToDisplay = maxDisplay;
			}

			//System.out.println("number to display = " + nbToDisplay);

			int[] sigOrder = ListUtils.order(ArrayUtils.toList(adjPValues));

			DoubleMatrix doubleMatrix = new DoubleMatrix(nbToDisplay, subDataset.getColumnDimension());

			List<Integer> sigRowIndexes = new ArrayList<Integer>();
			System.out.println("feature names for sig");
			for(int i=0 ; i<sigOrder.length ; i++) {
				if (i<nbToDisplay) {
					doubleMatrix.setRow(i, subDataset.getDoubleMatrix().getRow(sigOrder[i]));
					//System.out.println(subDataset.getFeatureNames().get(sigOrder[i]));
					sigRowIndexes.add(sigOrder[i]);
				}
			}

//			if (numberSigValues > maxDisplay) {
//				tool.getResult().addOutputItem(new Item("sig_results", "" + numberSigValues + " (" + maxDisplay + " most significative values will be displayed)", "Number of significative results (adj. p-value = " + pValue + ")", TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(2), "Significative results"));																				
//			} else {
//				tool.getResult().addOutputItem(new Item("sig_results", "" + numberSigValues, "Number of significative results (adj. p-value = " + pValue + ")", TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(2), "Significative results"));																									
//			}

			// adding dataset to results
			//
			File file = new File(tool.getOutdir() + "/" + test +"_significative_dataset.txt");
			Dataset sigDataset = new Dataset(subDataset.getSampleNames(), ListUtils.subList(subDataset.getFeatureNames(), ListUtils.toIntArray(sigRowIndexes)), doubleMatrix);
			sigDataset.setVariables(subDataset.getVariables());
			sigDataset.validate();
			sigDataset.save(file);
//			if (file.exists()) {
//				String tags = "datamatrix,expression";
//				tool.getResult().addOutputItem(new Item(test + "_sig_dataset", file.getName(), "Significative values dataset (adj. p-value = " + pValue + ")", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Significative results"));											
//
//				File redirectionFile = new File(tool.getOutdir() + "/clustering.redirection");
//				createClusteringRedirectionFile(redirectionFile, file);
//				if ( redirectionFile.exists() ) {
//					tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Clustering tool...)";
//					tool.getResult().addOutputItem(new Item(test + "_sig_dataset", file.getName(), "Significative values dataset (adj. p-value = " + pValue + ")", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Significative results"));											
//				}
//			}


			int[] rowOrder = null;
			if (isRegression || isCorrelation) {
				rowOrder = ListUtils.order(ListUtils.subList(ArrayUtils.toList(values1), ListUtils.toIntArray(sigRowIndexes)), true);				
			} else {
				rowOrder = ListUtils.order(ListUtils.subList(ArrayUtils.toList(statistics), ListUtils.toIntArray(sigRowIndexes)), true);
			}

			DataFrame dataFrame = new DataFrame(sigDataset.getFeatureNames().size(), 0);
			dataFrame.addColumn(statLabel, ListUtils.toStringList(ListUtils.ordered(ListUtils.subList(ArrayUtils.toList(statistics), ListUtils.toIntArray(sigRowIndexes)), rowOrder)));
			if (isCox || isCorrelation) {
				dataFrame.addColumn(label1, ListUtils.toStringList(ListUtils.ordered(ListUtils.subList(ArrayUtils.toList(values1), ListUtils.toIntArray(sigRowIndexes)), rowOrder)));				
			} else if (isRegression) {
				dataFrame.addColumn(label1, ListUtils.toStringList(ListUtils.ordered(ListUtils.subList(ArrayUtils.toList(values1), ListUtils.toIntArray(sigRowIndexes)), rowOrder)));				
				dataFrame.addColumn(label2, ListUtils.toStringList(ListUtils.ordered(ListUtils.subList(ArrayUtils.toList(values2), ListUtils.toIntArray(sigRowIndexes)), rowOrder)));								
			}
			dataFrame.addColumn(pValLabel, ListUtils.toStringList(ListUtils.ordered(ListUtils.subList(ArrayUtils.toList(pValues), ListUtils.toIntArray(sigRowIndexes)), rowOrder)));
			dataFrame.addColumn(adjPValLabel, ListUtils.toStringList(ListUtils.ordered(ListUtils.subList(ArrayUtils.toList(adjPValues), ListUtils.toIntArray(sigRowIndexes)), rowOrder)));
			dataFrame.setRowNames(ListUtils.ordered(ListUtils.subList(subDataset.getFeatureNames(), ListUtils.toIntArray(sigRowIndexes)), rowOrder));

			//System.out.println("names from data frame :\n" + ListUtils.toString(ListUtils.subList(subDataset.getFeatureNames(), ListUtils.toArray(sigRowIndexes)), "\n"));
			//System.out.println("names from data set   :\n" + ListUtils.toString(sigDataset.getFeatureNames(), "\n"));

			// adding table to results
			//
			file = new File(tool.getOutdir() + "/" + test + "_significative_table.txt");
			IOUtils.write(file, dataFrame.toString(true, true));
			if ( file.exists() ) {
				String table = "DIFF_EXPRESSION_TABLE";
				if (isCox) {
					table = "COX_TABLE";
				} else if (isCorrelation) {
					table = "CORRELATION_TABLE";
				} else if (isRegression) {
					table = "REGRESSION_TABLE";
				}
				//tool.getResult().addOutputItem(new Item(test + "_table", file.getName(), "Significative values table (adj, p-value = " + pValue + ")", TYPE.FILE, StringUtils.toList("TABLE," + table, ","), new HashMap<String, String>(2), "Significative results"));											
			}

			// adding heatmap to results
			//
//			Canvas sigHeatmap = null;
//			if (isCox || isCorrelation || isRegression) {
//				sigHeatmap = DiffExpressionUtils.generateHeatmap(sigDataset, className, columnOrder, rowOrder, label1, ListUtils.toDoubleArray(ListUtils.subList(ArrayUtils.toList(values1), ListUtils.toIntArray(sigRowIndexes))), adjPValLabel, ListUtils.toDoubleArray(ListUtils.subList(ArrayUtils.toList(adjPValues), ListUtils.toIntArray(sigRowIndexes))));
//			} else {
//				sigHeatmap = DiffExpressionUtils.generateHeatmap(sigDataset, className, columnOrder, rowOrder, statLabel, ListUtils.toDoubleArray(ListUtils.subList(ArrayUtils.toList(statistics), ListUtils.toIntArray(sigRowIndexes))), adjPValLabel, ListUtils.toDoubleArray(ListUtils.subList(ArrayUtils.toList(adjPValues), ListUtils.toIntArray(sigRowIndexes))));
//			}
//			if (sigHeatmap == null) {
//				tool.printError("ioexception_executet_classcomparison", "ERROR", "Error generating heatmap image");
//			} else {
//				try {
//					File sigHeatmapFile = new File(tool.getOutdir() + "/" + test + "_significative_heatmap.png");
//					sigHeatmap.save(sigHeatmapFile.getAbsolutePath());
//					if (sigHeatmapFile.exists()) {
//						tool.getResult().addOutputItem(new Item(test + "_significative_heatmap", sigHeatmapFile.getName(), test.toUpperCase() + " heatmap with significative values (adj. p-value = " + pValue + ")", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Significative results"));
//					}
//				} catch (IOException e) {
//					tool.printError("ioexception_executet_classcomparison", "ERROR", "Error saving heatmap image");
//				}
//			}
//			createFatiGoRedirection(ListUtils.ordered(subDataset.getFeatureNames(), sigOrder), ListUtils.toDoubleArray(ListUtils.ordered(ArrayUtils.toList(statistics), sigOrder)), test, tool.getResult(), tool.getOutdir());
		} else {
//			tool.getResult().addOutputItem(new Item("no_sig_results", "No significative results (p-value = " + pValue + ")", "Significative results", TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(2), "Significative results"));															
		}	
		
		String json = "{\\\"paramfilename\\\": \\\"input_params.txt\\\", \\\"testfilename\\\": \\\"" + test + ".txt\\\"}";
		tool.getResult().addOutputItem(new Item("diff_expr_" + StringUtils.randomString(8), json, "Significative results", TYPE.FILE, StringUtils.toList("DIFF_EXPRESSION_VIEWER"), new HashMap<String, String>(), "Significative results"));															
	}

	public static void createClusteringRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=clustering");
		redirectionInputs.add("jobname=clustering");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("method=upgma");
		redirectionInputs.add("gene_clustering=true");
		redirectionInputs.add("distance=euclidean");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}

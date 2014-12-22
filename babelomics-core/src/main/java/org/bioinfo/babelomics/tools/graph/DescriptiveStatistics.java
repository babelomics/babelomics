package org.bioinfo.babelomics.tools.graph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.methods.expression.clustering.ClusteringUtils;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.babelomics.utils.RCommand;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.chart.HistogramChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.format.io.exception.InvalidFileFormatException;
import org.bioinfo.data.format.io.parser.NewickParser;
import org.bioinfo.data.tree.multiway.MultipleTree;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;

public class DescriptiveStatistics extends BabelomicsTool {
	Dataset dataset;
	String test;
	String className;
	String axeTitle;
	List<String> values;
	List<Double> doubleVars;
	String correction;
	String type;

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("datalist", "the feature data containig the ranked list"));
		options.addOption(OptionFactory.createOption("tree", "",false));
		options.addOption(OptionFactory.createOption("histogram", "",false));
		options.addOption(OptionFactory.createOption("boxplot", "",false));
		options.addOption(OptionFactory.createOption("pcaplot", "",false));
		//options.addOption(OptionFactory.createOption("histogramboxplot", "",false,false));		
		options.addOption(OptionFactory.createOption("classname", "",false,true));
		options.addOption(OptionFactory.createOption("title", "",false,true,' '));
		options.addOption(OptionFactory.createOption("width", "",false,true));
		options.addOption(OptionFactory.createOption("height", "",false,true));
	}


	@Override
	protected void execute() {

		updateJobStatus("10", "init execution");
		dataset = null;
		className = commandLine.getOptionValue("classname", null);
		axeTitle = commandLine.getOptionValue("title", "");
		axeTitle = axeTitle.replaceAll("_____", " ");
		//axeTitle = ArrayUtils.toString(commandLine.getOptionValues("title")," ");

		try {
			
			if(commandLine.hasOption("histogram")){
				BoxPlotChart bp = (commandLine.getOptionValue("boxplot", null) != null) ? new BoxPlotChart("Box-plot","", axeTitle):null;
				HistogramChart hc = new HistogramChart("Histogram chart", axeTitle, "frec");
				preloadSeries(hc,bp);
				finishGraph(hc, bp);
			}else 
				if (commandLine.getOptionValue("boxplot")!=null){
					BoxPlotChart bp = new BoxPlotChart("Box-plot", "", "");
					preloadSeries(null,bp);
					finishGraph(null, bp);
				} else if(commandLine.getOptionValue("tree") != null){
					doTree();
				} else if(commandLine.hasOption("pcaplot")){
					doPcaPlot(babelomicsHomePath + "/bin/plots/pcaplot.r");
				}
		} catch(IOException e) {
			printError("DescriptiveStatistics", "An error occurred", e.toString(), e);
		} catch(Exception e) {

		}
	}


	private void doPcaPlot(String pcaPlotBinPath) throws InvalidParameterException, IOException  {

		String pcaplotFileName = "pca.png";

		List<String> env = new ArrayList<String>();
		env.add("datafile=" + commandLine.getOptionValue("datalist"));		
		env.add("outfile=" + this.getOutdir() + "/" + pcaplotFileName);
		//		Command cmd = new Command("/usr/local/R-292/bin/R CMD BATCH --verbose --no-save --no-restore " + pcaPlotBinPath + " " + this.getOutdir() + "/rLog.log", env);		
		//		System.out.println("command = " + cmd.getCommandLine());
		//		System.out.println("env = " + ListUtils.toString(env, " "));
		//		cmd.run();

		//		System.out.println("out std = " + cmd.getOutput());
		//		System.out.println("out std = " + cmd.getStatus());
		//		System.out.println("out std = " + cmd.getException());
		//		System.out.println("out error = " + cmd.getError());

		System.err.println("salidaPass1:");

		try {
			jobStatus.addStatusMessage("" + ("40"), "preparing execution");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
		String classField = "";
		updateJobStatus("40", "preparing execution");

		//dataset = new Dataset(new File(commandLine.getOptionValue("datalist")));

		File dataset = new File(commandLine.getOptionValue("datalist"));
		FileUtils.checkFile(commandLine.getOptionValue("datalist"));
		BufferedReader br = new BufferedReader(new FileReader(dataset.getAbsolutePath()));
		String line = "";
		String[] attr;
		boolean hasName = false;
		while((line = br.readLine()) != null) {
			line = line.trim();
			attr = line.split("\t");
			if(line.startsWith("#")) {
				System.out.println(attr[0]);
			}if(attr[0].indexOf("#NAMES") > -1) {
				hasName=true;
				RCommand rCommand = new RCommand(pcaPlotBinPath, this.getOutdir());

				if (className!=null){
					rCommand.addParam("grupos", className);	
				}

				rCommand.addParam("datafile", commandLine.getOptionValue("datalist"));
				rCommand.addParam("outfile", this.getOutdir() + "/" + pcaplotFileName);
				updateJobStatus("80", "exec");
				rCommand.exec();
				// saving results
				//
				updateJobStatus("90", "saving results");
				result.addOutputItem(new Item("pca_plot", pcaplotFileName,"pca plot (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","PCA_IMAGE"),new HashMap<String,String>(),"pca image"));
				break;
			}
		}
		br.close();

		if (!hasName)result.addOutputItem(new Item("pca_plot","PCA analisys can't be continue with you data. '#NAMES' variable doesn't found. This variable is used as 'id' to tag your samples in PCA analisys and plotting","PCA warning result",Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"PCA warning result"));

	}



	private void preloadSeries(HistogramChart hc, BoxPlotChart bp) throws IOException {
		dataset = new Dataset(new File(commandLine.getOptionValue("datalist")));
		if (className != null && !className.equalsIgnoreCase("all_classes") && "CATEGORICAL".equalsIgnoreCase(dataset.getVariables().getByName(className).getType()) ){
			System.out.println("---------------------"+dataset.getVariables().getByName(className).getType());
			values = dataset.getVariables().getByName(className).getLabels();
			for (String str: values){
				int[] colIndexByVariableValue = dataset.getColumnIndexesByVariableValue(className, str);
				doubleVars = new ArrayList<Double>(colIndexByVariableValue.length);
				DoubleMatrix matrixByVal =dataset.getSubMatrixByColumns(colIndexByVariableValue);
				addSeries(matrixByVal, hc, bp, str);
			}
		}
		else{
			addSeries(dataset.getDoubleMatrix(),hc,bp, null);
		}
	}


	private void addSeries(DoubleMatrix matrixByVal, HistogramChart hc, BoxPlotChart bp, String label) {
		String title = (label == null) ? "" : label;
		int[] variablesByval = null;
		List<String> samplesNames = dataset.getSampleNames();
		if(label != null){
			variablesByval = dataset.getColumnIndexesByVariableValue(className, label);	
		}
		if (hc != null){
			if (matrixByVal.getRowDimension() != 0 && matrixByVal.getColumnDimension() != 0  ) {
				for(int i=0; i < matrixByVal.getColumnDimension(); i++) {
					String sampleName ="";
					if(label != null){
						sampleName = samplesNames.get(variablesByval[i]);
						hc.addSeries(matrixByVal.getColumn(i), sampleName+"("+ label +")");
					}
					else{
						sampleName = dataset.getSampleNames().get(i);
						hc.addSeries(matrixByVal.getColumn(i), sampleName);	
					}

				}
			}
			else{
				hc.addSeries(matrixByVal.getColumn(0),title);	
			}

			//hc.removeLegend();		
		}
		if(bp != null ){
			if (matrixByVal.getRowDimension() != 0 && matrixByVal.getColumnDimension() != 0  ) {

				for(int i=0; i < matrixByVal.getColumnDimension(); i++) {
					String sampleName = "";
					if(label != null){
						sampleName = samplesNames.get(variablesByval[i]);
						bp.addSeries(matrixByVal.getColumn(i),sampleName+"("+ label +")", "");
					}else{
						sampleName = dataset.getSampleNames().get(i);
						bp.addSeries(matrixByVal.getColumn(i),sampleName, "");
					}	

				}
			}
			else{
				bp.addSeries(matrixByVal.getColumn(0),title, "");	
			}

			//bp.removeLegend();
		}
	}


	private void finishGraph(HistogramChart hc, BoxPlotChart bp){
		String imgFilename ="";
		int progress = 1;
		int finalProgress = 3;

		try {
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "reading ranked list");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		progress ++;
		updateJobStatus(""+progress,"reading ok");


		// Generate histogram
		if (hc != null){
			imgFilename = this.getOutdir() + "/histogram.png";
			try {
				ChartUtilities.saveChartAsPNG(new File(imgFilename), hc, 700, 500);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			result.addOutputItem(new Item("histogram_image","histogram.png","histogram image (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","HISTOGRAM_IMAGE"),new HashMap<String,String>(),"histogram image"));
		}

		if (bp != null){
			imgFilename = this.getOutdir() + "/boxplot.png";
			try {
				ChartUtilities.saveChartAsPNG(new File(imgFilename), bp, 700, 500);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			result.addOutputItem(new Item("boxplot_image","boxplot.png","boxplot image (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","BOXPLOT_IMAGE"),new HashMap<String,String>(),"boxplot image"));
		}


		// Save histogram
		progress++;
		updateJobStatus(""+(progress*100/finalProgress),"saving graph");

	}

	private void doTree() {
		try {
			jobStatus.addStatusMessage("" + ("40"), "preparing execution");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		updateJobStatus("40", "preparing execution");
		NewickParser nwParser = new NewickParser();
		MultipleTree tree;
		try {
			tree = nwParser.parse(new File(commandLine.getOptionValue("datalist")));
			String imgFilename = this.getOutdir() + "/tree.png";
			ClusteringUtils.saveImageTree(tree, "Tree",  imgFilename, true, false);
			try {
				jobStatus.addStatusMessage("" + ("80"), "saving");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}

			updateJobStatus("80", "saving");
			result.addOutputItem(new Item("tree_image","tree.png","tree image (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","TREE_IMAGE"),new HashMap<String,String>(),"tree image"));
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidFileFormatException e) {
			e.printStackTrace();
		}

	}

}


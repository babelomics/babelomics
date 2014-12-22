package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.methods.expression.differential.maSigPro;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class TimeSeries extends BabelomicsTool {

	public TimeSeries() {
		initOptions();
	}

	public void initOptions() {
		options.addOption(OptionFactory.createOption("test", "test: masigpro", false));
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("contin-class", "class variable"));
		options.addOption(OptionFactory.createOption("series-class", "class class variable"));
		options.addOption(OptionFactory.createOption("degree", "Polynomial degree", false));
		options.addOption(OptionFactory.createOption("q-value", "Q-value", false));
		options.addOption(OptionFactory.createOption("correction", "Multiple testing adjustment", false));
		options.addOption(OptionFactory.createOption("significance-level", "Significance level for model variable", false));
		options.addOption(OptionFactory.createOption("clustering-method", "Clustering method", false));
		options.addOption(OptionFactory.createOption("k-value", "Number of clusters (k-value) for K-means clustering", false));
	}


	@Override
	public void execute() {
		String test = commandLine.getOptionValue("test", "masigpro");
		executeMaSigPro();
	}

	public void executeMaSigPro() {

		// reading dataset
		//
		Dataset dataset = initDataset(new File(commandLine.getOptionValue("dataset")));

		String continClass = commandLine.getOptionValue("contin-class", null);
		String seriesClass = commandLine.getOptionValue("series-class", null);
		int degree = Integer.parseInt(commandLine.getOptionValue("degree", "2"));
		double q = Double.parseDouble(commandLine.getOptionValue("q-value", "0.05"));
		String correction = commandLine.getOptionValue("correction", "BH");
		double alfa = Double.parseDouble(commandLine.getOptionValue("significance-level", "0.05"));
		String clustering = commandLine.getOptionValue("clustering-method", "hclust");
		int kvalue = Integer.parseInt(commandLine.getOptionValue("k-value", "9"));

		List<String> continVars = dataset.getVariables().getByName(continClass).getValues();
		List<String> seriesVars = dataset.getVariables().getByName(seriesClass).getValues();

		// masigpro
		//
		updateJobStatus("40", "computing maSigPro test");
		List<String> lines = new ArrayList<String>();
		lines.add("#names\t" + ListUtils.toString(dataset.getSampleNames(), "\t"));
		lines.add("#contin\t" + ListUtils.toString(continVars, "\t"));
		lines.add("#series\t" + ListUtils.toString(seriesVars, "\t"));
		for(int i=0 ; i<dataset.getRowDimension() ; i++) {
			lines.add(dataset.getFeatureNames().get(i) + "\t" + ListUtils.toString(ArrayUtils.toStringList(dataset.getDoubleMatrix().getRow(i)), "\t"));
		}
		File inputFile = new File(outdir + "/input.maSigPro.txt");
		try {
			IOUtils.write(inputFile, lines);
		} catch (IOException e) {
			abort("ioexception_execute_masigpro", "error writting intermediate file", e.toString(), StringUtils.getStackTrace(e));
		}

		maSigPro masigpro = new maSigPro(babelomicsHomePath + "/bin/masigpro/masigpro.R");
		masigpro.setInputFilename(inputFile.getAbsolutePath());
		masigpro.setOutdir(outdir);

		masigpro.compute(degree, q, correction, alfa, clustering, kvalue);

		// saving data
		//
		updateJobStatus("80", "saving results");
		File outDirFile = new File(outdir);

		File outFile = null;
		File[] outFiles = null;

		// general info
		outFiles = FileUtils.listFiles(outDirFile, ".*summary.*");
		for (File f: outFiles) {
			File listFile = new File(f.getAbsolutePath() + ".id.list.txt");
			try {
				IOUtils.write(listFile, IOUtils.column(f, 0));
				if (listFile.exists()) {
					File redirectionFile = new File(f.getAbsoluteFile() + ".genoome.fatigo.redirection");		
					DiffExpressionUtils.createFatiGoRedirectionFile(redirectionFile, listFile);
					if ( redirectionFile.exists() ) {
						String tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiGO tool...)";
						result.addOutputItem(new Item("summaryfile", f.getName(), "Significant genes for '" + getCleanName(f) + "'", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "maSigPro output.List of significant genes"));
					}						
				} else {
					result.addOutputItem(new Item("summaryfile", f.getName(), "Significant genes for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "maSigPro output.List of significant genes"));
				}
			} catch (Exception e) {
				result.addOutputItem(new Item("summaryfile", f.getName(), "Significant genes for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "maSigPro output.List of significant genes"));				
			}
		}
		if ( (outFile = new File(outdir + "/pvalues.txt")).exists() ) {
			result.addOutputItem(new Item("pvaluesfile", outFile.getName(), "p-values file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.p-values and adjusted p-values of global model for all genes"));						
		}
		if ( (outFile = new File(outdir + "/influ_info.txt")).exists() ) {
			result.addOutputItem(new Item("influinfofile", outFile.getName(), "Influence data file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.Influence data (genes with influential values, possible outliers)"));						
		}
		//		if ( (outFile = new File(outdir + "/influ_data.png")).exists() ) {
		//			result.addOutputItem(new Item("infludataimg", outFile.getName(), "Influence data plot", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.Influence data (genes with influential values, possible outliers)"));						
		//		}

		// model data
		outFiles = FileUtils.listFiles(outDirFile, ".*coefficients\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigcoeffile", f.getName(), "Significant coefficients for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.Significant coefficients"));			
		}	
		outFiles = FileUtils.listFiles(outDirFile, ".*sig\\.pvalues\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigpvaluefile", f.getName(), "Significant p-values for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.Significant p-values"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*profiles\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigprofilefile", f.getName(), "Significant profiles for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.Significant profiles"));			
		}


		// visualization
		outFiles = FileUtils.listFiles(outDirFile, ".*heatmap\\.png");
		for (File f: outFiles) {
			result.addOutputItem(new Item("heatmapimg", f.getName(),  "'" + getCleanName(f) + "' heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.Heatmap plot"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*Groups\\.png");
		for (File f: outFiles) {
			result.addOutputItem(new Item("groupimg", f.getName(), "'" + getCleanName(f) + "' group", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "maSigPro output.Groups plot"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*Profiles\\.png");
		for (File f: outFiles) {
			result.addOutputItem(new Item("profileimg", f.getName(), "'" + getCleanName(f) + "' profile", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Additional results.Profiles plot"));			
		}				
	}


	private String getCleanName(File file) {
		String res = file.getName().replace(".txt", "").replace(".png", "");
		res = res.replace("groups_", "").replace("summary", "").replace("_sig.profiles", "").replace("_sig.pvalues", "").replace("_coefficients", "");
		res = res.replace("_Groups", "").replace("_heatmap", "").replace("_Profiles", "").replace("_sig.pvalues", "_heatmap");
		res = res.replace("_", " ");
		return res;
	}

}

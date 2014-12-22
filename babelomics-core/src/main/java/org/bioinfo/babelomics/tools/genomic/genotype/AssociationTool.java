package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.math.util.MathUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class AssociationTool extends GenotypeAnalysisTool {

	private static final String PLINK_OUTFILE_NAME = "plink";
	
	public AssociationTool() {

	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("test", "Valid values: assoc, fisher, linear, logistic, tdt", true, true));
		options.addOption(OptionFactory.createOption("maf", "Minor allele frequency, default value: 0.02", false));
		options.addOption(OptionFactory.createOption("pvalue-cutoff", "P-value cutoff, default value: 0.05", false));
		options.addOption(OptionFactory.createOption("log", "Wether to apply Odd ratio logarithm to the file result", false, false));
	}

	@Override
	public void execute() {
		logger.debug("executing association test");
		try {
			// update status
			jobStatus.addStatusMessage("10", "Parsing parameters");
			
			parseGenotypeCommonOptions();
			// specific options
			String test = commandLine.getOptionValue("test", "assoc");
			double maf = Double.parseDouble(commandLine.getOptionValue("maf", "0.02"));
			double pvalueCutoff = Double.parseDouble(commandLine.getOptionValue("pvalue-cutoff", "0.05"));

			// prepare the GenotypeAnalysis object for execution
			genotypeAnalysis = new GenotypeAnalysis(pedFilePath, mapFilePath);
			genotypeAnalysis.setPlinkPath(plinkPath);
			genotypeAnalysis.setOutdir(outdir);
			
			logger.debug("executing: "+plinkPath+" --ped "+pedFilePath+" --map "+mapFilePath+" --out "+outdir+"/"+PLINK_OUTFILE_NAME+" --maf "+maf + "--"+test);
			jobStatus.addStatusMessage("40", "Executing association test");
			
			genotypeAnalysis.association(test, maf);
			
			jobStatus.addStatusMessage("90", "Saving association results");
			// check if everything has gone well
			if(test != null) {
				saveResults(test, PLINK_OUTFILE_NAME);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
			result.addOutputItem(new Item("ioexcepcion_error", e.toString(), "An error occured", Item.TYPE.MESSAGE, Arrays.asList("ERROR"), new HashMap<String, String>(2), "Error message"));
		} catch (Exception e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
			result.addOutputItem(new Item("excepcion_error", e.toString(), "An error occured", Item.TYPE.MESSAGE, Arrays.asList("ERROR"), new HashMap<String, String>(2), "Error message"));
		} 				
	}

	private void saveResults(String test, String filename) throws IOException, InvalidIndexException {
		FeatureData featureData;
		DataFrame dataFrame = new DataFrame();
		String filePath = outdir+"/"+filename;
		if(test.equalsIgnoreCase("assoc")) {
			filePath += ".assoc";
			FileUtils.checkFile(filePath);
			result.addOutputItem(new Item(filename+".assoc_file", filename+".assoc", "Association result file form PLINK (test: chi-square)", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: chi-square"));
			result.addOutputItem(new Item(filename+".hh_file", filename+".hh", "List of heterozygous haploid genotypes (SNPs/individuals)", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: chi-square"));
			result.addOutputItem(new Item(filename+".log_file", filename+".log", "Log file from PLINK", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: chi-square"));
			// create dataframe
			dataFrame.addColumn("dbsnp", IOUtils.column(filePath, 1, "\\s+"));
			dataFrame.addColumn("chromosome", IOUtils.column(filePath, 0, "\\s+"));
			dataFrame.addColumn("position", IOUtils.column(filePath, 2, "\\s+"));
			List<String> pvalues = IOUtils.column(filePath, 7, "\\s+");
			double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
			minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
			dataFrame.addColumn("p_values", pvalues);
			dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
			dataFrame.addColumn("odd_ratio", IOUtils.column(filePath, 8, "\\s+"));
			dataFrame.removeRow(0);
		}else {
			if(test.equalsIgnoreCase("fisher")) {
				filePath += ".assoc."+test;
				FileUtils.checkFile(filePath);
				result.addOutputItem(new Item(filename+".assoc."+test+"_file", filename+".assoc."+test, "Association result file form PLINK (test: "+test+")", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: "+test));
				// create dataframe
				dataFrame.addColumn("dbsnp", IOUtils.column(filePath, 1, "\\s+"));
				dataFrame.addColumn("chromosome", IOUtils.column(filePath, 0, "\\s+"));
				dataFrame.addColumn("position", IOUtils.column(filePath, 2, "\\s+"));
				List<String> pvalues = IOUtils.column(filePath, 7, "\\s+");
				double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
				minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
				dataFrame.addColumn("p_values", pvalues);
				dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
				dataFrame.addColumn("odd_ratio", IOUtils.column(filePath, 8, "\\s+"));
				dataFrame.removeRow(0);
			}
			if(test.equalsIgnoreCase("linear") || test.equalsIgnoreCase("logistic")) {
				filePath += ".assoc.logistic";
				FileUtils.checkFile(filePath);
				result.addOutputItem(new Item(filename+".assoc.logistic_file", filename+".assoc.logistic", "Association result file form PLINK (test: "+test+")", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: "+test));
			}
			if(test.equalsIgnoreCase("tdt")) {
				FileUtils.checkFile(outdir+"/"+filename+".tdt");
				result.addOutputItem(new Item(filename+"."+test+"_file", filename+"."+test, "Association result file form PLINK (test: "+test+")", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: "+test));
			}
			result.addOutputItem(new Item(filename+".hh_file", filename+".hh", "List of heterozygous haploid genotypes (SNPs/individuals)", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: "+test));
			result.addOutputItem(new Item(filename+".log_file", filename+".log", "Log file from PLINK", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test: "+test));	
		}
		
		featureData = new FeatureData(dataFrame);
		featureData.save(new File(outdir+"/plink.featdata"));
	}
	
	private void saveFeatureData() {
		
	}
	
}

package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.IOException;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.tool.OptionFactory;

public class StratificationTool extends GenotypeAnalysisTool {
	
	public StratificationTool() {
		
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("test", "To perform complete linkage clustering of individuals on the basis of autosomal genome-wide SNP data, values: ibs-cluster", true, false));
		options.addOption(OptionFactory.createOption("ibs-test", "To test whether or not there are group differences in this metric, with respect to a binary phenotype", false, false));
		options.addOption(OptionFactory.createOption("ppc", "Pairwise population concordance (PPC). To only merge clusters that do not contain individuals differing at a certain p-value, by default 0.0001", false, true));
		options.addOption(OptionFactory.createOption("mc", "Maximum cluster size, not applied by default", false, true));
		options.addOption(OptionFactory.createOption("maf", "Minor allele frequency, default value: 0.02", false));
	}

	@Override
	public void execute() {
		logger.debug("executing Stratification analysis");
		logger.debug("executing association test");
		try {
			// update status
			jobStatus.addStatusMessage("10", "Reading options");
			
			parseGenotypeCommonOptions();
			// specific options
			String test = commandLine.getOptionValue("test", "ibs-cluster");
			double ppc = Double.parseDouble(commandLine.getOptionValue("ppc", "0.0001"));
			int mc = Integer.parseInt(commandLine.getOptionValue("mc", "0"));
			double maf = Double.parseDouble(commandLine.getOptionValue("maf", "0.02"));

			// prepare the GenotypeAnalysis object for execution
			genotypeAnalysis = new GenotypeAnalysis(pedFilePath, mapFilePath);
			genotypeAnalysis.setPlinkPath(plinkPath);
			genotypeAnalysis.setOutdir(outdir);
			
			logger.debug("executing: "+plinkPath+" --ped "+pedFilePath+" --map "+mapFilePath+" --out "+outdir+"/plink --maf "+maf + "--"+test+" --ppc "+ppc+" --mc "+mc);
			jobStatus.addStatusMessage("30", "Executing stratification job");
			
			genotypeAnalysis.stratification(test, commandLine.hasOption("ibs-test"), ppc, mc, maf);
			
		} catch (IOException e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} catch (Exception e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} 				
	}
}

package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;
import java.io.IOException;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.tool.OptionFactory;

public abstract class GenotypeAnalysisTool extends BabelomicsTool {

	protected GenotypeAnalysis genotypeAnalysis;
	
	protected String pedFilePath;
	protected String mapFilePath;
	protected String plinkPath;
	
	public GenotypeAnalysisTool() {
		genotypeAnalysis = new GenotypeAnalysis();
		initGenotypeCommonsOptions();
	}

	public void initGenotypeCommonsOptions() {
		// the typical ped and map files
		options.addOption(OptionFactory.createOption("ped-file", "PED file path", false));
		options.addOption(OptionFactory.createOption("map-file", "MAP file path", false));
		options.addOption(OptionFactory.createOption("zip-file", "ZIP file containing PED and MAP files", false));
		
		options.addOption(OptionFactory.createOption("plink", "PLINK file path", false));
	}

	protected void parseGenotypeCommonOptions() throws IOException {
		if(commandLine.hasOption("zip-file")) {
			// update status
			jobStatus.addStatusMessage("20", "Parsing parameters");
			
			FileUtils.checkFile(commandLine.getOptionValue("zip-file"));
			FileUtils.unzipFiles(commandLine.getOptionValue("zip-file"), outdir);
			// we expect just 1 PED file, if more than 1 the first one is selected
			File[] files = FileUtils.listFiles(new File(outdir), ".*.ped", true);
			if(files != null && files.length > 0) {
				pedFilePath = files[0].getAbsolutePath();
				logger.debug("PED file: "+pedFilePath);
			}
			// we expect just 1 MAP file, if more than 1 the first one is selected			
			files = FileUtils.listFiles(new File(outdir), ".*.map", true);
			if(files != null && files.length > 0) {
				mapFilePath = files[0].getAbsolutePath();
				logger.debug("MAP file: "+mapFilePath);
			}
		}
		if(commandLine.hasOption("ped-file")) {
			pedFilePath = commandLine.getOptionValue("ped-file");
			logger.debug("PED file: "+pedFilePath);
		}
		if(commandLine.hasOption("map-file")) {
			mapFilePath = commandLine.getOptionValue("map-file");
			logger.debug("MAP file: "+mapFilePath);
		}
		plinkPath = commandLine.getOptionValue("plink", babelomicsHomePath+"/bin/genotype/plink64");
		logger.debug("plink binary: "+plinkPath);
		
		// some obligatory checks
		FileUtils.checkFile(pedFilePath);
		FileUtils.checkFile(mapFilePath);
		FileUtils.checkFile(plinkPath);
	}
	
}

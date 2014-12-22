package org.bioinfo.babelomics.tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.GenericBioTool;
import org.bioinfo.tool.OptionFactory;

public abstract class BabelomicsTool extends GenericBioTool {
	
	protected String species;
	protected String babelomicsHomePath;
	
	
	public BabelomicsTool() {
//		babelomicsHomePath = System.getenv("BABELOMICS_HOME");
//		
//		try {
//			this.appendConfig(new File(babelomicsHomePath + "/conf/babelomics.properties"));
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		
		initCommonsOptions();
		initOptions();
	}
	
	public abstract void initOptions();
	
	private void initCommonsOptions() {
		getOptions().addOption(OptionFactory.createOption("tool", "to", "tool name", true));
		getOptions().addOption(OptionFactory.createOption("home", "home", "babelomics home path", true));
		getOptions().addOption(OptionFactory.createOption("species", "The specie of the ids", false));
	}

	
	/* (non-Javadoc)
	 * @see org.bioinfo.tool.GenericBioTool#parse(java.lang.String[])
	 */
	@Override
	public void parse(String[] args) throws ParseException, IOException {
		parse(args, false);
	}
	
	
	/* (non-Javadoc)
	 * @see org.bioinfo.tool.GenericBioTool#parse(java.lang.String[], boolean)
	 */
	@Override
	public void parse(String[] args, boolean stopAtNoOption) throws ParseException, IOException {
		super.parse(args, stopAtNoOption);
		// must be in commandLine, just in case we initialize...
		this.toolName = commandLine.getOptionValue("tool", "");
		this.species = commandLine.getOptionValue("species", "unknown");
		
		// invoke getCanonicalPath() to remove the '.' and '..' of the paths
		babelomicsHomePath = new File(commandLine.getOptionValue("home", "")).getCanonicalPath();
		logger.debug("BABELOMICS_HOME path: "+babelomicsHomePath);
		try {
			this.appendConfig(new File(babelomicsHomePath + "/conf/babelomics.properties"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		// default report in PDF is requiered in Babelomics
//		if(!commandLine.hasOption("report")) {
//			logger.debug("Adding PDF report option");
//			this.report = "pdf";
//		}
	}
	
	@Deprecated
	public void checkFile(String filePath) throws InvalidParameterException, IOException {
		checkFile(new File(filePath));
	}
	
	@Deprecated
	public void checkFile(File file) throws InvalidParameterException, IOException {
		if(file == null) {
			throw new InvalidParameterException("File is null");
		}
		if(!file.exists()) {
			throw new IOException("File '" + file.getAbsolutePath() + "' does not exist");
		}
		if(file.isDirectory()) {
			throw new IOException("File '" + file.getAbsolutePath() + "' is a directory");
		}
		if(!file.canRead()) {
			throw new IOException("File '" + file.getAbsolutePath() + "' cannot be read");	
		}
	}
	
	
	public Dataset initDataset(File file) {
		updateJobStatus("20", "reading dataset");
			
		try {
			// reading data
			//
			return new Dataset(file, true);
		} catch (Exception e) {
			abort("filenotfoundexception_initdataset_babelomicstool", "error reading dataset", e.toString(), StringUtils.getStackTrace(e));
		}
		return null;
	}
	

	public void updateJobStatus(String progress, String message) {
		logger.debug(message + "...\n");
		try {
			jobStatus.addStatusMessage(progress, message);
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_updatejobstatus_babelomicstools", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
	}
	
	
	
	/**
	 * @param species the species to set
	 */
	public void setSpecies(String species) {
		this.species = species;
	}

	/**
	 * @return the species
	 */
	public String getSpecies() {
		return species;
	}

	
	
}

package org.bioinfo.babelomics.tools.functional;

import java.io.IOException;

import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.tool.OptionFactory;

public class GesbapTool extends FunctionalProfilingTool {

	
	
	public GesbapTool() {
		
	}
	
	/* (non-Javadoc)
	 * @see org.bioinfo.babelomics.tools.functional.FunctionalProfilingTool#initOptions()
	 */
	@Override
	public void initOptions() {
		super.initOptions();
		OptionGroup inputaData = new OptionGroup();
		inputaData.setRequired(true);
		inputaData.addOption(OptionFactory.createOption("snp-file", "Tab file with two columns: SNP id and statistic", false, true));
		inputaData.addOption(OptionFactory.createOption("plink-assoc-file", "Plink association output (assoc, fisher, tdt, linera, logistic)", false, true));
//		inputaData.addOption(OptionFactory.createOption("ped-file", "PED file path", false));
//		inputaData.addOption(OptionFactory.createOption("map-file", "MAP file path", false));
//		inputaData.addOption(OptionFactory.createOption("zip-file", "ZIP file containing PED and MAP files", false));
		options.addOptionGroup(inputaData);
		
		options.addOption(OptionFactory.createOption("method", "Gene set analysis method values: fatiscan, logistic, default value fatiscan", false, true));
		options.addOption(OptionFactory.createOption("partitions", "Set the number of partitions, by default '30'", false, true));
		options.addOption(OptionFactory.createOption("output-format", "Values: short (just most significant partition) or long (term results for all partitions), by deafult 'short'", false, true));
	}

	@Override
	protected void execute() {
		try {
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");
			
			// parse options form parent CLI
			prepare();
			
			// parse command line specific options
			if(options.hasOption("snp-file")) {
				
			}else {
				if(options.hasOption("plink-assoc-file")) {
					
				}else {
					printError("gesbaptool_execute", "Error in GeSBAP tool", "No valid input data provided");
				}
			}
			
			
		} catch (IOException e) {
			e.printStackTrace();
			printError("gesbaptool_execute", "Error in GeSBAP tool", "An error ocurred during GeSBAP execution", e);
		} catch (ParseException e) {
			e.printStackTrace();
			printError("gesbaptool_execute", "Error in GeSBAP tool", "An error ocurred during GeSBAP execution", e);
		} catch (InvalidIndexException e) {
			e.printStackTrace();
			printError("gesbaptool_execute", "Error in GeSBAP tool", "An error ocurred during GeSBAP execution", e);
		} catch (Exception e) {
			e.printStackTrace();
			printError("gesbaptool_execute", "Error in GeSBAP tool", "An error ocurred during GeSBAP execution", e);
		}
	}

	private void createSnpFileFromPlinkAssocOutput() {
		
	}
}

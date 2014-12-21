package org.bioinfo.babelomics.tools.genomic.copynumber;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.methods.genomic.copynumber.CopyNumberAnalysisExecutor;
import org.bioinfo.babelomics.methods.genomic.copynumber.CopyNumberUtils;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class CopyNumberAnalysis extends BabelomicsTool {

	@Override
	public void initOptions() {
		// second input: the typical ped and map files
		options.addOption(OptionFactory.createOption("normalized-file", "Normalized file", true, true));

		options.addOption(OptionFactory.createOption("segmentation-method", "Segmentation method: dnacopy or glad", true, true));
		options.addOption(OptionFactory.createOption("alpha", "Significance levels for DNAcopy test to accept change-points (only for DNAcopy segmentation method)", false, true));

		options.addOption(OptionFactory.createOption("cgh-mcr", "Minimal common region", false, false));
		options.addOption(OptionFactory.createOption("gap-allowed", "Gap allowed (only for CGH-MCR)", false, true));
		options.addOption(OptionFactory.createOption("altered-low", "Gap altered low (only for CGH-MCR)", false, true));
		options.addOption(OptionFactory.createOption("altered-high", "Gap altered high (only for CGH-MCR)", false, true));
		options.addOption(OptionFactory.createOption("recurrence", "Recurrence (onlyfor CGH-MCR)", false, true));
	}

	@Override
	protected void execute() {

		// reading data
		//
		updateJobStatus("20", "reading input data");

		String normalizedFilename = commandLine.getOptionValue("normalized-file", null);

		String segmentation = commandLine.getOptionValue("segmentation-method", null);
		double alpha = Double.parseDouble(commandLine.getOptionValue("alpha", "0.01"));
		boolean cghMcr = commandLine.hasOption("cgh-mcr");
		int gapAllowed = Integer.parseInt(commandLine.getOptionValue("gap-allowed", "500"));
		double alteredLow = Double.parseDouble(commandLine.getOptionValue("altered-low", "0.2"));
		double alteredHigh = Double.parseDouble(commandLine.getOptionValue("altered-high", "0.8"));
		int recurrence = Integer.parseInt(commandLine.getOptionValue("recurrence", "50"));

		if ( normalizedFilename == null ) {
			abort("normalizedfilemissing_execute_copynumberanalysis", "input normalized file missing", "input normalized file missing", "input normalized file missing");
		}

		if ( ! new File(normalizedFilename).exists() ) {
			abort("normalizedfilenotexist_execute_copynumberanalysis", "input normalized file does not exist", "input normalized file  does not exist", "input normalized file  does not exist");
		}

//		Dataset dataset = null;
//		try {
//			dataset = new Dataset(new File(normalizedFilename), true);
//		} catch (IOException e) {
//			abort("normalizedfilenotvalid_execute_copynumberanalysis", "input dataset is not valid", "input dataset is not valid", e.getMessage());
//		}
//		
//		System.out.println("samples names = " + ListUtils.toString(dataset.getSampleNames(), ", "));
//		System.out.println("sample attributes = " + ListUtils.toString(MapUtils.getKeys(dataset.getAttributes()), ", "));
//		System.out.println("features column names = " + ListUtils.toString(dataset.getFeatureData().getDataFrame().getColumnNames(), ", "));
//		System.out.println("feature attributes = " + ListUtils.toString(MapUtils.getKeys(dataset.getFeatureData().getAttributes()), ", "));
//		
//		try {
//			dataset.getFeatureData().save(new File("/tmp/features.txt"));
//			dataset.save("/tmp/dataset.test");
//		} catch (IOException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
		
		String inputCghFileName = outdir + "/copynumber.in";
		try {
			CopyNumberUtils.saveAsCGHInput(new Dataset(new File(normalizedFilename), true), inputCghFileName);
		} catch (Exception e) {
			abort("normalizedfilenotvalid_execute_copynumberanalysis", "input dataset is not valid", "input dataset is not valid", e.getMessage());
		}
		if ( ! (new File(inputCghFileName).exists()) ) {
			abort("datasetnotvalid_execute_copynumberanalysis", "cannot convert dataset to input copy number file", "cannot convert dataset to input copy number file", "cannot convert dataset to input copy number file");			
		}
		
		CopyNumberAnalysisExecutor cpExecutor = null;

		if ( "dnacopy".equalsIgnoreCase(segmentation) ) {
			cpExecutor = new CopyNumberAnalysisExecutor(babelomicsHomePath + "/bin/copynumber/DNAcopy.r");
		} else if ( "glad".equalsIgnoreCase(segmentation) ) {
			cpExecutor = new CopyNumberAnalysisExecutor(babelomicsHomePath + "/bin/copynumber/GLAD.r");
		} else {
			abort("unknownanalysis_execute_copynumberanalysis", "copy number analysis unknown", "copy number analysis unknown, valid values are dnacopy and glad", "copy number analysis unknown, valid values are dnacopy and glad");
		}

		File cghFile = new File(outdir + "/cgh.txt");;
		File segmentedFile = new File(outdir + "/segmented.txt");
		
		//cpExecutor.setNormalizedFile(normalizedFilename);
		cpExecutor.setNormalizedFile(inputCghFileName);
		cpExecutor.setSegmentatedFilename(segmentedFile.getAbsolutePath());
		cpExecutor.setCghFilename(cghFile.getAbsolutePath());
		cpExecutor.setCghMcrBinPath(null);
		
		if ( cghMcr ) {
			cpExecutor.setCghMcrBinPath(babelomicsHomePath + "/bin/copynumber/cghMCR.r");
			cpExecutor.setAlteredLow(alteredLow);
			cpExecutor.setAlteredHigh(alteredHigh);
			cpExecutor.setGapAllowed(gapAllowed);
			cpExecutor.setRecurrence(recurrence);
		}

		// executing segmentation
		//
		updateJobStatus("40", "executing segmentation");

		try {
			cpExecutor.run();			
		} catch (InvalidParameterException e) {
			printError("invalidparameteexception_execute_copynumberanalysis", "error executing segmentation", e.getMessage(), e);
		} catch (IOException e) {
			printError("ioexception_execute_copynumberanalysis", "error executing segmentation", e.getMessage(), e);
		}

		// saving results
		//
		updateJobStatus("90", "saving results");

		if ( segmentedFile.exists() ) {
			result.addOutputItem(new Item("segmentedfile", segmentedFile.getName(), "Segmented " + segmentation.toUpperCase() + " file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), segmentation.toUpperCase() + " segmentation"));						
		} else {
			printError("segmentedfileerror_execute_copynumberanalyis", "error creating segmented file", "error creating segmented file");
		}
		if ( cghMcr ) {
			if ( cghFile.exists() ) {			
				List<String> lines = null;
				try {
					lines = IOUtils.grep(cghFile, "#ERROR.*");
					if ( lines != null && lines.size() > 0 ) {
						String error = lines.get(0).replace("#ERROR\t", ""); 
						printError("cghfileerror_execute_copynumberanalyis", "Error", error);			
					} else {
						result.addOutputItem(new Item("cghfile", cghFile.getName(), "CGH file for " + segmentation.toUpperCase() + " segmentation", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "CGH file"));						
					}				
				} catch (IOException e) {
					printError("cghfileerror_execute_copynumberanalyis", "Error", "error creating cgh file");			
				}
			} else {
				printError("cghfileerror_execute_copynumberanalyis", "Error", "error creating cgh file");			
			}
		}
	}
}

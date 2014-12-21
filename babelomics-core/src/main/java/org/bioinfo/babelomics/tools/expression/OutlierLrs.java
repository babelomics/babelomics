package org.bioinfo.babelomics.tools.expression;

import java.io.File;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class OutlierLrs extends BabelomicsTool {


	public OutlierLrs() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("dataset", "the data"));
		getOptions().addOption(OptionFactory.createOption("low-end", "Low-end of permutations", false));
		getOptions().addOption(OptionFactory.createOption("up-end", "Up-end variable", false));		
	}

	@Override
	public void execute() {
		try {
			Dataset expresionDataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
			
			String permutation = commandLine.getOptionValue("low-end");
			String percentile = commandLine.getOptionValue("up-end");
			executeOutLRS(expresionDataset,permutation,percentile);
		} catch (Exception e) {
			logger.error("Error opening the dataset", e.toString());
		} 
	}


	private void executeOutLRS(Dataset dataset, String permutation,String percentile ) {
		logger.info("executing outLRS, not implemented yet");
	}
	
}

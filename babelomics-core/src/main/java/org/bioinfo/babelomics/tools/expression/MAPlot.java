package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.microarray.ExpressionUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class MAPlot extends BabelomicsTool {

	public MAPlot() {
		initOptions();
	}

	@Override
	public void initOptions() {

		getOptions().addOption(OptionFactory.createOption("dataset", "Dataset (data matrix)", false));
		//		getOptions().addOption(OptionFactory.createOption("technology", "Techonology type: agilent, genepix", false));
		//		getOptions().addOption(OptionFactory.createOption("channels", "Number of channels: 1 or 2", false));

		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	@Override
	public void execute() {
		outdir = outdir + "/";
		String datasetFileName = commandLine.getOptionValue("dataset", null);

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( datasetFileName == null ) {
			abort("missingdata_execute_maplot", "missing input data", "missing input data", "missing input data");						
		}

		File fFile = new File(datasetFileName);

		try {
			jobStatus.addStatusMessage("10", "reading dataset");

			Dataset ds = new Dataset(fFile);
			
			if ( ds == null || ds.getFeatureNames() == null || ds.getFeatureNames().size() == 0 ) {
				printError("invaliddataset_execute_maplot", "Invalid dataset", "Dataset " + fFile.getName() + " is not valid", "");
			}
			
			File inFile = new File(outdir + "infile.txt");
			
			List<String> lines = new ArrayList<String>();
			for(int row=0 ; row<ds.getRowDimension() ; row++) {
				lines.add(ds.getFeatureNames().get(row) + "\t" + ArrayUtils.toString(ds.getDoubleMatrix().getRow(row), "\t"));
			}
			IOUtils.write(inFile, lines);

			// ma plot
			//
			jobStatus.addStatusMessage("50", "creating MA plot");

			String maPlotBinPath = babelomicsHomePath + "/bin/plots/plotMA_from_single_matrix.r";


			// ma plots
			//
			if ( inFile.exists() ) {
				ExpressionUtils.createMAPlot(maPlotBinPath, inFile.getAbsolutePath(), "MA_", true, "ma_plot.Rout", outdir);
				addOutputItems(outdir, "MA_", "png", "ma_plot", "MA plots");				
			}

		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_maplot", "File not found error", "Impossible to read input file: " + fFile.getName() , "");
		} catch (IOException e) {
			printError("ioexception_execute_maplot", "IO error", "Error ocurred when accessing input file: " + fFile.getName(), "");
		}

	}


	private void addOutputItems(String dir, String prefix, String suffix, String id, String groupName) {
		String sampleName;
		String pattern = prefix + ".*" + suffix;
		File [] files = FileUtils.listFiles(new File(dir), pattern);

		System.out.println("addOutputItems...");
		for (int i=0 ; i<files.length ; i++) {
			//System.out.println("file " + i + " -> " + files[i].getName());
			
			sampleName = files[i].getName().replace(prefix, "").replace(".png", "");
			result.addOutputItem(new Item(id + "_" + i, files[i].getName(), groupName + " for " + sampleName, TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), groupName));
		}
	}
}

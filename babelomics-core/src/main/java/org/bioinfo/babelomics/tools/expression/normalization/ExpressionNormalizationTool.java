package org.bioinfo.babelomics.tools.expression.normalization;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.io.file.compress.CompressFactory;
import org.bioinfo.io.file.compress.GenericCompressManager;
import org.bioinfo.microarray.AffymetrixExpressionUtils;
import org.bioinfo.microarray.ExpressionUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;

public class ExpressionNormalizationTool extends BabelomicsTool {

	int nbChannels = 0;
	String technology = "";
	List<String> rawFileNames = null;
	File tmpDir = null;
	String sampleInfoPath = "";

	public ExpressionNormalizationTool() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("compressed-file", "Compressed file containning the raw files", false));
		getOptions().addOption(OptionFactory.createOption("raw-dir", "Directory where the raw files are located", false));
		getOptions().addOption(OptionFactory.createOption("compressed-file-tags", "Data tags separated by commas, valid values are: microarray, expression, affymetrix, agilent, genepix, one-channel, two-channels", false));

		// affy normalization
		//
		getOptions().addOption(OptionFactory.createOption("no-cel-convert", "Disable to convert CEL files to GCOS text file format (by default, the conversion is enabled). Only for Affymetrix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("rma", "RMA analysis. Only for Affymetrix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("plier", "Plier analysis. Only for Affymetrix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("calls", "Present-absent calls analysis. Only for Affymetrix normalization.", false, false));


		// agilent or genepix normalization
		//
		getOptions().addOption(OptionFactory.createOption("sample-names", "Sample names. Only for agilent or genepix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("bg-correction", "Background correction: agilent, normexp, rma, half, subtract, minimum, movingmin, edwards, none. Default background correction: rma. Only for Agilent or GenePix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("wa-normalization", "Within arrays normalization: none, loess, printtiploess, median, composite, control, robustspline. Default within arrays normalization: loess. Only for two-colors Agilent or GenePix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("ba-normalization", "Between arrays normalization: none, quantiles, scale, vsn. Default between arrays normalization: scale. Only for Agilent or GenePix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("flags-no-fitted", "If this option is set then spots will not be used in the fitting of the parameters of the normalization steps. Only for Agilent or GenePix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("flags-as-missing", "If this option is set then spots will have a missing (NA) normalized value and A-value as well. Only for Agilent or GenePix normalization.", false, false));		
	}

	/**
	 * 
	 * @throws IOException
	 */
	public void prepare() throws IOException {
		outdir = outdir + "/";
		tmpDir = new File(outdir + "tmp");
		sampleInfoPath = tmpDir + "/sampleinfo.txt";
		String compressedFileName = commandLine.getOptionValue("compressed-file", null);
		String rawDirName = commandLine.getOptionValue("raw-dir", null);
		List<String> tags = StringUtils.toList(commandLine.getOptionValue("compressed-file-tags", ""), ",");


		System.out.println("----------> " + ListUtils.toString(tags, ","));

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( compressedFileName == null && rawDirName == null ) {
			abort("missingdata_execute_expressionnormalization", "Error", "Missing input data", "Missing input data");						
		}

		if ( tags == null || tags.size() == 0 ) {
			abort("unidentifieddata_execute_expressionnormalization", "Error", "Unidentified input data", "Unidentified input data");						
		}


		//		String technology = commandLine.getOptionValue("technology", null);
		//		String channels = commandLine.getOptionValue("channels", null);
		//		int nbChannels = 0;

		if ( tags.contains("one-channel") ) {
			nbChannels = 1;
		} else if ( tags.contains("two-channels") ) {
			nbChannels = 2;
		}

		if ( tags.contains("agilent") ) {
			technology = "agilent";
		} else if ( tags.contains("genepix") ) {
			technology = "genepix";
		} else if ( tags.contains("affymetrix") ) {
			technology = "affy";
			nbChannels = 1;
		}

		System.out.println("nbChannels = " + nbChannels + ", technology = " + technology);

		if ( technology == null ) {
			abort("missingtechnology_execute_expressionnormalization", "Error", "Missing technology tag", "Missing technology tag");						
		}

		if ( !technology.equalsIgnoreCase("agilent") && !technology.equalsIgnoreCase("genepix") && !technology.equalsIgnoreCase("affy") ) {
			abort("invalidtechnology_execute_expressionnormalization", "Error", "Invalid technology type " + technology + ". Valid values are: agilent, genepix", "Invalid technology type " + technology + ". Valid values are: agilent, genepix, affymetrix.");									
		}

		if ( nbChannels == 0 ) {
			abort("missingnbchannels_execute_expressionnormalization", "Error", "Missing number of channels", "Missing number of channels");						
		}


		jobStatus.addStatusMessage("10", "reading dataset");

		if ( compressedFileName != null ) {
			// getting raw files from compressed file
			//
			jobStatus.addStatusMessage("30", "decomprising dataset");

			System.out.println("input dataset = " + compressedFileName);
			System.out.println("tmp dir = " + tmpDir.getAbsolutePath());

			GenericCompressManager compresor = CompressFactory.getCompressManager(new File(compressedFileName));
			if (compresor != null) {
				try {
					rawFileNames = compresor.decompress(compressedFileName, tmpDir.getAbsolutePath());
				} catch (Exception e) {
					System.out.println("***** exception when decompressing input file, maybe it is not a zip-file like");
					abort("missingrawfiles_execute_expressionnormalization", "Error", "Missing raw files", "Missing raw files");
				}
			} else {
				tmpDir.mkdir();
				File outfile = new File(tmpDir.getAbsolutePath() + "/" + new File(compressedFileName).getName());

				System.out.println("***** compresor is null, copy " + compressedFileName + " ---> " + outfile.getAbsolutePath());

				InputStream inputStream = new FileInputStream(new File(compressedFileName));
				OutputStream outputStream = new FileOutputStream(outfile);

				long copiedBytes = FileUtils.copy(inputStream, outputStream);
				System.out.println("copied bytes = " + copiedBytes);

				inputStream.close();
				outputStream.close();

				//				FileSystem
				//				tmpDir.setWritable(true, false);
				//				tmpDir.mkdir();
				//				//outfile.setWritable(true);
				//				System.out.println("copy file to: " + outfile + ", is writable? " + outfile.canWrite());
				//				FileUtils.copy(new File(compressedFileName), outfile);
			}
		} else {
			// getting raw files from directory
			//
			tmpDir = new File(rawDirName); 
		}

		File[] rawFiles = FileUtils.listFiles(tmpDir, "affy".equalsIgnoreCase(technology) ? ".+.cel" : ".+", true);
		rawFileNames = new ArrayList<String>(rawFiles.length);
		for(File file: rawFiles) {
			if ( !"sampleinfo.txt".equalsIgnoreCase(file.getName())) {
				rawFileNames.add(file.getAbsolutePath());
			}
		}

		// sanity check
		//
		if ( rawFileNames == null || rawFileNames.size() == 0 ) {
			abort("missingrawfiles_execute_expressionnormalization", "Error", "Missing raw files", "Missing raw files");
		}		
	}

	/**
	 * 
	 */
	@Override
	protected void execute() {
		try {
			prepare();

			if ( "affy".equalsIgnoreCase(technology) ) {
				affyNormalization();
			} else {
				normalization();
			}
		} catch (Exception e) {
			printError("exception_execute_genepixexpression2cnormalization", e.toString(), e.getMessage(), e);
		}
	}

	private void normalization() throws IOException, InvalidIndexException {

		String sampleNames = commandLine.getOptionValue("sample-names", null);
		String bgCorrection = commandLine.getOptionValue("bg-correction", "minimum");
		String waNormalization = commandLine.getOptionValue("wa-normalization", "median");
		String baNormalization = commandLine.getOptionValue("ba-normalization", "none");
		boolean flagsNotFitted = commandLine.hasOption("flags-no-fitted");
		boolean flagsAsMissing = commandLine.hasOption("flags-as-missing");

		if (rawFileNames.size() == 1) {
			if (!"none".equalsIgnoreCase(baNormalization)) {
				printWarning("ba_normalization_warning", "Warnig", "Between arrays normalization has been set to 'none' since your dataset contains one single array");
			}
			baNormalization = "none";
		}
		
		String intensityPlotBinPath = babelomicsHomePath + "/bin/plots/plot_image_" + technology.toLowerCase() + ".r";

		// input parameters
		//
		result.addOutputItem(new Item("tech_input_param", technology + " (" + (nbChannels == 1 ? "one-channel" : "two-channels") + ")", "Technology", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("bg_input_param", bgCorrection, "Background correction", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		if ( nbChannels == 2 ) {
			result.addOutputItem(new Item("wa_input_param", waNormalization, "Within arrays normalization", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));			
		}
		result.addOutputItem(new Item("ba_input_param", baNormalization, "Between arrays normalization", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		if ( flagsNotFitted ) {
			result.addOutputItem(new Item("notfitted_input_param", "Flagged spots are not be used in the fitting of the parameters of the normalization steps", "Flag", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		}
		if ( flagsAsMissing ) {
			result.addOutputItem(new Item("missing_input_param", "Flagged spots have a missing (NA) normalized value and A-value as well", "Flag", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		}


		String readingScript = babelomicsHomePath + "/bin/normalizexp/" + (nbChannels == 1 ? "onecolor" : "twocolor");
		String normalizationScript = babelomicsHomePath + "/bin/normalizexp/" + (nbChannels == 1 ? "onecolor" : "twocolor");

		if ( "genepix".equalsIgnoreCase(technology) ) {
			readingScript += "_genepix_reading.r";
			normalizationScript += "_genepix_normalizing.r";
		} else if ( "agilent".equalsIgnoreCase(technology) ) {
			readingScript += "_agilent_reading.r";
			normalizationScript += "_agilent_normalizing.r";			
		} else {
			abort("invalidtechnology_execute_expressionnormalization", "Error", "Invalid technology type " + technology + ". Valid values are: agilent, genepix", "Invalid technology type " + technology + ". Valid values are: agilent, genepix, affymetrix");
		}

		// normalizing data
		//
		jobStatus.addStatusMessage("50", "normalizing data");
		String maPlotBinPath = "";

		if ( nbChannels == 1 ) {
			ExpressionUtils.OneColorNormalization(readingScript, normalizationScript, rawFileNames, (sampleNames != null ? StringUtils.toList(sampleNames, ","): getSamples(rawFileNames)), bgCorrection, baNormalization, flagsNotFitted, flagsAsMissing, outdir);
			maPlotBinPath = babelomicsHomePath + "/bin/plots/plotMA_from_single_matrix.r";
		} else if ( nbChannels == 2 ){
			ExpressionUtils.TwoColorsNormalization(readingScript, normalizationScript, rawFileNames, (sampleNames != null ? StringUtils.toList(sampleNames, ","): getSamples(rawFileNames)), bgCorrection, waNormalization, baNormalization, flagsNotFitted, flagsAsMissing, outdir);			
			maPlotBinPath = babelomicsHomePath + "/bin/plots/plotMA_from_MA_matrices.r";
		} else {
			abort("invalidchannels_execute_expressionnormalization", "Error", "Invalid number of channels (" + nbChannels + ")", "Invalid number of channels (" + nbChannels + ")");			
		}

		// saving normalization results
		//
		jobStatus.addStatusMessage("90", "saving normalization results");


		File file;
		if ( new File(outdir + "/" + ExpressionUtils.getNormalizedFileName()).exists() && 
				new File(outdir + "/" + ExpressionUtils.getFeatureDataFileName()).exists() ) {

			file = new File(outdir + "/normalized_dataset.txt"); 
			ExpressionUtils.createDataset(outdir + "/" + ExpressionUtils.getNormalizedFileName(), outdir + "/" + ExpressionUtils.getFeatureDataFileName(), 1, file.getAbsolutePath(), sampleInfoPath);

			if ( file.exists() ) {				
				String tags = "datamatrix,expression";
				result.addOutputItem(new Item("normalized", file.getName(), "Normalized dataset ", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Normalization output files"));

				File redirectionFile = new File(outdir + "/normalized.redirection");
				createPreprocessingRedirectionFile(redirectionFile, file);
				if ( redirectionFile.exists() ) {
					tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
					result.addOutputItem(new Item("normalized", file.getName(), "Normalized dataset ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Normalization output files"));
				}
				saveBoxPlot(file, false, "Box-plot", "boxplot", "Box-plots");				
			}

			file = new File(outdir + "/normalized_dataset.featdata"); 			
			if ( file.exists() ) {				
				result.addOutputItem(new Item("normalized", file.getName(), "Feature data ", TYPE.FILE, StringUtils.toList("idlist", ","), new HashMap<String, String>(2), "Normalization output files"));

				file = new File(outdir + "/normalized_dataset.gff3"); 
				//	DatasetUtils.dataset2Gff3(outdir + "/normalized_dataset.txt", file.getAbsolutePath());
				if ( file.exists() ) {
					result.addOutputItem(new Item("gff3", file.getName(), "GFF3 format", TYPE.FILE, StringUtils.toList("gff", ","), new HashMap<String, String>(2), "Normalization output files"));
				}
			}

		} else {
			abort("expressionnormalization", "Error", "Could not normalize your data", "Could not normalize your data");						
		}


		file = new File(outdir + "/" + ExpressionUtils.getaValuesFileName()); 
		if ( file.exists() ) {
			File aFile = new File(outdir + "/Avalues.txt");
			ExpressionUtils.updateAValuesFile(file.getAbsolutePath(), outdir + "/" + ExpressionUtils.getFeatureDataFileName(), 1, aFile.getAbsolutePath());
			result.addOutputItem(new Item("avalues", aFile.getName(), "A-values", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Normalization output files"));
		}


		// ma plots
		//
		if ( nbChannels == 1 ) {
			ExpressionUtils.createMAPlot(maPlotBinPath, outdir + "/" + ExpressionUtils.getNormalizedFileName(), "MA_", false, "ma_plot.Rout", outdir);
			addOutputItemImgs(outdir, "MA_", "png", "ma_plot", "MA plot", "MA plots");				
		} else if ( nbChannels == 2 ) {
			ExpressionUtils.createMAPlot(maPlotBinPath, outdir + "/" + ExpressionUtils.getNormalizedFileName(), outdir + "/" + ExpressionUtils.getaValuesFileName(), "MA_", false, "ma_plot.Rout", outdir);
			addOutputItemImgs(outdir, "MA_", "png", "ma_plot", "MA plot", "MA plots");				
		}

		// intensity images
		//
		if ( new File(outdir + "/" + ExpressionUtils.getNormalizedFileName()).exists() && 
				new File(outdir + "/" + ExpressionUtils.getFeatureDataFileName()).exists() ) {
			ExpressionUtils.createIntensityPlot(intensityPlotBinPath, outdir + "/" + ExpressionUtils.getNormalizedFileName(), outdir + "/" + ExpressionUtils.getFeatureDataFileName(), "norm_", false, "norm_intensity_plot.Rout", outdir);
			addOutputItemImgs(outdir, "norm_", "png", "norm_intensity_image", "Intensity image", "Intensity images");
		}

	}


	private void addOutputItemImgs(String dir, String prefix, String suffix, String id, String label, String groupName) {
		String sampleName;
		String pattern = prefix + ".*" + suffix;
		File [] files = FileUtils.listFiles(new File(dir), pattern);

		System.out.println("addOutputItems...");
		for (int i=0 ; i<files.length ; i++) {
			System.out.println("file " + i + " -> " + files[i].getName());

			sampleName = files[i].getName().replace(prefix, "").replace(".png", "");
			result.addOutputItem(new Item(id + "_" + i, files[i].getName(), label + " for " + sampleName, TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), groupName));
		}
	}

	/**
	 * 
	 * @throws IOException
	 * @throws InvalidParameterException
	 * @throws InvalidIndexException
	 */
	private void affyNormalization() throws IOException, InvalidParameterException, InvalidIndexException {

		String aptBinPath = babelomicsHomePath + "/bin/apt";

		boolean celConvert = !commandLine.hasOption("no-cel-convert");
		boolean rma = commandLine.hasOption("rma");
		boolean plier = commandLine.hasOption("plier");
		boolean calls = commandLine.hasOption("calls");

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( !rma && !plier && !calls ) {
			abort("missinganalysis_execute_affynormalization", "Error", "Missing analysis, valid values are rma, plier and calls", "Missing analysis, valid values are rma, plier and calls");			
		}

		// input parameters
		//
		result.addOutputItem(new Item("affy_input_param", "affymetrix (one-channel)", "Technology", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		if ( rma ) {
			result.addOutputItem(new Item("rma_input_param", "RMA", "Analysis", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		}
		if ( plier ) {
			result.addOutputItem(new Item("plier_input_param", "Plier", "Analysis", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		}
		if ( calls ) {
			result.addOutputItem(new Item("calls_input_param", "Present-absent calls", "Analysis", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		}	

		// creating the cel_files containning the cel files to normalize
		//
		File celFiles = new File(outdir + "/cel_files.txt");
		IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFileNames, "\n"));

		// converting to CEL text
		//
		if ( celConvert ) {			
			jobStatus.addStatusMessage("40", "converting CEL to GCOS text file format");

			File convertDir = new File(tmpDir.getAbsolutePath() + "/convert");
			convertDir.mkdir();
			AffymetrixExpressionUtils.aptCelConvert(aptBinPath + "/apt-cel-convert", celFiles.getAbsolutePath(), convertDir.getAbsolutePath());
			File[] rawFiles = FileUtils.listFiles(convertDir, ".+.CEL", true);
			rawFileNames = ArrayUtils.toStringList(rawFiles);
			//System.out.println("-----------> converting to gcos text file format, raw file names = " + ListUtils.toString(rawFileNames, ","));
			IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFileNames, "\n"));

			//			AffymetrixExpresionUtils.aptCelConvert(aptBinPath + "/apt-cel-convert", celFiles.getAbsolutePath(), tmpDir.getAbsolutePath());
			//
			//			File[] rawFiles = FileUtils.listFiles(tmpDir, ".+.CEL", true);
			//			rawFileNames = ArrayUtils.toStringList(rawFiles);
			//			System.out.println("-----------> converting to gcos text file format, raw file names = " + ListUtils.toString(rawFileNames, ","));
			//
			//			IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFileNames, "\n"));
		}		

		//Config config = new Config();
		config.append(new File(babelomicsHomePath + "/conf/apt.conf"));
		String chipName = getChipName(rawFileNames, config.getKeys());

		if ( chipName == null ) {
			abort("exception_execute_affynormalization", "Error", "Array type not supported", "Array type not supported");			
		}

		String infoStr = config.getProperty(chipName);
		Map<String, String> chipInfo = MapUtils.stringToMap(infoStr);
		chipInfo.put("name", chipName);

		System.out.println(" chip info = " + chipInfo.toString());


		String chipType = chipInfo.get("type");
		if ( chipType == null ) {
			abort("exception_execute_affynormalization", "Error", "Could not find out the chip type", "Could not find out the chip type");			
		}

		if ( !"3-prime".equalsIgnoreCase(chipType) && !"wt".equalsIgnoreCase(chipType) ) {
			abort("exception_execute_affynormalization", "Error", "Array type (" + chipType + ") not supported", "Array type (" + chipType + ") not supported");			
		}
		// normalizing data
		//
		jobStatus.addStatusMessage("50", "normalizing data");

		List<String> analysis = new ArrayList<String>();
		if ( "3-prime".equalsIgnoreCase(chipInfo.get("type")) ) {

			if ( rma ) analysis.add("rma");
			if ( calls ) analysis.add("pm-mm,mas5-detect.calls=1.pairs=1");
			if ( plier ) analysis.add("plier-mm");

			String cdfFile = config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("cdf");
			AffymetrixExpressionUtils.aptProbesetSummarize(aptBinPath + "/apt-probeset-summarize", analysis, cdfFile, celFiles.getAbsolutePath(), outdir);
		} else {

			if ( rma ) analysis.add("rma");
			if ( calls ) analysis.add("dabg");
			if ( plier ) analysis.add("plier-gcbg");

			String clfFile = config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("clf");
			String pgfFile = config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("pgf");
			String bgpFile = config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("bgp");
			String qccFile = config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("qcc");
			String mpsFile =  config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("mps");
			AffymetrixExpressionUtils.aptProbesetSummarize(aptBinPath + "/apt-probeset-summarize", analysis, clfFile, pgfFile, bgpFile, qccFile, mpsFile, celFiles.getAbsolutePath(), outdir);
		}

		//		System.err.println("cmd output: " + sp.getRunnableProcess().getOutput());
		//		System.err.println("cmd error: " + sp.getRunnableProcess().getError());

		// saving normalization results
		//
		jobStatus.addStatusMessage("90", "saving normalization results");

		File file;
		//List<String> tags = StringUtils.toList("data,datamatrix,expression", ",");


		String maPlotBinPath = babelomicsHomePath + "/bin/plots/plotMA_from_single_matrix.r";

		boolean results = false;
		
		file = new File(outdir + "/rma.summary.txt"); 
		if ( file.exists() ) {
			results = true;
			
			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			ExpressionUtils.createMAPlot(maPlotBinPath, file.getAbsolutePath(), "MA_RMA_", false, "ma_plot.Rout", outdir);

			saveAsDataset(file);

			String tags = "datamatrix,expression";
			result.addOutputItem(new Item("rma.summary", file.getName(), "RMA summary ", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "RMA.Summary"));
			File redirectionFile = new File(outdir + "/rma.summary.redirection");
			createPreprocessingRedirectionFile(redirectionFile, file);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
				result.addOutputItem(new Item("rma.summary", file.getName(), "RMA summary ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "RMA.Summary"));
			}

			saveBoxPlot(file, false, "RMA box-plot", "rmaimg", "RMA.Box-plot");

			addOutputItemImgs(outdir, "MA_RMA_", "png", "ma_plot", "RMA plot", "RMA.MA plot");							
		}

		file = new File(outdir + "/plier-mm.summary.txt"); 
		if ( file.exists() ) {
			results = true;
			
			IOUtils.write(file, cleanLines(IOUtils.readLines(file), true));
			ExpressionUtils.createMAPlot(maPlotBinPath, file.getAbsolutePath(), "MA_PLIER_MM_", false, "ma_plot.Rout", outdir);

			saveAsDataset(file);

			String tags = "datamatrix,expression";
			result.addOutputItem(new Item("plier-mm.summary", file.getName(), "Plier MM summary ", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Plier MM.Summary"));								
			File redirectionFile = new File(outdir + "/plier_mm.summary.redirection");
			createPreprocessingRedirectionFile(redirectionFile, file);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
				result.addOutputItem(new Item("plier-mm.summary", file.getName(), "Plier MM summary ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Plier MM.Summary"));								
			}

			saveBoxPlot(file, false, "Plier MM box-plot", "plierimg", "Plier MM.Box-plot");

			addOutputItemImgs(outdir, "MA_PLIER_MM", "png", "ma_plot", "MA plot", "Plier MM.MA plot");							
		}

		file = new File(outdir + "/plier-gcbg.summary.txt"); 
		if ( file.exists() ) {
			results = true;

			IOUtils.write(file, cleanLines(IOUtils.readLines(file), true));
			ExpressionUtils.createMAPlot(maPlotBinPath, file.getAbsolutePath(), "MA_PLIER_GCBG_", false, "ma_plot.Rout", outdir);

			saveAsDataset(file);

			String tags = "datamatrix,expression";
			result.addOutputItem(new Item("plier-gcbg.summary", file.getName(), "Plier GCBG summary ", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Plier GCBG.Summary"));								
			File redirectionFile = new File(outdir + "/plier_gcbg.summary.redirection");
			createPreprocessingRedirectionFile(redirectionFile, file);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
				result.addOutputItem(new Item("plier-gcbg.summary", file.getName(), "Plier GCBG summary ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Plier GCBG.Summary"));								
			}

			saveBoxPlot(file, false, "Plier GCBG box-plot", "plierimg", "Plier GCBG.Box-plot");

			addOutputItemImgs(outdir, "MA_PLIER_GCBG_", "png", "ma_plot", "MA plot", "Plier GCBG.MA plot");							
		}

		file = new File(outdir + "/pm-mm.mas5-detect.summary.txt"); 
		if ( file.exists() ) {
			results = true;

			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			saveAsDataset(file);

			String tags = "datamatrix,expression";
			result.addOutputItem(new Item("pm-mm.summary", file.getName(), "PM-MM summary ", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Present-absent calls"));								
			File redirectionFile = new File(outdir + "/pm_mm.summary.redirection");
			createPreprocessingRedirectionFile(redirectionFile, file);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
				result.addOutputItem(new Item("pm-mm.summary", file.getName(), "PM-MM summary ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Present-absent calls"));								
			}

			//saveBoxPlot(file, "PM-MM box-plot", "pmmmimg", "Present-absent calls");				
		}

		file = new File(outdir + "/pm-mm.mas5-detect.calls.txt"); 
		if ( file.exists() ) {
			results = true;

			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			result.addOutputItem(new Item("pm-mm.calls", file.getName(), "Calls ", TYPE.FILE, new ArrayList<String>(1), new HashMap<String, String>(1), "Present-absent calls"));								
		}

		file = new File(outdir + "/dabg.summary.txt"); 
		if ( file.exists() ) {
			results = true;

			List<String> lines = cleanLines(IOUtils.readLines(file)); 
			IOUtils.write(file, lines);
			saveAsDataset(file);

			String tags = "datamatrix,expression";
			result.addOutputItem(new Item("dabg.summary", file.getName(), "DABG summary ", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Present-absent calls"));
			File redirectionFile = new File(outdir + "/dabg.summary.redirection");
			createPreprocessingRedirectionFile(redirectionFile, file);
			if ( redirectionFile.exists() ) {
				tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
				result.addOutputItem(new Item("dabg.summary", file.getName(), "DABG summary ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(1), "Present-absent calls"));
			}



			File callFile = new File(outdir + "/dabg.calls.txt");
			List<String> newLines = new ArrayList<String>();
			String [] values;
			double value;
			List<String> valueList = new ArrayList<String>();
			for(String line: lines) {

				if ( line != null && line.trim() != null && line.trim().length() > 0 ) {
					line = line.trim();
					if ( line.startsWith("#NAMES") ) {
						newLines.add(line);
					} if ( !line.startsWith("#") ) {
						values = line.split("\t");
						valueList.clear();
						valueList.add(values[0]);
						for(int i=1 ; i<values.length ; i++) {
							value = Double.parseDouble(values[i]);
							if ( value < 0.05 ) {
								valueList.add("P");
							} else if ( value > 0.065 ) {
								valueList.add("A");
							} else {
								valueList.add("M");
							}
						}
						newLines.add(ListUtils.toString(valueList, "\t"));
					}
				}
			}
			IOUtils.write(callFile, newLines);
			result.addOutputItem(new Item("dabg.summary.calls", callFile.getName(), "Calls ", TYPE.FILE, new ArrayList<String>(1), new HashMap<String, String>(1), "Present-absent calls"));


			//saveBoxPlot(file, "DABG box-plot", "dabgimg", "Present-absent calls");				
		}

		if (!results) {
			abort("expressionnormalization", "Error", "Could not normalize your data", "Could not normalize your data");						
		}

		//		String maPlotBinPath = babelomicsHomePath + "/bin/plots/plotMA_from_single_matrix.r";
		//		ExpressionUtils.createMAPlot(maPlotBinPath, ExpressionUtils.getNormalizedFileName(), "MA_", "ma_plot.Rout", outdir);
		//		addOutputItemImgs(outdir, "MA_", "png", "ma_plot", "MA plot");				

	}


	/**
	 * 
	 * @param rawFilenames
	 * @param chipNames
	 * @return
	 * @throws InvalidParameterException 
	 * @throws IOException 
	 * @throws Exception
	 */
	private String getChipName(List<String> rawFilenames, List<String> chipNames) throws InvalidParameterException, IOException {
		String chipName = null;
		List<String> lines = null;
		List<String> foundChips = new ArrayList<String>();
		List<List<String>> results = new ArrayList<List<String>>(rawFilenames.size());

		//System.out.println("----------> getChipName, chip names = " + ListUtils.toString(chipNames, ", "));

		for(int i=0 ; i<rawFilenames.size() ; i++) {

			//System.out.println("----------> getChipName, file name = " + rawFilenames.get(i));

			results.add(i, new ArrayList<String>());
			lines = IOUtils.head(new File(rawFilenames.get(i)), 20);
			for(String name: chipNames) {
				for(String line: lines) {
					if ( line.contains(name) ) {
						System.out.println("found " + name + " in file " + rawFilenames.get(i));
						results.get(i).add(name);
						foundChips.add(name);
					}
				}
				//if ( results.get(i) != null ) break;
			}
		}

		if ( foundChips.size() == 0 ) {
			throw new InvalidParameterException("array type not supported");
		}


		int maxLength = 0;
		for(int i=0 ; i<foundChips.size() ; i++) {
			if ( foundChips.get(i).length() > maxLength ) {
				maxLength = foundChips.get(i).length();
				chipName = foundChips.get(i);
			}
		}

		//System.out.println("***** results = " + ListUtils.toString(results));
		//chipName = results.get(0);
		List<String> errors = new ArrayList<String>();
		for (int i=0 ; i<rawFilenames.size() ; i++) {
			if ( ! results.get(i).contains(chipName) ) {
				String msg = "mismatch CEL files corresponding to different chips, please, check your cel files.\n";
				//				for(int j=0 ; j<rawFilenames.size() ; j++) {
				//					if ( ! results.get(i).contains(chipName) ) {
				//						msg = msg + "Raw file '" + new File(rawFilenames.get(j)).getName() + "' is a '" + results.get(j) + "' array\n";
				//					}
				//				}
				throw new InvalidParameterException(msg);
			}
		}
		return chipName;
	}

	/**
	 * 
	 * @param lines
	 * @return
	 */
	private List<String> cleanLines(List<String> lines) {
		return cleanLines(lines, false);
	}

	private List<String> cleanLines(List<String> lines, boolean log2) {
		double value;
		String []values = null;
		List<String> newLine = new ArrayList<String>();
		List<String> result = new ArrayList<String>();
		for(String line: lines) {
			if ( line.startsWith("#") ) {
			} else if ( line.startsWith("probeset_id") ) {
				//line = line.replace("\t", ",");
				//result.add(line.replace("probeset_id,", "#NAMES\t"));
				result.add(line.replace("probeset_id", "#NAMES"));
			} else {
				if ( log2 ) {
					double log2Value = Math.log(2);

					values = line.split("\t");
					if ( values != null && values.length > 0 ) {
						newLine.clear();
						newLine.add(values[0]);
						for(int i=1 ; i<values.length ; i++) {
							try {
								value = Math.log(Double.parseDouble(values[i])) / log2Value;
							} catch (Exception e) {
								value = Double.NaN;
							}
							newLine.add(String.valueOf(value));
						}
						result.add(ListUtils.toString(newLine, "\t"));
					}
				} else {
					result.add(line);
				}
			}
		}
		return result;
	}

	/**
	 * 
	 * @param file
	 * @throws IOException
	 * @throws InvalidIndexException
	 */
	private void saveAsDataset(File file) throws IOException, InvalidIndexException {
		Dataset dataset = new Dataset(file);
		dataset.load();
		if ( new File(sampleInfoPath).exists() ) {
			String sampleInfo = ListUtils.toString(IOUtils.readLines(new File(sampleInfoPath)), "\n");
			dataset.setVariables(sampleInfo);
		}
		dataset.save();
	}

	public void saveBoxPlot(File file, boolean log2, String title, String resultId, String group) throws IOException, InvalidIndexException {
		File imgFile = new File(file.getAbsolutePath().replace(".txt", ".png"));
		Dataset dataset = new Dataset(file, true);
		if (!dataset.load() && !dataset.validate()) {
			abort("exception_execute_clustering", "Error", "Error loading dataset " + file.getName() + ": " + ListUtils.toString(dataset.getMessages().getErrorMessages(),". "), "");				
		}
		//Dataset dataset = new Dataset(file, true);
		BoxPlotChart bpc = new BoxPlotChart(title, "", "");
		bpc.getLegend().setVisible(false);

		if ( dataset.getColumnDimension() != 0 && dataset.getRowDimension() != 0 ) {
			if ( log2 ) {
				double log2Value = Math.log(2);
				double[]srcValues = null;
				double [] values = new double[dataset.getRowDimension()];

				for(int i=0; i<dataset.getColumnDimension(); i++) {
					srcValues = dataset.getDoubleMatrix().getColumn(i);
					for(int j=0; j<dataset.getRowDimension(); j++) {
						values[j] = Math.log(srcValues[j]) / log2Value;		
					}
					bpc.addSeries(values, "samples", dataset.getSampleNames().get(i));
				}
			} else {
				for(int i=0; i<dataset.getColumnDimension(); i++) {
					try {
						bpc.addSeries(dataset.getDoubleMatrix().getColumn(i), "samples", dataset.getSampleNames().get(i));
					} catch (Exception e) {
						System.out.println("***** file name: " + file.getAbsolutePath());
						System.out.println("***** is dataset.getDoubleMarix is null ? " + (dataset.getDoubleMatrix()==null));
						System.out.println("***** is column " + i + " null ? " + (dataset.getDoubleMatrix().getColumn(i) == null));
						System.out.println("***** is sample " + i + " null ? " + (dataset.getSampleNames() == null || dataset.getSampleNames().get(i) == null));
					}
				}
			}
		}

		try {
			ChartUtilities.saveChartAsPNG(imgFile, bpc, 400, 256);
			if ( imgFile.exists() ) {
				result.addOutputItem(new Item(resultId, imgFile.getName(), title, TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), group));							
			} else {
				printError("error saving " + title, "error saving " + title, "error saving " + title);
			}
		} catch (IOException e) {
			e.printStackTrace();
			printError("error generating " + title, "error generating " + title, "error generating " + title);
		}
	}

	/**
	 * 
	 * @param rawFileNames
	 * @return
	 */
	private List<String> getSamples(List<String> rawFileNames) {
		int index;
		File file;
		List<String> samples = new ArrayList<String>(rawFileNames.size());
		for(int i=0 ; i<rawFileNames.size() ; i++) {
			file = new File(rawFileNames.get(i));
			index = file.getName().lastIndexOf('.');
			if ( (index > 0) && (index <= file.getName().length() - 2) ) {
				samples.add(file.getName().substring(0, index));
			} else {
				samples.add("Sample_" + i);
			}
		}
		return samples;
	}

	public void createPreprocessingRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=preprocessing");
		redirectionInputs.add("jobname=preprocessing");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("merge_replicates=mean");
		redirectionInputs.add("impute_missing=mean");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}	
}

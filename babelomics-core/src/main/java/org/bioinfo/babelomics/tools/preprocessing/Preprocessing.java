package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Preprocessing extends BabelomicsTool {

	public Preprocessing() {
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "The dataset"));
		options.addOption(OptionFactory.createOption("logarithm-base", "The logarithm base to apply transformation, possible values: e, 2 10 for log base 2, log base 2 and log base 10", false));
		options.addOption(OptionFactory.createOption("exp", "The exponential function, possible values: e, 2, 10", false));
		options.addOption(OptionFactory.createOption("merge-replicates", "Method to merge replicates, valid values are: mean or median", false));
		options.addOption(OptionFactory.createOption("filter-missing", "Minimum percentage of existing values, from 0 to 100", false));
		options.addOption(OptionFactory.createOption("impute-missing", "Method to impute missing values, valid values are: zero, mean, median, knn", false));
		options.addOption(OptionFactory.createOption("kvalue", "K-value for knn impute method, default 15", false));
		options.addOption(OptionFactory.createOption("extract-ids", "This option will extract the IDs (first column) in the dataset", false));
		options.addOption(OptionFactory.createOption("gene-file-filter", "This option will remove all the patterns of the genes that are not present in this gene file", false));
		options.addOption(OptionFactory.createOption("convert-ids", "This option will convert IDs. User can convert IDs to the following terms: ensembl_gene, ensembl_transcript", false));
		options.addOption(OptionFactory.createOption("convert-ids-merge", "Method to merge the converted IDs if replicated, valid values are mean or median", false));
		
		//		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		//		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	public void execute2() {

	}

	@Override
	public void execute() {
		boolean preprocessed = false;
		
		Dataset dataset = null; //, auxDataset = null;
		double filterPercentage = 0;
		int kvalue = 15;
		String imputeMethodMsg = "";
		String logBase = commandLine.getOptionValue("logarithm-base", "none");
		String exp = commandLine.getOptionValue("exp",  "none");
		String mergeMethod = commandLine.getOptionValue("merge-replicates",  "none");
		try {
			filterPercentage = Double.parseDouble(commandLine.getOptionValue("filter-missing",  "0.0"));
			if (filterPercentage>100 || filterPercentage<0) {
				filterPercentage = 0;
			}
		} catch(Exception e) {
			filterPercentage = 0;
		}
//		String minFilterPercentage = commandLine.getOptionValue("min_perc_int","0"); 
		String imputeMethod = commandLine.getOptionValue("impute-missing",  "none");
		String extractIds = commandLine.getOptionValue("extract-ids",  "none");
		String filterFilename = commandLine.getOptionValue("gene-file-filter",  "none");
		String convertIds = commandLine.getOptionValue("convert-ids",  "none");
		String convertIdsMerge = commandLine.getOptionValue("convert-ids-merge",  "mean");

		int progress = 1;
		int finalProgress = 3;
		if ( logBase != null && !("none".equalsIgnoreCase(logBase)) ) finalProgress++;
		if ( exp != null && !("none".equalsIgnoreCase(exp)) ) finalProgress++;
		if ( mergeMethod != null && !("none".equalsIgnoreCase(mergeMethod)) ) finalProgress++;
		if ( filterPercentage > 0 ) finalProgress++;
		if ( imputeMethod != null && !("none".equalsIgnoreCase(imputeMethod)) ) finalProgress++; 
		if ( extractIds != null && !("none".equalsIgnoreCase(extractIds)) ) finalProgress++; 
		if ( filterFilename != null && !("none".equalsIgnoreCase(filterFilename)) && new File(filterFilename).exists() ) finalProgress++; 
		if ( convertIds != null && !("none".equalsIgnoreCase(convertIds)) ) finalProgress++; 

		
		// input parameters
		//
		result.addOutputItem(new Item("log_input_param", logBase, "Logarithmic transformation", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("exp_input_param", exp, "Exponential function", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("merge_input_param", mergeMethod, "Merge replicates", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("filter_input_param", filterPercentage + " %", "Minimum percentage of existing values", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("impute_input_param", ("knn".equalsIgnoreCase(imputeMethod) ? ("knn, k-value = " + commandLine.getOptionValue("k-value", "15")) : imputeMethod) , "Impute missing values", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("extractid_input_param", "" + (extractIds != null && !extractIds.equalsIgnoreCase("none")), "Extract IDs", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("filenamefilter_input_param", (filterFilename != null && !"none".equalsIgnoreCase(filterFilename) ? "none" : new File(filterFilename).getName()), "ID-file filter", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("convertid_input_param", (convertIds != null && !convertIds.equalsIgnoreCase("none") ? (convertIds + ", merge method: " + convertIdsMerge) : "none"), "Convert IDs", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		
		imputeMethodMsg = imputeMethod;
		if ( "knn".equalsIgnoreCase(imputeMethod) ) {
			try {
				kvalue = Integer.parseInt(commandLine.getOptionValue("k-value", "15"));
			} catch (NumberFormatException e) {
				abort("numberformatexception_execute_preprocessing", "Invalid k-value '" + commandLine.getOptionValue("kvalue", "15") + "' for knn imputation", e.toString(), StringUtils.getStackTrace(e));
			}
			imputeMethodMsg = imputeMethodMsg + ", k-value = " + kvalue;
		}

		try {
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "reading dataset");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		String datasetPath = commandLine.getOptionValue("dataset");
		if ( datasetPath == null ) {
			abort("missingdataset_execute_preprocessing", "Missing dataset", "Missing dataset", "Missing dataset");
		}
		File datasetFile = new File(commandLine.getOptionValue("dataset"));
		try {
			dataset = new Dataset(datasetFile);
		} catch (Exception e) {
			abort("exception_execute_preprocessing", "error reading dataset '" + datasetFile.getName() + "'", e.toString(), StringUtils.getStackTrace(e));
		}		

		//		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
		//			dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), "");
		//		}

		if ( dataset == null ) {
			abort("datasetisnull_execute_preprocessing", "dataset is null", "Dataset is null after reading file '" + datasetFile.getName() + "'", "Dataset is null after reading file " + datasetFile.getAbsolutePath());
		}

		progress++;



		// apply logarithm
		//
		if ( logBase != null && !("none".equalsIgnoreCase(logBase)) ) {

			try {
				jobStatus.addStatusMessage(StringUtils.decimalFormat((double)progress*100/finalProgress, "##.00"), "applying logarithm base " + logBase);
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("executing logarithm base " + logBase + "...\n");

			try {
				if ( dataset.getDoubleMatrix() == null ) { 
					try {
						dataset.load();
					} catch (Exception e) {
						abort("exception_logbase_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before applying logarithm base " + logBase, e.toString(), StringUtils.getStackTrace(e));
					}
					dataset.validate();
				}

				dataset.setDoubleMatrix(dataset.getDoubleMatrix().applyLogarithm(logBase));				
				dataset.validate();					
			} 
			catch (CloneNotSupportedException e) {
				abort("cloneexecption_logbase_execute_preprocessing", "logarithm base preprocessing error", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("end of executing logarithm base " + logBase + "\n");
			progress++;
			preprocessed = true;
		}

		// apply exponential
		//
		if ( exp != null && !("none".equalsIgnoreCase(exp)) ) {

			try {
				jobStatus.addStatusMessage(StringUtils.decimalFormat((double)progress*100/finalProgress, "##.00"), "applying exponential " + exp);
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("executing exponential " + exp + "...\n");

			try {
				if ( dataset.getDoubleMatrix() == null ) { 
					try {
						dataset.load();
					} catch (Exception e) {
						abort("exception_logbase_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before applying logarithm base " + logBase, e.toString(), StringUtils.getStackTrace(e));
					}
					dataset.validate();
				}

				dataset.setDoubleMatrix(dataset.getDoubleMatrix().applyExponential(exp));				
				dataset.validate();					
			} 
			catch (CloneNotSupportedException e) {
				abort("cloneexecption_exponential_execute_preprocessing", "exponential preprocessing error", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("end of executing exponential " + exp + "\n");
			progress++;
			preprocessed = true;
		}

		// merge replicated rows
		//
		if ( mergeMethod != null && !("none".equalsIgnoreCase(mergeMethod)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "merging replicated rows (" + mergeMethod + ")");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("merging replicated rows (" + mergeMethod + ")...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					abort("exception_mergereplicated_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before merging replicated rows with method " + mergeMethod, e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();
			}

			try {
				dataset = dataset.mergeReplicatedFeatures(mergeMethod);
//				auxDataset = dataset.mergeReplicatedFeatures(mergeMethod);
//				auxDataset.save("/tmp/PreprocessingTest/newdataset.txt");
//				dataset = auxDataset;
			} catch (Exception e) {
				abort("exception_mergereplicated_execute_preprocessing", "Error merging replicated rows with method " + mergeMethod, e.toString(), StringUtils.getStackTrace(e));
			}

			logger.debug("end of merging replicated rows (" + mergeMethod + ")\n");
			progress++;
			preprocessed = true;
		}


		// filter missing values according to the given percentage
		//
		if ( filterPercentage > 0 ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "filtering missing values by percentage " + filterPercentage);
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("filtering missing values by percentage " + filterPercentage + "...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					abort("exception_filterpermissingvalues_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before filtering by missing values, percentage " + filterPercentage + "%", e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();
			}

			try {
				dataset = dataset.filterRowsByPercOfMissingValues(filterPercentage);
				dataset.validate();
			} catch (Exception e) {
				abort("exception_filterbymissingvalues_execute_preprocessing", "Error filtering by missing values, percentage value of " + filterPercentage + "%", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("end of filtering missing values by percentage " + filterPercentage + "\n");
			progress++;
			preprocessed = true;
		}


		// imputing missing values
		//
		if ( imputeMethod != null && !("none".equalsIgnoreCase(imputeMethod)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "imputing missing values (" + imputeMethodMsg + ")");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("imputing missing values (" + imputeMethodMsg + ")...\n");
			try {
				if ( dataset.getDoubleMatrix() == null ) { 
					try {
						dataset.load();
					} catch (Exception e) {
						abort("exception_imputingmissingvalues_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before imputing missing values with method " + imputeMethodMsg, e.toString(), StringUtils.getStackTrace(e));
					}
					dataset.validate();
				}

				if ( "knn".equalsIgnoreCase(imputeMethod) ) {
					Dataset newDataset = knnImpute(dataset, kvalue);
					if ( newDataset != null ) {
						dataset = newDataset;
					}
				} else {
					dataset.setDoubleMatrix(dataset.getDoubleMatrix().imputeMissingValuesInRows(imputeMethod));
				}
				dataset.validate();

			} catch (CloneNotSupportedException e) {
				abort("cloneexeption_imputemissingvalues_execute_preprocessing", "Error imputing missing values with method " + imputeMethodMsg, e.toString(), StringUtils.getStackTrace(e));
			}	
			logger.debug("end of imputing missing values (" + imputeMethodMsg + ")\n");
			
			progress++;
			preprocessed = true;
		}

		// extract ids
		//
		File idsFile = new File(outdir + "/id_list.txt");
		if ( extractIds != null && !("none".equalsIgnoreCase(extractIds)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "extracting ids...");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("extracting ids...\n");

			try {
				//System.out.println("*************** features names size = " + (dataset.getFeatureNames().size()));
				if ( dataset.getFeatureNames() == null || dataset.getFeatureNames().size() == 0) {
					printWarning("invaliddataset_execute_preprocessing", "Invalid dataset", "Dataset " + datasetFile.getName() + " is not valid");
				} else {
					IOUtils.write(idsFile, dataset.getFeatureNames());
				}
			} catch (IOException e1) {
				printWarning("ioexception_execute_preprocessing", "IO error", "Error ocurred when accessing input file: " + datasetFile.getName(), "");
			}				

			logger.debug("end of extracting ids\n");
			progress++;
		}


		// filter by gene names
		//
		if ( filterFilename != null && !("none".equalsIgnoreCase(filterFilename)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "filtering by names");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("filtering by names...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					abort("exception_filterbynames_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before filtering by names", e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();
			}

			List<String> validGenes = null;
			try {
				//System.out.println("---> " + filterFilename + ", content = " + IOUtils.toString(filterFilename));

				validGenes = IOUtils.grep(filterFilename, "[^#].+", false);
				//				validGenes = IOUtils.grep(filterFilename, "f");

				//System.out.println("size of valid genes = " + validGenes.size());
				if ( validGenes.size() > 0 ) {
					//System.out.println("list = " + ListUtils.toString(validGenes));
				}


			} catch (IOException e) {
				abort("ioexecption_filterbynames_execute_preprocessing", "Error reading names from file '" + new File(filterFilename).getName() + "' when filtering by names", e.toString(), StringUtils.getStackTrace(e));
			}

			//System.out.println("valid genes = " + ListUtils.toString(validGenes));

			if ( validGenes != null && validGenes.size() > 0 ) {
				List<String> genes = dataset.getFeatureNames();					
				List<Integer> rows = new ArrayList<Integer>();
				for (int row=0 ; row<genes.size() ; row++) {
					if ( !validGenes.contains(genes.get(row))) {
						rows.add(row);
					}
				}

				try {
					dataset = dataset.filterRows(rows);
				} catch (Exception e) {
					abort("exception_filterbynames_execute_preprocessing", "Error filtering rows by names", e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();	


				//System.out.println("------------> feature data is null ? " + (dataset.getFeatureData() == null));
			}

			logger.debug("end of filtering by names\n");
			progress++;
			preprocessed = true;
		}

		// convert ids
		//
		//File idsFile = new File(outdir + "/id_list.txt");
		if ( convertIds != null && !("none".equalsIgnoreCase(convertIds)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "converting ids...");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("converting ids...\n");

			try {
				List<Integer> rows = new ArrayList<Integer>();
				List<String> ids = dataset.getFeatureNames();
				List<String> outIds = new ArrayList<String>(ids.size());
				
				DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));
								
				List<String> dbNames = new ArrayList<String>();
				//dbNames.add("ensembl_gene");
				dbNames.add(convertIds);
				
//				List<String> outIds = new ArrayList<String>(); 
//				Map<String, List<String>> idsMap = new HashMap<String, List<String>>(); 
				
				System.out.println("-> infrared config = " + babelomicsHomePath + "/conf/infrared.conf");
				
				System.out.println("-> db connector = " + dbConnector.toString());
				System.out.println("-> db names = " + ListUtils.toString(dbNames, ","));
				List<Map<String, FeatureList<XRef>>> list = new XRefDBManager(dbConnector).getListByDBNames(ids, dbNames);

				List<String> noFound = new ArrayList<String>();
				if ( list != null && list.size() > 0 ) {
					//System.out.println("list size = " + list.size());
					String key;
					for(int i=0 ; i<list.size() ; i++) {
						//System.out.println("Converting " + ids.get(i));
						key = convertIds;
						//for(String key: MapUtils.getKeys(list.get(i))) {
							//System.out.print("\tto " + key + " ---> ");
							if (list.get(i).get(key)!=null && list.get(i).get(key).size()>0) {
								//for(XRef xref: list.get(i).get(key)) {
								//	System.out.print(xref.getId() + "\t");
								//}
								//System.out.print(list.get(i).get(key).get(0).getId());
								outIds.add(list.get(i).get(key).get(0).getId());
							} else {
								rows.add(i);
								outIds.add(ids.get(i));
								noFound.add(ids.get(i));
							}
							//System.out.println("");
						//}
					}
					
					dataset.setFeatureNames(outIds);
					if (noFound.size()==outIds.size()) {
						printWarning("preprocessing", "Warning", "None input ID could be converted (" + convertIds + "). The result dataset keeps the original IDs.");
					} else {
						if (noFound.size()>0) {
							System.out.println("rows to remove: " + rows.size() + ", indices: " + ListUtils.toString(rows, ","));
							System.out.println("before, number of rows = " + dataset.getRowDimension());
							//dataset.save(new File(this.getOutdir() + "/before_remove_rows.txt"));
							dataset = dataset.filterRows(rows);
							//dataset.save(new File(this.getOutdir() + "/after_remove_rows.txt"));
							System.out.println("after, number of rows = " + dataset.getRowDimension());
						
							String msg = ListUtils.toString(noFound, ", ");
							printWarning("preprocessing", "Warning", "None (" + convertIds + ") ID found for " + noFound.size() + " input IDs: " + (msg.length()>36 ? (msg.substring(0,36)+"..."): msg) + ". These IDs were removed from the result dataset.");
						}
						dataset = dataset.mergeReplicatedFeatures(convertIdsMerge);
						preprocessed = true;
					}
				} else {
					printWarning("preprocessing", "Warning", "None input ID could be converted (" + convertIds + "). The result dataset keeps the original IDs.");
				}
				
				
				//System.out.println("*************** features names size = " + (dataset.getFeatureNames().size()));
//				if ( dataset.getFeatureNames() == null || dataset.getFeatureNames().size() == 0) {
//					printWarning("invaliddataset_execute_preprocessing", "Invalid dataset", "Dataset " + datasetFile.getName() + " is not valid");
//				} else {
//					IOUtils.write(idsFile, dataset.getFeatureNames());
//				}
			} catch (Exception e) {
				e.printStackTrace();
				printWarning("exception_execute_preprocessing", "IO error", "Error ocurred when accessing input file: " + datasetFile.getName(), "");
			}				

			logger.debug("end of converting ids\n");
			progress++;
		}


		logger.debug("saving dataset...\n");
		try {
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "saving results");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}



		try {
			if ( preprocessed ) {

				File file = new File(this.getOutdir() + "/preprocessed.txt");
				dataset.save(file);


				if ( file.exists() ) {
					
					String tags = "datamatrix,expression";
					result.addOutputItem(new Item("prepocessed_file", file.getName(), "Preprocessed file", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Preprocessed data"));			

					List<String> tagList = new ArrayList<String>();
					File redirectionFile = new File(outdir + "/clustering.redirection");
					createClusteringRedirectionFile(redirectionFile, file);
					if ( redirectionFile.exists() ) {
						tagList.add("REDIRECTION(" + redirectionFile.getName() + ":Send to Clustering tool...)");
					}

					if (dataset.getVariables() != null && dataset.getVariables().size() > 0 ) {
						redirectionFile = new File(outdir + "/classcomparison.redirection");
						createClassComparisonRedirectionFile(redirectionFile, file);
						if ( redirectionFile.exists() ) {
							tagList.add("REDIRECTION(" + redirectionFile.getName() + ":Send to Class-comparison tool...)");
						}

						redirectionFile = new File(outdir + "/correlation.redirection");
						createCorrelationRedirectionFile(redirectionFile, file);
						if ( redirectionFile.exists() ) {
							tagList.add("REDIRECTION(" + redirectionFile.getName() + ":Send to Correlation tool...)");
						}

						redirectionFile = new File(outdir + "/classprediction.redirection");
						createClassPredictionRedirectionFile(redirectionFile, file);
						if ( redirectionFile.exists() ) {
							tagList.add("REDIRECTION(" + redirectionFile.getName() + ":Send to Class-prediction tool...)");
						}
					}
					if (tagList.size() > 0) {
						result.addOutputItem(new Item("prepocessed_file", file.getName(), "Preprocessed file", TYPE.FILE, tagList, new HashMap<String, String>(2), "Preprocessed data"));
					}
				}
			}

			if ( idsFile.exists() ) {				
				String tags = "idlist";											
				result.addOutputItem(new Item("ids_file", idsFile.getName(), "File containing IDs", TYPE.DATA, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "ID list"));
				File redirectionFile = new File(outdir + "/idconverter.redirection");
				createIDConverterRedirectionFile(redirectionFile, idsFile);
				if ( redirectionFile.exists() ) {
					tags = "REDIRECTION(" + redirectionFile.getName() + ":Send to ID converter tool...)";
					result.addOutputItem(new Item("ids_file", idsFile.getName(), "File containing IDs", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "ID list"));
				}
			}

		} catch (IOException e) {
			abort("ioexception_savingresults_execute_preprocessing", "error saving output file", e.toString(), StringUtils.getStackTrace(e));
		}		
	}


	private void createClusteringRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=clustering");
		redirectionInputs.add("jobname=clustering");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("method=upgma");
		redirectionInputs.add("distance=euclidean");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void createClassComparisonRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=class-comparison");
		redirectionInputs.add("jobname=class-comparison");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("correction=fdr");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void createCorrelationRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=correlation");
		redirectionInputs.add("jobname=correlation");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("test=pearson");
		redirectionInputs.add("correction=fdr");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void createClassPredictionRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=class-prediction");
		redirectionInputs.add("jobname=class-prediction");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("svm=svm");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void createIDConverterRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=id-converter");
		redirectionInputs.add("jobname=id converter");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("listfile_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("listfile=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("listfile_wum_data=true");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * 
	 * @param dataset
	 * @param kvalue
	 * @return
	 */
	private Dataset knnImpute(Dataset dataset, int kvalue) {
		Dataset newDataset = null;
		String wd = outdir + "/knn";
		try {
			File inputFile = new File(wd + "/in.txt");
			File outputFile = new File(wd + "/out.txt");
			if ( new File(wd).isDirectory() || FileUtils.createDirectory(wd) ) {
				dataset.save(inputFile);

				String cmdStr = babelomicsHomePath + "/bin/KNNimpute -K=" + kvalue + " " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath();
				Command cmd = new Command(cmdStr); 
				SingleProcess sp = new SingleProcess(cmd);
				sp.runSync();

				System.out.println("cmd = " + cmdStr);				
				if ( outputFile.exists() ) {
					List<String> values;
					List<String> lines = IOUtils.grep(outputFile, "[^#].+");
					DoubleMatrix matrix = new DoubleMatrix(dataset.getRowDimension(), dataset.getColumnDimension());
					for(int row=0 ; row<lines.size() ; row++) {
						values = StringUtils.toList(lines.get(row), "\t");
						values.remove(0);
						matrix.setRow(row, ArrayUtils.toDoubleArray(ListUtils.toStringArray(values)));
					}
					newDataset = dataset;
					newDataset.setDoubleMatrix(matrix);
				}
			}
		} catch (IOException e) {
			newDataset = null;
		}
		//FileUtils.deleteDirectory(new File(wd));
		return newDataset;
	}
}

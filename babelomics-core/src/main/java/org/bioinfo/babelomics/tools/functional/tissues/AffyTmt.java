package org.bioinfo.babelomics.tools.functional.tissues;

import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.methods.functional.InfraredUtils;
import org.bioinfo.babelomics.tools.functional.FunctionalUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.db.DBConnection;
import org.bioinfo.db.api.PreparedQuery;
import org.bioinfo.db.api.Query;
import org.bioinfo.db.handler.BeanArrayListHandler;
import org.bioinfo.db.handler.MatrixHandler;
import org.bioinfo.db.handler.ResultSetHandler;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.MultipleTestCorrection;
import org.bioinfo.math.util.MathUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class AffyTmt extends Tmt {

	protected String dbName = "gnf_human";
	DBConnection dbConn; 
	
	//	private Map<String, String> ensemblIdtoGene1 = new HashMap<String, String>(), ensemblIdtoGene2 = new HashMap<String, String>();
	//	private int replicated1 = 0, replicated2 = 0;
	//	private List<String> genesWithoutEnsemblIds1 = new ArrayList<String>(), genesWithoutEnsemblIds2 = new ArrayList<String>();
	//	private Map<String, List<String>> genesWithMultipleEnsemblIds1 = new HashMap<String, List<String>>(), genesWithMultipleEnsemblIds2 = new HashMap<String, List<String>>();
	//	
	//	private Map<String, List<String>> null, probesMap2 = null;config.getProperty("INFRARED.HOST")
	//	private List<String> genesWithoutProbes1 = null, genesWithoutProbes2 = null;

	public AffyTmt() {
		
	}

	public void initOptions() {
		super.initOptions();
		options.addOption(OptionFactory.createOption("normalization-method", "Normalization method, valid values are 'mas5' and 'gcrma'. Defalut value: 'mas5'", false));
		options.addOption(OptionFactory.createOption("multiple-probes", "Multiple probes expression value, valid values are 'mean', 'greatest', 'lowest', 'perc25', 'perc50' and 'perc75'. Default value: 'mean'", false));
	}

	public void execute() {
		File f1 = new File(commandLine.getOptionValue("list1"));
		File f2 = commandLine.hasOption("list2") ? new File(commandLine.getOptionValue("list2")) :  null;

		String normMethod = commandLine.getOptionValue("normalization-method", "mas5");
		String multipleProbes = commandLine.getOptionValue("multiple-probes", "mean");

		List<String> tissues = StringUtils.toList(commandLine.getOptionValue("tissues"), ",");
		
		System.out.println("tissues from getOptionValue: >>> " + commandLine.getOptionValue("tissues") + " <<<");
		System.out.println("tissues: >>>> " + ListUtils.toString(tissues) + " <<<<");		

		try {		
			dbName = "mmu".equalsIgnoreCase(species) ? "gnf_mouse" : "gnf_human";
			config.append(new File(babelomicsHomePath + "/conf/infrared.properties"));
			dbConn = new DBConnection("mysql", config.getProperty("INFRARED.HOST"), config.getProperty("INFRARED.PORT"), dbName, config.getProperty("INFRARED.USER"), config.getProperty("INFRARED.PASSWORD"));

			if ( tissues.contains("alltissues") ) {
				tissues = getAllTissues(species);
			}

			System.out.println("tissues:\n" + ListUtils.toString(tissues));

			DBConnector dbConnector = new DBConnector("mouse".equalsIgnoreCase(species) ? "mmu" : "hsa", new File(babelomicsHomePath + "/conf/infrared.properties"));
			//DBConnector dbConnector = new DBConnector(species, "mysqlweb", "3306", "biouser", "biopass");
			System.out.println("species = " + species + ", db name = " +  dbName + ", db connector = " + dbConnector.toString());

			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			// handling list #1
			//
			List<String> geneList1, uniqueGeneList1, dupGeneList1 = null;
			List<String> genesToConvert1, ensemblList1, noConverted1 = null;
			Map<String, List<String>> ensemblMap1 = null;
			geneList1 = IOUtils.readLines(f1);
			if ( geneList1 == null || geneList1.size() == 0 ) {
				throw new Exception("No genes found in list #1");
			}
			uniqueGeneList1 = ListUtils.unique(geneList1);
			if ( geneList1.size() != uniqueGeneList1.size() ) {
				dupGeneList1 = ListUtils.duplicated(geneList1);
				logger.debug("removing " + dupGeneList1.size() + " duplicated genes from List #1: " + ListUtils.toString(dupGeneList1, ","));
				System.out.println("removing " + dupGeneList1.size() + " duplicated genes from List #1: " + ListUtils.toString(dupGeneList1, ","));
			}
			
			ensemblList1 = new ArrayList<String> ();
			genesToConvert1 = new ArrayList<String> ();
			for(String id: uniqueGeneList1) {
				if ( FunctionalUtils.isEnsemblID(id) ) {
					ensemblList1.add(id);
				} else {
					genesToConvert1.add(id);
				}
			}
			if ( genesToConvert1.size() > 0 ) {
				logger.debug("found " + genesToConvert1.size() + " IDs non-Ensembl ID in List #1: " + ListUtils.toString(genesToConvert1, ","));
				System.out.println("found " + genesToConvert1.size() + " IDs non-Ensembl ID in List #1: " + ListUtils.toString(genesToConvert1, ","));
				ensemblMap1 = InfraredUtils.getEnsemblMap(dbConnector, genesToConvert1);
				if ( ensemblMap1 == null || ensemblMap1.size() == 0 ) {
					logger.debug("No Ensembl IDs found for your input genes in List #1 when converting your input genes to Ensembl ID");
					System.out.println("No Ensembl IDs found for your input genes in List #1 when converting your input genes to Ensembl ID");
				}
				noConverted1 = new ArrayList<String> ();
				for(String key: ensemblMap1.keySet()) {
					if ( ensemblMap1.get(key) != null && ensemblMap1.get(key).size() > 0 ) {
						ensemblList1.addAll(ensemblMap1.get(key));
					} else {
						noConverted1.add(key);
					}
				}
			}
			uniqueGeneList1 = ListUtils.unique(ensemblList1);
			if ( uniqueGeneList1.size() == 0 ) {
				throw new Exception("No Ensembl IDs found for list #1 when converting to Ensembl ID");
			}
			
			Map<String, List<String>> probesMap1 = getProbes(species, uniqueGeneList1);
			System.out.println(probesMap1);
			if ( probesMap1 == null || probesMap1.size() == 0 ) {
				throw new Exception("No Affymetrix probes found for list #1");
			}


			// handling list #2
			//
			List<String> geneList2, uniqueGeneList2, dupGeneList2 = null;
			List<String> genesToConvert2, ensemblList2, noConverted2 = null;
			Map<String, List<String>> ensemblMap2 = null;

			if ( f2 == null ) {
				List<String> genes = getAllGenes(species);
				
				ensemblList2 = new ArrayList<String>();
				for (String gene: genes) {
					if ( !uniqueGeneList1.contains(gene) ) {
						ensemblList2.add(gene);
					}
				}
			} else {
				geneList2 = IOUtils.readLines(f2);
				if ( geneList2 == null || geneList2.size() == 0 ) {
					throw new Exception("No genes found in list #2");
				}
				uniqueGeneList2 = ListUtils.unique(geneList2);
				if ( geneList2.size() != uniqueGeneList2.size() ) {
					dupGeneList2 = ListUtils.duplicated(geneList2);
					logger.debug("removing " + dupGeneList2.size() + " duplicated genes from List #2: " + ListUtils.toString(dupGeneList2, ","));
					System.out.println("removing " + dupGeneList2.size() + " duplicated genes from List #2: " + ListUtils.toString(dupGeneList2, ","));
				}
				
				ensemblList2 = new ArrayList<String> ();
				genesToConvert2 = new ArrayList<String> ();
				for(String id: uniqueGeneList2) {
					if ( FunctionalUtils.isEnsemblID(id) ) {
						ensemblList2.add(id);
					} else {
						genesToConvert2.add(id);
					}
				}
				if ( genesToConvert2.size() > 0 ) {
					logger.debug("found " + genesToConvert1.size() + " IDs non-Ensembl ID in List #2: " + ListUtils.toString(genesToConvert2, ","));
					ensemblMap2 = InfraredUtils.getEnsemblMap(dbConnector, genesToConvert2);
					if ( ensemblMap2 == null || ensemblMap2.size() == 0 ) {
						logger.debug("No Ensembl IDs found for your input genes in List #2 when converting your input genes to Ensembl ID");
						System.out.println("No Ensembl IDs found for your input genes in List #2 when converting your input genes to Ensembl ID");
					}
					noConverted2 = new ArrayList<String> ();
					for(String key: ensemblMap2.keySet()) {
						if ( ensemblMap2.get(key) != null && ensemblMap2.get(key).size() > 0 ) {
							ensemblList2.addAll(ensemblMap2.get(key));
						} else {
							noConverted2.add(key);
						}
					}
				}			
			}
			
			uniqueGeneList2 = ListUtils.unique(ensemblList2);
			if ( ensemblList2.size() == 0 ) {
				throw new Exception("No Ensembl IDs found for list #2 when converting to Ensembl ID");
			}				
			
			logger.debug("unique genes: " + uniqueGeneList2);
			System.out.println("unique genes: " + uniqueGeneList2);
			Map<String, List<String>> probesMap2 = getProbes(species, uniqueGeneList2);
			if ( probesMap2 == null || probesMap2.size() == 0 ) {
//				abort("exception_execute_tmt", "No Affymetrix probes found for list #2", "1", "2");
//				return;
				throw new Exception("No Affymetrix probes found for list #2");
			}

			// getting libraries
			//
			System.out.println("tissues = " + tissues);
			Map<String, String> libraryMap = getTissueLibraries(species, tissues);
			List<String> libraryNames = MapUtils.getKeys(libraryMap);

			if ( libraryMap == null || libraryMap.size() ==  0 ) {
				throw new Exception("No tissue libraries found for your input parameters. Please, change your parameters and try again");
			}
			
			// getting frequency matrixes
			//
			jobStatus.addStatusMessage("40", "getting frequency matrixes");
			logger.debug("getting frequency matrixes...\n");

			DoubleMatrix matrix1 = getFreqMatrix(uniqueGeneList1, libraryNames, libraryMap, probesMap1, normMethod, multipleProbes);
			if ( matrix1 == null || matrix1.getColumnDimension() ==  0 || matrix1.getRowDimension() == 0) {
				throw new Exception("Impossible to create the frequency matrix for list #1 for your input parameters. Please, change your parameters and try again");
			}
						
			DoubleMatrix matrix2 = getFreqMatrix(uniqueGeneList2, libraryNames, libraryMap, probesMap2, normMethod, multipleProbes);
			if ( matrix1 == null || matrix1.getColumnDimension() ==  0 || matrix1.getRowDimension() == 0) {
				throw new Exception("Impossible to create the frequency matrix for list #2 for your input parameters. Please, change your parameters and try again");
			}

			// computing t-test
			//
			jobStatus.addStatusMessage("60", "computing t-test");
			logger.debug("computing t-test...\n");
			
			TestResultList<TTestResult> res = runTtest(matrix1, matrix2);				
//			TTest tTest = new TTest();
//			TestResultList<TTestResult> res = tTest.tTest(matrix1, matrix2);				
			MultipleTestCorrection.BHCorrection(res);

			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(libraryNames.size(), 0);

			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(libraryNames, rowOrder));
						
			if ( dupGeneList1 != null && dupGeneList1.size() > 0 ) {
				IOUtils.write(new File(outdir + "/duplicated_list1.txt"), dupGeneList1);
				result.addOutputItem(new Item("duplicated_list1_file", "duplicated_list1.txt", "Duplicated genes from list #1 (they were not taken into accout)", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "List management"));
			}
			if ( dupGeneList2 != null && dupGeneList2.size() > 0 ) {
				IOUtils.write(new File(outdir + "/duplicated_list2.txt"), dupGeneList2);
				result.addOutputItem(new Item("duplicated_list2_file", "duplicated_list2.txt", "Duplicated genes from list #2 (they were not taken into accout)", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "List management"));
			}

			if ( noConverted1 != null && noConverted1.size() > 0 ) {
				IOUtils.write(new File(outdir + "/no_converted_list1.txt"), noConverted1);
				result.addOutputItem(new Item("no_converted_list1_file", "no_converted_list1.txt", "Not found Ensembl IDs for the following genes of list #1", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "List management"));
			}
			if ( noConverted2 != null && noConverted2.size() > 0 ) {
				IOUtils.write(new File(outdir + "/no_converted_list2.txt"), noConverted2);
				result.addOutputItem(new Item("no_converted_list2_file", "no_converted_list2.txt", "Not found Ensembl IDs for the following genes of list #2", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "List management"));
			}
			
			FeatureData featureData = new FeatureData(dataFrame);
			featureData.save(new File(outdir + "/affy-tmt.txt"));
			//result.addOutputItem(new Item("tmt_file", "affy-tmt.txt", "Output file:", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Output file"));
			result.addOutputItem(new Item("tmt_file", "affy-tmt.txt", "Output file:", TYPE.FILE, StringUtils.toList("TABLE,DIFF_EXPRESSION_TABLE"), new HashMap<String, String>(), "Output file"));
			
//			FileUtils.writeStringToFile(new File(getOutdir() + "/replicates.txt"), "Removed " + replicated1 + " genes from list 1\nRemoved " + replicated2 + " genes from list 2");

//			writeProbesMap(new File(getOutdir() + "/gene_probes_1.txt"), probesMap1, ensemblIdtoGene1);
//			writeProbesMap(new File(getOutdir() + "/gene_probes_2.txt"), probesMap2, ensemblIdtoGene1);
//
//			if ( genesWithoutEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_without_ensemblid_1.txt"), genesWithoutEnsemblIds1, "#Genes from list 1 without ensembl ID\n");
//			if ( genesWithoutEnsemblIds2 != null ) writeGeneList(new File(getOutdir() + "/genes_without_ensemblid_2.txt"), genesWithoutEnsemblIds2, "#Genes from list 2 without ensembl ID\n");
//
//			if ( genesWithoutProbes1 != null ) writeGeneList(new File(getOutdir() + "/genes_without_probes_1.txt"), genesWithoutProbes1, "#Genes from list 1 without probes\n");
//			if ( genesWithoutProbes2 != null ) writeGeneList(new File(getOutdir() + "/genes_without_probes_2.txt"), genesWithoutProbes2, "#Genes from list 2 without probes\n");
//
//			if ( genesWithMultipleEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_with_multiple_ensemblid_1.txt"), genesWithMultipleEnsemblIds1, "#Genes from list 1 with multiple Ensembl IDs, they have been removed from the analysis\n");
//			if ( genesWithMultipleEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_with_multiple_ensemblid_2.txt"), genesWithMultipleEnsemblIds2, "#Genes from list 2 with multiple Ensembl IDs, they have been removed from the analysis\n");


		} catch (Exception e) {
			abort("exception_execute_tmt", "Affymetrix profiling tissue error", e.getMessage(), StringUtils.getStackTrace(e));
		}
	}


//	private void writeProbesMap(File file, Map<String, List<String>> probes, Map<String, String> genes) throws IOException {
//		StringBuilder sb = new StringBuilder();
//		sb.append("#gene").append("\tensembl ID\tprobes\n");
//		for(String key: probes.keySet()) {
//			sb.append(genes.get(key)).append("\t").append(key).append("\t").append(ListUtils.toString(probes.get(key), " ")).append("\n");
//		}
//		//		FileUtils.writeStringToFile(file, sb.toString());
//	}

//	private void writeGeneList(File file, List<String> list, String msg) throws IOException {
//		if ( list != null && list.size() > 0 ) {
//			StringBuilder sb = new StringBuilder();
//			sb.append(msg);
//			for(String item: list) {
//				sb.append(item).append("\n");
//			}
//			//			FileUtils.writeStringToFile(file, sb.toString());
//		}
//	}

//	private void writeGeneList(File file, Map<String, List<String>> map, String msg) throws IOException {
//		if ( map != null && map.size() > 0 ) {
//			StringBuilder sb = new StringBuilder();
//			sb.append(msg);
//			for(String key: map.keySet()) {
//				sb.append(key).append("\t").append(ListUtils.toString(map.get(key), " ")).append("\n");
//			}
//			//			FileUtils.writeStringToFile(file, sb.toString());
//		}
//	}


//	private List<String> createListFromMap(Map<String, List<String>> map) {
//		List<String> list = new ArrayList<String> ();
//		for(String key: map.keySet()) {
//			list.add(map.get(key).get(0));
//		}
//		return list;
//	}

//	private Map<String, List<String>> cleanGeneMap(Map<String, List<String>> geneMap, boolean firstList) {
//		List<String> list;
//		Map<String, List<String>> map = new HashMap<String, List<String>> ();
//		for(String key: geneMap.keySet()) {
//			list = geneMap.get(key);
//			if ( list == null ) {
//				if ( firstList ) {
//					genesWithoutEnsemblIds1.add(key);
//				} else {
//					genesWithoutEnsemblIds2.add(key);					
//				}
//			} else if ( list.size() > 1 ) {
//				if ( firstList ) {
//					genesWithMultipleEnsemblIds1.put(key, list);					
//				} else {
//					genesWithMultipleEnsemblIds2.put(key, list);					
//				}
//			} else {
//				map.put(key, geneMap.get(key));
//			}
//		}
//		return map;
//	}

	private DoubleMatrix getFreqMatrix(List<String> geneList, List<String> libraryNames, Map<String, String> libraryMap, Map<String, List<String>> probeMap, String normMethod, String freqMethod) {

		List<String> genes = new ArrayList<String>();
		for(String gene: geneList) {
			if ( probeMap.get(gene) != null ) {
				genes.add(gene);
			}
		}

		DoubleMatrix matrix = new DoubleMatrix(libraryNames.size(), genes.size());

		double freq;
		List<Double> values = new ArrayList<Double>();

		Query query;
		Object[][] expressions;
		ResultSetHandler rsh = new MatrixHandler();
		String q, prefix = "select expression from " + ("gcrma".equalsIgnoreCase(normMethod) ? "expression_gcRMA" : "expression_MAS5") + " where";

		int total = 0, row = 0, column = 0;
		for(String libraryKey: libraryNames) {
			column = 0;
			for(String gene: genes) {

				//				System.out.println("(row, column) = (" + row + ", " + column + "), (library, gene) = (" + libraryKey + ", " + gene + "), freqMethod = " + freqMethod); 

				values.clear();
				for(String probe: probeMap.get(gene)) {
					q = prefix + " probeset_id = \"" + probe + "\" and tissue = \"" + libraryMap.get(libraryKey) + "\""; 
					//						System.out.println("dbName = " + dbName + ", query = " + q);
					try {
						query = dbConn.createSQLQuery(q);
						expressions = (Object[][]) query.execute(rsh);
						if ( expressions != null && expressions.length > 0 ) {
							total = 0;
							for(int i=0 ; i<expressions.length ; i++) {
								//									System.out.print(expressions[i][0] + " ");
								total += (Integer) expressions[i][0];
							}
							//								System.out.println("mean = " + (1.0 * total / expressions.length));
							values.add(1.0 * total / expressions.length);
						}
					} catch (Exception e) {
						probeMap = null;
						printError("exception_get_probes_tmt", "get probes tmt error", e.toString(), e);
					}
				}
				freq = getFrequency(values, freqMethod); 
				//System.out.println("probes for gene " +  gene + " values = " + ListUtils.toString(values) + ", freq = " + freq);					

				matrix.set(row, column, freq);
				column++;
			}
			row++;
		}			
		return matrix;
	}

	private double getFrequency(List<Double> values, String freqMethod) {

		double[] array = ArrayUtils.toDoubleArray(values);

		if ( "mean".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.mean(array);
		} else if ( "max".equalsIgnoreCase(freqMethod) ) {
			return ArrayUtils.max(array);
		} else if ( "min".equalsIgnoreCase(freqMethod) ) {
			return ArrayUtils.min(array);
		} else if ( "25".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.percentile(array, 25);
		} else if ( "50".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.percentile(array, 50);
		} else if ( "75".equalsIgnoreCase(freqMethod) ) {			
			return MathUtils.percentile(array, 75);
		}

		return Double.NaN;
	}

	
	public List<String> getAllTissues(String organism) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		List<String> tissues;

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct tissue_name from tissue_info order by tissue_name asc"; 
		Query query;
		query = dbConn.createSQLQuery(q);
		tissues = (ArrayList<String>) query.execute(rsh);

		return tissues;
	}

	public Map<String, List<String>> getProbes(String organism, List<String> geneList) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		Map<String, List<String>> probeMap = new HashMap<String, List<String>>();

		System.out.println("getProbes, database = " + dbConn.getDatabase());
		
		String q;
		PreparedQuery prepQuery;
		List<String> probes;
		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		q = "select probeset_id from ensemblGene where ensemblGene_id = ?"; 
		//			System.out.println("dbName = " + dbName + ", query = " + q);
		prepQuery = dbConn.createSQLPrepQuery(q);
		for(String gene: geneList) {
			prepQuery.setParams(gene);
			probes = (ArrayList<String>) prepQuery.execute(rsh);
			System.out.println("getProbes in AffyTmt: gene: "+gene+"  probes: "+probes);
			if ( probes != null && probes.size() > 0 ) {
				System.out.println("getProbes in AffyTmt: gene: "+gene);
				probeMap.put(gene, probes);
			}
		}
		prepQuery.close();
		return probeMap;
	}

	public List<String> getAllGenes(String organism) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		List<String> genes;

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct ensemblGene_id from ensemblGene"; 
		Query query = dbConn.createSQLQuery(q);
		genes = (ArrayList<String>) query.execute(rsh);

		return genes;
	}

	public Map<String, String> getTissueLibraries(String organism, List<String> tissues) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		Map<String, String> libraryMap = new HashMap<String, String>();

		String q;
		PreparedQuery query;
		List<Library> libraries;
		ResultSetHandler rsh = new BeanArrayListHandler(Library.class);
		for(String tissue: tissues) {
			q = "select distinct tissue_id, tissue_name from tissue_info where tissue_name = ?"; 
			query = dbConn.createSQLPrepQuery(q);
			query.setParams(tissue);
			libraries = (ArrayList<Library>) query.execute(rsh);
			for(Library library: libraries) {
				libraryMap.put(library.name, library.id);
			}
		}

		return libraryMap;
	}

}

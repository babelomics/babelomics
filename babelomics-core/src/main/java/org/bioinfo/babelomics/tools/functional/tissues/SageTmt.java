package org.bioinfo.babelomics.tools.functional.tissues;

import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.methods.functional.InfraredUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
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
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class SageTmt extends Tmt {

	protected String dbName = "sage_long_human_07_11_2007";
	DBConnection dbConn = null; 

	//	private Map<String, String> ensemblIdtoGene1 = new HashMap<String, String>(), ensemblIdtoGene2 = new HashMap<String, String>();
	//	private int replicated1 = 0, replicated2 = 0;
	//	private List<String> genesWithoutEnsemblIds1 = new ArrayList<String>(), genesWithoutEnsemblIds2 = new ArrayList<String>();
	//	private Map<String, List<String>> genesWithMultipleEnsemblIds1 = new HashMap<String, List<String>>(), genesWithMultipleEnsemblIds2 = new HashMap<String, List<String>>();
	//	
	//	private Map<String, List<String>> null, probesMap2 = null;
	//	private List<String> genesWithoutProbes1 = null, genesWithoutProbes2 = null;

	public SageTmt() {
		
	}

	public void initOptions() {
		super.initOptions();
		options.addOption(OptionFactory.createOption("histologies", "the list of tissues separated by commas. Enter 'all histologies' to take into account all available histologies"));
		options.addOption(OptionFactory.createOption("tag-type", "Type of SAGE tags, valid values are 'short' and 'long'", false));
		options.addOption(OptionFactory.createOption("exclude-cell-lines", "Exclude cell lines", false,false));
		options.addOption(OptionFactory.createOption("min-tags", "Minimum number of tags. Defalut value: 5000", false));
		options.addOption(OptionFactory.createOption("perc-null-libraries", "Percentage for null values accepted in libraries (ranging from 0 to 100). Default value: 80", false));
		options.addOption(OptionFactory.createOption("perc-null-genes", "Percentage for null values accepted in genes (ranging from 0 to 100). Default value: 80", false));
	}

	public void execute() {

		File f1 = new File(commandLine.getOptionValue("list1"));
		File f2 = commandLine.hasOption("list2") ? new File(commandLine.getOptionValue("list2")) :  null;
		//		
		List<String> tissues = StringUtils.toList(commandLine.getOptionValue("tissues"), ",");
		List<String> histologies = StringUtils.toList(commandLine.getOptionValue("histologies"), ",");

		String tagType = commandLine.getOptionValue("tag-type", "short");
		boolean excludeCellLines = commandLine.hasOption("exclude-cell-lines");
		int minTags = Integer.parseInt(commandLine.getOptionValue("minTags", "5000"));
		double filterGenes = Double.parseDouble(commandLine.getOptionValue("perc-null-genes", "80"));
		double filterLibraries = Double.parseDouble(commandLine.getOptionValue("perc-null-libraries", "80"));
		
		
		
		//String multipleProbes = commandLine.getOptionValue("multiple-probes", "mean");
		//
		//		
		try {		
			
			dbName = "sage_" + ("long".equalsIgnoreCase(tagType) ? "long" : "short") + "_" + ("mouse".equalsIgnoreCase(species) ? "mouse" : "human") + "_07_11_2007";
			config.append(new File(babelomicsHomePath + "/conf/infrared.properties"));
			dbConn = new DBConnection("mysql", config.getProperty("INFRARED.HOST"), config.getProperty("INFRARED.PORT"), dbName, config.getProperty("INFRARED.USER"), config.getProperty("INFRARED.PASSWORD"));
			System.out.println("db name = " + dbName);
			
			
			if ( tissues.contains("alltissues") ) {
				tissues = getAllTissues(species);
			}

//			if ( histologies.contains("all histologies") ) {
//				histologies = getAllHistologies(species);
//			}

			System.out.println("tissues:\n" + ListUtils.toString(tissues));
			
			DBConnector dbConnector = new DBConnector("mouse".equalsIgnoreCase(species) ? "mmu" : "hsa", new File(babelomicsHomePath + "/conf/infrared.properties"));
			System.out.println("db connector = " + dbConnector.toString());

			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			// handling list #1
			//
			List<String> geneList1, uniqueGeneList1, dupGeneList1 = null;
			List<String> unigeneList1, noConverted1 = null;
			Map<String, List<String>> unigeneMap1 = null;
			geneList1 = IOUtils.readLines(f1);
			uniqueGeneList1 = ListUtils.unique(geneList1);
			if ( geneList1.size() != uniqueGeneList1.size() ) {
				dupGeneList1 = ListUtils.duplicated(geneList1);
				logger.debug("removing " + dupGeneList1.size() + " duplicated genes from List #1: " + ListUtils.toString(dupGeneList1, ","));
			}
			unigeneMap1 = InfraredUtils.getUnigeneMap(dbConnector, uniqueGeneList1);
			if ( unigeneMap1 == null || unigeneMap1.size() == 0 ) {
				throw new Exception("No Unigene IDs for list #1 when convertion to Unigene IDs");
			}
			noConverted1 = new ArrayList<String> ();
			unigeneList1 = new ArrayList<String> ();
			for(String key: unigeneMap1.keySet()) {
				if ( unigeneMap1.get(key) != null && unigeneMap1.get(key).size() > 0 ) {
					unigeneList1.addAll(unigeneMap1.get(key));
				} else {
					noConverted1.add(key);
				}
			}
			uniqueGeneList1 = ListUtils.unique(unigeneList1);

			System.out.println("initial list1 = " + ListUtils.toString(geneList1, ","));
			System.out.println("unigene list1 = " + ListUtils.toString(unigeneList1, ","));
			System.out.println("genes no converted to unigene = " + ListUtils.toString(noConverted1, ","));
			System.out.println("unique gene list1 = " + ListUtils.toString(uniqueGeneList1, ","));

//			Map<String, List<String>> probesMap1 = getProbes(species, uniqueGeneList1);
//			if ( probesMap1 == null || probesMap1.size() == 0 ) {
//				throw new Exception("No Affymetrix probes found for gene list #1");
//			}

			// handling list #2
			//
			List<String> geneList2 = null, uniqueGeneList2 = null, dupGeneList2 = null;
			List<String> unigeneList2 = null, noConverted2 = null;
			Map<String, List<String>> unigeneMap2 = null;

			if ( f2 == null ) {
				List<String> genes = getAllGenes(species);

				uniqueGeneList2 = new ArrayList<String>();
				for (String gene: genes) {
					if ( !uniqueGeneList1.contains(gene) ) {
						uniqueGeneList2.add(gene);
					}
				}
			} else {
				geneList2 = IOUtils.readLines(f2);
				uniqueGeneList2 = ListUtils.unique(geneList2);
				if ( geneList2.size() != uniqueGeneList2.size() ) {
					dupGeneList2 = ListUtils.duplicated(geneList2);
					logger.debug("removing " + dupGeneList1.size() + " duplicated genes from List #2: " + ListUtils.toString(dupGeneList2, ","));
				}
					unigeneMap2 = InfraredUtils.getUnigeneMap(dbConnector, uniqueGeneList2);
					if ( unigeneMap2 == null || unigeneMap2.size() == 0 ) {
						throw new Exception("No Unigene IDs found for list #2 when converting to Unigene IDs");
					}
					noConverted2 = new ArrayList<String> ();
					unigeneList2 = new ArrayList<String> ();
					for(String key: unigeneMap2.keySet()) {
						if ( unigeneMap2.get(key) != null && unigeneMap2.get(key).size() > 0 ) {
							unigeneList2.addAll(unigeneMap2.get(key));
						} else {
							noConverted2.add(key);
						}
					}
					uniqueGeneList2 = ListUtils.unique(unigeneList2);
				}

//			Map<String, List<String>> probesMap2 = getProbes(species, uniqueGeneList2);
//			if ( probesMap2 == null || probesMap2.size() == 0 ) {
//				throw new Exception("No Affymetrix probes found for gene list #2");
//			}

			System.out.println("initial list2 = " + ListUtils.toString(geneList2, ","));
			System.out.println("unigene list2 = " + ListUtils.toString(unigeneList2, ","));
			System.out.println("genes no converted to unigene = " + ListUtils.toString(noConverted2, ","));
			System.out.println("unique gene list2 = " + ListUtils.toString(uniqueGeneList2, ","));
			
			

			// getting libraries: first list: ids, second list: names
			//
			List<List<String>> libraries = getLibraries(tissues, histologies, excludeCellLines, minTags);
//			Map<String, String> libraryMap = getTissueLibraries(species, tissues);
//			List<String> libraryNames = MapUtils.getKeys(libraryMap);

			List<String> libraryIDs = libraries.get(0);
			List<String> libraryNames = libraries.get(1);
			System.out.println("libraries ids: " + ListUtils.toString(libraries.get(0), ","));
			System.out.println("libraries names: " + ListUtils.toString(libraries.get(1), ","));
			
			if ( libraryIDs == null || libraryIDs.size() ==  0 ) {
				throw new Exception("No tissue libraries found for your input parameters. Please, change your parameters and try again");
			}
			//System.out.println("frequencies for list1");
			//getFrequencies(uniqueGeneList1);
			
			//System.out.println("frequencies for list2");
			//getFrequencies(uniqueGeneList2);

			List<String> unigenes = new ArrayList<String>();
			unigenes.addAll(uniqueGeneList1);
			unigenes.addAll(uniqueGeneList2);
			unigenes = ListUtils.unique(unigenes);
			
			Map<String, Double> frequencyMap = getFrequencyMap(unigenes);

			if ( frequencyMap == null || frequencyMap.size() ==  0 ) {
				throw new Exception("Impossible to create the frequency map for your input parameters. Please, change your parameters and try again");
			}

			// getting frequency matrixes
			//
			jobStatus.addStatusMessage("40", "getting frequency matrixes");
			logger.debug("getting frequency matrixes...\n");

			System.out.println("matrix 1");
			DoubleMatrix matrix1 = getFreqMatrix(uniqueGeneList1, libraryIDs, frequencyMap);
			System.out.println("matrix1 :\n" + matrix1);
			if ( matrix1 == null || matrix1.getColumnDimension() ==  0 || matrix1.getRowDimension() == 0) {
				throw new Exception("Impossible to create the frequency matrix for list #1 for your input parameters. Please, change your parameters and try again");
			}
			
			System.out.println("matrix 2");
			DoubleMatrix matrix2 = getFreqMatrix(uniqueGeneList2, libraryIDs, frequencyMap);
			System.out.println("matrix2 :\n" + matrix2);
			if ( matrix2 == null || matrix2.getColumnDimension() ==  0 || matrix2.getRowDimension() == 0) {
				throw new Exception("Impossible to create the frequency matrix for list #2 for your input parameters. Please, change your parameters and try again");
			}
			
			System.out.println("filtering....\n");
			
			// filtering...
			//
			List<Integer> columnIndexes1 = getColumns(matrix1, filterLibraries);			
			System.out.println("matrix1, columns: " + ListUtils.toString(columnIndexes1, ","));
			List<Integer> rowIndexes1 = getRows(matrix1, filterGenes);
			System.out.println("matrix1, row: " + ListUtils.toString(rowIndexes1, ","));

			List<Integer> columnIndexes2 = getColumns(matrix2, filterLibraries);
			System.out.println("matrix2, columns: " + ListUtils.toString(columnIndexes2, ","));
			List<Integer> rowIndexes2 = getRows(matrix2, filterGenes);
			System.out.println("matrix2, row: " + ListUtils.toString(rowIndexes2, ","));
			
			List<Integer> rowIndexes = new ArrayList<Integer>();
//			rowIndexes.addAll(rowIndexes1);
//			rowIndexes.addAll(rowIndexes2);
//			rowIndexes = ListUtils.unique(rowIndexes);
			rowIndexes = ListUtils.intersection(rowIndexes1, rowIndexes2);
			
			if ( columnIndexes1 == null || columnIndexes1.size() ==  0 || columnIndexes2 == null || columnIndexes2.size() ==  0 || rowIndexes == null || rowIndexes.size() == 0) {
				throw new Exception("Impossible to create the frequency matrixes for your filtering parameters. Please, change your parameters and try again");
			}
			
			matrix1 = new DoubleMatrix(matrix1.getSubMatrix(ListUtils.toIntArray(rowIndexes), ListUtils.toIntArray(columnIndexes1)).getData());
			System.out.println("matrix1 :\n" + matrix1);
			if ( matrix1 == null || matrix1.getColumnDimension() ==  0 || matrix1.getRowDimension() == 0) {
				throw new Exception("Impossible to create the frequency matrix for list #1 for your filtering parameters. Please, change your parameters and try again");
			}

			matrix2 = new DoubleMatrix(matrix2.getSubMatrix(ListUtils.toIntArray(rowIndexes), ListUtils.toIntArray(columnIndexes2)).getData());
			System.out.println("matrix2 :\n" + matrix2);
			if ( matrix2 == null || matrix2.getColumnDimension() ==  0 || matrix2.getRowDimension() == 0) {
				throw new Exception("Impossible to create the frequency matrix for list #2 for your filtering parameters. Please, change your parameters and try again");
			}

			// computing t-test (but removing NaN before running the t-test)
			//
			jobStatus.addStatusMessage("60", "computing t-test");
			logger.debug("computing t-test...\n");
			

			List<String> names = new ArrayList<String>(rowIndexes.size());
			for(int index: rowIndexes) {
				names.add(libraryNames.get(index));
			}

			TestResultList<TTestResult> res = runTtest(matrix1, matrix2);

//			TTest tTest = new TTest();
//			TestResultList<TTestResult> res = tTest.tTest(matrix1, matrix2);				
			MultipleTestCorrection.BHCorrection(res);

			int[] rowOrder = ArrayUtils.order(res.getStatistics(), true);

			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(names.size(), 0);

			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(names, rowOrder));

			System.out.println("results :\n" + dataFrame.toString(true, true));

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
			featureData.save(new File(outdir + "/sage-tmt.txt"));
			//result.addOutputItem(new Item("sage_file", "sage-tmt.txt", "Output file:", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Output file"));
			result.addOutputItem(new Item("sage_file", "sage-tmt.txt", "Output file:", TYPE.FILE, StringUtils.toList("TABLE,DIFF_EXPRESSION_TABLE"), new HashMap<String, String>(), "Output file"));
			
		} catch (Exception e) {
			abort("exception_execute_tmt", "SAGE tissue profiling error", e.getMessage(), StringUtils.getStackTrace(e));
		}
	}

	private DoubleMatrix getFreqMatrix(List<String> geneList, List<String> libraryIDs, Map<String, Double> frequencyMap) {
		String key;
		Double value;
		DoubleMatrix matrix = new DoubleMatrix(libraryIDs.size(), geneList.size());
		for(int row=0 ; row<libraryIDs.size() ; row++) {
			for(int column=0 ; column<geneList.size() ; column++) {
				value = Double.NaN;
				key = libraryIDs.get(row) + "|" + geneList.get(column);
				if ( frequencyMap.containsKey(key) ) {
					value = frequencyMap.get(key);
				} else {
					//System.out.println("not found: " + key);
				}
				matrix.set(row, column, value);
			}
		}
		return matrix;
	}

	private Map<String, Double> getFrequencyMap(List<String> unigenes) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {

		Map<String, Double> res = new HashMap<String, Double>();

		String q;
		PreparedQuery query;
		Object[][] info;
		ResultSetHandler rsh = new MatrixHandler();

		for(int i=0 ; i<unigenes.size() ; i++) {
			q = "select f.numeric_library_id, f.frequency from frequencies f, best_tag b where f.tag = b.tag and b.uniGene_cluster_number = ?";
			query = dbConn.createSQLPrepQuery(q);
			query.setParams(unigenes.get(i));
				info = (Object[][]) query.execute(rsh);
				if ( info != null && info.length > 0 ) {
					for(int k=0; k<info.length ; k++) {
						res.put(info[k][0] + "|" + unigenes.get(i), Double.parseDouble(""+info[k][1]));
					}
	        }
		}
		
		return res;
	}
	
	public List<List<String>> getLibraries(List<String> tissues, List<String> histologies, boolean cellLine, int minTags) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {

		List<List<String>> res = new ArrayList<List<String>>(2);

	    String endStatement;

	    // did user choose exclude cell lines?
	    if ( cellLine ) {
	        endStatement = " and t.preparation != 'cell line';";
	    }
	    else{
	        endStatement = ";";
	    }

		String q;
		PreparedQuery query;
		Object[][] info;
		ResultSetHandler rsh = new MatrixHandler();

		List<String> names = new ArrayList<String> ();
		List<String> ids = new ArrayList<String> ();
		
		for(int i=0 ; i<tissues.size() ; i++) {
	        for(int j=0 ; j<histologies.size() ; j++) {
				q = "select t.numeric_library_id, t.library_name from tissue_info t, libraries l where t.numeric_library_id=l.numeric_library_id and t.tissue = ? and t.histology = ? and l.unique_tags >= " + minTags + endStatement; 
				query = dbConn.createSQLPrepQuery(q);
				query.setParams(tissues.get(i), histologies.get(j));
				info = (Object[][]) query.execute(rsh);
				if ( info != null && info.length > 0 ) {
					for(int k=0; k<info.length ; k++) {
						ids.add("" + info[k][0]);
						names.add("" + info[k][1]);
					}
				}
	        }
		}
		res.add(ids);
		res.add(names);
		
		return res;
	}

	
	public List<String> getAllTissues(String organism) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		List<String> tissues;

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct tissue from tissue_info order by tissue asc"; 
		Query query;
		query = dbConn.createSQLQuery(q);
		tissues = (ArrayList<String>) query.execute(rsh);

		return tissues;
	}
	
	public List<String> getAllHistologies(String organism) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		List<String> histologies;

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct histology from tissue_info order by histology asc"; 
		Query query;
		query = dbConn.createSQLQuery(q);
		histologies = (ArrayList<String>) query.execute(rsh);

		return histologies;
	}

	public List<String> getAllGenes(String organism) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		List<String> genes;

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct uniGene_cluster_number from best_tag"; 
		Query query = dbConn.createSQLQuery(q);
		genes = (ArrayList<String>) query.execute(rsh);

		return genes;
	}
}

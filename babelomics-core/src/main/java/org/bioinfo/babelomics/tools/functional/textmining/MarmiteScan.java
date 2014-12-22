package org.bioinfo.babelomics.tools.functional.textmining;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.methods.functional.textmining.MarmiteTest;
import org.bioinfo.babelomics.methods.functional.textmining.MarmiteUtils;
import org.bioinfo.babelomics.methods.functional.textmining.Score;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;

public class MarmiteScan extends BabelomicsTool {


	public MarmiteScan() {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("list", "ranked gene list"));
		options.addOption(OptionFactory.createOption("bioentity-name", "Valid values: disease, chemical, drug, sympton"));
		options.addOption(OptionFactory.createOption("bioentity-score-filter", "Minimum number of genes with a score (0-10000)", false));
		options.addOption(OptionFactory.createOption("bioentity-number-filter", "Number of bio-entities in results (0-10000)", false));
		options.addOption(OptionFactory.createOption("sort", "Sort list", false));

		//		options.addOption(OptionFactory.createOption("gene-name-list", "gene name in list", false));
		//		options.addOption(OptionFactory.createOption("partition-number", "Number of partitions "));
		//		options.addOption(OptionFactory.createOption("significance", "p-value for statistical significance"));
		//		options.addOption(OptionFactory.createOption("sort", "Sort list"));

	}

	@Override
	public void execute() {
		File listFile = new File(commandLine.getOptionValue("list"));

		String bioentity = commandLine.getOptionValue("bioentity-name");
		int scoreFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-score-filter", "5"));
		int numberFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-number-filter", "50"));
		String sort = commandLine.getOptionValue("sort", "desc");
		
		int order;
		if ( "asc".equalsIgnoreCase(sort) ) {
			order = MarmiteUtils.ASCENDING_SORT;
		} else {
			order = MarmiteUtils.DESCENDING_SORT;
		}
		

		executeMarmiteScan(listFile, bioentity, scoreFilter, numberFilter, order);			
		//executeMarmiteScan(dataset, bioEntityName, bioentityScoreFilter, bioentityNumberFilter,geneNameList, partitionNumber, significance,  sort);			
	}

	//private void executeMarmiteScan(Dataset dataset, String bioEntityName, String bioentityScoreFilter, String bioentityNumberFilter,String geneNameList, String partitionNumber, String significance, String sort) {
	private void executeMarmiteScan(File listFile, String bioEntityName, int bioentityScoreFilter, int bioentityNumberFilter, int order) {

		try {
			int numberOfPartitions = 30;
			List<String> idList = IOUtils.column(listFile, 0);
			List<Double> statistic = ListUtils.toDoubleList(IOUtils.column(listFile, 1));
			int[] indexes = ListUtils.order(statistic, (order == MarmiteUtils.ASCENDING_SORT ? false : true));

			statistic = ListUtils.ordered(statistic, indexes);
			idList = ListUtils.ordered(idList, indexes);

			double inc = -(double)(statistic.get(0)-statistic.get(statistic.size()-1))/(numberOfPartitions+1);
			double acum = statistic.get(0) + inc;

			String entity;
			int thresholdPosition;
			List<String> list1,list2;
			Map<String, Double> statisticMap = new HashMap<String, Double>();
			Map<String, String> sideMap = new HashMap<String, String>();
			//Map<String, Double> adjPvalueMap = new HashMap<String, Double>();
			Map<String, Double> adjPvalueMap = new HashMap<String, Double>();
			Map<String, Double> pvalueMap = new HashMap<String, Double>();

			//			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)).subList(0, nEntities));
			//			dataFrame.addColumn("side", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getSides()), rowOrder)).subList(0, nEntities));
			//			dataFrame.addColumn("box-plot", boxplotList);
			//			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)).subList(0, nEntities));
			//			dataFrame.addColumn("adj. p-value", adjPvalueList);


			Map<String, MarmiteTest> inputMap = new HashMap<String, MarmiteTest>();
			Map<String, TestResultList<KolmogorovSmirnovTestResult>> outputMap = new HashMap<String, TestResultList<KolmogorovSmirnovTestResult>>();

			System.out.println("idList size = " + idList.size());

			for(int i=0; i<numberOfPartitions; i++){

				thresholdPosition = MarmiteUtils.getThresholdPosition(acum, statistic, order);

				System.err.println("** partition = " + i + ": threshold = " + acum + " (" + thresholdPosition + ") ");
				System.err.println("**** list #1: " + 0 + ".." + thresholdPosition + ", list #2: " + (thresholdPosition + 1) + ".." + (idList.size()-1));

				if ( thresholdPosition < idList.size()-1 ) {
					
					list1 = idList.subList(0, thresholdPosition);
					list2 = idList.subList(thresholdPosition + 1, idList.size()-1);

					System.err.println("l1.size: " + list1.size() + "l2.size: " + list2.size());

					try {
						MarmiteTest marmiteTest = new MarmiteTest(babelomicsHomePath + "/conf/infrared.properties");

						TestResultList<KolmogorovSmirnovTestResult> res = marmiteTest.run(list1, list2, bioEntityName, bioentityScoreFilter);

						//auxEntities = marmiteTest.getValidEntities();
						for(int j=0 ; j<res.getAdjPValues().length ; j++) {
							entity = marmiteTest.getValidEntities().get(j);
							if ( !adjPvalueMap.containsKey(entity) || (res.getAdjPValues()[j] < adjPvalueMap.get(entity)) ) {
								statisticMap.put(entity, res.getStatistics()[j]);
								sideMap.put(entity, res.getSides()[j]);
								pvalueMap.put(entity, res.getPValues()[j]);
								adjPvalueMap.put(entity, res.getAdjPValues()[j]);

								inputMap.put(entity, marmiteTest);
								outputMap.put(entity, res);						
							}
						}
					} catch (InvalidParameterException e) {
						System.err.println("---------> Exception executing partition #" + i + ", sizes=(" + list1.size() + ", " + list2.size() + "): " + e.getMessage());
					}
				}

				acum += inc;
			}

			// convert maps to lists, in order to guarantee the same order
			List<Double> statistics = new ArrayList<Double>();
			List<String> sides = new ArrayList<String>();
			List<Double> pvalues = new ArrayList<Double>();
			List<Double> adjPvalues = new ArrayList<Double>();
			List<String> entities = MapUtils.getKeys(adjPvalueMap);
			for(String e: entities) {
				statistics.add(statisticMap.get(e));
				sides.add(sideMap.get(e));
				pvalues.add(pvalueMap.get(e));
				adjPvalues.add(adjPvalueMap.get(e));
			}

			// sort lists according to de adj p-values
			indexes = ListUtils.order(adjPvalues);
			entities = ListUtils.ordered(entities, indexes);
			statistics = ListUtils.ordered(statistics, indexes);
			sides = ListUtils.ordered(sides, indexes);
			pvalues = ListUtils.ordered(pvalues, indexes);
			adjPvalues = ListUtils.ordered(adjPvalues, indexes);




			// entities, statistics, sides, pvalues, adjpvalues,

			DataFrame dataFrame = new DataFrame(entities.size(), 0);
			List<String> resGeneList1 = new ArrayList<String>();
			List<String> resGeneList2 = new ArrayList<String>();
			List<String> boxplotList = new ArrayList<String>();
			BoxPlotChart boxplot;
			List<Double> doubleList1,doubleList2;
			List<Score> scoreList;
			String cleanName, filename, aux;

			int nEntities = Math.min(bioentityNumberFilter, indexes.length); 
			for(int i=0 ; i<nEntities ; i++) {
				entity = entities.get(i);
				cleanName = entity.toLowerCase().replace(" ", "_");

				// list 1
				filename = cleanName + "_list1.txt";
				aux = MarmiteUtils.getEntityGenes(inputMap.get(entity).getEntityMap1().get(entity));
				IOUtils.write(new File(outdir + "/" + filename), aux.replace(",", "\n"));
				if ( aux.length() > 16 ) {
					resGeneList1.add(aux.substring(0, 15) + "...::::" + filename);						
				} else {
					resGeneList1.add(aux + "::::" + filename);	
				}

				// list 2
				filename = cleanName + "_list2.txt";
				aux = MarmiteUtils.getEntityGenes(inputMap.get(entity).getEntityMap2().get(entity));
				IOUtils.write(new File(outdir + "/" + filename), aux.replace(",", "\n"));
				if ( aux.length() > 16 ) {
					resGeneList2.add(aux.substring(0, 15) + "...::::" + filename);						
				} else {
					resGeneList2.add(aux + "::::" + filename);	
				}

				// boxplot
				scoreList = inputMap.get(entity).getEntityMap1().get(entity);
				doubleList1 = new ArrayList<Double> (scoreList.size());
				for (int j=0 ; j<scoreList.size(); j++) {
					doubleList1.add((double) scoreList.get(j).getValue());
				}						

				scoreList = inputMap.get(entity).getEntityMap2().get(entity); 
				doubleList2 = new ArrayList<Double> (scoreList.size());
				for (int j=0 ; j<scoreList.size(); j++) {
					doubleList2.add((double) scoreList.get(j).getValue());
				}						

				//System.out.println("list1, entity = " + entity + ": " + ListUtils.toString(list1));
				//System.out.println("list2, entity = " + entity + ": " + ListUtils.toString(list2));

				boxplot = MarmiteUtils.createBoxplot(entity, adjPvalues.get(i), doubleList1, doubleList2);
				filename = cleanName + ".png";
				boxplotList.add(filename);
				ChartUtilities.saveChartAsPNG(new File(getOutdir() + "/" + filename), boxplot, 200, 100);
				//result.addOutputItem(new Item(cleanName + "_boxplot", filename, "Boxplot for " + entity, TYPE.IMAGE));
			}


			dataFrame.addColumn("statistic", ListUtils.toStringList(statistics.subList(0, nEntities)));
			dataFrame.addColumn("side", ListUtils.toStringList(sides.subList(0, nEntities)));
			dataFrame.addColumn("box-plot", boxplotList);
			dataFrame.addColumn("p-value", ListUtils.toStringList(pvalues.subList(0, nEntities)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(adjPvalues.subList(0, nEntities)));

			dataFrame.setRowNames(entities.subList(0, nEntities));

//			// output file
//			filename = "marmitescan_output.txt";
//			FeatureData featureData = new FeatureData(dataFrame);
//			featureData.save(new File(getOutdir() + "/" + filename));
//			result.addOutputItem(new Item("marmitescan_file", filename, "MARMITE SCAN output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Results"));

			// output table
			filename = "marmitescan_table.txt";
			dataFrame.addColumn("Genes (list #1)", resGeneList1);
			dataFrame.addColumn("Genes (list #2)", resGeneList2);
			IOUtils.write(getOutdir() + "/" + filename, dataFrame.toString(true, true));
			result.addOutputItem(new Item("marmite_table", filename, "MARMITE SCAN output table", TYPE.FILE, StringUtils.toList("TABLE,MARMITE_TABLE", ","), new HashMap<String, String>(), "Results"));			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidIndexException e) {
			e.printStackTrace();
		}	

	}

}

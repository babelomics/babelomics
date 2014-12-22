package org.bioinfo.babelomics.tools.functional.textmining;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.methods.functional.textmining.MarmiteTest;
import org.bioinfo.babelomics.methods.functional.textmining.MarmiteUtils;
import org.bioinfo.babelomics.methods.functional.textmining.Score;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;

public class Marmite extends BabelomicsTool {

	public Marmite() {
	}

	@Override
	public void initOptions() {

		options.addOption(OptionFactory.createOption("list1", "gene list #1"));
		options.addOption(OptionFactory.createOption("list2", "gene list #2"));
		options.addOption(OptionFactory.createOption("bioentity-name", "Valid values: disease, chemical, drug, sympton"));
		options.addOption(OptionFactory.createOption("bioentity-score-filter", "Minimum number of genes with a score (0-10000)", false));
		options.addOption(OptionFactory.createOption("bioentity-number-filter", "Number of bio-entities in results (0-10000)", false));
		options.addOption(OptionFactory.createOption("duplicates", "to remove duplicated IDs, valid values: never, each, ref. By default, never", false));
	}

	@Override
	public void execute() {
		File f1 = new File(commandLine.getOptionValue("list1"));
		File f2 = new File(commandLine.getOptionValue("list2"));

		String bioentity = commandLine.getOptionValue("bioentity-name");
		int scoreFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-score-filter", "5"));
		int numberFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-number-filter", "50"));

		int duplicatesMode = MarmiteUtils.REMOVE_NEVER;
		String duplicates = commandLine.getOptionValue("remove-duplicates", "never");
		if(duplicates.equalsIgnoreCase("each")) duplicatesMode = MarmiteUtils.REMOVE_EACH;
		if(duplicates.equalsIgnoreCase("ref")) duplicatesMode = MarmiteUtils.REMOVE_REF;
		if(duplicates.equalsIgnoreCase("all")) duplicatesMode = MarmiteUtils.REMOVE_ALL;		

		executeMarmite(f1, f2, bioentity, scoreFilter, numberFilter, duplicatesMode);
	}



	private void executeMarmite(File f1, File f2, String bioentity, int scoreFilter, int numberFilter, int duplicatesMode) {

		try {
			String filename;
			String cleanName;
			String aux;
			MarmiteTest marmiteTest = new MarmiteTest(babelomicsHomePath + "/conf/infrared.properties");
			
			// reading data
			//
			jobStatus.addStatusMessage("15", "reading data");
			logger.debug("reading data...\n");
			
			List<String> geneList1 = IOUtils.column(f1, 0);
			List<String> geneList2 = IOUtils.column(f2, 0);

			// managing duplicates
			//
			List<List<String>> cleanLists = MarmiteUtils.removeDuplicates(geneList1, geneList2, duplicatesMode);

			geneList1 = cleanLists.get(0);
			geneList2 = cleanLists.get(1);
			
			// saving clean list
			//
			filename = "clean_list1.txt";
			IOUtils.write(new File(getOutdir() + "/" + filename), geneList1);
			result.addOutputItem(new Item("list1_file", filename, "List #1 (after managing duplicates mode: " + MarmiteUtils.getDuplicatesModeString(duplicatesMode) + ")", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Input"));

			filename = "clean_list2.txt";
			IOUtils.write(new File(getOutdir() + "/" + filename), geneList2);
			result.addOutputItem(new Item("list2_file", filename, "List #2 (after managing duplicates mode: " + MarmiteUtils.getDuplicatesModeString(duplicatesMode) + ")", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Input"));			

			// run test
			//
			TestResultList<KolmogorovSmirnovTestResult> res = marmiteTest.run(geneList1, geneList2, bioentity, scoreFilter);

			// saving exclude entities from list
			//
			filename = "excludedentities_list1.txt";
			MarmiteUtils.writeExcludedEntitiesFile(new File(getOutdir() + "/" + filename), marmiteTest.getExcludedEntities1(), marmiteTest.getEntityMap1());
			result.addOutputItem(new Item("excludedentities_list1_file", filename, "Entities exclude from List #1", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Input"));
			
			filename = "excludedentities_list2.txt";
			MarmiteUtils.writeExcludedEntitiesFile(new File(getOutdir() + "/" + filename), marmiteTest.getExcludedEntities2(), marmiteTest.getEntityMap2());
			result.addOutputItem(new Item("excludedentities_list2_file", filename, "Entities exclude from List #2", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Input"));
			
			//System.out.println("kolmogorov result :\n" +  res.toString());
			//jobStatus.addStatusMessage("#", "ks test, elapsed time = " + ((System.currentTimeMillis()-time)/1000F) + " s");

			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getAdjPValues()));
			
			int nEntities = Math.min(numberFilter, rowOrder.length); 
			List<String> adjPvalueList = ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)).subList(0, nEntities);		
			
			DataFrame dataFrame = new DataFrame(marmiteTest.getValidEntities().size(), 0);
			List<String> rowNames = ListUtils.ordered(marmiteTest.getValidEntities(), rowOrder); 
			List<String> resGeneList1 = new ArrayList<String>();
			List<String> resGeneList2 = new ArrayList<String>();
			List<String> boxplotList = new ArrayList<String>();
			BoxPlotChart boxplot;
			List<Score> scoreList;
			List<Double> list1, list2;
			for(int i=0 ; i<nEntities ; i++) {
				cleanName = rowNames.get(i).toLowerCase().replace(" ", "_");
				
				// list 1
				filename = cleanName + "_list1.txt";
				aux = MarmiteUtils.getEntityGenes(marmiteTest.getEntityMap1().get(rowNames.get(i)));
				IOUtils.write(new File(outdir + "/" + filename), aux.replace(",", "\n"));
				if ( aux.length() > 16 ) {
					resGeneList1.add(aux.substring(0, 15) + "...::::" + filename);						
				} else {
					resGeneList1.add(aux + "::::" + filename);	
				}
				
				// list 2
				filename = cleanName + "_list2.txt";
				aux = MarmiteUtils.getEntityGenes(marmiteTest.getEntityMap2().get(rowNames.get(i)));
				IOUtils.write(new File(outdir + "/" + filename), aux.replace(",", "\n"));
				if ( aux.length() > 16 ) {
					resGeneList2.add(aux.substring(0, 15) + "...::::" + filename);						
				} else {
					resGeneList2.add(aux + "::::" + filename);	
				}
				
				// boxplot
				scoreList = marmiteTest.getEntityMap1().get(rowNames.get(i));
				list1 = new ArrayList<Double> (scoreList.size());
				for (int j=0 ; j<scoreList.size(); j++) {
					list1.add((double) scoreList.get(j).getValue());
				}						

				scoreList = marmiteTest.getEntityMap2().get(rowNames.get(i)); 
				list2 = new ArrayList<Double> (scoreList.size());
				for (int j=0 ; j<scoreList.size(); j++) {
					list2.add((double) scoreList.get(j).getValue());
				}						

				//System.out.println("list1, entity = " + entity + ": " + ListUtils.toString(list1));
				//System.out.println("list2, entity = " + entity + ": " + ListUtils.toString(list2));

				boxplot = MarmiteUtils.createBoxplot(rowNames.get(i), Double.parseDouble(adjPvalueList.get(i)), list1, list2);
				filename = cleanName + ".png";
				boxplotList.add(filename);
				ChartUtilities.saveChartAsPNG(new File(getOutdir() + "/" + filename), boxplot, 200, 100);
				//result.addOutputItem(new Item(cleanName + "_boxplot", filename, "Boxplot for " + rowNames.get(i), TYPE.IMAGE));
			}


			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)).subList(0, nEntities));
			dataFrame.addColumn("side", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getSides()), rowOrder)).subList(0, nEntities));
			dataFrame.addColumn("box-plot", boxplotList);
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)).subList(0, nEntities));
			dataFrame.addColumn("adj. p-value", adjPvalueList);
			
			dataFrame.setRowNames(rowNames.subList(0, nEntities));
			
			// output file
//			filename = "marmite_output.txt";
//			FeatureData featureData = new FeatureData(dataFrame);
//			featureData.save(new File(getOutdir() + "/" + filename));
//			result.addOutputItem(new Item("marmite_file", filename, "MARMITE output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Results"));

			// output table
			filename = "marmite_table.txt";
			dataFrame.addColumn("Genes (list #1)", resGeneList1);
			dataFrame.addColumn("Genes (list #2)", resGeneList2);
			IOUtils.write(getOutdir() + "/" + filename, dataFrame.toString(true, true));
			result.addOutputItem(new Item("marmite_table", filename, "MARMITE output table", TYPE.FILE, StringUtils.toList("TABLE,MARMITE_TABLE", ","), new HashMap<String, String>(), "Results"));

			// saving data
			//
			jobStatus.addStatusMessage("90", "saving results");

		} catch (Exception e) {
			e.printStackTrace();
			printError("exception_executemarmite_marmite", e.toString(), e.getMessage(), e);
		}
	}



}

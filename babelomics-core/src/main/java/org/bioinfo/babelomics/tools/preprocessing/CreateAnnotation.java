package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.InfraredUtils;
import org.bioinfo.babelomics.tools.functional.FunctionalDbDescriptor;
import org.bioinfo.babelomics.tools.functional.FunctionalProfilingTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.GeneDBManager;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class CreateAnnotation extends FunctionalProfilingTool {

	FeatureData list = null;
	List<String> ids = null;
	boolean allGenome;
	String outputFormat;

	public CreateAnnotation() {
		initOptions();
	}

	@Override
	public void initOptions() {
		super.initOptions();
		
		getOptions().addOption(OptionFactory.createOption("listfile", "File containning the IDs to convert", false));
		getOptions().addOption(OptionFactory.createOption("list", "IDs to convert separated by commas", false));
		getOptions().addOption(OptionFactory.createOption("all-genome", "All genome",false));

		// output format
		getOptions().addOption(OptionFactory.createOption("output-format", "Output format: compact or extended. Default: compact", false));
	}

	@Override
	public void prepare() throws IOException, ParseException, InvalidIndexException {
		super.prepare();
		
		// check input IDs
		//
		String inputIds = commandLine.getOptionValue("list", null);			
		if ( inputIds == null ) {
			if ( commandLine.hasOption("listfile") ) {
				// list 1		
				list = new FeatureData(new File(commandLine.getOptionValue("listfile")), true);
			} else if ( commandLine.hasOption("all-genome") ) {
				allGenome = Boolean.parseBoolean(commandLine.getOptionValue("all-genome"));
			} else {
				abort("error_execute_createannotation", "Missing input data", "Missing input data", "Missing input data");
			}
		} else {
			ids = StringUtils.toList(inputIds, ",");
		}

		outputFormat = commandLine.getOptionValue("output-format", "compact");			
		
	}

	@Override
	protected void execute() {
		try {	

			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));	
			System.out.println("species = " + species + ", db connector = " + dbConnector.toString());

			prepare();
			
			if ( allGenome ) {
				GeneDBManager geneDBManager = new GeneDBManager(dbConnector);
				ids = geneDBManager.getAllEnsemblIds();				
			} else if ( list != null ) {
				ids = list.getDataFrame().getRowNames();
			}

			if ( ids == null || ids.size() == 0) {
				abort("error_execute_createannotation", "No IDs found", "No IDs found", "No IDs found");
			}

			FeatureList<AnnotationItem> al = null;

			//System.out.println("dbConnector = " + dbConnector.toString());
			//System.out.println("ids = " + ListUtils.toString(ids, ","));
			//System.out.println("number of filters = " + filterList.size());

			String name = null;

			for(FunctionalFilter filter: filterList) {
				
				// db attributes
				FunctionalDbDescriptor filterInfo = new FunctionalDbDescriptor(filter);
				
				al = InfraredUtils.getAnnotations(dbConnector, ids, filter);

				if ( al != null ) {
					name = filterInfo.getName();
					
					if ( outputFormat.equalsIgnoreCase("extended") ) {
						IOUtils.write(new File(outdir + "/" + name + ".txt"), al.toString());
					} else {
						Map<String, List<String>> map = new HashMap<String, List<String>>();
						for(int i=0 ; i<al.size() ; i++) {
							if ( ! map.containsKey(al.get(i).getId()) ) {
								map.put(al.get(i).getId(), new ArrayList<String>()); 
							}
							map.get(al.get(i).getId()).add(al.get(i).getFunctionalTermId());
						}
						StringBuilder sb = new StringBuilder();
						for(String key: MapUtils.getKeys(map)) {
							sb.append(key).append("\t").append(ListUtils.toString(map.get(key), ",")).append("\n");
						}
						IOUtils.write(new File(outdir + "/" + name + ".txt"), sb.toString());

						//System.out.println(sb.toString());
					}

					result.addOutputItem(new Item(name, name + ".txt", filterInfo.getTitle() + " annotation table", Item.TYPE.FILE, StringUtils.toList("TABLE,ANNOTATION_TABLE", ","), new HashMap<String, String>(), "Annotation tables"));
					result.addOutputItem(new Item(name + "_data", name + ".txt", filterInfo.getTitle() + " annotation table", Item.TYPE.DATA, StringUtils.toList("annotation", ","), new HashMap<String, String>(), ""));
					//result.addOutputItem(new Item(name, name + ".txt", getDBTitle(filter) + " annotation: ", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Results"));					
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.infrared.core.feature.XRef.XRefItem;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class IDConverter extends BabelomicsTool {

	public IDConverter() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("id-file", "File containning the IDs to convert", false));
		getOptions().addOption(OptionFactory.createOption("id-list", "IDs to convert separated by commas", false));

		getOptions().addOption(OptionFactory.createOption("db-names", "DB names to convert", false));
		getOptions().addOption(OptionFactory.createOption("output-format", "Output format: compact or extended. Default: extended", false));		
	}

	@Override
	protected void execute() {
		
		List<String> ids = null;
		List<String> dbNames = null;

		// check input IDs
		if(commandLine.hasOption("id-list") || commandLine.hasOption("id-file")) {
			if(commandLine.hasOption("id-list")) {
				logger.debug("Reading IDs from CLI list: "+commandLine.getOptionValue("id-list"));
				ids = StringUtils.toList(commandLine.getOptionValue("id-list"), ",");	
			}else {
				try {
					logger.debug("Reading IDs from file: "+commandLine.getOptionValue("id-file"));
					ids = IOUtils.column(new File(commandLine.getOptionValue("id-file")), 0);
				} catch (IOException e) {
					abort("ioexception_execute_idconverter", e.getMessage(), e.getMessage(), e.getMessage());
				}
			}
			if(ids == null || ids.size() == 0) {
				abort("error_execute_idconverter", "No IDs found", "No IDs found", "No IDs found");
			}
		}else {
			abort("error_execute_idconverter", "No IDs found", "No IDs found", "No IDs found");
		}
		
		// check output dbnames
		if(commandLine.hasOption("db-names") ) {
			dbNames = StringUtils.toList(commandLine.getOptionValue("db-names"), ","); 
		}
		if(dbNames == null || dbNames.size() == 0) {
			abort("error_execute_idconverter", "Missing DB names", "Missing DB names", "Missing DB names");
		}
		
		// checking rest of parameters
		String outputFormat = commandLine.getOptionValue("output-format", "extended");
		if(species == null || species.equalsIgnoreCase("unknown")) {
			abort("error_execute_idconverter", "Missing species", "Missing species", "Missing species");
		}
		
		try {
			logger.debug("Infrared properties: "+babelomicsHomePath + "/conf/infrared.properties");
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));
			
			StringBuilder line = new StringBuilder();
			List<String> alls = new ArrayList<String>();
			List<String> outIds = new ArrayList<String>(); 
			Map<String, List<String>> idsMap = new HashMap<String, List<String>>(); 
			
			logger.debug("DB names = " + ListUtils.toString(dbNames, ","));
			List<XRef> xrefsList = new XRefDBManager(dbConnector).getByDBName(ids, dbNames);
			logger.debug(xrefsList.toString());
			
			alls.add("#NAMES"+"\t"+ListUtils.toString(xrefsList.get(0).getDbNames(), "\t"));
			if(xrefsList != null && xrefsList.size() > 0) {
				
				for(int i=0 ; i<xrefsList.size() ; i++) {
					line.delete(0, line.length());
					line.append(ids.get(i));
					for(String key: xrefsList.get(i).getDbNames()) {
						outIds.clear();
						for(XRefItem xrefItem: xrefsList.get(i).getXrefItems().get(key)) {
							outIds.add(xrefItem.getDisplayName());
						}
						
						if(!idsMap.containsKey(key) ) {
							idsMap.put(key, new ArrayList<String>());
							idsMap.get(key).add("#NAMES\t" + key);
						}
						idsMap.get(key).add(ids.get(i) + "\t" + ListUtils.toString(outIds,","));
						line.append("\t").append(ListUtils.toString(outIds,","));
					}
					alls.add(line.toString());
				}
				// save results
				String fileName = "ids_summary.txt";
				IOUtils.write(new File(outdir + "/" + fileName), alls);
				result.addOutputItem(new Item("all", fileName, "Conversion (summary file)", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Output files"));
				
				
				StringBuilder conversionResult = null;
				StringBuilder identifiersResult = null;
				logger.debug("DBnames: "+MapUtils.getKeys(idsMap).toString());
				for(String dbname: MapUtils.getKeys(idsMap)) {
					conversionResult = new StringBuilder();
					identifiersResult = new StringBuilder();
					
					conversionResult.append("#NAMES").append("\t").append(dbname).append("\n");
					for(XRef xref: xrefsList) {
						if(xref.getXrefItems().get(dbname) != null && xref.getXrefItems().get(dbname).size() > 0) {
							if(outputFormat.equalsIgnoreCase("extended")) {
								for(XRefItem xrefItem: xref.getXrefItems().get(dbname)) {
									conversionResult.append(xref.getId()).append("\t").append(xrefItem.getDisplayName()).append("\n");
									identifiersResult.append(xrefItem.getDisplayName()).append("\n");
								}
							}else {
								conversionResult.append(xref.getId()).append("\t");
								for(XRefItem xrefItem: xref.getXrefItems().get(dbname)) {
									conversionResult.append(xrefItem.getDisplayName()).append(",");
									identifiersResult.append(xrefItem.getDisplayName()).append("\n");
								}
								// no need of checking if ',' exist, it must exist, if not we are not here.
								conversionResult.replace(conversionResult.length()-1, conversionResult.length(), "\n");
							}
						}else {
							// no IDs found 
							conversionResult.append(xref.getId()).append("\n");
						}
					}
					
					// write files
					IOUtils.write(new File(outdir + "/" + dbname + "_conversion.txt"), conversionResult.toString().trim());
					result.addOutputItem(new Item(dbname, dbname + "_conversion.txt", "ID conversion table", Item.TYPE.FILE, StringUtils.toList("TABLE,IDCONVERTER_TABLE", ","), new HashMap<String, String>(), "ID conversion tables (" + outputFormat + " format)"));
					
					IOUtils.write(new File(outdir + "/" + dbname +".txt"), identifiersResult.toString().trim());
					result.addOutputItem(new Item(dbname + "_ids", dbname +".txt", dbname + " IDs", Item.TYPE.DATA, StringUtils.toList("idlist", ","), new HashMap<String,String>(), "Output files"));
					result.addOutputItem(new Item(dbname + "_file_ids", dbname +".txt", dbname + " IDs", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Output files"));
				}
			}else {		
				// no results found
				IOUtils.write(new File(outdir + "/" + "ids_summary.txt"), alls);			
			}
			
		} catch (Exception e) {
			this.printError("exception_execute_idconverter", e.toString(), e.toString(), e);
		}
	}

}

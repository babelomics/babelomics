package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.FileInputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Properties;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

import es.blast2go.prog.B2GPipe;

public class Blast2GoTool extends BabelomicsTool {


	// config params
	private Properties config;
	public final static String PREFIX = "b2g";
	B2GPipe b2Pipe;
	
	
	public Blast2GoTool() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("xml-file", "Blast file (xml format)"));
		getOptions().addOption(OptionFactory.createOption("xml-version", "Blast file xml version [joint,separate]", false));
		
		getOptions().addOption(OptionFactory.createOption("blast-hits", "Number of blast hits (0-50)", false));		
		getOptions().addOption(OptionFactory.createOption("bh-min-region-length", "Minimal length (in amino acids) of the matching region for a blast hit to be considered (0-300)", false));		
		getOptions().addOption(OptionFactory.createOption("bh-description-filter", "Blast hit description filter", false));	
		getOptions().addOption(OptionFactory.createOption("bh-description-position", "Blast hit description position (0-10)", false));
		getOptions().addOption(OptionFactory.createOption("add-id", "Add ID to Blast definition", false, false));
		getOptions().addOption(OptionFactory.createOption("use-descriptor-annotator", "Use blast descriptor annotator?", false, false));

		getOptions().addOption(OptionFactory.createOption("e-value", "E-value hit filter ((1e-)", false));
		getOptions().addOption(OptionFactory.createOption("annotation-cutoff", "Annotation cut-off (0-100)", false));
		getOptions().addOption(OptionFactory.createOption("go-weight", "GO weight (0-30)", false));
		getOptions().addOption(OptionFactory.createOption("hit-coverage", "Hsp hit coverage cut-off (0-100)", false));

		//EC codes
		getOptions().addOption(OptionFactory.createOption("ida", "ida", false));
		getOptions().addOption(OptionFactory.createOption("ipi", "ipi", false));
		getOptions().addOption(OptionFactory.createOption("imp", "imp", false));
		getOptions().addOption(OptionFactory.createOption("igi", "igi", false));
		getOptions().addOption(OptionFactory.createOption("iep", "iep", false));
		getOptions().addOption(OptionFactory.createOption("exp", "exp", false));
		
		getOptions().addOption(OptionFactory.createOption("iss", "iss", false));
		getOptions().addOption(OptionFactory.createOption("iso", "iss", false));
		getOptions().addOption(OptionFactory.createOption("isa", "iss", false));
		getOptions().addOption(OptionFactory.createOption("ism", "ism", false));
		getOptions().addOption(OptionFactory.createOption("igc", "igc", false));
		getOptions().addOption(OptionFactory.createOption("rca", "rca", false));
		
		getOptions().addOption(OptionFactory.createOption("tas", "tas", false));
		getOptions().addOption(OptionFactory.createOption("nas", "nas", false));
		getOptions().addOption(OptionFactory.createOption("ic", "ic", false));	
		getOptions().addOption(OptionFactory.createOption("nd", "nd", false));
		getOptions().addOption(OptionFactory.createOption("iea", "iea", false));		
		getOptions().addOption(OptionFactory.createOption("nr", "nr", false));




	}

	@Override
	public void execute() {
		
		// load config
		loadConfig();
		
		try{
			
			// xml file
			String xmlFileName = commandLine.getOptionValue("xml-file");
			File xmlFile = new File(xmlFileName);
			FileUtils.checkFile(xmlFile);
			
			// xml version
			String xmlVersion = "separate";
			if(isValidParam("xml-version","separate","joint")) xmlVersion = commandLine.getOptionValue("xml-version");			
			
			// blast Hits
			int blastHits = 20;
			if(isValidParam("blast-hits")) blastHits = Integer.parseInt(commandLine.getOptionValue("blast-hits"));
			
			// bh-min-region-length
			int bhMinRegionLength = 33;
			if(isValidParam("bh-min-region-length")) bhMinRegionLength = Integer.parseInt(commandLine.getOptionValue("bh-min-region-length"));
							
			// description filter
			String bhDescriptionFilter = "";
			if(isValidParam("bh-description-filter")) bhDescriptionFilter = commandLine.getOptionValue("bh-description-filter");
			
			// bh-min-region-length
			int bhDescriptionPosition = 5;
			if(isValidParam("bh-description-position")) bhDescriptionPosition = Integer.parseInt(commandLine.getOptionValue("bh-description-position"));
			
			// add id
			boolean addId = false;
			if(commandLine.hasOption("add-id")) addId = true;
			
			// use descriptor annotation
			boolean useDescriptorAnnotator = false;
			if(commandLine.hasOption("use-descriptor-annotator")) useDescriptorAnnotator = true;
			
			// e-value
			int eValue = 6;
			if(isValidParam("e-value")) eValue = Integer.parseInt(commandLine.getOptionValue("e-value"));
						
			// annotation cutoff
			int annotationCutOff = 55;
			if(isValidParam("annotation-cutoff")) annotationCutOff = Integer.parseInt(commandLine.getOptionValue("annotation-cutoff"));
			
			// e-value
			int goWeight = 5;
			if(isValidParam("go-weight")) goWeight = Integer.parseInt(commandLine.getOptionValue("go-weight"));
			
			// e-value
			int hitCoverage = 0;
			if(isValidParam("hit-coverage")) hitCoverage = Integer.parseInt(commandLine.getOptionValue("hit-coverage"));
			
			// ec codes
			Hashtable<String, String> ecs = new Hashtable<String, String>();
			ecs.put("ida", commandLine.getOptionValue("ida","1"));
			ecs.put("ipi", commandLine.getOptionValue("ipi","1"));
			ecs.put("imp", commandLine.getOptionValue("imp","1"));
			ecs.put("igi", commandLine.getOptionValue("igi","1"));			
			ecs.put("iep", commandLine.getOptionValue("iep","1"));
			ecs.put("exp", commandLine.getOptionValue("exp","1"));
			
			ecs.put("iss", commandLine.getOptionValue("iss","0.8"));
			ecs.put("iso", commandLine.getOptionValue("iso","0.8"));
			ecs.put("isa", commandLine.getOptionValue("isa","0.8"));
			ecs.put("ism", commandLine.getOptionValue("ism","0.8"));
			ecs.put("igc", commandLine.getOptionValue("igc","0.7"));
			ecs.put("rca", commandLine.getOptionValue("rca","0.9"));
			
			ecs.put("tas", commandLine.getOptionValue("tas","0.9"));
			ecs.put("nas", commandLine.getOptionValue("nas","0.8"));
			ecs.put("ic", commandLine.getOptionValue("ic","0.9"));
			ecs.put("nd", commandLine.getOptionValue("nd","0.5"));
			ecs.put("iea", commandLine.getOptionValue("iea","0.7"));
			ecs.put("nr", commandLine.getOptionValue("nr","0"));
	
			jobStatus.addStatusMessage("10", "Starting tool");
			
			if(config!=null){
				
				b2Pipe = new B2GPipe(outdir + "/",PREFIX,xmlFileName,blastHits,bhMinRegionLength,bhDescriptionFilter,bhDescriptionPosition,addId,useDescriptorAnnotator,eValue,annotationCutOff,goWeight,hitCoverage);
				b2Pipe.setEvidenceCodes(ecs);				
				b2Pipe.setDataBase(config.getProperty("BLAST2GO_HOST_NAME"),config.getProperty("BLAST2GO_DB_NAME"),config.getProperty("BLAST2GO_DB_USER"), config.getProperty("BLAST2GO_DB_PASSWORD"));
				
				jobStatus.addStatusMessage(""+StringUtils.decimalFormat((Math.random()*20 + 10), "##"), "Executing tool");
				
				b2Pipe.runPipe();
			
				jobStatus.addStatusMessage(""+StringUtils.decimalFormat((Math.random()*10 + 80), "##"), "Saving results");
				
				fillResultXML();
				
				
			} else {
				throw new Exception("Properties file $BABELOMICS_HOME/conf/blast2go.properties cannot be read");
			}

			
		}catch(ParseException pe){
			logger.error(pe.getMessage());
		}catch(NumberFormatException nfe){
			logger.error(nfe.getMessage());
		}catch(Exception e){
			e.printStackTrace();
		}
	
	}

	private boolean isValidParam(String paramName) throws ParseException{
		return isValidParam(paramName,null);
	}
	
	private boolean isValidParam(String paramName,String...validOptions) throws ParseException{
		if(commandLine.hasOption(paramName) && !(commandLine.getOptionValue(paramName).trim().length()==0) && !commandLine.getOptionValue(paramName).trim().equalsIgnoreCase("null")) {
			if(validOptions!=null){
				boolean contains = false;
				for(String validOption: validOptions){
					if(validOption.equals(commandLine.getOptionValue(paramName).trim())){
						contains = true;
					}
				}
				if(!contains) throw new ParseException("invalid value for parameter " + paramName);
				return true;
			} else {
				return true;	
			}			
		}
		return false;
	}
	
	
	private void fillResultXML() {
		
		// RESULT FILES 
		result.addOutputItem(new Item("annotation", PREFIX + "_output.annot", "Annotations (.txt)", Item.TYPE.FILE, Arrays.asList("DATA","annotation","GO"), new HashMap<String, String>(), "Result files"));
		result.addOutputItem(new Item("resultTable", PREFIX + "_resultTable.xml", "Result table", Item.TYPE.FILE, "Result files"));
		// BLAST2GO FILES
//		result.addOutputItem(new Item("property_file", PREFIX + "_b2gPipe.properties", "Blast2GO properties file", Item.TYPE.FILE, "Blast2GO files"));
		result.addOutputItem(new Item("project",PREFIX + "_output.dat", "Blast2GO project file (.dat)", Item.TYPE.FILE, "Blast2GO files"));
		
		// IMAGES
		  // Blast results
		result.addOutputItem(new Item("DbDist", PREFIX + "_DbDist.png", "Source database distribution", Item.TYPE.IMAGE, "Images.Blast result"));
		result.addOutputItem(new Item("simDist", PREFIX + "_simDist.png", "Sequence similarity distribution", Item.TYPE.IMAGE, "Images.Blast result"));
		result.addOutputItem(new Item("eValueDist", PREFIX + "_eValueDist.png", "e-Value distribution", Item.TYPE.IMAGE, "Images.Blast result"));
		  // Mapping results 
		result.addOutputItem(new Item("ecDistSeq", PREFIX + "_ecDistSeq.png", "Sequence evidence code distribution", Item.TYPE.IMAGE, "Images.Mapping results"));
		result.addOutputItem(new Item("ecDistHit", PREFIX + "_ecDistHit.png", "Hit evidence code distribution", Item.TYPE.IMAGE, "Images.Mapping results"));
		result.addOutputItem(new Item("speciesDist", PREFIX + "_speciesDist.png", "Species distribution", Item.TYPE.IMAGE, "Images.Mapping results"));
		  // Annotation results
		result.addOutputItem(new Item("dataStat", PREFIX + "_dataStat.png", "General performance", Item.TYPE.IMAGE, "Images.Annotation results"));
		result.addOutputItem(new Item("annotStat", PREFIX + "_annotStat.png", "Number of annotations per sequence", Item.TYPE.IMAGE, "Images.Annotation results"));
		result.addOutputItem(new Item("levelStat", PREFIX + "_levelStat.png", "Gene Ontology term level distribution", Item.TYPE.IMAGE, "Images.Annotation results"));
		result.addOutputItem(new Item("goCountBp", PREFIX + "_goCountBp.png", "Number of sequences per gene ontology term (molecular function)", Item.TYPE.IMAGE, "Images.Annotation results"));
		result.addOutputItem(new Item("goCountMf", PREFIX + "_goCountMf.png", "Number of sequences per gene ontology term (biological process)", Item.TYPE.IMAGE, "Images.Annotation results"));		
		result.addOutputItem(new Item("goCountCc", PREFIX + "_goCountCc.png", "Number of sequences per gene ontology term (cellular component)", Item.TYPE.IMAGE, "Images.Annotation results"));

	}
	

	private Properties loadConfig(){
		config = new Properties();
		try {
			config.load(new FileInputStream(new File(babelomicsHomePath + "/conf/blast2go.properties")));
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		} 
		return config;
	}
	
}

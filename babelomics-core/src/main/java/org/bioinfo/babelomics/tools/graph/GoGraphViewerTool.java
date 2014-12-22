package org.bioinfo.babelomics.tools.graph;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.security.InvalidParameterException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Properties;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

import es.blast2go.prog.graph.GetGraphApi;
import es.blast2go.prog.graph.GoGraphException;


public class GoGraphViewerTool  extends BabelomicsTool{

	// config params
	private Properties config;
	public final static String PREFIX = "b2g";
	
	// tool params
	private String associationFile;
	private String goDomain;
	private int idsPerNodeFilter;
	private int annotScoreNodeFilter;
	private double annotScoreParameter;
	private String graphColoring;
		
	@Override
	public void initOptions() {				
		options.addOption(OptionFactory.createOption("association-file", "the association file containing db ids"));		
		// optional
		options.addOption(OptionFactory.createOption("go-domain", "[b,m,c] biological process, molecular function, cellular component (b by default)",false,true));
		options.addOption(OptionFactory.createOption("ids-per-node-filter", "[0,1000] Number of IDs per node filter (0 by default)",false,true));
		options.addOption(OptionFactory.createOption("annot-score-node-filter", "[0,1000] Thin out graphs by the annot-score (0 by default)",false,true));
		options.addOption(OptionFactory.createOption("annot-score-parameter", "[0,1] Alpha of the annot-score formula (0.6 by default)",false,true));
		options.addOption(OptionFactory.createOption("graph-coloring", "[byScore,bySeqCount,byDesc] Mode of graph coloring",false,true));		
	}

	@Override
	protected void execute() {		

		// load config
		loadConfig();
		
		try {
			
			// association file		
			associationFile = commandLine.getOptionValue("association-file");
			if(!new File(associationFile).exists()) {
				throw new FileNotFoundException(associationFile + " not found");
			} else {
				File input = new File(associationFile);
				File output = new File(outdir + "/" +  new File(associationFile).getName());				
				FileUtils.touch(output);
				FileUtils.copy(input,output);
				associationFile = new File(associationFile).getName();
			}
			
			// go domain
			goDomain = "b";
			if(commandLine.hasOption("go-domain")) goDomain = commandLine.getOptionValue("go-domain");
			if(!goDomain.equals("b") && !goDomain.equals("m") && !goDomain.equals("c")) throw new InvalidParameterException("value " + goDomain + " for go-domain parameter is not valid");
	
			// ids per node filter
			idsPerNodeFilter = 0;		
			if(commandLine.hasOption("ids-per-node-filter") && !commandLine.getOptionValue("ids-per-node-filter").trim().equals("")) idsPerNodeFilter = Integer.parseInt(commandLine.getOptionValue("ids-per-node-filter"));
			if(idsPerNodeFilter<0 | idsPerNodeFilter>1000) throw new InvalidParameterException("value " + idsPerNodeFilter + " for ids-per-node-filter parameter is out of range");
			
			// annots-score node filter
			annotScoreNodeFilter = 0;		
			if(commandLine.hasOption("annot-score-node-filter") && !commandLine.getOptionValue("annot-score-node-filter").trim().equals("")) annotScoreNodeFilter = Integer.parseInt(commandLine.getOptionValue("annot-score-node-filter"));
			if(annotScoreNodeFilter<0 | annotScoreNodeFilter>1000) throw new InvalidParameterException("value " + annotScoreNodeFilter + " for annot-score-node-filter parameter is out of range");
			
			// annots-score node filter
			annotScoreParameter = 0.6;
			if(commandLine.hasOption("annot-score-parameter") && !commandLine.getOptionValue("annot-score-parameter").trim().equals("")) annotScoreParameter = Double.parseDouble(commandLine.getOptionValue("annot-score-parameter"));
			if(annotScoreParameter<0 | annotScoreParameter>1) throw new InvalidParameterException("value " + annotScoreParameter + " for annot-score-parameter parameter is out of range");
			
			//
			graphColoring = "byScore";
			if(commandLine.hasOption("graph-coloring")) graphColoring = commandLine.getOptionValue("graph-coloring");
			if(!graphColoring.equals("byScore") && !graphColoring.equals("bySeqCount") &&  !graphColoring.equals("byDesc")) throw new InvalidParameterException("value " + graphColoring + " for graph-coloring parameter is invalid");
			
			
			if(config!=null){
				
				// init graph api
				GetGraphApi graph = new GetGraphApi(outdir + "/",PREFIX,new File(associationFile).getName(),goDomain,idsPerNodeFilter,graphColoring,annotScoreParameter,annotScoreNodeFilter,"orange",12,"");
				graph.setDownloader(config.getProperty("JNLP_DOWNLOADER_HOST_NAME"));				
				graph.setDataBase(config.getProperty("BLAST2GO_HOST_NAME"),config.getProperty("BLAST2GO_DB_NAME"),config.getProperty("BLAST2GO_DB_USER"), config.getProperty("BLAST2GO_DB_PASSWORD"));
				
				try {
					// run
					graph.run();
					
					// fill output
					fillResultXML();

				  // continuar normal
				} catch (GoGraphException e) {
					// fill output
					fillFailedResultXML(e.getMessage());
				}
					
				
			} else {
				throw new Exception("Properties file $BABELOMICS_HOME/conf/blast2go.properties cannot be read");
			}
					
		} catch (Exception e){			
			e.printStackTrace();
		}
			
		
	}
		
	private void fillResultXML() {
		// images
		result.addOutputItem(new Item("imageJpgsmall", PREFIX + "_graphimagesmall.jpg", "Graph as JPG image (low resolution)", Item.TYPE.IMAGE, "Images"));
		result.addOutputItem(new Item("imageJpg", PREFIX + "_graphimage.jpg", "Graph as JPG image (high resolution)", Item.TYPE.FILE, "Images"));
		result.addOutputItem(new Item("imagePng", PREFIX + "_graphimage.png", "Graph as PNG", Item.TYPE.FILE, "Images"));
		result.addOutputItem(new Item("imageSVG", PREFIX + "_graphimage.svg", "Graph as SVG", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(), "Images"));		
		// jnlp link
		result.addOutputItem(new Item("jnlp_link", PREFIX + "_graphimage.svg", "Click here to open the graph with GOGraphViz (Java Web-Start Application)", Item.TYPE.FILE, Arrays.asList("GO_GRAPH_VIZ_JNLP"), new HashMap<String, String>(), "Interactive Graph Visualisation"));
		// graph as text		
		result.addOutputItem(new Item("graphTxt", PREFIX + "_graph.txt", "Graph as textual representation",  Item.TYPE.FILE, "Images"));
	}

	private void fillFailedResultXML(String message) {
		result.addOutputItem(new Item("error",message,"Error",Item.TYPE.MESSAGE,Arrays.asList("ERROR"), new HashMap<String,String>(),"Images"));
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

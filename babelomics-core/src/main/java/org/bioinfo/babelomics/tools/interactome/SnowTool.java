package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.Node;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.stats.inference.KolmogorovSmirnovTest;
import org.bioinfo.math.util.MathUtils;
import org.bioinfo.networks.protein.InteractomeParser;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinNetworkToFile;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public abstract class SnowTool extends BabelomicsTool{

	protected String interactome;
	protected String outputFileName;
	protected String type;
	protected String group;
	protected int randoms;
	protected boolean components;
	protected boolean bicomponents;
	protected boolean intermediate;
	protected ProteinNetwork proteinNetwork;
	protected List<ProteinNetwork> subProteinNetworkRandoms;
	
	
	public SnowTool(){
		initSnowToolOptions();
	}
	
	public void initSnowToolOptions() {
		options.addOption(OptionFactory.createOption("interactome", "Interactome: hsa, sce or own (for you own interactome)", true, true));
		//options.addOption(OptionFactory.createOption("type", "An argument saying if you want genes, proteins(default) or transcripts", false, true));
		///** List tags lo manda la web, por lo tanto no hace falta que use type, sin embargo si se usa **/
		options.addOption(OptionFactory.createOption("list-tags", "An argument saying if you want genes, proteins(default) or transcripts", false, true));
		options.addOption(OptionFactory.createOption("components", "If we want the number of components we put 1, otherwise we put 0", false, true));
		options.addOption(OptionFactory.createOption("bicomponents", "If we want the number of bicomponents we put 1, otherwise we put 0", false, true));
		options.addOption(OptionFactory.createOption("o-name", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		options.addOption(OptionFactory.createOption("group", "Values all(by default) or curated. It is An argument saying whether you want a curated interactome or the whole interactome", true, true));
		options.addOption(OptionFactory.createOption("intermediate", "If we want the intermediate we put 1, otherwise we put 0", false, true));
		options.addOption(OptionFactory.createOption("randoms", "Number of randoms", false, true));
		//options.addOption(OptionFactory.createOption("side", "side for kolmogorov and wilkoxon test. Can be less or greater", false, true));

	}
	
	public void initExecute(){
		
		interactome = commandLine.getOptionValue("interactome");
		if(commandLine.hasOption("randoms"))
			randoms = Integer.parseInt(commandLine.getOptionValue("randoms"));
		
		//type = !commandLine.hasOption("type") ? "proteins":commandLine.getOptionValue("type");
		String listTags = !commandLine.hasOption("list-tags") ? "proteins":commandLine.getOptionValue("list-tags");
		String listTagsArray[] = listTags.split(",");
		for (int i = 0; i < listTagsArray.length; i++) {
			if(listTagsArray[i].equalsIgnoreCase("protein") || listTagsArray[i].equalsIgnoreCase("transcript"))
				type = "proteins";
			if(listTagsArray[i].equalsIgnoreCase("gene"))
				type = "genes";
			if(listTagsArray[i].equalsIgnoreCase("ranked"))
				type = "genes";
			if(listTagsArray[i].contains("vcf"))
				type = "vcf";
		}
		
		//type = !commandLine.hasOption("type") ? "proteins":commandLine.getOptionValue("type");
		group = !commandLine.hasOption("group") ? "all": commandLine.getOptionValue("group");
		components = (commandLine.hasOption("components") && !"0".equalsIgnoreCase(commandLine.getOptionValue("components")));
		bicomponents = (commandLine.hasOption("bicomponents") && !"0".equalsIgnoreCase(commandLine.getOptionValue("bicomponents")));
		intermediate = (commandLine.hasOption("intermediate") && !"0".equalsIgnoreCase(commandLine.getOptionValue("intermediate")));
		//side = !commandLine.hasOption("side") ? "less" : commandLine.getOptionValue("side");
		outputFileName = outdir + "/" + commandLine.getOptionValue("o-name", "result");
		loadInteractome();
	}
	public void loadInteractome(){
		ProteinNetworkToFile file = new ProteinNetworkToFile();
		String folderInteractions = this.babelomicsHomePath + "/conf/interactions/";
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = null;
		try {
			File inputFile;
			if (interactome.equals("own")){
				if(commandLine.hasOption("sif-file") ){
					inputFile = new File(commandLine.getOptionValue("sif-file"));
					logger.debug("checking if inputFile exist...");
					FileUtils.checkFile(inputFile);
					interactomeGraph = InteractomeParser.parseFromSifFile(commandLine.getOptionValue("sif-file") );
					proteinNetwork = new ProteinNetwork(interactomeGraph);
	
					if(commandLine.hasOption("topo-file"))
						proteinNetwork.loadTopologicalValues(commandLine.getOptionValue("topo-file"));
					else
						proteinNetwork.calcTopologicalValues();
					if(commandLine.hasOption("o-sif-topo-file")) {
						file.toTopologicalFile(outputFileName+"_topo.txt", proteinNetwork);
						System.out.println("File created...");
						return;
					}
				} else {
					System.err.println("Missing custom interactome");
					return;					
				}
			} else {
				//you have to be very strict with the name of the files, otherwise the program won't work fine
				String localType = "";
				if(type.equalsIgnoreCase("genes") || type.equalsIgnoreCase("vcf"))
					localType = "genes";
				if(type.equalsIgnoreCase("proteins") || type.equalsIgnoreCase("transcripts"))
					localType = "proteins";
				String parseFromSifFile = folderInteractions+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr.sif";
				interactomeGraph = InteractomeParser.parseFromSifFile(parseFromSifFile);
				System.out.println("Interactome read: "+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr.sif");
				proteinNetwork = new ProteinNetwork(interactomeGraph);
				proteinNetwork.loadTopologicalValues(folderInteractions+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr_topo.txt");
				System.out.println("Interactome topo values read: "+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr_topo.txt");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public String getInteractomeMsg() {
		String interactomeMsg = "Homo sapiens";
		if ("hsa".equalsIgnoreCase(interactome) ) {
			interactomeMsg = "Homo sapiens";
		} else if ("sce".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Saccharomyce cerevisiae";
		} else if ("bta".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Bos taurus";
		} else if ("dme".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Drosophila melanogaster";
		} else if ("mmu".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Mus musculus";
		} else if ("ath".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Arabidopsis thaliana";
		} else if ("eco".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Escherichia coli";
		} else if ("own".equalsIgnoreCase(interactome)) {
			interactomeMsg = new File(commandLine.getOptionValue("sif-file")).getName();
		} else {
			interactomeMsg = "unknown";
		}
			return interactomeMsg;
	}
	public StringBuilder deleteLastLine(StringBuilder sb){
		return this.deleteLastCh(sb, System.getProperty("line.separator"));
//		if(!sb.toString().equals(""))
//			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
//		return sb.toString();
	}
	public StringBuilder deleteLastCh(StringBuilder sb, String ch){
		StringBuilder sbReturn = new StringBuilder();
		if(!sb.toString().equals("")){
			//System.out.println(sb);
			sb.deleteCharAt(sb.lastIndexOf(ch));
		}
		return sbReturn.append(sb.toString());
	}
	
	public List<ProteinVertex> nodes2ProteinVertices(List<Node> nodes){
		List<ProteinVertex> list = new ArrayList<ProteinVertex>();
		for(Node node : nodes)
			list.add(new ProteinVertex(node.getId()));
		return list;
	}
	public void addOutputAppletItem(File xmlFile, int index) {
		if (xmlFile.exists()) {
			String url = "SnowViewer?filename=" + xmlFile.getName() + "&width=600&height=600";
		//	result.addOutputItem(new Item("viewer" + index + "_param", url, "Viewer for network #" + index, TYPE.HTML, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Minimun Connected Network description"));
			result.addOutputItem(new Item("viewer" + index + "_param", url, "Viewer for network", TYPE.HTML, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Minimun Connected Network description"));

			
			url = "SnowViewer?filename=" + xmlFile.getName();
//			result.addOutputItem(new Item("viewer" + index + "_param_new_window", url, "Open applet for network #" + index + " in a new window", TYPE.LINK, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Network viewer"));
			result.addOutputItem(new Item("viewer" + index + "_param_new_window", url, "Open applet for network  in a new window", TYPE.LINK, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Minimun Connected Network description"));

		}
	}
	public String ksTest(List<Double> relBetList1, List<Double> relBetList2, List<Double> connList1, List<Double> connList2, List<Double> clustList1, List<Double> clustList2/*, String side*/) throws IOException{
		if( relBetList1.size() == 0 || relBetList2.size() == 0 || connList1.size() == 0 || connList2.size() == 0 || clustList1.size() == 0 || clustList2.size() == 0)
			return "";
		StringBuilder sb = new StringBuilder();
		sb.append("#parameter\tpval\tside\n");
		KolmogorovSmirnovTestResult result;
		
		result = getPValue(relBetList1, relBetList2/*, side*/);
		sb.append("betweenness\t"+result.getPValue()+"\t"+result.getSide()+"\n");
		result = getPValue(connList1, connList2/*, side*/);
		sb.append("connections\t"+result.getPValue()+"\t"+result.getSide()+"\n");
		result = getPValue(clustList1,clustList2/*, side*/);
		sb.append("coefficient\t"+result.getPValue()+"\t"+result.getSide());
		return sb.toString();
	}
	public KolmogorovSmirnovTestResult getPValue(List<Double> list1, List<Double> list2/*, String side*/) throws IOException{
//		KSTest kstest = new KSTest();
//		return kstest.resultKolmogorovSmirnov(ListUtils.toDoubleArray(list1), ListUtils.toDoubleArray(list2), side).getPValue();
		KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();
		return  ksTest.compute(ListUtils.toDoubleArray(list1), ListUtils.toDoubleArray(list2));
	}
//	public String getSymbol(double n1, double n2){
//		String symbol = "<";
//		if(n1 > n2)
//			symbol = ">";
//		return symbol;
//	}
	public String getSymbol(String side){
		String symbol = "<";
		if(side.equals("greater"))
			symbol = ">";
		return symbol;
	}
	public void createImages(String fileName, List<Double> list1, String legend1, List<Double> list2, String legend2, String itemName, String itemLabel, String itemGroup, String xAxis) throws IOException {
		File f = new File(fileName);
		BoxPlotChart bpc = new BoxPlotChart("", "", "");
		bpc.getLegend().setVisible(true);
		if(!list2.isEmpty())
			bpc.addSeries(list2, legend2, legend2);
		if(!list1.isEmpty())
			bpc.addSeries(list1, legend1, legend1);
		bpc.save(f.getAbsolutePath()+".png", 300, 250, "png");
		result.addOutputItem(new Item(itemName, f.getName()+".png", itemLabel, TYPE.IMAGE, new ArrayList<String>(), new HashMap<String, String>(2), itemGroup));		
	
	}
	public int componentsMoreThanOneNode(List<List<ProteinVertex>> components){
		int compMoreThan1Node = 0;
		for(List<ProteinVertex> comp : components){
			if(comp.size() > 1)
				compMoreThan1Node++;
		}
		return compMoreThan1Node;
	}
	public int[] getRange(double []componentsRandoms){
		int[] resultRange = new int[2]; 
		resultRange[0] = (int)MathUtils.percentile(componentsRandoms, 2.5);
		resultRange[1] = (int)MathUtils.percentile(componentsRandoms, 97.5);
		return resultRange;
//		range = "["+Double.toString(minRange)+","+Double.toString(maxRange)+"]";
//		result.addOutputItem(new Item("components_number",  range, "Number of components [95% confidence interval]", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network description. Topology description"));
		
	}
	
//	protected void createJson(ProteinNetwork proteinNetwork, String sourceDotFile, List<List<ProteinVertex>> componentsListSub, Set<String> intermediatesSub, int node, Map<String, String> mapList) throws IOException{
//
//		Dot<ProteinVertex, DefaultEdge> dot = new Dot<ProteinVertex, DefaultEdge>();
//		Svg svg = new Svg();
//		Json<ProteinVertex, DefaultEdge> json = new Json<ProteinVertex, DefaultEdge>();
//		List<File> fileList = new ArrayList<File>();
//		List<String> layoutsName = new ArrayList<String>();
//
//		if(componentsListSub == null)
//			logger.error("not components calculated, please, set the --components option");
//
//		IOUtils.write(sourceDotFile, dot.toDot(proteinNetwork.getInteractomeGraph()));
//
//		File neatoFile = new File(outputFileName+"_list"+node+"_dot.svg");
//		IOUtils.write(neatoFile, svg.toSvg(sourceDotFile, "neato"));
//		fileList.add(neatoFile);
//		layoutsName.add("neato");
//
//		File f = new File(outputFileName+"_list"+node+".json");
//		IOUtils.write(f.getAbsoluteFile(), json.toJson(this.interactome, fileList, layoutsName, proteinNetwork.getInteractomeGraph(), intermediatesSub, componentsListSub, mapList));
//	}
	protected void addOutputSvgViewer(File file, int index, Set<String> intermediates, Set<String> seedList) {
		List<String> tags = new ArrayList<String>();
		String intermediateTag = "";
		String seedTag = "";
		for(String intermediate : intermediates)
			intermediateTag+=intermediate+"|";
		for(String seed : seedList)
			seedTag+=seed+"|";
		
		/** Borro la última coma **/
		if(!intermediateTag.equals("")){
			intermediateTag = intermediateTag.substring(0,intermediateTag.length()-1);
		}
		/** Borro la última coma **/
		if(!seedTag.equals("")){
			seedTag = seedTag.substring(0,seedTag.length()-1);
		}
		
		tags.add("INTERACTOME_VIEWER");
		//tags.add(this.type);
		//tags.add("REDIRECT_TOOL");
		if (file.exists()) {
			String list = "list";
			if(index==1)
				list+="1";
			if(index==2)
				list+="2";
			tags.add(list);
			tags.add(this.interactome);
			tags.add(intermediateTag);
			tags.add(seedTag);
			result.addOutputItem(new Item("svg_viewer" + index + "_param", file.getName(), "Minimun Connected Network interactions", TYPE.FILE, tags, new HashMap<String, String>(2), "Results: Minimum Connected Network selected.Viewer"));

//			url = "SnowViewer2?filename=" + jsonFile.getName();
//			result.addOutputItem(new Item("svg_viewer" + index + "_param_new_window", url, "Open svg for network #" + index + " in a new window", TYPE.LINK, StringUtils.toList("SERVER,INCLUDE_REFS,SNOW", ","), new HashMap<String, String>(2), "Network viewer 2"));
		}
	}
	
	
}


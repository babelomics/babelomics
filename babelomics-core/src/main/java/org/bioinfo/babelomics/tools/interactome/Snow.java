package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.ListInfo;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.Node;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.Subgraph;
import org.bioinfo.data.graph.alg.Calc;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.feature.DBName;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.infrared.core.feature.XRef.XRefItem;
import org.bioinfo.infrared.core.funcannot.GO;
import org.bioinfo.infrared.funcannot.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.GODBManager;
import org.bioinfo.math.util.MathUtils;
import org.bioinfo.networks.protein.InteractomeParser;
import org.bioinfo.networks.protein.KSTest;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinNetworkToFile;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.networks.protein.files.Dot;
import org.bioinfo.networks.protein.files.Json;
import org.bioinfo.networks.protein.files.Svg;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Snow  extends BabelomicsTool{

	private File inputFile;
	private String interactome;
	private String outputFileName;
	private String type;
	private String group;
	private ProteinNetwork proteinNetwork;
	private ProteinNetwork subProteinNetwork1, subProteinNetwork2;
	private List<ProteinNetwork> subProteinNetworkRandoms;
	private boolean components;
	private boolean bicomponents;
	private boolean intermediate;
	private boolean images;
	private boolean json;
	//private boolean xml;
	private boolean sif;
	private Set<String> intermediatesSub1, intermediatesSub2; 
	private List<List<ProteinVertex>> componentsListSub1, componentsListSub2;
	private int bicomponentsNumberList1, bicomponentsNumberList2;
	private int randomSize;
	/** ensembl_gene, input_id**/
	private Map<String, String> mapList1, mapList2;
	
	private DBConnector dbConnector;
	private XRefDBManager xrefDBMan;
	private String decimalFormat;
	private int listMaxSize;
	
	private DBName dbName;

	//private String wBinPath;

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("interactome", "Interactome: hsa, sce or own (for you own interactome)", true));
		options.addOption(OptionFactory.createOption("sif-file", "An input file containing a SIF interactome", false));
		options.addOption(OptionFactory.createOption("topo-file","t", "An input file containing topological values", false, true));
		options.addOption(OptionFactory.createOption("list1", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("list2", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("type", "An argument saying if you want genes, proteins(default) or transcripts", false, true));
		

		options.addOption(OptionFactory.createOption("randoms", "Number of randoms", false, true));
		options.addOption(OptionFactory.createOption("intermediate", "If we want the intermediate we put 1, otherwise we put 0", false, true));
		options.addOption(OptionFactory.createOption("components", "If we want the number of components we put 1, otherwise we put 0", false, false));
		options.addOption(OptionFactory.createOption("bicomponents", "If we want the number of bicomponents we put 1, otherwise we put 0", false, false));
		options.addOption(OptionFactory.createOption("json", "It will create an output .json file", false, false));
		options.addOption(OptionFactory.createOption("o-sif-topo-file", "Create a full topological file from a sif file", false, false));
		options.addOption(OptionFactory.createOption("o-name", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		options.addOption(OptionFactory.createOption("side", "side for kolmogorov test. Can be less or greater", false, true));
		options.addOption(OptionFactory.createOption("images", "Print the images for the statistics", false, false));
		
		//options.addOption(OptionFactory.createOption("xml", "Output xml file with the representation of the graph", false, false));
		//options.addOption(OptionFactory.createOption("json", "Output json file with the representation of the graph", false, false));

		options.addOption(OptionFactory.createOption("sif", "Output sif file with the representation of the graph", false, false));
		options.addOption(OptionFactory.createOption("group", "Values all(by default) or curated. It is An argument saying whether you want a curated interactome or the whole interactome", false, true));


	}
	@Override
	protected void execute() {
		File f = null;
		this.decimalFormat= "#.####";
		this.listMaxSize = 500;
		List<String> listToVertex1 = new ArrayList<String>();
		List<String> listToVertex2 = new ArrayList<String>();
		mapList1 = new HashMap<String, String>();
		mapList2 = new HashMap<String, String>();
		
//		wBinPath = babelomicsHomePath + "/bin/snow/wilcoxtest.r";
		subProteinNetwork1 = null;
		subProteinNetwork2 = null;

		subProteinNetworkRandoms = new ArrayList<ProteinNetwork>();
		components = commandLine.hasOption("components");
		bicomponents = commandLine.hasOption("bicomponents");
		images = commandLine.hasOption("images");
		json = commandLine.hasOption("json");
	//	xml = commandLine.hasOption("xml");
		sif = commandLine.hasOption("sif");

		//xml = true;
//		json = commandLine.hasOption("json");
		json = true;
		outputFileName = outdir + "/" + commandLine.getOptionValue("o-name", "result");
		
		
		interactome = commandLine.getOptionValue("interactome");
		String interactomeMsg = getInteractomeMsg();
		dbConnector = new DBConnector(interactome, new File(babelomicsHomePath + "/conf/infrared.properties"));
		xrefDBMan = new XRefDBManager(dbConnector);
		
		
		
		result.addOutputItem(new Item("interactome_param", interactomeMsg, "Species", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		String side = !commandLine.hasOption("side") ? "less" : commandLine.getOptionValue("side");
		result.addOutputItem(new Item("side_param", side, "Side for kolmogorov test", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		type = !commandLine.hasOption("type") ? "proteins":commandLine.getOptionValue("type");
		result.addOutputItem(new Item("type_param", type, "ID type", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		group = !commandLine.hasOption("group") ? "all": commandLine.getOptionValue("group");
		String groupMsg = group.equalsIgnoreCase("all") ? "all ppis" : "ppis detected by at least two methods (curated)";
		result.addOutputItem(new Item("group_param", groupMsg, "Group interaction", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		
		intermediate = (commandLine.hasOption("intermediate") && !"0".equalsIgnoreCase(commandLine.getOptionValue("intermediate")));

		ProteinNetworkToFile file = new ProteinNetworkToFile();

		intermediatesSub1 = new HashSet<String>();
		intermediatesSub2 = new HashSet<String>();
		try {

			String folderInteractions = this.babelomicsHomePath + "/conf/interactions/";
			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = null;

			
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
				String localType = "genes";
				if(!type.equalsIgnoreCase("genes"))
					localType = "proteins";
				String parseFromSifFile = folderInteractions+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr.sif";
				interactomeGraph = InteractomeParser.parseFromSifFile(parseFromSifFile);
				System.out.println("Interactome read: "+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr.sif");
				proteinNetwork = new ProteinNetwork(interactomeGraph);
				proteinNetwork.loadTopologicalValues(folderInteractions+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr_topo.txt");
				System.out.println("Interactome topo values read: "+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr_topo.txt");
//				System.err.println("Unknown interactome");
//				return;
			}
			
			
			//StringBuilder sbMeans = createMeansHeader();
			StringBuilder sbTopo = createTopoHeader();
			StringBuilder sbExternalTopo = createExternalTopoHeader();
			StringBuilder sbComponents = createComponentsHeader();

			if(commandLine.hasOption("list1")) {
				logger.debug("Starting list1.........");
				String nodeFile = commandLine.getOptionValue("list1");
				int node=1;
				FileUtils.checkFile(new File(nodeFile));
				List<String> list1 = IOUtils.readLines(nodeFile);
				
				
				listToVertex1 = new ArrayList<String>();
				if(type.equalsIgnoreCase("transcripts") || type.equalsIgnoreCase("proteins")){
					this.mapList1 = transcriptToUniprot(list1);
					listToVertex1.addAll(mapList1.keySet());
				}
				else if(type.equalsIgnoreCase("genes")){
					this.mapList1 = getGenEnsemble(list1);
					listToVertex1.addAll(mapList1.keySet());
				}
				if(this.mapList1.isEmpty())
					listToVertex1 = list1;
				
				if(listToVertex1.size() > listMaxSize){
					result.addOutputItem(new Item("list1_too_big", "The list is too big, only considered the 500 first on the interactome ", "List", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
					listToVertex1 = listToVertex1.subList(0, listMaxSize);
				}
				String categoryOutput = "";
				if(!commandLine.hasOption("list2") || "none".equalsIgnoreCase(commandLine.getOptionValue("list2"))) {
					categoryOutput = "Network parameters evaluation.Minimal Connected Network topological evaluation";
				}
				else if(commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2"))) {
					categoryOutput = "Network parameters evaluation.Minimal Connected Network topological evaluation for List 1";
				}
				logger.debug("nodes read: " + list1.toString());
				//System.out.println("nodes read: " + listToVertex1.toString());

				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toVertex(listToVertex1));
				System.out.println("Before intermediate: "+subgraph.getVertices().size()+" nodes");
				this.randomSize = subgraph.getVertices().size();
				dbName = getOriginalDbName(mapList1.values());
				
				if(intermediate){
					intermediatesSub1 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
					
					for (String id : intermediatesSub1) {
						mapList1.put(id, getIdToDbName(id));
					}
//					if(intermediatesSub1.size()>0)
//						result.addOutputItem(new Item("external_nodes_list_"+node, intermediatesSub1.toString().substring(1, intermediatesSub1.toString().length()-1), "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));
//					else
//						result.addOutputItem(new Item("external_nodes_list_"+node, "No added", "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
					System.out.println("With intermediate: "+subgraph.getVertices().size()+" nodes");
				}
				
				subProteinNetwork1 = createSubnet(subgraph);

//				sbTopo.append(getTopologicalValues(subProteinNetwork1, node, false));
//				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
//				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topological values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));
//
//				sbExternalTopo.append(getExternalTopologicalValues(subProteinNetwork1, 1, this.intermediatesSub1));
//				f = new File(outputFileName+"_sn_external"+node+"_topo.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbExternalTopo.toString());
//				result.addOutputItem(new Item("sn_external"+node+"_topo_param", f.getName(), "External proteins introduced to the network", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));

//				sbMeans.append(getTopologicalMeanValues(subProteinNetwork1, node));
//				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
//				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));


				if(components) {
					logger.debug("Starting list1 components.........");
					System.out.println("Starting list1 components.........");
					componentsListSub1 = subProteinNetwork1.getInteractomeGraph().getAllInformationComponents(true);
						
					if(commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2"))) 
						result.addOutputItem(new Item("components_number_1", Integer.toString(componentsListSub1.size()), "Number of components", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));

					int compMoreThan1Node = componentsMoreThanOneNode(componentsListSub1);
					result.addOutputItem(new Item("comp_more_than_1_node", Integer.toString(compMoreThan1Node), "Number of components with more than 1 node", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));

				}
				else{
					//We assume the subgraph as a component
					componentsListSub1 = new ArrayList<List<ProteinVertex>>();
					componentsListSub1.add(subProteinNetwork1.getInteractomeGraph().getVertices());
				}
				if(bicomponents){
					bicomponentsNumberList1 = subProteinNetwork1.getInteractomeGraph().getNumberOfBicomponents();
					//	result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
					result.addOutputItem(new Item("bicomponents_list1", Integer.toString(bicomponentsNumberList1), "Number of Bicomponents", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));

				}
				if(commandLine.hasOption("randoms") && (!commandLine.hasOption("list2") || "none".equalsIgnoreCase(commandLine.getOptionValue("list2")))) {
					logger.debug("Starting randoms.........");
					int randoms = Integer.parseInt(commandLine.getOptionValue("randoms"));
					if(randomSize > 0)
						createRandoms(randoms, randomSize);
					result.addOutputItem(new Item("randoms_param", ""+randoms, "Number of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
					//result.addOutputItem(new Item("randoms_size_param", ""+randomSize, "Size of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
					statsOneListAnalysis(side);
				}
			
				sbComponents.append(getComponentsValues(subProteinNetwork1, node, componentsListSub1));
				sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_comp.txt");
				IOUtils.write(f.getAbsoluteFile(), sbComponents.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_components_param", f.getName(), "Component values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));
			
				sbTopo.append(getTopologicalValues(subProteinNetwork1, node, false));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topological values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));

				if(this.intermediatesSub1 != null){
					sbExternalTopo.append(getExternalTopologicalValues(subProteinNetwork1, 1, this.intermediatesSub1));
					f = new File(outputFileName+"_sn_external"+node+"_topo.txt");
					IOUtils.write(f.getAbsoluteFile(), sbExternalTopo.toString());
					result.addOutputItem(new Item("sn_external"+node+"_topo_param", f.getName(), "External proteins introduced to the network", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));
				}
//				sbMeans.append(getTopologicalMeanValues(subProteinNetwork1, node, componentsListSub1, bicomponentsNumberList1));
//				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
//				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));
				

			}

			if(commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2"))) {
				//sbMeans = createMeansHeader();
				sbTopo = createTopoHeader();
				sbComponents = createComponentsHeader();
				logger.debug("Starting list2.........");
				String nodeFile = commandLine.getOptionValue("list2");
				int node=2;
				FileUtils.checkFile(new File(nodeFile));
				List<String> list2 = IOUtils.readLines(nodeFile);
				
				listToVertex2 = new ArrayList<String>();
				if(type.equalsIgnoreCase("transcripts") || type.equalsIgnoreCase("proteins")){
					this.mapList2 = transcriptToUniprot(list2);
					listToVertex2.addAll(mapList2.keySet());
				}
				else if(type.equalsIgnoreCase("genes")){
					this.mapList2 = getGenEnsemble(list2);
					listToVertex2.addAll(mapList2.keySet());
				}
				if(this.mapList2.isEmpty())
					listToVertex2 = list2;
				
				if(listToVertex2.size() > listMaxSize){
					result.addOutputItem(new Item("list2_too_big", "The list 2 is too big, only considered the 500 first on the interactome ", "List", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
					listToVertex2 = listToVertex2.subList(0, listMaxSize);
				}
				
				logger.debug("nodes read: " + list2.toString());
				String categoryOutput = "Network parameters evaluation.Minimal Connected Network topological evaluation for List 2";
			
				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toVertex(listToVertex2));

				if(intermediate){
					intermediatesSub2 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
					for (String id : intermediatesSub2) {
						mapList2.put(id, getIdToDbName(id));
					}
//					if(intermediatesSub2.size()>0)
//						result.addOutputItem(new Item("external_nodes_list_"+node, intermediatesSub2.toString().substring(1, intermediatesSub2.toString().length()-1), "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));
//					else
//						result.addOutputItem(new Item("external_nodes_list_"+node, "No added", "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
				}
				subProteinNetwork2 = createSubnet(subgraph);

				
//				sbTopo.append(getTopologicalValues(subProteinNetwork2, node, false));
//				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
//				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topografical values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));
//				
//				sbExternalTopo = createExternalTopoHeader();
//				sbExternalTopo.append(getExternalTopologicalValues(subProteinNetwork2, 2, this.intermediatesSub2));
//				f = new File(outputFileName+"_sn_external"+2+"_topo.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbExternalTopo.toString());
//				result.addOutputItem(new Item("sn_external"+2+"_topo_param", f.getName(), "External proteins introduced to the network", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));

				
//				sbMeans.append(getTopologicalMeanValues(subProteinNetwork2, node));
//				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
//				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				if(components) {
					logger.debug("Starting list2 components.........");
					componentsListSub2 = subProteinNetwork2.getInteractomeGraph().getAllInformationComponents(true);
						
					result.addOutputItem(new Item("components_number_2", Integer.toString(componentsListSub2.size()), "Number of components", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));
					int compMoreThan1Node = componentsMoreThanOneNode(componentsListSub2);
					result.addOutputItem(new Item("comp_more_than_1_node_list2", Integer.toString(compMoreThan1Node), "Number of components with more than 1 node", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));
				
				}
				else{
					//We assume the subgraph as a component
					componentsListSub2 = new ArrayList<List<ProteinVertex>>();
					componentsListSub2.add(subProteinNetwork2.getInteractomeGraph().getVertices());
				}
				if(bicomponents){
						bicomponentsNumberList2 = subProteinNetwork2.getInteractomeGraph().getNumberOfBicomponents();
						result.addOutputItem(new Item("bicomponents_list2", Integer.toString(bicomponentsNumberList2), "Number of Bicomponents", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));

				}
				if(!commandLine.hasOption("randoms") && commandLine.hasOption("list1") && commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2"))){
					statsTwoListsAnalisys(side);
				}


				sbComponents.append(getComponentsValues(subProteinNetwork2, node, componentsListSub2));
				sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_comp.txt");
				IOUtils.write(f.getAbsoluteFile(), sbComponents.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_components_param", f.getName(), "Component values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));
			
				sbTopo.append(getTopologicalValues(subProteinNetwork2, node, false));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topografical values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));
				
				if(this.intermediatesSub2 != null){
					sbExternalTopo = createExternalTopoHeader();
					sbExternalTopo.append(getExternalTopologicalValues(subProteinNetwork2, 2, this.intermediatesSub2));
					f = new File(outputFileName+"_sn_external"+2+"_topo.txt");
					IOUtils.write(f.getAbsoluteFile(), sbExternalTopo.toString());
					result.addOutputItem(new Item("sn_external"+2+"_topo_param", f.getName(), "External proteins introduced to the network", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput+".Files"));
				}
				
//				sbMeans.append(getTopologicalMeanValues(subProteinNetwork2, node, componentsListSub2, bicomponentsNumberList2));
//				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
//				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
//				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),categoryOutput));

			}

//			if(commandLine.hasOption("randoms")){
//				logger.debug("Starting randoms.........");
//				int randoms = Integer.parseInt(commandLine.getOptionValue("randoms"));
//				if(randomSize > 0)
//					createRandoms(randoms, randomSize);
//				result.addOutputItem(new Item("randoms_param", ""+randoms, "Number of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//				result.addOutputItem(new Item("randoms_size_param", ""+randomSize, "Size of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//			}

			//result.addOutputItem(new Item("components_param", (commandLine.hasOption("components") ? "yes" : "no"), "Calculate the number of components", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			//result.addOutputItem(new Item("bicomponents_param", (commandLine.hasOption("bicomponents") ? "yes" : "no"), "Calculate the number of bicomponents", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("intermediate_param", (intermediate ? "1" : "0"), "Max. number of external proteins introduced", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

//			if(commandLine.hasOption("randoms") && (!commandLine.hasOption("list2") || "none".equalsIgnoreCase(commandLine.getOptionValue("list2")))) {
//				statsOneListAnalysis(side);
//				System.out.println("adas");
//			} else if(!commandLine.hasOption("randoms") && commandLine.hasOption("list1") && commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2")))
//				statsTwoListsAnalisys(side);
//			else{
//				logger.error("Not correct arguments for statistic test");
//				return;
//			}
				
			if(subProteinNetwork1 != null){
				Set<String> convertedIntermediates = new HashSet<String>();
				for(String key : intermediatesSub1){
					
					convertedIntermediates.add(mapList1.get(key));
				}
				File fSif;
				fSif = new File(outputFileName+"_subnetwork1.sif");
				IOUtils.write(fSif.getAbsoluteFile(), graphToSif(subProteinNetwork1.getInteractomeGraph(), mapList1));
				addOutputSvgViewer(fSif, 1, convertedIntermediates);
				try {
					String mcnInteractors =	getMcnInteractors(subProteinNetwork1, intermediatesSub1, mapList1);
					String googleTable="NETWORKMINERNOTRANKED_TABLE"; 
					f = new File(outputFileName+"_mcn_interactors_list1.txt");
					IOUtils.write(f.getAbsoluteFile(), mcnInteractors);
					result.addOutputItem(new Item("mcn_interactors_list1", f.getName(), "Minimun Connected Network interactors", Item.TYPE.FILE,StringUtils.toList("TABLE,"+googleTable, ",") ,new HashMap<String,String>(),"Interactors list 1"));
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			if(subProteinNetwork2 != null){
				Set<String> convertedIntermediates = new HashSet<String>();
				for(String key : intermediatesSub2){
					convertedIntermediates.add(mapList2.get(key));
				}
				File fSif;
				fSif = new File(outputFileName+"_subnetwork2.sif");
				IOUtils.write(fSif.getAbsoluteFile(), graphToSif(subProteinNetwork2.getInteractomeGraph(), mapList2));
				addOutputSvgViewer(fSif, 2, convertedIntermediates);
			}

			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		}
			
	}

	private Map<String, String> transcriptToUniprot(List<String> list) {
		Map<String, String> map = new HashMap<String, String>();
//		FeatureList<XRef> xrefsEns = new FeatureList<XRef>();
		for(String id : list){
			try {
//				xrefsEns  = xrefDBMan.getByDBName(id, "uniprot_swissprot_accession");
				XRef xrefEns  = xrefDBMan.getByDBName(id, "uniprot_swissprot_accession");
//				if(xrefsEns != null && !xrefsEns.isEmpty() && !xrefsEns.get(0).getId().equals(id))
//					map.put(xrefsEns.get(0).getId(), id);
				if(xrefEns != null && !xrefEns.getXrefItems().get("uniprot_swissprot_accession").isEmpty() /*&& !xrefEns.getXrefItems().get("uniprot_swissprot_accession").get(0).equals(id)*/)
					map.put(xrefEns.getXrefItems().get("uniprot_swissprot_accession").get(0).getDisplayName(), id);
			} catch (Exception e) {
				map.put(id,id);
				e.printStackTrace();
			}
		}
		return map;
	}

	private String getInteractomeMsg() {
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
		} else if ("own".equalsIgnoreCase(interactome)) {
			interactomeMsg = new File(commandLine.getOptionValue("sif-file")).getName();
		} else {
			interactomeMsg = "unknown";
		}
			return interactomeMsg;
	}

	private void createJson(ProteinNetwork proteinNetwork, String sourceDotFile, List<List<ProteinVertex>> componentsListSub, Set<String> intermediatesSub, int node, Map<String, String> mapList) throws IOException{

		Dot<ProteinVertex, DefaultEdge> dot = new Dot<ProteinVertex, DefaultEdge>();
		Svg svg = new Svg();
		Json<ProteinVertex, DefaultEdge> json = new Json<ProteinVertex, DefaultEdge>();
		List<File> fileList = new ArrayList<File>();
		List<String> layoutsName = new ArrayList<String>();

		if(componentsListSub == null)
			logger.error("not components calculated, please, set the --components option");

		IOUtils.write(sourceDotFile, dot.toDot(proteinNetwork.getInteractomeGraph()));

		File neatoFile = new File(outputFileName+"_list"+node+"_dot.svg");
		IOUtils.write(neatoFile, svg.toSvg(sourceDotFile, "neato"));
		fileList.add(neatoFile);
		layoutsName.add("neato");
		
//		File dotFile = new File(outputFileName+"_list"+node+"_dot.svg");
//		IOUtils.write(dotFile, svg.toSvg(sourceDotFile, "dot"));
//		fileList.add(dotFile);
//		layoutsName.add("dot");

//		File twopiFile = new File(outputFileName+"_list"+node+"_twopi.svg");
//		IOUtils.write(twopiFile, svg.toSvg(sourceDotFile, "twopi"));
//		fileList.add(twopiFile);
//		layoutsName.add("twopi");

		File f = new File(outputFileName+"_list"+node+".json");
		IOUtils.write(f.getAbsoluteFile(), json.toJson(this.interactome, fileList, layoutsName, proteinNetwork.getInteractomeGraph(), intermediatesSub, componentsListSub,mapList));
	}

	private void createImages(String fileName, List<Double> list1, String legend1, List<Double> list2, String legend2, String itemName, String itemLabel, String itemGroup, String xAxis) throws IOException {
		File f = new File(fileName);
		BoxPlotChart bpc = new BoxPlotChart("", "", "");
		
		if(!list2.isEmpty())
			bpc.addSeries(list2, legend2, legend2);
		if(!list1.isEmpty())
			bpc.addSeries(list1, legend1, legend1);
		bpc.save(f.getAbsolutePath()+"_bp.png", 300, 250, "png");
		result.addOutputItem(new Item(itemName+"bp", f.getName()+"_bp.png", itemLabel, TYPE.IMAGE, new ArrayList<String>(), new HashMap<String, String>(2), itemGroup));
		
		

	}


	private String getSymbol(double n1, double n2){
		String symbol = "<";
		if(n1 > n2)
			symbol = ">";
		return symbol;
	}
	private void statsOneListAnalysis(String side) throws IOException{
		// 1st Analysis
		File f = null;
		logger.debug("Starting 1st Analysis..................");
		List<ProteinVertex> list = this.subProteinNetwork1.getInteractomeGraph().getVertices();


		List<Double> relBetInter = proteinNetwork.getBetRelList();
		List<Double> relBetList1 = new ArrayList<Double>();

		List<Double> connInter = proteinNetwork.getConnList();
		List<Double> connList1 = new ArrayList<Double>();

		List<Double> clustInter = proteinNetwork.getClustList();
		List<Double> clustList1 = new ArrayList<Double>();

		createTopoFilterListNoIntermediates(list, relBetList1, connList1, clustList1);
		//createOutputTopoList(list,this.mapList1,"1");
		String toWrite = ksTest(relBetList1, relBetInter, connList1, connInter, clustList1, clustInter, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_list_inter_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			String symbol = "";
			double mean = 0.0;
			
			//result.addOutputItem(new Item("list_inter_kol_param", f.getName(), "File", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List - Interactome"));

			mean = MathUtils.mean(ListUtils.toDoubleArray(relBetList1));
			symbol = getSymbol(mean, proteinNetwork.getMeanRelBet());
			//result.addOutputItem(new Item("list_inter_kol_param_bet",  Double.toString(getPValue(relBetList1, relBetInter, side)), "P-value betweenness: List "+symbol+" Interactome", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List - Interactome"));
			result.addOutputItem(new Item("list_inter_kol_param_bet",  StringUtils.decimalFormat(getPValue(relBetList1, relBetInter, side), decimalFormat), "Relative betweenness: List "+symbol+" Interactome pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Role of list within interactome of reference.Betweenness"));
			
			mean = MathUtils.mean(ListUtils.toDoubleArray(connList1));
			symbol = getSymbol(mean, proteinNetwork.getMeanConnections());
			//result.addOutputItem(new Item("list_inter_kol_param_conn",  Double.toString(getPValue(connList1, connInter, side)), "P-value connections: List "+symbol+" Interactome", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List - Interactome"));
			result.addOutputItem(new Item("list_inter_kol_param_conn",  StringUtils.decimalFormat(getPValue(connList1, connInter, side), decimalFormat), "Connections: List "+symbol+" Interactome pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Role of list within interactome of reference.Connections"));

			mean = MathUtils.mean(ListUtils.toDoubleArray(clustList1));
			symbol = getSymbol(mean, proteinNetwork.getMeanClust());
			//result.addOutputItem(new Item("list_inter_kol_param_clu",  Double.toString(getPValue(clustList1, clustInter, side)), "P-value clustering: List "+symbol+" Interactome", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List - Interactome"));
			result.addOutputItem(new Item("list_inter_kol_param_clu",  StringUtils.decimalFormat(getPValue(clustList1, clustInter, side), decimalFormat), "Clustering: List "+symbol+" Interactome pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Role of list within interactome of reference.Clustering coeff"));

		}
//		else
//			result.addOutputItem(new Item("list_inter_kol_param", "Empty results", "List - Interactome", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Statistic results.Kolmogorov-Smirnov test"));

		if(images){
			createImages(outputFileName+"_list_inter_relBet", relBetList1, "list1", relBetInter, "inter", "list_inter_relBet", "Plot", "Network parameters evaluation.Role of list within interactome of reference.Betweenness", "Relative Betweenness");
			createImages(outputFileName+"_list_inter_conn", connList1, "list1", connInter, "inter", "list_inter_conn", "Plot", "Network parameters evaluation.Role of list within interactome of reference.Connections", "Connections");
			createImages(outputFileName+"_list_inter_clust", clustList1, "list1", clustInter, "inter", "list_inter_clust", "Plot", "Network parameters evaluation.Role of list within interactome of reference.Clustering coeff","Clustering coeff");
		}
		createOutputTopoList(list,this.mapList1,"1");	
		logger.debug("Finished 1st Analysis..................");
		// Starting 2nd analysis
		logger.debug("Starting 2nd Analysis..................");
		List<Double> relBetSubnet1 = subProteinNetwork1.getBetRelList();
		List<Double> relBetRandoms = new ArrayList<Double>();

		List<Double> connSubnet1 = subProteinNetwork1.getConnList();
		List<Double> connRandoms = new ArrayList<Double>();

		List<Double> clustSubnet1 = subProteinNetwork1.getClustList();
		List<Double> clustRandoms = new ArrayList<Double>();
		for(ProteinNetwork proteinNetwork : subProteinNetworkRandoms){
			int  randomVertex = (int) (Math.random() * proteinNetwork.getInteractomeGraph().getVertices().size());
			ProteinVertex v = proteinNetwork.getInteractomeGraph().getVertices().get(randomVertex);
			relBetRandoms.add(v.getRelativeBetweenness());
			connRandoms.add(Double.parseDouble(Integer.toString(proteinNetwork.getInteractomeGraph().getDegreeOf(v))));
			clustRandoms.add(v.getClusteringCoefficient());
		}

		toWrite = ksTest(relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms, side);
		if(!toWrite.equals("")){
			String symbol;
			double mean1, mean2;
			f = new File(outputFileName+"_sn_random_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			//result.addOutputItem(new Item("sn_random_kol_param", f.getName(), "File", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.Subnet - Random"));
			
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(relBetSubnet1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(relBetRandoms));
			symbol = getSymbol(mean1, mean2);
			result.addOutputItem(new Item("sn_random_kol_param_bet", StringUtils.decimalFormat(getPValue(relBetSubnet1, relBetRandoms, side), decimalFormat), "Relative betweenness: List1 "+symbol+" Random pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Minimal Connected Network topological evaluation.Betweenness"));
			
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(connSubnet1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(connRandoms));
			symbol = getSymbol(mean1, mean2);
			result.addOutputItem(new Item("sn_random_kol_param_conn",  StringUtils.decimalFormat(getPValue(connSubnet1, connRandoms, side), decimalFormat), "Connections: List1 "+symbol+" Random pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Minimal Connected Network topological evaluation.Connections"));
			
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(clustSubnet1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(clustRandoms));
			symbol = getSymbol(mean1, mean2);
			result.addOutputItem(new Item("sn_random_kol_param_clu",  StringUtils.decimalFormat(getPValue(clustSubnet1, clustRandoms, side), decimalFormat), "Clustering: List1 "+symbol+" Random pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Minimal Connected Network topological evaluation.Clustering Coeff"));

		}
//		else
//			result.addOutputItem(new Item("sn_random_kol_param", "Empty results", "Subnet - Random", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Statistic results.Kolmogorov-Smirnov test"));

		if(images){
			createImages(outputFileName+"_sn_random_relBet", relBetSubnet1, "subnet1", relBetRandoms, "randoms", "sn_random_relBet", "Plot", "Network parameters evaluation.Minimal Connected Network topological evaluation.Betweenness","Relative Betweenness");
			createImages(outputFileName+"_sn_random_conn", connSubnet1, "subnet1", connRandoms, "randoms", "sn_random_conn", "Plot", "Network parameters evaluation.Minimal Connected Network topological evaluation.Connections","Connections");
			createImages(outputFileName+"_sn_random_clust", clustSubnet1, "subnet1", clustRandoms, "randoms", "sn_random_clust", "Plot", "Network parameters evaluation.Minimal Connected Network topological evaluation.Clustering Coeff","Clustering Coeff");
		}
		logger.debug("Finished 2nd Analysis..................");
	}
	private void statsTwoListsAnalisys(String side) throws IOException{
		// 1st Analysis
		logger.debug("Starting 1st Analysis..................");
		File f = null;

		List<Double> relBetList1 = new ArrayList<Double>();
		List<Double> connList1 = new ArrayList<Double>();
		List<Double> clustList1 = new ArrayList<Double>();
		List<Double> relBetList2 = new ArrayList<Double>();
		List<Double> connList2 = new ArrayList<Double>();
		List<Double> clustList2 = new ArrayList<Double>();

		List<ProteinVertex> list1 = this.subProteinNetwork1.getInteractomeGraph().getVertices();
		List<ProteinVertex> list2 = this.subProteinNetwork2.getInteractomeGraph().getVertices();


		createTopoFilterList(list1, relBetList1, connList1, clustList1, 1);
		createTopoFilterList(list2, relBetList2, connList2, clustList2, 2);
		createOutputTopoList(list1,this.mapList1,"1");
		createOutputTopoList(list2,this.mapList2,"2");
		
		String toWrite = ksTest(relBetList1, relBetList2, connList1, connList2, clustList1, clustList2, side);
		if(!toWrite.equals("")){
			double mean1, mean2;
			String symbol;
			f = new File(outputFileName+"_list1_list2_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			//result.addOutputItem(new Item("list1_list2_kol", f.getName(), "File", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List1 - List2"));
			
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(relBetList1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(relBetList2));
			symbol = getSymbol(mean1, mean2);
//			result.addOutputItem(new Item("list1_list2_kol_bet",  Double.toString(getPValue(relBetList1, relBetList2, side)), "P-value betweenness: List1 "+symbol+" List2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List1 - List2"));
			result.addOutputItem(new Item("list1_list2_kol_bet",  StringUtils.decimalFormat(getPValue(relBetList1, relBetList2, side),decimalFormat), "Relative betweenness: List1 "+symbol+" List2 pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.List1 against List2.Betweenness"));

			mean1 = MathUtils.mean(ListUtils.toDoubleArray(connList1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(connList2));
			symbol = getSymbol(mean1, mean2);
//			result.addOutputItem(new Item("list1_list2_kol_con",  Double.toString(getPValue(connList1, connList2, side)), "P-value connections: List1 "+symbol+" List2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List1 - List2"));
			result.addOutputItem(new Item("list1_list2_kol_con",  StringUtils.decimalFormat(getPValue(connList1, connList2, side),decimalFormat), "Connections: List1 "+symbol+" List2 pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.List1 against List2.Connections"));
	
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(clustList1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(clustList2));
			symbol = getSymbol(mean1, mean2);
//			result.addOutputItem(new Item("list1_list2_kol_clu",  Double.toString(getPValue(clustList1, clustList2, side)), "P-value clustering: List1 "+symbol+" List2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.List1 - List2"));
			result.addOutputItem(new Item("list1_list2_kol_clu",  StringUtils.decimalFormat(getPValue(clustList1, clustList2, side),decimalFormat), "Clustering: List1 "+symbol+" List2 pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.List1 against List2.Clustering Coeff"));

		}
//		else
//			result.addOutputItem(new Item("list1_list2_kol", "Empty results", "List1 - List2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));

		if(images){
			createImages(outputFileName+"_list1_list2_relBet", relBetList1, "list1", relBetList2, "list2", "list1_list2_relBet", "Plot", "Network parameters evaluation.List1 against List2.Betweenness","Relative Betweenness");
			createImages(outputFileName+"_list1_list2_conn", connList1, "list1", connList2, "list2", "list1_list2_conn", "Plot", "Network parameters evaluation.List1 against List2.Connections","Connections");
			createImages(outputFileName+"_list1_list2_clust", clustList1, "list1", clustList2, "list2", "list1_list2_clust", "Plot", "Network parameters evaluation.List1 against List2.Clustering Coeff","Clustering Coeff");
		}

		//2nd Analysis
		logger.debug("Starting 2nd Analysis..................");
		List<Double> relBetSubnet1 = subProteinNetwork1.getBetRelList();
		List<Double> relBetSubnet2 = subProteinNetwork2.getBetRelList();

		List<Double> connSubnet1 = subProteinNetwork1.getConnList();
		List<Double> connSubnet2 = subProteinNetwork2.getConnList();

		List<Double> clustSubnet1 = subProteinNetwork1.getClustList();
		List<Double> clustSubnet2 = subProteinNetwork2.getClustList();

		toWrite = ksTest(relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_sn1_sn2_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);		
			//result.addOutputItem(new Item("sn1_sn2_kol", f.getName(), "File", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test.Subnet1 - Subnet2"));
			double mean1, mean2;
			String symbol;
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(relBetSubnet1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(relBetSubnet2));
			symbol = getSymbol(mean1, mean2);
			result.addOutputItem(new Item("sn1_sn2_kol_bet",  StringUtils.decimalFormat(getPValue(relBetSubnet1, relBetSubnet2, side),decimalFormat), "Relative betweenness: Subnet1 "+symbol+" Subnet2 pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.List 1 Minimal Connected Network against List 2 Minimal Connected Network.Betweenness"));
			
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(connSubnet1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(connSubnet2));
			symbol = getSymbol(mean1, mean2);
			result.addOutputItem(new Item("sn1_sn2_kol_con",  StringUtils.decimalFormat(getPValue(connSubnet1, connSubnet2, side),decimalFormat), "Connections: Subnet1 "+symbol+" Subnet2 pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.List 1 Minimal Connected Network against List 2 Minimal Connected Network.Connections"));
			
			mean1 = MathUtils.mean(ListUtils.toDoubleArray(clustSubnet1));
			mean2 = MathUtils.mean(ListUtils.toDoubleArray(clustSubnet2));
			symbol = getSymbol(mean1, mean2);
			result.addOutputItem(new Item("sn1_sn2_kol_clu",  StringUtils.decimalFormat(getPValue(clustSubnet1, clustSubnet2, side),decimalFormat), "Clustering: Subnet1 "+symbol+" Subnet2 pval", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.List 1 Minimal Connected Network against List 2 Minimal Connected Network.Clustering Coeff"));

		}
//		else
//			result.addOutputItem(new Item("sn1_sn2_kol", "Empty results", "Subnet1 - Subnet2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));

		if(images){
			createImages(outputFileName+"_sn1_sn2_relBet", relBetSubnet1, "subnet1", relBetSubnet2, "subnet2", "sn1_sn2_relBet", "Plot", "Network parameters evaluation.List 1 Minimal Connected Network against List 2 Minimal Connected Network.Betweenness","Relative Betweenness");
			createImages(outputFileName+"_sn1_sn2_conn", connSubnet1, "subnet1", connSubnet2, "subnet2", "sn1_sn2_conn", "Plot", "Network parameters evaluation.List 1 Minimal Connected Network against List 2 Minimal Connected Network.Connections", "Connections");
			createImages(outputFileName+"_sn1_sn2_clust", clustSubnet1, "subnet1", clustSubnet2, "subnet2", "sn1_sn2_clust", "Plot", "Network parameters evaluation.List 1 Minimal Connected Network against List 2 Minimal Connected Network.Clustering Coeff","Clustering Coeff");
		}

	}
	private String ksTest(List<Double> relBetList1, List<Double> relBetList2, List<Double> connList1, List<Double> connList2, List<Double> clustList1, List<Double> clustList2, String side) throws IOException{
		if( relBetList1.size() == 0 || relBetList2.size() == 0 || connList1.size() == 0 || connList2.size() == 0 || clustList1.size() == 0 || clustList2.size() == 0)
			return "";
		StringBuilder sb = new StringBuilder();
		sb.append("#parameter\tpval\tside\n");
		sb.append("betweenness\t"+getPValue(relBetList1, relBetList2, side)+"\t"+side+"\n");
		sb.append("connections\t"+getPValue(connList1, connList2, side)+"\t"+side+"\n");
		sb.append("coefficient\t"+getPValue(clustList1,clustList2, side)+"\t"+side);
		return sb.toString();
	}
	
	private double getPValue(List<Double> list1, List<Double> list2, String side) throws IOException{
		KSTest kstest = new KSTest();
		return kstest.resultKolmogorovSmirnov(ListUtils.toDoubleArray(list1), ListUtils.toDoubleArray(list2), side).getPValue();
	}
	private void createOutputTopoList(List<ProteinVertex> list, Map<String, String> mapList, String index) throws IOException{
		
		StringBuilder sb = new StringBuilder();
		double relBete, conne, cluste;
		sb.append("#inputId\tid\tBetweenness\tConnections\tClustering").append(System.getProperty("line.separator"));
		for(ProteinVertex protein : list){
			if(mapList.get(protein.getId()) == null)
				continue;
			sb.append(mapList.get(protein.getId())).append("\t");
			sb.append(protein.getId()).append("\t");
			
			relBete = proteinNetwork.getRelBetweennessVertex(protein.getId());
			sb.append(relBete).append("\t");
			
			conne = proteinNetwork.getDegreeVertex(protein.getId());
			sb.append(conne).append("\t");

			cluste = proteinNetwork.getClusterinVertex(protein.getId());
			sb.append(cluste);
			
			sb.append(System.getProperty("line.separator"));
			
		}
		
		File f = new File(outputFileName+"_list_"+index+"_topo.txt");
		IOUtils.write(f, sb.toString());
		result.addOutputItem(new Item("list_"+index+"_topo", f.getName(), "Topological values of list "+index+" within interactome", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Role of list within interactome of reference.Files"));
	
	}
	private void createTopoFilterListNoIntermediates(List<ProteinVertex> list, List<Double> relBetList, List<Double> connList, List<Double> clustList) throws IOException{
		double relBete, conne, cluste;
		for(ProteinVertex protein : list){
			if( (this.intermediatesSub1 != null && this.intermediatesSub1.contains(protein.getId())) || proteinNetwork.getInteractomeGraph().getVertex(protein.getId()) == null)
				continue;
			relBete = proteinNetwork.getRelBetweennessVertex(protein.getId());
			relBetList.add(relBete);
			conne = proteinNetwork.getDegreeVertex(protein.getId());
			connList.add(conne);
			cluste = proteinNetwork.getClusterinVertex(protein.getId());
			clustList.add(cluste);
		}

	}
	private void createTopoFilterList(List<ProteinVertex> list, List<Double> relBetList, List<Double> connList, List<Double> clustList, int index) throws IOException{
		double relBete, conne, cluste;
		for(ProteinVertex protein : list){
			if( (this.intermediatesSub1 != null && this.intermediatesSub1.contains(protein.getId())) || proteinNetwork.getInteractomeGraph().getVertex(protein.getId()) == null)
				continue;
			relBete = proteinNetwork.getRelBetweennessVertex(protein.getId());
			relBetList.add(relBete);
			conne = proteinNetwork.getDegreeVertex(protein.getId());
			connList.add(conne);
			cluste = proteinNetwork.getClusterinVertex(protein.getId());
			clustList.add(cluste);
			
		}
	}
	
	
	
	private void createRandoms(int randoms, int randomSize) throws IOException{

//		StringBuilder sbMeans = createMeansHeader();
//		StringBuilder sbTopo = createTopoHeader();
//		StringBuilder sbComponents =createComponentsHeader();
		double []componentsRandoms = new double[randoms];
		for(int i=1; i<=randoms; i++){

			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), randomSize);
//			logger.debug("Randoms["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
			if(intermediate) {
				Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
				logger.debug("Randoms intermediate["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getAllEdges().size());
			}
			ProteinNetwork subProteinNetwork = createSubnet(subgraph);
			logger.debug("Subnet created");
//			sbTopo.append(getTopologicalValues(subProteinNetwork, i, true)).append(System.getProperty("line.separator"));
//			sbMeans.append(getTopologicalMeanValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));

			if(components) {
				List<List<ProteinVertex>> componentsList =  subProteinNetwork.getInteractomeGraph().getAllInformationComponents(true);
				componentsRandoms[i-1] = (double)componentsList.size();
			//	sbComponents.append(getComponentsValues(subProteinNetwork, i, componentsList));
			}
			logger.debug("Components created");
			subProteinNetworkRandoms.add(subProteinNetwork);
			
			
			
		}
		//Here we add the info for components in result.xml
		addOutputComponentConfidenceInterval(componentsRandoms);

		//File f = new File(outputFileName+"_sn_1-"+(randoms)+"_topo.txt");
		//randomFilesToString(sbTopo, f.getAbsolutePath());
		//result.addOutputItem(new Item("randoms_topo_param", f.getName(), "Topografical values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Randoms results"));

		//f = new File(outputFileName+"_sn_1-"+(randoms)+"_means.txt");
		//randomFilesToString(sbMeans, f.getAbsolutePath());
		//result.addOutputItem(new Item("randoms_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Randoms results"));

		//f = new File(outputFileName+"_sn_1-"+(randoms)+"_comp.txt");
		//randomFilesToString(sbComponents, f.getAbsolutePath());
		//result.addOutputItem(new Item("randoms_means_param", f.getName(), "Components values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Randoms results"));
	}
	private void addOutputComponentConfidenceInterval(double []componentsRandoms){
		//Here we calculate 95% confidence interval
		int[] range = getRange(componentsRandoms);
		String rangeString = "["+Integer.toString(range[0])+", "+Integer.toString(range[1])+"]";
		result.addOutputItem(new Item("components_number", componentsListSub1.size()+" "+rangeString, "Number of components [95% confidence interval]", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Network parameters evaluation.Minimal Connected Network topological evaluation"));
	}
//	private void addOutputComponentMore1Node(List<List<ProteinVertex>> componentsListSub){
//		int compMoreThan1Node = componentsMoreThanOneNode(componentsListSub);
//		result.addOutputItem(new Item("comp_more_than_1_node", Integer.toString(compMoreThan1Node), "Number of components with more than 1 node", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimal Connected Network topological evaluation"));
//
//	}
	private int componentsMoreThanOneNode(List<List<ProteinVertex>> components){
		int compMoreThan1Node = 0;
		for(List<ProteinVertex> comp : components){
			if(comp.size() > 1)
				compMoreThan1Node++;
		}
		return compMoreThan1Node;
	}
	private int[] getRange(double []componentsRandoms){
		int[] resultRange = new int[2]; 
		resultRange[0] = (int)MathUtils.percentile(componentsRandoms, 2.5);
		resultRange[1] = (int)MathUtils.percentile(componentsRandoms, 97.5);
		return resultRange;
		
	}
	
	
	
	private ProteinNetwork createSubnet(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph) {
		ProteinNetwork subProteinNetwork = new ProteinNetwork(subgraph);
		subProteinNetwork.calcTopologicalValues();
		subProteinNetwork.calcTopologicalMeanValues();
		return subProteinNetwork;
	}

	private StringBuilder createMeansHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Std: Standard deviation").append(System.getProperty("line.separator"));
		sb.append("#Subnet\tMeanBetweenness\tStdBetweenness\t").append("MeanConnections\tStdConnections\t").append("MeanClustering\tStdClustering");
		if(components){
			sb.append("\tComp\t1Comp");
		}
		sb.append("\tNodes");
		if(bicomponents){
			sb.append("\tBiComp\t");
		}
		sb.append(System.getProperty("line.separator"));
		return sb;
	}

	private StringBuilder createTopoHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tId\tBetweenness\tConnections\tClustering").append(System.getProperty("line.separator"));
		return sb;
	}
	
	private StringBuilder createExternalTopoHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tId\tBetweenness\tConnections\tClustering").append(System.getProperty("line.separator"));
		return sb;
	}

	private StringBuilder createComponentsHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#CN: Component number").append(System.getProperty("line.separator"));
		sb.append("#Subnet\tCN\tDiameter\tSize\tNodes").append(System.getProperty("line.separator"));
		return sb;
	}


	private String getTopologicalValues(ProteinNetwork subProteinNetwork, int subnet, boolean random){
		StringBuilder sb = new StringBuilder();
		
		for (ProteinVertex proteinVertex : subProteinNetwork.getInteractomeGraph().getVertices()) {
			if(proteinVertex == null)
				continue;
			String inputId= "";
			sb.append("sn"+(subnet)).append("\t");
			if( (type.equals("transcripts") || type.equals("genes"))  && random == false){
				inputId = this.getValue(proteinVertex.getId());
			}
			
			sb.append(inputId+proteinVertex.getId()).append("\t");
			sb.append(proteinVertex.getRelativeBetweenness()).append("\t").append(subProteinNetwork.getInteractomeGraph().getDegreeOf(proteinVertex)).append("\t").append(proteinVertex.getClusteringCoefficient());
			sb.append(System.getProperty("line.separator"));
		}
		if(!sb.toString().equals(""))
			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
		return sb.toString();
	}
	private String getExternalTopologicalValues(ProteinNetwork subProteinNetwork, int subnet, Set<String> intermediates){
		StringBuilder sb = new StringBuilder();
		for (ProteinVertex proteinVertex : subProteinNetwork.getInteractomeGraph().getVertices()) {
			if(proteinVertex == null || !intermediates.contains(proteinVertex.getId()))
				continue;
			sb.append("sn"+(subnet)).append("\t");
			sb.append(proteinVertex.getId()).append("\t");
			sb.append(proteinVertex.getRelativeBetweenness()).append("\t").append(subProteinNetwork.getInteractomeGraph().getDegreeOf(proteinVertex)).append("\t").append(proteinVertex.getClusteringCoefficient());
			sb.append(System.getProperty("line.separator"));
		}
		if(!sb.toString().equals(""))
			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
		return sb.toString();
	}

//	private String getTopologicalMeanValues(ProteinNetwork subProteinNetwork, int subnet, List<List<ProteinVertex>> listComponents, int bicomponentsNumber){
//		StringBuilder sb = new StringBuilder();
//		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = subProteinNetwork.getInteractomeGraph();
//		sb.append("sn"+(subnet)).append("\t");
//		sb.append(subProteinNetwork.getTopologicalMeanValuesToString());
//
//		if(components){
//			int oneComponent = 0;
//			int moreOneComponent = 0;
//			//List<List<ProteinVertex>> listComponents = interactomeGraph.getAllInformationComponents(true);
//			for(List<ProteinVertex> listVertex : listComponents){
//				if(listVertex.size() == 1)
//					oneComponent++;
//				else if(listVertex.size() > 1)
//					moreOneComponent++;
//			}
//			sb.append("\t"+moreOneComponent+"\t"+oneComponent);
//		}
//		sb.append("\t").append(interactomeGraph.getVertices().size());
//		if(bicomponents){
//			
//			sb.append("\t").append(bicomponentsNumber);
//			//sb.append("\t").append(interactomeGraph.getNumberOfBicomponents());
//			
//		}
//
//		return sb.toString();
//
//	}
	private String getComponentsValues(ProteinNetwork subProteinNetwork, int subnet, List<List<ProteinVertex>> componentsList){
		StringBuilder sb = new StringBuilder();
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = subProteinNetwork.getInteractomeGraph();

		if(components){
			List<Double> componentsDiameter = Calc.calcDiameter(interactomeGraph, componentsList);
			for(int i=0; i < componentsList.size(); i++){
				sb.append("sn"+(subnet)).append("\t");
				StringBuilder sbNodes = new StringBuilder();
				for (int k = 0; k < componentsList.get(i).size(); k++) {
					String nodeId = "";
					if(type.equals("transcripts") || type.equals("genes")){
						nodeId = getValue(componentsList.get(i).get(k).getId());
					}
					sbNodes.append(nodeId+componentsList.get(i).get(k).getId());
					if(k!=componentsList.get(i).size()-1)
						sbNodes.append(",");
				}
				sb.append((i)+"\t"+componentsDiameter.get(i)+"\t"+componentsList.get(i).size()+"\t"+sbNodes.toString()+ System.getProperty("line.separator"));
			}
		}
		return sb.toString();
	}
	private String getValue(String idNode){
		String transcriptId="";
		if(this.mapList1.get(idNode) != null)
			transcriptId = this.mapList1.get(idNode)+":";
		else if(this.mapList2.get(idNode) != null)
			transcriptId = this.mapList2.get(idNode)+":";
		return transcriptId;
	}
//	private void randomFilesToString(StringBuilder sb, String file) throws IOException{
//		if(!sb.toString().equals("")){
//			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
//			IOUtils.write(file, sb.toString());
//		}
//	}
	private List<ProteinVertex> toVertex(List<String> vertices){
		List<ProteinVertex> verticesList = new ArrayList<ProteinVertex>(vertices.size());
		for(String proteinName : vertices){
			if(!proteinName.equals(""))
				verticesList.add(new ProteinVertex(proteinName));
		}
		return verticesList;
	}

	private Map<String,String> getGenEnsemble(List<String> vertices){
		Map<String, String> listGenEnsembl = new HashMap<String, String>();
//		FeatureList<XRef> xrefsEns = new FeatureList<XRef>();
		for(String proteinName : vertices){
			try {
//				xrefsEns  = xrefDBMan.getByDBName(proteinName, "ensembl_gene");
				XRef xrefsEns  = xrefDBMan.getByDBName(proteinName, "ensembl_gene");
				if(xrefsEns != null && !xrefsEns.getXrefItems().get("ensembl_gene").isEmpty() /*&& !xrefsEns.getXrefItems().get("ensembl_gene").get(0).getDisplayName().equals(proteinName)*/)
					listGenEnsembl.put(xrefsEns.getXrefItems().get("ensembl_gene").get(0).getDisplayName(), proteinName);
			} catch (Exception e) {
				listGenEnsembl.put(proteinName,proteinName);
				e.printStackTrace();
			}
		}
		return listGenEnsembl;
	}
//	public void addOutputAppletItem(File xmlFile, int index) {
//		if (xmlFile.exists()) {
//			String url = "SnowViewer?filename=" + xmlFile.getName() + "&width=600&height=600";
//			result.addOutputItem(new Item("viewer" + index + "_param", url, "Viewer for network #" + index, TYPE.HTML, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Network viewer"));
//
//			url = "SnowViewer?filename=" + xmlFile.getName();
//			result.addOutputItem(new Item("viewer" + index + "_param_new_window", url, "Open applet for network #" + index + " in a new window", TYPE.LINK, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Network viewer"));
//		}
//	}
	private String graphToSif(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> graph, Map<String, String> mapList){
		StringBuilder sb = new StringBuilder();
		String source = "";
		String target = "";
		/**Non Isolated vertices **/
		for (DefaultEdge e : graph.getAllEdges()) {
			source = graph.getVertex(e.getSource()).getId();
			target = graph.getVertex(e.getTarget()).getId();
			if(mapList.containsKey(source))
				source = mapList.get(source);
			if(mapList.containsKey(target))
				target = mapList.get(target);
			sb.append(source).append("\t").append("pp").append("\t").append(target);
			sb.append(System.getProperty("line.separator"));
		}
		/** Isolated vertices **/
		for(ProteinVertex vertex : graph.getAllVertices()){
			if(graph.getDegreeOf(vertex) == 0){
				if(mapList.containsKey(vertex.getId()))
					sb.append(mapList.get(vertex.getId())).append(System.getProperty("line.separator"));
			}
		}
		if(!sb.toString().equals(""))
			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
		return sb.toString();
	}
	private void addOutputSvgViewer(File file, int index, Set<String> intermediates) {
		List<String> tags = new ArrayList<String>();
		String intermediateTag = "";
		for(String intermediate : intermediates)
			intermediateTag+=intermediate+"|";
		
		if (file.exists()) {
			String list = "list";
			if(index==1)
				list+="1";
			if(index==2)
				list+="2";
			tags.add("INTERACTOME_VIEWER");
			tags.add(list);
			tags.add(this.interactome);
			tags.add(intermediateTag);
			result.addOutputItem(new Item("svg_viewer" + index + "_param", file.getName(), "Minimun Connected Network interactions", TYPE.FILE, tags, new HashMap<String, String>(2), "Network viewer "+list));

//			url = "SnowViewer2?filename=" + jsonFile.getName();
//			result.addOutputItem(new Item("svg_viewer" + index + "_param_new_window", url, "Open svg for network #" + index + " in a new window", TYPE.LINK, StringUtils.toList("SERVER,INCLUDE_REFS,SNOW", ","), new HashMap<String, String>(2), "Network viewer 2"));
		}
	}
	
	private String getMcnInteractors(ProteinNetwork subProteinNetwork, Set<String> intermediates, Map<String, String> mapList) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException{
		
		StringBuilder sb = new StringBuilder();
		String tab = "\t";
		sb.append("#input_id").append(tab).append("id").append(tab).append("type").append(tab).append("rank").append(tab).append("bet").append(tab).append("clust").append(tab);
		sb.append("conn").append(tab).append("go").append(tab).append("kegg").append(System.getProperty("line.separator"));
		
		/** First significant no external nodes **/
		for(ProteinVertex pr : subProteinNetwork.getInteractomeGraph().getAllVertices()){
			
			String id = "-";
			if(mapList.get(pr.getId()) != null)
				id = mapList.get(pr.getId());
			sb.append(id).append(tab);
			sb.append(pr.getId()).append(tab);
			String listType = "list";
			if(intermediates.contains(pr.getId()))
				listType = "external";
			sb.append(listType).append(tab);
			sb.append(pr.getRelativeBetweenness()).append(tab);
			sb.append(getSigleMcnInteractors(subProteinNetwork,pr.getId()));
			sb.append(System.getProperty("line.separator"));
		}
//		/** Second significant external nodes **/
//		for(String id : intermediates){
//			String convertedId = "-";
//			if(mapNames.get(id) != null)
//				convertedId = mapNames.get(id);
//			sb.append(convertedId).append(tab);
//			sb.append(id).append(tab);
//			sb.append("external").append(tab);
//			/** Input value **/
//			sb.append("-").append(tab);
//			sb.append(getSigleMcnInteractors(subProteinNetwork, id));
//			sb.append(System.getProperty("line.separator"));
//		}
		return sb.toString();
	}
	private String getSigleMcnInteractors(ProteinNetwork subProteinNetwork, String id) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException{
		GODBManager gof = new GODBManager(dbConnector);
		StringBuilder sb = new StringBuilder();
		String tab = "\t";
		ProteinVertex prVertex = subProteinNetwork.getInteractomeGraph().getVertex(id);
		sb.append(prVertex.getRelativeBetweenness()).append(tab);
		sb.append(prVertex.getClusteringCoefficient()).append(tab);
		sb.append(subProteinNetwork.getInteractomeGraph().getDegreeOf(prVertex)).append(tab);
		
		if(xrefDBMan.getDBConnector().getDbConnection().getDatabase() != null){
			List<String> dbNames = new ArrayList<String>();
			dbNames.add("go");
			dbNames.add("kegg");
			AnnotationDBManager annotationMng = new AnnotationDBManager(dbConnector);
			Map<String, String> keggItems  = annotationMng.getAnnotationTermNames("kegg");
			XRef xrefEns  = xrefDBMan.getByDBName(id, dbNames);
			List<XRefItem> xrefitemListGO = xrefEns.getXrefItems().get("go");
			if(xrefitemListGO != null && !xrefitemListGO.isEmpty()){
				for(XRefItem xrefitem : xrefitemListGO){
					GO go = gof.getByAccesion(xrefitem.getDisplayName());
					sb.append(go.getName()).append(",");
				}
				sb = this.deleteLastCh(sb, ",");
			}
			
			List<XRefItem> xrefitemListKegg = xrefEns.getXrefItems().get("kegg");
			if(xrefitemListKegg != null && !xrefitemListKegg.isEmpty()){
				sb.append(tab);
				for(XRefItem xrefitem : xrefitemListKegg){
					if(keggItems.containsKey(xrefitem.getDisplayName())){
						sb.append(keggItems.get(xrefitem.getDisplayName())).append(",");
					}
				}
				sb = this.deleteLastCh(sb, ",");
			}
		}
		return sb.toString();
	}
	/** Here we discover the type of input id, if it is uniprot, ensembl_gene, hgnc_symbol 
	 **/
	private DBName getOriginalDbName(Collection<String> list)  {
		Map<String, Integer> dbNameNumber = new HashMap<String, Integer>();
		int recognized = 0;
		if(xrefDBMan.getDBConnector().getDbConnection().getDatabase() == null){
			return null;
		}
		int max = Integer.MIN_VALUE;
		DBName dbNameReturn = null;
		for(String node : list){
			List<DBName> dbnames;
			try {
				dbnames = xrefDBMan.getAllDBNamesById(node);
				if(dbnames.size() > 0){
					for (DBName dbName : dbnames) {
						String name = dbName.getDbname();
						int number = 1;
						if( !dbNameNumber.containsKey(name) ){
								dbNameNumber.put(name, number);
						}
						else{
							number = dbNameNumber.get(name);
							number++;
							dbNameNumber.put(name, number);
						}
						if(number >= max){
							max = number;
							dbNameReturn = dbName;
						}
					}
				}
			} catch (Exception e) {
				recognized++;
				/** Esto ocurre cuando no tenemos la especie en infrared, por ejemplo ECO **/
				//System.err.println("Error GsnowPreprocessin getOriginalDbName: "+node.getOriginalId()+" error: "+e.getLocalizedMessage());
			}
		}
		System.out.println("no recognized: "+recognized);
		System.out.println("DB Matched: "+dbNameReturn);
		return dbNameReturn;
	}
	private String getIdToDbName(String id) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		String nodeId = "-";
		String dbName = this.dbName.getDbname();
		if(xrefDBMan.getDBConnector().getDbConnection().getDatabase() != null){
			List<String> dbNames = new ArrayList<String>();
			dbNames.add(dbName);
			XRef xrefEns  = xrefDBMan.getByDBName(id, dbNames);
			List<XRefItem> xrefitemList = xrefEns.getXrefItems().get(dbName);
			if(xrefitemList != null && !xrefitemList.isEmpty()){
				for(XRefItem xrefitem : xrefitemList){
					nodeId = xrefitem.getDisplayName();
				}
			}
		}
		return nodeId;
	}
	private StringBuilder deleteLastCh(StringBuilder sb, String ch){
		StringBuilder sbReturn = new StringBuilder();
		if(!sb.toString().equals("")){
			//System.out.println(sb);
			sb.deleteCharAt(sb.lastIndexOf(ch));
		}
		return sbReturn.append(sb.toString());
	}
}


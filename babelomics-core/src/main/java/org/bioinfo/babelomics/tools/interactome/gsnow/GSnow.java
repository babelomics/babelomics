package org.bioinfo.babelomics.tools.interactome.gsnow;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bioinfo.babelomics.tools.interactome.RandomsSnowManager;
import org.bioinfo.babelomics.tools.interactome.SnowPrinter;
import org.bioinfo.babelomics.tools.interactome.SnowTool;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.ListInfo;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.Node;
import org.bioinfo.babelomics.tools.interactome.gsnow.annovar.AnnovarInterface;
import org.bioinfo.babelomics.tools.interactome.gsnow.annovar.Vcf;
import org.bioinfo.commons.io.TextFileReader;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.Subgraph;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.funcannot.GODBManager;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.util.MathUtils;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.bioinfo.babelomics.utils.XrefManager;

public class GSnow extends SnowTool {

    private ListInfo listInfo;
    private GSnowItem significantItem;
    private Map<Integer, GSnowItem> gsnowItems;
//    private DBConnector dbConnector;
//    private XRefDBManager xrefDBMan;
//    private GODBManager gof;
    private int numberOfMinNodes = 0;
    /**
     * Here we indicate the number of nodes from where we will start the GSnow, in case of no seed list it will get calculate automatically *
     */
    private int numberOfMaxNodes = 200;
    /**
     * Here we indicate the max number of nodes that we will analyse,  in case of no seed list it will get calculate automatically *
     */
    private double significantValue;
    /**
     * Here we indicate which is the minimum significant value that we are interested in. *
     */
    private double defaultSignificantValue = 0.05;
    /**
     * Here we indicate the default significant value that we are interested in. *
     */
    //private int maxNumberOfRandomsComponents; // Here we indicate the maximum size of the randoms from the components
    private String order;
    private Integer numberItems;
    private double cutOff;
    //	private boolean snow;
    private Map<String, String> mapNames;
    /**
     * Important info: mapNames<String,String> = <Ensenbl_gen, normal_name>  *
     */
    private String decimalFormat;

    /** goterms **/
    Map<String, List<String>> goTerms;

    @Override
    public void initOptions() {

        options.addOption(OptionFactory.createOption("list", "An input file containing a node per line", false, true));
        options.addOption(OptionFactory.createOption("seedlist", "An input file containing a node per line. This is a seed list", false, true));
        options.addOption(OptionFactory.createOption("order", "Here we indicate wether we want to order the list ascending(default) or descending", false, true));
        options.addOption(OptionFactory.createOption("number-items", "Here we indicate how many nodes we want to process", false, true));
        options.addOption(OptionFactory.createOption("cut-off", "Here we indicate where we cut the list from (considering the order)", false, true));
        options.addOption(OptionFactory.createOption("significant-value", "Here we indicate the significant value", false, true));
//		options.addOption(OptionFactory.createOption("snow", "Here if there is this option we calculate a SNOW", false, true));

        /** This next two options are only used when we want to generate the randoms for the test, only from time to time(maybe twice a year). **/
        options.addOption(OptionFactory.createOption("size-min", "Minimum size for randoms", false, true));
        options.addOption(OptionFactory.createOption("size-max", "Maximum size for randoms", false, true));
    }

    @Override
    protected void execute() {
        this.decimalFormat = "#.####";
        initExecute();
//		if(commandLine.hasOption("list") && !(commandLine.hasOption("size-min") && commandLine.hasOption("size-max")))
//			executeGSnow();
        if (commandLine.hasOption("list"))
            executeGSnow();
        else if (options.hasOption("size-min") && options.hasOption("size-max"))
            executeGenerateRandoms();
        else
            System.err.println("[Error] Some parameters missing.");
    }


    private void executeGSnow() {
        try {

            /** Get goterms**/
            org.bioinfo.babelomics.utils.GOManager goManager = new org.bioinfo.babelomics.utils.GOManager();
            goTerms = goManager.getTerms();


            SnowPrinter snowPrinter = new SnowPrinter();
            String interactomeMsg = getInteractomeMsg();

            order = commandLine.hasOption("order") ? commandLine.getOptionValue("order") : "ascending";
            cutOff = commandLine.hasOption("cut-off") ? Double.parseDouble(commandLine.getOptionValue("cut-off")) : Double.NaN;
            significantValue = commandLine.hasOption("significant-value") ? Double.parseDouble(commandLine.getOptionValue("significant-value")) : defaultSignificantValue;
//			snow = commandLine.hasOption("snow") ? true : false;
//			snow = (commandLine.hasOption("snow") && !(commandLine.getOptionValue("snow").equalsIgnoreCase("0")));

//            dbConnector = new DBConnector(interactome, new File(babelomicsHomePath + "/conf/infrared.properties"));
//            xrefDBMan = new XRefDBManager(dbConnector);
//            gof = new GODBManager(dbConnector);
            logger.debug("Starting list.........");

            result.addOutputItem(new Item("interactome_param", interactomeMsg, "Interactome", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String, String>(), "Input parameters"));
            result.addOutputItem(new Item("group_param", group, "Interactome group (curated or all)", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String, String>(), "Input parameters"));
            result.addOutputItem(new Item("intermediate_param", (intermediate ? "1" : "0"), "Number of external proteins introduced", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String, String>(), "Input parameters"));

            String folder = loadFile();
            System.out.println(folder);

            String seedListFile = "";
            if (commandLine.hasOption("seedlist") && !commandLine.getOptionValue("seedlist").equalsIgnoreCase("none")) {
                seedListFile = commandLine.getOptionValue("seedlist");
                logger.info("Working with seedlist: " + seedListFile);
                FileUtils.checkFile(new File(seedListFile));

            } else {
                logger.info("Not working with seedlist" + commandLine.getOptionValue("seedlist"));
            }

            String nodeFile = commandLine.getOptionValue("list");
            /** Si no hay vcf seguimos el flujo normal del programa**/
            if (type.equalsIgnoreCase("vcf")) {
                Vcf vcf = new Vcf(nodeFile);
                AnnovarInterface annovar = new AnnovarInterface("/home/ralonso/annovar", "/home/ralonso/annovar/humandb");
                List<String> list = annovar.annotate(vcf, outdir);
                nodeFile = outputFileName + "_list.txt";
                File f = new File(nodeFile);
                IOUtils.write(f.getAbsoluteFile(), list);
                FileUtils.checkFile(f);
            }

            Map<Integer, List<Double>> dataMatrixRandoms = null;

            /** 1º Here the map is (size of list, value of this list of each random) **/
            if (seedListFile.equalsIgnoreCase("")) {
                dataMatrixRandoms = getDataMatrixRandoms(folder);
            }
            numberItems = commandLine.hasOption("number-items") ? Integer.parseInt(commandLine.getOptionValue("number-items")) : numberOfMaxNodes;

            GSnowPreprocessing preprocessing = new GSnowPreprocessing(proteinNetwork, type, interactome, numberOfMaxNodes, order, numberItems, cutOff);

//			String nodeFile = commandLine.getOptionValue("list");
            FileUtils.checkFile(new File(nodeFile));
            listInfo = preprocessing.preprocess(seedListFile, nodeFile);

            mapNames = new HashMap<String, String>();

            for (Node n : listInfo.getNodes()) {
//				if(n.getOriginalId().equals("NOD2"))
//					System.out.println("hola");
                mapNames.put(n.getId(), n.getOriginalId());
            }
            updateJobStatus("5", "List preprocessed");
            String inputDataTitle = "Input data";
            String seedInputDataTitle = inputDataTitle + ". Seed List";
            if (!seedListFile.equalsIgnoreCase("")) {
                inputDataTitle = inputDataTitle + ". List";
            }
            File f = new File(outputFileName + "_nodes.txt");
            IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodesList(listInfo.getNodes()));
            result.addOutputItem(new Item("nodes_file", f.getName(), "Final list file", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), inputDataTitle));
            result.addOutputItem(new Item("nodes_file_number", String.valueOf(listInfo.getNodes().size()), "Final list number", Item.TYPE.TEXT, new ArrayList<String>(), new HashMap<String, String>(), inputDataTitle));

            f = new File(outputFileName + "_not_matched_nodes.txt");
            IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodes(listInfo.getNotMatchNodes()));
            result.addOutputItem(new Item("list_info_not_matched_nodes", f.getName(), "Not matched nodes file", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), inputDataTitle));
            result.addOutputItem(new Item("list_info_not_matched_nodes_number", String.valueOf(listInfo.getNotMatchNodes().size()), "Not matched nodes number", Item.TYPE.TEXT, new ArrayList<String>(), new HashMap<String, String>(), inputDataTitle));

            f = new File(outputFileName + "_repeated_nodes.txt");
            IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodes(listInfo.getRepeatedNodes()));
            result.addOutputItem(new Item("list_info_repeated_nodes_file", f.getName(), "Duplicated nodes file", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), inputDataTitle));
            result.addOutputItem(new Item("list_info_repeated_nodes", (listInfo.getRepeatedNodes().size() == 0) ? "0" : String.valueOf(listInfo.getRepeatedNodes().size()), "Duplicated nodes number", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), inputDataTitle));

            if (!seedListFile.equalsIgnoreCase("")) {
                f = new File(outputFileName + "_seed_nodes.txt");
                List<Node> seedNodes = new ArrayList<Node>();
                for (Node node : listInfo.getNodes()) {
                    if (node.isSeed())
                        seedNodes.add(node);
                }
                IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodesList(seedNodes));
                result.addOutputItem(new Item("seed_nodes_file", f.getName(), "Final list file", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), seedInputDataTitle));
                result.addOutputItem(new Item("seed_nodes_file_number", String.valueOf(listInfo.getSeedNodes().size()), "Final list number", Item.TYPE.TEXT, new ArrayList<String>(), new HashMap<String, String>(), seedInputDataTitle));

                f = new File(outputFileName + "_not_matched_seed_nodes.txt");
                IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodes(listInfo.getSeedNodesNotMatched()));
                result.addOutputItem(new Item("list_info_not_matched_seed_nodes", f.getName(), "Not matched nodes file", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), seedInputDataTitle));
                result.addOutputItem(new Item("list_info_not_matched_seed_nodes_number", String.valueOf(listInfo.getSeedNodesNotMatched().size()), "Not matched nodes number", Item.TYPE.TEXT, new ArrayList<String>(), new HashMap<String, String>(), seedInputDataTitle));

                f = new File(outputFileName + "_repeated_seed_nodes.txt");
                IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodes(listInfo.getSeedNodesDuplicates()));
                result.addOutputItem(new Item("list_info_repeated_seed_nodes_file", f.getName(), "Duplicated nodes file", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), seedInputDataTitle));
                result.addOutputItem(new Item("list_info_repeated_seed_nodes", (listInfo.getSeedNodesDuplicates().size() == 0) ? "0" : String.valueOf(listInfo.getSeedNodesDuplicates().size()), "Duplicated nodes number", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), seedInputDataTitle));

            }

            if (listInfo.getNodes().size() == 0 || listInfo.getNodes().size() < this.numberOfMinNodes) {
                logger.info("List to small after preprocessing it. The size should have at least " + numberOfMinNodes + " after preprocessing. ");
                result.addOutputItem(new Item("small_size_list_error", "List to small after preprocessing it. The size should have at least " + numberOfMinNodes + " after preprocessing. ", "An error occurred", Item.TYPE.MESSAGE, Arrays.asList("ERROR"), new HashMap<String, String>(), "Error"));
                return;
            }
            /** 2º Here the map is (size of list, nodes of list) **/
            gsnowItems = getGSnowItems();
            updateJobStatus("10", "List preprocessed");

            /** 3º Here the map is (size of list, value of this list) **/
            getDataMatrixList();
            updateJobStatus("15", "List preprocessed");

            /** 3.1 In case of seed list we have to get another  dataMatrixRandoms**/
            if (!seedListFile.equalsIgnoreCase("")) {
                significantItem = getSeedPvalue();
            } else {
                /** 4º Here the map is (size of list, value of this list compared with randoms) **/
                compareDataMatrixListRandoms(dataMatrixRandoms);
                updateJobStatus("25", "List preprocessed");

                /** 5º Here we calculate the significant item from gsnowItems the **/
                significantItem = getSignificatValue();
            }


            if (significantItem != null)
            /** En el siguiente método se dibuja la tabla y el network visualizator, dentro se comprueba si se ha pedido un snow o no para dibujar las estadísticas**/
                this.drawResults(significantItem);

//			if(snow){
//				System.out.println("Executing SNOW");
//				this.executeSnow();
//			}
//			else
//				System.out.println("Not SNOW executed");

        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    /** 1º Here the map is (size of list, value of this list of each random) **/
    /**
     * Here we will load the data of the randoms in Map<Integer, List<Double>> *
     */
    public Map<Integer, List<Double>> getDataMatrixRandoms(String folder) throws IOException {

        Map<Integer, List<Double>> values1 = new HashMap<Integer, List<Double>>();
        List<Double> values2 = null;
        TextFileReader tfr = new TextFileReader(folder);
        String line = null;
        String[] fields;
        int size = 0;
        boolean minimunSize = false;
        while ((line = tfr.readLine()) != null) {
            if (!line.startsWith("#")) {

                values2 = new ArrayList<Double>();
                fields = line.split("\t");
                size = Integer.parseInt(fields[0]);
                for (int i = 1; i < fields.length; i++) {
                    values2.add(Double.parseDouble(fields[i]));
                }
                values1.put(size, values2);
                if (!minimunSize) {
                    //here we take the minimun size of the randoms
                    this.numberOfMinNodes = size;
                    minimunSize = true;
                }
            }
        }
        tfr.close();
        /** here we take the maximum0 size of the randoms **/
        this.numberOfMaxNodes = size;
        return values1;
    }

    /** 2º Here the map is (size of list, nodes of list) **/
    /**
     * i.e: 3={[ENSG00000198373	ENSG00000198373	-1.5, ENSG00000139946	ENSG00000139946	-1.5, ENSG00000132109	ENSG00000132109	-1.5] - rawValue = 0.0 - comparedValue = 0.0} *
     */
    private Map<Integer, GSnowItem> getGSnowItems() {

        Map<Integer, GSnowItem> gsnowItems = new HashMap<Integer, GSnowItem>();
        List<Node> nodes = new ArrayList<Node>();
        GSnowItem gsnowItem;
        List<Node> auxNodes;

        /** Here we indicate the numberOfMinNodes, that was really calculated before with the randoms **/
        if (this.listInfo.getSeedNodes().size() > this.numberOfMinNodes)
            this.numberOfMinNodes = this.listInfo.getSeedNodes().size();

        /** For the seed list, we start the gsnowItems from the size of the seed list **/
        if (!this.listInfo.getSeedNodes().isEmpty()) {
            for (int i = 0; i < this.listInfo.getSeedNodes().size(); i++)
                nodes.add(this.listInfo.getNodes().get(i));
            auxNodes = new ArrayList<Node>();
            auxNodes.addAll(nodes);
            gsnowItem = new GSnowItem();
            gsnowItem.setNodes(auxNodes);
            gsnowItems.put(this.listInfo.getSeedNodes().size(), gsnowItem);

        }
        /** This code works when there is a seed list and when no, because it has to start in 0 or in the seed list size, that can be 0 **/
        for (int i = this.listInfo.getSeedNodes().size(); i < this.listInfo.getNodes().size(); i++) {
            nodes.add(this.listInfo.getNodes().get(i));
            auxNodes = new ArrayList<Node>();
            auxNodes.addAll(nodes);
            gsnowItem = new GSnowItem();
            gsnowItem.setNodes(auxNodes);
            gsnowItems.put(i + 1, gsnowItem);
        }

        return gsnowItems;
    }

    /** 3º Here the map is (size of list, value of this list) **/
    /** This method will fill Map<Integer, GSnowItem> gsnowItems with its raw value, score (we don't fill the comparedValue) **/
    /**
     * Este es el que más tarda, sin duda *
     */
    private void getDataMatrixList() {
        List<ProteinVertex> nodes = new ArrayList<ProteinVertex>();
        List<List<ProteinVertex>> componentsList;
        List<Double> componentsSize;
        double rawValue;

        /** For the non-seed list **/
        /** for(int i=0; i < numberOfMinNodes; i++){
         nodes.add(new ProteinVertex(listInfo.getNodes().get(i).getId()));
         } **/

        for (int i = 0; i < this.numberOfMinNodes; i++) {
            nodes.add(new ProteinVertex(listInfo.getNodes().get(i).getId()));
        }
        SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph;

        /** this.listInfo.getNodes().size() is always <= numberOfMaxNodes, because I cut the list (if it is longer than numberOfMaxNodes) in the preprocessing **/
        /** For the non-seed list **/
        /**    for(int i = numberOfMinNodes; i <= this.listInfo.getNodes().size(); i++){ **/
        /** For the seed list **/
        for (int i = this.numberOfMinNodes; i <= this.listInfo.getNodes().size(); i++) {
            subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), nodes);
            if (intermediate) {
                Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
            }

            componentsSize = new ArrayList<Double>();
            componentsList = subgraph.getAllInformationComponents(true);
            for (List<ProteinVertex> component : componentsList)
                componentsSize.add((double) component.size());
            rawValue = getRawValue(subgraph, componentsList);

            gsnowItems.get(i).setRawValue(rawValue);
            gsnowItems.get(i).setComponentsSize(componentsSize);

            /** Calculating pendiente for the score **/
            double score = 0.0;
            score = getScore(rawValue, gsnowItems.get(i).getNodes().size());
//			score = (rawValue-1)/(gsnowItems.get(i).getNodes().size()-1);
            gsnowItems.get(i).setScore(score);

            if (i < listInfo.getNodes().size())
                nodes.add(new ProteinVertex(listInfo.getNodes().get(i).getId()));
        }
    }

    private double getRawValue(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph, List<List<ProteinVertex>> componentsList) {
        return (double) subgraph.getVertices().size() / (double) componentsList.size();
    }

    private double getScore(double rawValue, int size) {
        return (rawValue - 1) / (size - 1);
    }

    private GSnowItem getSeedPvalue() throws IOException {
        GSnowItem significantItemLocal = null;
//		Map<Integer, GSnowItem> gsnowItemsMap = new HashMap<Integer, GSnowItem>();
        List<ProteinVertex> seedList = new ArrayList<ProteinVertex>();
        for (String seed : listInfo.getSeedNodes()) {
            seedList.add(new ProteinVertex(seed));
        }
        /** This variable is used for knowing the size of the random**/
        int sizeRandomCutoffValue = 50;
        /** This variable is used for knowing the amount of executions that we are runnning **/
        int maxExecutions = 5;

        List<GSnowItem> gsnowItemsList = getMaxRelatives();
        Collections.sort(gsnowItemsList);

        //boolean found = false;
        for (int i = 0; i < maxExecutions; i++) {
            GSnowItem gSnowItem = gsnowItemsList.get(i);
//			gsnowItemsMap.put(gSnowItem.getNodes().size(), gSnowItem);
//			if(found)
//				continue;
            int randomsNumber = 500 + 500/** le sumo 500 luego lo bajaremos despues de los referees, hablar Luz **/;
            if (gSnowItem.getNodes().size() < sizeRandomCutoffValue || gSnowItem.getNodes().size() == sizeRandomCutoffValue) {
                randomsNumber = 1000;
            }

            int vertexToAdd = gSnowItem.getNodes().size() - seedList.size();
            List<Double> valueList = new ArrayList<Double>();
            for (int j = 0; j < randomsNumber; j++) {
                SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), seedList, vertexToAdd);
                if (intermediate) {
                    Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
                }
                List<Double> componentsSize = new ArrayList<Double>();
                List<List<ProteinVertex>> componentsList = subgraph.getAllInformationComponents(true);
                for (List<ProteinVertex> component : componentsList)
                    componentsSize.add((double) component.size());
                double rawValue = getRawValue(subgraph, componentsList);
                valueList.add(rawValue);
            }
            double comparedValue = getComparedValue(valueList, gSnowItem.getRawValue());
            if (comparedValue <= this.significantValue) {
                significantItemLocal = gSnowItem;
                significantItemLocal.setComparedValue(comparedValue);
//				found = true;
                break;
            }
        }
        significantItemLocal = siginificantItemOutputs(significantItemLocal);
        if (significantItemLocal != null)
            printGSnowItems(gsnowItemsList);
        return significantItemLocal;
    }
    /** 4º Here the map is (size of list, value of this list comparing with randoms) **/
    /**
     * Here we will get the compared value against the randoms *
     */
    private void compareDataMatrixListRandoms(Map<Integer, List<Double>> dataMatrixRandoms) {
        /** For the non-seed list **/
        /** for(int i = numberOfMinNodes; i <= this.listInfo.getNodes().size(); i++){ **/
        /** For the seed list **/
        for (int i = this.numberOfMinNodes; i <= this.listInfo.getNodes().size(); i++) {
            if (dataMatrixRandoms.get(i) == null) {
                logger.info("Warning: component with size " + i + " has not been found.");
                continue;
            }
            double rawValue = gsnowItems.get(i).getRawValue();
            List<Double> valueList = dataMatrixRandoms.get(i);
            double comparedValue = getComparedValue(valueList, rawValue);
            //System.out.println(comparedValue);
            gsnowItems.get(i).setComparedValue(comparedValue);
        }

    }

    private double getComparedValue(List<Double> valueList, double rawValue) {
        Double comparedValue = null;
        List<Double> intValue = new ArrayList<Double>();
        int rawBiggerValue = 0;
        int valueBiggerRaw = 0;
        for (double value : valueList) {
            if (rawValue > value) {
                intValue.add(1.0);
                rawBiggerValue++;
            } else {
                intValue.add(0.0);
                valueBiggerRaw++;
            }
        }
        comparedValue = (1 - MathUtils.mean(ListUtils.toDoubleArray(intValue)));
        return comparedValue;
    }

    /**
     * 5º Here we get the significant item in significantItem from gsnowItems *
     */
    private GSnowItem getSignificatValue() throws IOException {
        GSnowItem significantItemLocal = null;
        List<GSnowItem> gsnowItemsList = getMaxRelatives();
        Collections.sort(gsnowItemsList);

        if (gsnowItemsList.isEmpty())
            return null;

        for (GSnowItem gSnowItem : gsnowItemsList) {
            if (gSnowItem.getComparedValue() <= this.significantValue) {
                significantItemLocal = gSnowItem;
                break;
            }
        }
        significantItemLocal = siginificantItemOutputs(significantItemLocal);

        if (significantItemLocal != null)
            printGSnowItems(gsnowItemsList);

        return significantItemLocal;
    }

    private List<GSnowItem> getMaxRelatives() {
        double currentValue;
        double nextCurrentValue;
        double previousCurrentValue;
        List<GSnowItem> gsnowItemsList = new ArrayList<GSnowItem>();
        /** Aquí cogemos los picos, ahora los no significativos y los significativos, ¿sería necesario unicamente los significativos? **/
        for (int i = this.numberOfMinNodes + 1; i < this.listInfo.getNodes().size(); i++) {
            /** Cojemos el rawValue de la lista con size i-1, i y i+1 **/
            previousCurrentValue = gsnowItems.get(i - 1).getRawValue();
            currentValue = gsnowItems.get(i).getRawValue();
            nextCurrentValue = gsnowItems.get(i + 1).getRawValue();
            if (previousCurrentValue < currentValue && currentValue > nextCurrentValue) {
                gsnowItemsList.add(gsnowItems.get(i));
            }
        }
        return gsnowItemsList;
    }

    private GSnowItem siginificantItemOutputs(GSnowItem significantItemLocal) throws IOException {
        if (significantItemLocal == null) {
            logger.info("There is no significant minimal connected network");
            result.addOutputItem(new Item("significant_value", " There is no significant minimal connected network", "Significant value", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Results: Minimum Connected Network selected"));
            return null;
        }
        significantItemLocal.setSignificant(true);
        /** We write on web page the results **/
        logger.info("Significant value:" + significantItemLocal.getComparedValue() + ", Significant size: " + significantItemLocal.getNodes().size());
        String significantStringItem = StringUtils.decimalFormat(significantItemLocal.getComparedValue(), decimalFormat);

        if (significantItemLocal.getComparedValue() == 0)
            significantStringItem = "0.0001";
        result.addOutputItem(new Item("significant_value", "<" + significantStringItem, "MCN selected pval", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Results: Minimum Connected Network selected.All results"));
        result.addOutputItem(new Item("significant_size", Integer.toString(significantItemLocal.getNodes().size()), "MCN selected size", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Results: Minimum Connected Network selected.All results"));

        /** We draw the plot **/
        SnowPrinter snowPrinter = new SnowPrinter();
        File f = new File(outputFileName + "_size_pvalue.json");
        IOUtils.write(f.getAbsoluteFile(), snowPrinter.gsnowItemToJson(gsnowItems, numberOfMinNodes));
        List<String> tags = new ArrayList<String>();
        tags.add("NETWORKMINER_JSON");
        result.addOutputItem(new Item("plot_size_pvalue", f.getName(), "Plot", TYPE.IMAGE, tags, new HashMap<String, String>(2), "Results: Minimum Connected Network selected.All results"));
        return significantItemLocal;
    }

    //	private void printGSnowItems(List<GSnowItem> gsnowItemsList){
//		printGSnowItems(gsnowItemsList, true);
//	}
    private void printGSnowItems(List<GSnowItem> gsnowItemsList) {
        try {
            SnowPrinter snowPrinter = new SnowPrinter();
            File f = new File(outputFileName + "_all.txt");
            IOUtils.write(f.getAbsoluteFile(), snowPrinter.printGsnowItems(gsnowItemsList));
            result.addOutputItem(new Item("all_results", f.getName(), "All results", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Results: Minimum Connected Network selected.All results"));

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void drawResults(GSnowItem significantItem) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {

        List<Node> significantList = gsnowItems.get(significantItem.getNodes().size()).getNodes();
//		System.out.println("significantItem.getNodes(): "+significantItem.getNodes().size());
//		for(Node n : significantItem.getNodes())
//			System.out.println(n.getOriginalId());
//		
//		System.out.println("gsnowItems.get(significantItem.getNodes().size()).getNodes(): "+gsnowItems.get(significantItem.getNodes().size()).getNodes().size());
//		for(Node n : gsnowItems.get(significantItem.getNodes().size()).getNodes())
//			System.out.println(n.getOriginalId());


        List<ProteinVertex> nodes = nodes2ProteinVertices(significantList);
        Set<String> intermediates = new HashSet<String>();
//		List<List<ProteinVertex>> components = new ArrayList<List<ProteinVertex>>();
        SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), nodes);
        if (intermediate) {
            System.out.println("Antes intermediario [" + subgraph.getVertices().size() + "]: " + subgraph.getVertices().size());
            //boolean found = false;
            intermediates = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
            for (ProteinVertex vertex : nodes) {
                if (intermediates.contains(vertex.getId())) {
                    intermediates.remove(vertex.getId());
                }
            }
            for (String id : intermediates) {
                /** Get the input_id identifier for the intermediates **/
                String convertedId = id;
                if (this.listInfo.getOriginalDbName() == null)
                    continue;
//                if (xrefDBMan.getDBConnector().getDbConnection().getDatabase() != null) {
                if (this.interactome != null) {
                    List<String> dbNames = new ArrayList<String>();
                    String originalDbName = this.listInfo.getOriginalDbName().getDbname();
//					dbNames.add(originalDbName);
//					XRef xrefEns  = xrefDBMan.getByDBName(id, dbNames);
//					List<XRefItem> xrefitemList = xrefEns.getXrefItems().get(originalDbName);
//
                    List<String> list = new ArrayList<String>();
                    list.add(id);
//					xrefList = xrefDBMan.getAllIdentifiersByIds(list);
                    XrefManager xrefManager = new XrefManager(list, this.interactome);
                    Map<String, List<String>> xrefList = xrefManager.getXrefs(originalDbName);
                    if (!xrefList.get(id).isEmpty())
                        convertedId = xrefList.get(id).get(0);

//					if(xrefitemList != null && !xrefitemList.isEmpty()){
//						for(XRefItem xrefitem : xrefitemList){
//							convertedId = xrefitem.getDisplayName();
//						}
//					}
                }
                mapNames.put(id, convertedId);
            }
            System.out.println("Despues intermediario [" + subgraph.getVertices().size() + "]: " + subgraph.getVertices().size());
        }

//		components = subgraph.getAllInformationComponents(true);
        updateJobStatus("35", "List preprocessed");


//		if(snow){
//			/** Number of components with more than 1 node **/
//			int compMoreThan1Node = componentsMoreThanOneNode(components);
//			result.addOutputItem(new Item("comp_more_than_1_node", Integer.toString(compMoreThan1Node), "Number of components with more than 1 node", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));
//			int bicomponents = subgraph.getNumberOfBicomponents();
//			result.addOutputItem(new Item("bicomponents", Integer.toString(bicomponents), "Number of Bicomponents", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));
//		}

        ProteinNetwork subProteinNetwork = new ProteinNetwork(subgraph);
        subProteinNetwork.calcTopologicalValues();
        subProteinNetwork.calcTopologicalMeanValues();
        File f;

        if (subProteinNetwork != null) {
            Set<String> convertedIntermediates = new HashSet<String>();
            for (String key : intermediates) {
                convertedIntermediates.add(mapNames.get(key));
            }
            Set<String> convertedSeed = new HashSet<String>();
            for (String key : listInfo.getSeedNodes()) {
                convertedSeed.add(mapNames.get(key));
            }

            f = new File(outputFileName + "_mcn.sif");
            IOUtils.write(f.getAbsoluteFile(), graphToSif(significantItem, subgraph));
            addOutputSvgViewer(f, 1, convertedIntermediates, convertedSeed);
        }

        String googleTable = "NETWORKMINERRANKED_TABLE";
        if (!this.listInfo.isRanked())
            googleTable = "NETWORKMINERNOTRANKED_TABLE";

        if (intermediate && intermediates != null) {
            SnowPrinter snowPrinter = new SnowPrinter();
            f = new File(outputFileName + "_external_nodes_added.txt");
            IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodes(intermediates));
            result.addOutputItem(new Item("external_nodes_list", f.getName(), "Whole External nodes added", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Results: Minimum Connected Network selected.Interactors"));
        }

        f = new File(outputFileName + "_mcn_interactors.txt");
        IOUtils.write(f.getAbsoluteFile(), this.getMcnInteractors(subProteinNetwork, intermediates));                          //StringUtils.toList("TABLE,"+googleTable+",GO_TERM,KEGG_TERM", ",")
        result.addOutputItem(new Item("mcn_interactors", f.getName(), "Minimun Connected Network interactors", Item.TYPE.FILE, StringUtils.toList("TABLE," + googleTable, ","), new HashMap<String, String>(), "Results: Minimum Connected Network selected.Interactors"));


        /** Here we will put the redirections**/

        /**    Redirect to SNOW **/
        String listOfInteractors = "";
        listOfInteractors = significantItem.getNodesOriginalIds();
        redirect("REDIRECT_TOOL_SNOW", outputFileName + "_interactorsid_list_snow.txt", listOfInteractors);

        listOfInteractors = ListUtils.toString(getListOfInteractors(subgraph), System.getProperty("line.separator"));
        redirect("REDIRECT_TOOL_FATIGO", outputFileName + "_interactorsid_list_fatigo.txt", listOfInteractors);

        //		if(snow){
//			System.out.println("Excuting statsTest");
//			statsTest(subProteinNetwork, components);
//		}
//		else{
//			System.out.println("Not Excuting statsTest");
//		}

    }

    private void redirect(String tool, String fileName, String fileContent) throws IOException {
        String genomeFileName = writeGenomeInServer();
        File f = new File(fileName);
        IOUtils.write(f.getAbsoluteFile(), fileContent);

        List<String> tags = new ArrayList<String>();
        tags.add(tool);
        tags.add(this.interactome);
        tags.add(this.type);
        tags.add(this.group);
        tags.add(genomeFileName);
        String intermediate = this.intermediate ? "1" : "0";
        tags.add(intermediate);
        result.addOutputItem(new Item("mcn_id_interactors", f.getName(), "Interactors", Item.TYPE.FILE, tags, new HashMap<String, String>(), "Results: Minimum Connected Network selected.Continue processing"));

    }

    private String writeGenomeInServer() throws IOException {
        /** This method gets the list of the INTERCTORS if the genome and write them in a file
         * tener mucho cuidado xq es exactamente el mismo nombre que luego va a leer el newtworkwidget.js, ahí tb hay que cambiarlo**/
        String versionGenome = "ene11";
        String genomeFileName = "interactome_interactors_" + this.interactome + "_" + versionGenome + "_" + this.type + "_" + this.group + ".txt";
        File f = new File(outdir + "/" + genomeFileName);
        List<String> genomeInteractors = getListOfInteractors(proteinNetwork.getInteractomeGraph());
        IOUtils.write(f, genomeInteractors);
        return genomeFileName;
    }

    /**
     * This method gets the list of the inputs id + external ids
     */
    private List<String> getListOfInteractors(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph) {

        List<String> graph = new ArrayList<String>();
        for (DefaultEdge e : subgraph.getAllEdges()) {
            String source = subgraph.getVertex(e.getSource()).getId();
            String target = subgraph.getVertex(e.getTarget()).getId();
            if (mapNames.get(source) != null)
                source = mapNames.get(source);
            if (mapNames.get(target) != null)
                target = mapNames.get(target);
            graph.add(source);
            graph.add(target);
        }
        return graph;
    }

    private String getMcnInteractors(ProteinNetwork subProteinNetwork, Set<String> intermediates) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {

        StringBuilder sb = new StringBuilder();
        String tab = "\t";
        sb.append("#input_id").append(tab).append("id").append(tab).append("type").append(tab).append("rank").append(tab).append("bet").append(tab).append("clust").append(tab);
        sb.append("conn").append(tab).append("go_name").append(tab).append("go_id").append(System.getProperty("line.separator"));


        /** First significant no external nodes **/
        for (Node node : this.significantItem.getNodes()) {
            sb.append(node.getOriginalId()).append(tab);
            sb.append(node.getId()).append(tab);
            String listType = node.isSeed() ? "seedlist" : "list";
            sb.append(listType).append(tab);
            sb.append(node.getValue()).append(tab);
            sb.append(getSingleMcnInteractors(subProteinNetwork, node.getId()));
            sb.append(System.getProperty("line.separator"));
        }
        /** Second significant external nodes **/
        for (String id : intermediates) {
            String convertedId = "-";
            if (mapNames.get(id) != null)
                convertedId = mapNames.get(id);
            sb.append(convertedId).append(tab);
            sb.append(id).append(tab);
            sb.append("external").append(tab);
            /** Input value **/
            sb.append("-").append(tab);
            sb.append(getSingleMcnInteractors(subProteinNetwork, id));
            sb.append(System.getProperty("line.separator"));
        }
        return sb.toString();
    }

    private String getSingleMcnInteractors(ProteinNetwork subProteinNetwork, String id) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
//        GODBManager gof = new GODBManager(dbConnector);
        StringBuilder sb = new StringBuilder();
        String tab = "\t";
        ProteinVertex prVertex = subProteinNetwork.getInteractomeGraph().getVertex(id);
        sb.append(prVertex.getRelativeBetweenness()).append(tab);
        sb.append(prVertex.getClusteringCoefficient()).append(tab);
        sb.append(subProteinNetwork.getInteractomeGraph().getDegreeOf(prVertex)).append(tab);

        if (this.interactome != null) {


            List<String> ids = new ArrayList<String>();
            ids.add(id);
            XrefManager xrefManager = new XrefManager(ids, this.interactome);
            Map<String, List<String>> xrefs = xrefManager.getXrefs("go");
            String goNames = "";
            String goIds = "";
            if (!xrefs.isEmpty()) {
                for (String goId : xrefs.get(id)) {
//                System.out.println("key = " + key);
//                System.out.println("xrefs.get(key) = " + xrefs.get(key));
                    goIds += goId + "|";
                    String goTerm = "";
                    if(goTerms.get(goId)!=null)
                        goTerm = goTerms.get(goId).get(0);

                    goNames +=  goTerm + "|";
                }
            }
            sb.append(goNames);
            sb = this.deleteLastCh(sb, "|");
            sb.append("\t");
            sb.append(goIds);
            sb = this.deleteLastCh(sb, "|");


//			List<String> dbNames = new ArrayList<String>();
//			dbNames.add("go");
//			dbNames.add("kegg");
//			AnnotationDBManager annotationMng = new AnnotationDBManager(dbConnector);
//			Map<String, String> keggItems  = annotationMng.getAnnotationTermNames("kegg");
//			XRef xrefEns  = xrefDBMan.getByDBName(id, dbNames);
//			List<XRefItem> xrefitemListGO = xrefEns.getXrefItems().get("go");
//			if(xrefitemListGO != null && !xrefitemListGO.isEmpty()){
//				for(XRefItem xrefitem : xrefitemListGO){
//					GO go = gof.getByAccesion(xrefitem.getDisplayName());
//					sb.append(go.getName()).append(",");
//				}
//				sb = this.deleteLastCh(sb, ",");
//			}
//
//			List<XRefItem> xrefitemListKegg = xrefEns.getXrefItems().get("kegg");
//			if(xrefitemListKegg != null && !xrefitemListKegg.isEmpty()){
//				sb.append(tab);
//				for(XRefItem xrefitem : xrefitemListKegg){
//					if(keggItems.containsKey(xrefitem.getDisplayName())){
//						sb.append(keggItems.get(xrefitem.getDisplayName())).append(",");
//					}
//				}
//				sb = this.deleteLastCh(sb, ",");
//			}
        }
        return sb.toString();
    }

    /**
     * Esta es la parte del Snow que faltaría por ejecutar, en el caso de que se quisiera seguir el análisis de NM con un SNOW *
     */
    private void statsTest(ProteinNetwork subProteinNetwork, List<List<ProteinVertex>> components) throws IOException {
        //2nd Analysis
        logger.debug("Starting 2nd Analysis..................");
        List<Double> relBetSubnet1 = subProteinNetwork.getBetRelList();
        List<Double> relBetRandoms = new ArrayList<Double>();

        List<Double> connSubnet1 = subProteinNetwork.getConnList();
        List<Double> connRandoms = new ArrayList<Double>();

        List<Double> clustSubnet1 = subProteinNetwork.getClustList();
        List<Double> clustRandoms = new ArrayList<Double>();
        RandomsSnowManager randomsManager = new RandomsSnowManager(this.proteinNetwork, intermediate, true);

//		RandomsSnowManager randomsManager = new RandomsSnowManager(subProteinNetwork, intermediate, true);
        List<ProteinNetwork> subNetworks = randomsManager.createRandoms(randoms, this.significantItem.getNodes().size());
        double[] componentsRandoms = new double[randoms];
        int randonNumber = 0;
        for (ProteinNetwork proteinNetwork : subNetworks) {
            int randomVertex = (int) (Math.random() * proteinNetwork.getInteractomeGraph().getVertices().size());
            ProteinVertex v = proteinNetwork.getInteractomeGraph().getVertices().get(randomVertex);
            relBetRandoms.add(v.getRelativeBetweenness());
            connRandoms.add(Double.parseDouble(Integer.toString(proteinNetwork.getInteractomeGraph().getDegreeOf(v))));
            clustRandoms.add(v.getClusteringCoefficient());
            componentsRandoms[randonNumber] = (double) proteinNetwork.getInteractomeGraph().getAllInformationComponents(true).size();
            randonNumber++;
        }

        int[] range = getRange(componentsRandoms);
        String rangeString = "[" + Integer.toString(range[0]) + ", " + Integer.toString(range[1]) + "]";
        result.addOutputItem(new Item("components_number", components.size() + " " + rangeString, "Number of components [95% confidence interval]", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Minimun Connected Network selected. Topology description"));

        String toWrite = ksTest(relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms/*, side*/);
        if (!toWrite.equals("")) {
            String symbol;

            KolmogorovSmirnovTestResult resultKs;

            resultKs = getPValue(relBetSubnet1, relBetRandoms/*, side*/);
            symbol = getSymbol(resultKs.getSide());
            String stringResult = StringUtils.decimalFormat(resultKs.getPValue(), decimalFormat);
            //if(resultKs.getPValue() == 0)
            //	stringResult = "0.0001";
            result.addOutputItem(new Item("sn_random_kol_param_bet", stringResult, "Relative betweenness: Subnet " + symbol + " Random pval ", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Minimun Connected Network selected. Topology description"));
            createImages(outputFileName + "_sn_random_relBet", relBetSubnet1, "subnet1", relBetRandoms, "randoms", "sn_random_relBet", "Plot", "Minimun Connected Network selected. Topology description", "Relative Betweenness");

            resultKs = getPValue(connSubnet1, connRandoms);
            symbol = getSymbol(resultKs.getSide());
            stringResult = StringUtils.decimalFormat(resultKs.getPValue(), decimalFormat);
            //if(resultKs.getPValue() == 0)
            //	stringResult = "0.0001";
            result.addOutputItem(new Item("sn_random_kol_param_conn", stringResult, "Connection degree: Subnet " + symbol + " Random pval ", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Minimun Connected Network selected. Topology description"));
            createImages(outputFileName + "_sn_random_conn", connSubnet1, "subnet1", connRandoms, "randoms", "sn_random_conn", "Plot", "Minimun Connected Network selected. Topology description", "Connections");

            resultKs = getPValue(clustSubnet1, clustRandoms);
            symbol = getSymbol(resultKs.getSide());
            stringResult = StringUtils.decimalFormat(resultKs.getPValue(), decimalFormat);
            //if(resultKs.getPValue() == 0)
            //	stringResult = "0.0001";
            result.addOutputItem(new Item("sn_random_kol_param_clu", stringResult, "Clustering coefficient: Subnet " + symbol + " Random pval ", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Minimun Connected Network selected. Topology description"));
            createImages(outputFileName + "_sn_random_clust", clustSubnet1, "subnet1", clustRandoms, "randoms", "sn_random_clust", "Plot", "Minimun Connected Network selected. Topology description", "Clustering coeff");

            File f = new File(outputFileName + "_sn_random_kol.txt");
            IOUtils.write(f.getAbsolutePath(), toWrite);
            result.addOutputItem(new Item("sn_random_kol_param", f.getName(), "File", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Minimun Connected Network selected. Topology description"));
        } else
            result.addOutputItem(new Item("sn_random_kol_param", "Empty results", "Subnet - Random", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String, String>(), "Minimun Connected Network selected. Topology description"));
        logger.debug("Finished 2nd Analysis..................");
    }

    private String graphToSif(GSnowItem significantItem, SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph) {
        StringBuilder sb = new StringBuilder();
        /** Vertices with edges **/
        for (DefaultEdge e : subgraph.getAllEdges()) {
            String source = subgraph.getVertex(e.getSource()).getId();
            String target = subgraph.getVertex(e.getTarget()).getId();
            if (mapNames.get(source) != null)
                source = mapNames.get(source);
            if (mapNames.get(target) != null)
                target = mapNames.get(target);
            sb.append(source).append("\t").append("pp").append("\t").append(target);
            sb.append(System.getProperty("line.separator"));

        }
        /** Isolated vertices **/
        for (ProteinVertex vertex : subgraph.getAllVertices()) {
            if (subgraph.getDegreeOf(vertex) == 0) {
                for (Node n : significantItem.getNodes()) {
                    if (n.getId().equalsIgnoreCase(vertex.getId()))
                        sb.append(n.getOriginalId()).append(System.getProperty("line.separator"));
                }
            }
        }
        if (!sb.toString().equals(""))
            sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
        return sb.toString();
    }

    public void executeGenerateRandoms() {
        try {
            SnowPrinter snowPrinter = new SnowPrinter();
            int sizeMin = Integer.parseInt(commandLine.getOptionValue("size-min"));
            int sizeMax = Integer.parseInt(commandLine.getOptionValue("size-max"));
            RandomsSnowManager rndManager = new RandomsSnowManager(proteinNetwork, intermediate, false);

//			List<List<Double>> allGsnowValues = new ArrayList<List<Double>>();
            List<List<List<ProteinVertex>>> components;
            List<Double> gsnowValues;
            List<ProteinNetwork> list;

            double tInicio = 0.0;
            double tFinal = 0.0;
            double tTotal = 0.0;
            IOUtils.append(outputFileName, snowPrinter.printGSnowValuesHeader(randoms));
            for (int i = sizeMin; i <= sizeMax; i++) {
                tInicio = System.currentTimeMillis();
                System.out.println("Creating randoms with size: " + i);
                list = rndManager.createRandoms(randoms, i);
                components = rndManager.getComponents(list);
                gsnowValues = rndManager.getGSnowValues(list, components);
                IOUtils.append(outputFileName, i + "\t" + snowPrinter.printGSnowValues(gsnowValues));
                System.gc();
//				allGsnowValues.add(rndManager.getGSnowValues(list, components));
                tFinal = System.currentTimeMillis();
                tTotal = (tFinal - tInicio) / 1000;
                System.out.println("\t spend time: " + tTotal);
            }
//			gsnowValues2File(allGsnowValues, sizeMin, randoms);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
//	private void gsnowValues2File(List<List<Double>> allGsnowValues, int sizeMin, int randoms) throws IOException{
//		SnowPrinter snowPrinter = new SnowPrinter();
//		String values = snowPrinter.printAllGSnowValues(allGsnowValues, sizeMin, randoms);
//		IOUtils.write(outputFileName, values);
//	}

    private String loadFile() {
        String folder = this.babelomicsHomePath + "/conf/interactions/gsnow/";
        String intermediateString = intermediate ? "intermediate" : "nointermediate";
        String localType = type;
        if (type.equals("transcripts")) {
            localType = "proteins";
            logger.debug("Transcripts recognised, loading " + localType + " randoms");
        }
        if (type.equals("vcf")) {
            localType = "genes";
            logger.debug("vcf recognised, loading " + localType + " randoms");
        }

        folder += interactome + "_" + localType + "_" + group + "db_" + intermediateString + ".txt";
        return folder;

    }

    public class GSnowItem implements Comparable<GSnowItem> {

        /**
         * List<Node> nodes; *
         */
        List<Node> nodes;

        /**
         * This is the value(vertices/components) before compared with randoms *
         */
        Double rawValue;

        /**
         * This is the value(mean of significant values) after compared with randoms *
         */
        Double comparedValue;

        /**
         * This is the value :(1-rawValue)/(1-List<Node>.size) *
         */
        Double score;

        /**
         * This says if this is a significant item or not *
         */
        boolean significant;

        List<Double> componentsSize;

        public GSnowItem() {
            nodes = new ArrayList<Node>();
            componentsSize = new ArrayList<Double>();
            significant = false;
            comparedValue = null;
            rawValue = null;
            score = null;

        }

        public List<Node> getNodes() {
            return nodes;
        }

        public void setNodes(List<Node> nodes) {
            this.nodes = nodes;
        }

        public Double getRawValue() {
            return rawValue;
        }

        public void setRawValue(Double rawValue) {
            this.rawValue = rawValue;
        }

        public Double getComparedValue() {
            return comparedValue;
        }

        public void setComparedValue(Double comparedValue) {
            this.comparedValue = comparedValue;
        }

        public Double getScore() {
            return score;
        }

        public void setScore(Double score) {
            this.score = score;
        }

        public List<Double> getComponentsSize() {
            return componentsSize;
        }

        public void setComponentsSize(List<Double> componentsSize) {
            this.componentsSize = componentsSize;
        }

        public boolean isSignificant() {
            return significant;
        }

        public void setSignificant(boolean significant) {
            this.significant = significant;
        }

        @Override
        public String toString() {
            return "{" + nodes.toString() + " - rawValue = " + rawValue + " - comparedValue = " + comparedValue + "}";
        }

        public String getNodesOriginalIds() {
            String separator = System.getProperty("line.separator");
            return getNodesOriginalIds(separator);
        }

        public String getNodesOriginalIds(String separator) {
            StringBuilder sb = new StringBuilder();
            for (Node n : nodes) {
                sb.append(n.getOriginalId()).append(separator);
            }
            if (sb.toString().equals(""))
                return sb.toString();
            return sb.substring(0, sb.toString().length() - 1);
        }

        public String getNodesIds() {
            String separator = System.getProperty("line.separator");
            return getNodesIds(separator);
        }

        public String getNodesIds(String separator) {
            StringBuilder sb = new StringBuilder();
            for (Node n : nodes) {
                sb.append(n.getId()).append(separator);
            }
            if (sb.toString().equals(""))
                return sb.toString();
            return sb.substring(0, sb.toString().length() - 1);
        }

        @Override
        public int compareTo(GSnowItem o) {
            if (score < o.score)
                return 1;
            return 0;
        }
    }


}


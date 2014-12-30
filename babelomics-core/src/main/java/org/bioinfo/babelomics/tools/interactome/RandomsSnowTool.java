package org.bioinfo.babelomics.tools.interactome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.tool.OptionFactory;

public class RandomsSnowTool extends SnowTool{
	
//	private int randoms;
	private int randomSize;
	private boolean topoValues;
	private SnowPrinter snowPrinter;
	private RandomsSnowManager rndManager;
	//List<Double> gsnowValues;
	
	@Override
	public void initOptions() {
		
		options.addOption(OptionFactory.createOption("random-size", "Randoms minimal connection network size", true, true));
		options.addOption(OptionFactory.createOption("topo-values", "Calculate also the topological values for each node", false, true));
		
	}

	@Override
	protected void execute() {
		initExecute();
		//gsnowValues = new ArrayList<Double>(randoms);
		
		
		snowPrinter = new SnowPrinter();
		randoms = Integer.parseInt(commandLine.getOptionValue("randoms"));	
		randomSize = Integer.parseInt(commandLine.getOptionValue("random-size"));
		topoValues = (commandLine.hasOption("topo-values") && !"0".equalsIgnoreCase(commandLine.getOptionValue("topo-values")));
		
		
		rndManager = new RandomsSnowManager(proteinNetwork, intermediate, topoValues);
		subProteinNetworkRandoms = new ArrayList<ProteinNetwork>(randoms);
		subProteinNetworkRandoms = rndManager.createRandoms(randoms, randomSize);
		try {
			if(topoValues){
				this.meansValuesList2File(outputFileName+"_sn_1-"+(randoms)+"_means.txt");
				this.topologicalValuesList2File(outputFileName+"_sn_1-"+(randoms)+"_topo.txt");
			}
			if(components)
				this.componentsList2File(outputFileName+"_sn_1-"+(randoms)+"_comp.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	public void meansValuesList2File(String fileName) throws IOException{
		IOUtils.write(fileName,printMeansValuesList());
	}
	public String printMeansValuesList(){
		StringBuilder sb = new StringBuilder();
		sb.append(snowPrinter.createMeansHeader());
		int subnet = 1;
		for(ProteinNetwork proteinNetwork : subProteinNetworkRandoms){
			sb.append(snowPrinter.printMeanValues(proteinNetwork, subnet, false)).append(System.getProperty("line.separator"));
			subnet++;
		}
		sb = this.deleteLastLine(sb);
		return sb.toString();
	}
	public void topologicalValuesList2File(String fileName) throws IOException{
		IOUtils.write(fileName, printTopologicalValuesList());
	}
	public String printTopologicalValuesList(){
		StringBuilder sb = new StringBuilder();
		sb.append(snowPrinter.createTopoHeader());
		int subnet = 1;
		for(ProteinNetwork proteinNetwork : subProteinNetworkRandoms){
			sb.append(snowPrinter.printTopologicalValues(proteinNetwork, subnet, false)).append(System.getProperty("line.separator"));
			subnet++;
		}
		sb = this.deleteLastLine(sb);
		return sb.toString();
	}
	public void componentsList2File(String fileName) throws IOException{
		IOUtils.write(fileName, this.printComponentsList());
	}
	private String printComponentsList(){
		List<List<List<ProteinVertex>>> componentsList = rndManager.getComponents(subProteinNetworkRandoms);
		List<List<Double>> diameters = rndManager.getDiameters(subProteinNetworkRandoms, componentsList);
		return snowPrinter.printComponents(componentsList, diameters);
	}
}


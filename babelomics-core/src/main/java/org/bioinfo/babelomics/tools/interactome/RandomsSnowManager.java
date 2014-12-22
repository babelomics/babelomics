package org.bioinfo.babelomics.tools.interactome;

import java.util.ArrayList;
import java.util.List;

import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.Subgraph;
import org.bioinfo.data.graph.alg.Calc;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinVertex;

public class RandomsSnowManager {

	private ProteinNetwork proteinNetwork;
	private boolean intermediate;
	private boolean topoValues;
	
	public RandomsSnowManager(ProteinNetwork proteinNetwork, boolean intermediate, boolean topoValues) {
		this.proteinNetwork = proteinNetwork;
		this.intermediate = intermediate;
		this.topoValues = topoValues;
	}
	
	public List<ProteinNetwork> createRandoms(int randoms, int randomSize){
		List<ProteinNetwork> subProteinNetworkRandoms = new ArrayList<ProteinNetwork>();
		for(int i=0; i<randoms; i++){
			ProteinNetwork subProteinNetwork = createRandom(randomSize);
			subProteinNetworkRandoms.add(subProteinNetwork);
		}
		return subProteinNetworkRandoms;
	}
	public ProteinNetwork createRandom(int randomSize){
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), randomSize);
		if(intermediate) {
			Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
		}
		ProteinNetwork subProteinNetwork = createSubnet(subgraph);
		return subProteinNetwork;
	}
	private ProteinNetwork createSubnet(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph) {
		ProteinNetwork subProteinNetwork = new ProteinNetwork(subgraph);
		if(topoValues){
			subProteinNetwork.calcTopologicalValues();
			subProteinNetwork.calcTopologicalMeanValues();
		}
		return subProteinNetwork;
	}
	
	public List<List<List<ProteinVertex>>> getComponents(List<ProteinNetwork> subProteinNetworkRandoms){
		
		List<List<ProteinVertex>> componentsList;
		List<List<List<ProteinVertex>>> componentsProteinNetwork = new ArrayList<List<List<ProteinVertex>>>(subProteinNetworkRandoms.size());
		for(ProteinNetwork proteinNetwork : subProteinNetworkRandoms){
			 componentsList  = getComponents(proteinNetwork);
			 componentsProteinNetwork.add(componentsList);
		}
		return componentsProteinNetwork;
	}
	public List<List<ProteinVertex>> getComponents(ProteinNetwork subProteinNetwork){
		List<List<ProteinVertex>> componentsList;
		componentsList  = subProteinNetwork.getInteractomeGraph().getAllInformationComponents(true);
		return componentsList;
	}
	public List<List<Double>> getDiameters(List<ProteinNetwork> subProteinNetworkRandoms, List<List<List<ProteinVertex>>> componentsList){
		List<List<Double>> diameters = new ArrayList<List<Double>>();
		for(int i=0; i < subProteinNetworkRandoms.size(); i++){
			diameters.add(getDiameters(subProteinNetworkRandoms.get(i), componentsList.get(i)));
		}
		return diameters;
	}
	public List<Double> getDiameters(ProteinNetwork subProteinNetwork, List<List<ProteinVertex>> componentsList){
		return Calc.calcDiameter(subProteinNetwork.getInteractomeGraph(), componentsList);
	}
	
	
	public List<Double> getGSnowValues(List<ProteinNetwork> subProteinNetworkRandoms, List<List<List<ProteinVertex>>> componentsProteinNetwork){
		List<Double> gsnowValues = new ArrayList<Double>(componentsProteinNetwork.size());
		double gsnowValue = 0.0;
		
		for(int i=0; i < subProteinNetworkRandoms.size(); i++){
			gsnowValue = (double)subProteinNetworkRandoms.get(i).getInteractomeGraph().getVertices().size() / (double)componentsProteinNetwork.get(i).size();
			gsnowValues.add(gsnowValue);
		}
		return gsnowValues;
	}
}


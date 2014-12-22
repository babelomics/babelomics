package org.bioinfo.babelomics.tools.interactome;

import java.util.ArrayList;
import java.util.List;

import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.Subgraph;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinVertex;

public class SnowRandomUtils {
	protected boolean intermediate;
	protected boolean topoValues;
	protected ProteinNetwork proteinNetwork;
	
	public SnowRandomUtils(ProteinNetwork proteinNetwork, boolean intermediate, boolean topoValues){
		this.proteinNetwork = proteinNetwork;
		this.intermediate = intermediate;
		this.topoValues = topoValues;
	}
	public List<ProteinNetwork> createRandoms(int randoms, int randomSize){
		List<ProteinNetwork> subProteinNetworkRandoms = new ArrayList<ProteinNetwork>();
		for(int i=1; i<=randoms; i++){
			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), randomSize);
			if(intermediate) {
				Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
			}
			//System.out.println("Randoms["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
			ProteinNetwork subProteinNetwork = createSubnet(subgraph);
			subProteinNetworkRandoms.add(subProteinNetwork);
			//System.out.println(subProteinNetwork.toString());
		}
		return subProteinNetworkRandoms;

	}
	private ProteinNetwork createSubnet(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph) {
		ProteinNetwork subProteinNetwork = new ProteinNetwork(subgraph);
		if(topoValues){
			subProteinNetwork.calcTopologicalValues();
			subProteinNetwork.calcTopologicalMeanValues();
		}
		return subProteinNetwork;
	}
}


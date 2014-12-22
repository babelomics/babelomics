package org.bioinfo.babelomics.tools.interactome;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bioinfo.babelomics.tools.interactome.gsnow.GSnow.GSnowItem;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.Node;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinVertex;

import com.google.gson.Gson;

public class SnowPrinter {

	private boolean components;
	private boolean bicomponents;
	private String lineSeparator;
	private String tab;
	
	public SnowPrinter(){
		lineSeparator = System.getProperty("line.separator");
		tab = "\t";
	}
	
	
	public String printComponents(List<List<List<ProteinVertex>>> compList, List<List<Double>> compsDiameter){
		StringBuilder sb = new StringBuilder();
		sb.append(this.createComponentsHeader());
		int subnet = 0;
		for(int i=0; i < compList.size(); i++){
			subnet++;
			sb.append(printComponents(subnet, compList.get(i), compsDiameter.get(i), false));
		}
		return sb.toString();
	}
	
	public String printComponents(int subnet, List<List<ProteinVertex>> componentsList,List<Double> componentsDiameter,  boolean header){
		StringBuilder sb = new StringBuilder();
		if(header)
			sb.append(this.createComponentsHeader());
		for(int i=0; i < componentsList.size(); i++){
			sb.append("sn"+(subnet)).append("\t");
			StringBuilder sbNodes = new StringBuilder();
			for (int k = 0; k < componentsList.get(i).size(); k++) {
				String nodeId = "";
				/****
				 *  Hay que tener cuidado con esto si al final acoplo Snow.java a estas nuevas clases, por la conversión de genes/proteínas/transcritos
				 *	La variable randoms se pasa por parámetro.
				 *
				 * if(type.equals("transcripts") || type.equals("genes")){
				 *	nodeId = getValue(componentsList.get(i).get(k).getId());
				 * }
				*****/
				sbNodes.append(nodeId+componentsList.get(i).get(k).getId());
				if(k!=componentsList.get(i).size()-1)
					sbNodes.append(",");
			}
			sb.append((i)+"\t"+componentsDiameter.get(i)+"\t"+componentsList.get(i).size()+"\t"+sbNodes.toString()+ lineSeparator);
		}
		return sb.toString();
	}
	
	
	
	public String printTopologicalValues(ProteinNetwork subProteinNetwork, int subnet, boolean header){
		StringBuilder sb = new StringBuilder();
		if(header)
			sb.append(this.createTopoHeader());
		for (ProteinVertex proteinVertex : subProteinNetwork.getInteractomeGraph().getVertices()) {
			if(proteinVertex == null) 
				continue;
			String inputId= "";
			sb.append("sn"+(subnet)).append("\t");
			/****
			 *  Hay que tener cuidado con esto si al final acoplo Snow.java a estas nuevas clases, por la conversión de genes/proteínas/transcritos
			 *	La variable randoms se pasa por parámetro.
			 *
			 *	if( (type.equals("transcripts") || type.equals("genes"))  && random == false){
			 *	inputId = this.getValue(proteinVertex.getId());
			 *  }
			 ***/
			sb.append(inputId+proteinVertex.getId()).append("\t");
			sb.append(proteinVertex.getRelativeBetweenness()).append("\t").append(subProteinNetwork.getInteractomeGraph().getDegreeOf(proteinVertex)).append("\t").append(proteinVertex.getClusteringCoefficient());
			sb.append(lineSeparator);
		}
		if(!sb.toString().equals(""))
			sb.deleteCharAt(sb.lastIndexOf(lineSeparator));
		return sb.toString();
	}
	public String printMeanValues(ProteinNetwork subProteinNetwork, int subnet, boolean header){
		StringBuilder sb = new StringBuilder();
		if(header)
			sb.append(this.createMeansHeader());
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = subProteinNetwork.getInteractomeGraph();
		sb.append("sn"+(subnet)).append("\t");
		sb.append(subProteinNetwork.getTopologicalMeanValuesToString());

		if(components){
			int oneComponent = 0;
			int moreOneComponent = 0;
			List<List<ProteinVertex>> listComponents = interactomeGraph.getAllInformationComponents(true);
			for(List<ProteinVertex> listVertex : listComponents){
				if(listVertex.size() == 1)
					oneComponent++;
				else if(listVertex.size() > 1)
					moreOneComponent++;
			}
			sb.append("\t"+moreOneComponent+"\t"+oneComponent);
		}
		sb.append("\t").append(subProteinNetwork.getInteractomeGraph().getVertices().size());
		if(bicomponents){
			sb.append("\t").append(interactomeGraph.getNumberOfBicomponents());
		}
		return sb.toString();
	}
	public StringBuilder createTopoHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tId\tBet\tConn\tClust").append(lineSeparator);
		return sb;
	}
	
	public StringBuilder createComponentsHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tComp\tDia\tSize\tNod").append(lineSeparator);
		return sb;
	}
	
	public StringBuilder createMeansHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tMeanBet\tStdBet\t").append("MeanCon\tStdCon\t").append("MeanCls\tStdCls");
		if(components){
			sb.append("\tComp\t1Comp");
		}
		sb.append("\tNodes");
		if(bicomponents){
			sb.append("\tBiComp\t");
		}
		sb.append(lineSeparator);
		return sb;
	}
//	public String printAllGSnowValues(List<List<Double>> values, int init, int randoms){
//		StringBuilder sb = new StringBuilder();
//		sb.append(this.printGSnowValuesHeader(randoms));
//		
//		int i=init;
//		for(List<Double> list : values){
//			sb.append(i+"\t");
//			sb.append(printGSnowValues(list));
////			for(double value : list){
////				sb.append(value).append("\t");
////			}
//			sb = deleteLastCh(sb);
//			sb.append(lineSeparator);
//			i++;
//		}
//		sb = deleteLastCh(sb);
//		return sb.toString();
//	}
	public String printGSnowValues(List<Double> values){
		StringBuilder sb = new StringBuilder();
			
		for(double value : values){
			sb.append(value).append("\t");
		}
		sb = deleteLastCh(sb);
		return sb.toString();
	}
	
	public String printGSnowValuesHeader(int randoms){
		StringBuilder sb = new StringBuilder();
		sb.append("#size\t");
		for(int z=1; z<=randoms; z++)
			sb.append("rnd"+z+"\t");
		return sb.toString();
	}
	public StringBuilder deleteLastCh(StringBuilder sb){
		if(sb.toString().equals(""))
			return sb;
		return sb.deleteCharAt(sb.length()-1);
	}
	
//	public String printGsnowItems(Map<Integer, GSnowItem> gsnowItems){
//		List<Double> orderedList1 = new ArrayList<Double>();
//		StringBuilder sb1 = new StringBuilder();
//		for(int i : gsnowItems.keySet()){
//			orderedList1.add(gsnowItems.get(i).getScore());
//		}
//		GSnowItem gsnowItem1;
//		Collections.sort(orderedList1);
//		for(double d : orderedList1){
//			for(int i : gsnowItems.keySet()){
//				if(d == gsnowItems.get(i).getScore()){
//					gsnowItem1 = gsnowItems.get(i);
//					sb1.insert(0, printGsnowItem(gsnowItem1)+lineSeparator);
//					//sb1.append(printGsnowItem(gsnowItem1)).append(lineSeparator);
//				}
//			}
//		}
//		sb1.insert(0,"#Size"+tab+"Nodes"+tab+"P-value"+tab+"NodesNumber/ComponentsNumber"+tab+"Score"+lineSeparator);
//		sb1 = deleteLastCh(sb1);
//		return sb1.toString();
//		
//	}
	public String printGsnowItems(List<GSnowItem> gsnowItems){
		StringBuilder sb = new StringBuilder();
		sb.append("#Size"+tab+"Nodes"+tab+"P-value"+tab+"NodesNumber/ComponentsNumber"+tab+"Score"+lineSeparator);
		for (GSnowItem gSnowItem : gsnowItems) {
			sb.append(printGsnowItem(gSnowItem)).append(lineSeparator);
		}
		sb = deleteLastCh(sb);
		return sb.toString();
	}
//	public String printGsnowItem(GSnowItem gsnowItem){
//		return printGsnowItem(gsnowItem, true);
//	}
	public String printGsnowItem(GSnowItem gsnowItem){
		/** Cuando printPvalue es false, si el pvalue es 0, pondremos un - (se hace para la seed list)**/
		String pValue = "-";
		String rawValue = "-";
		String score="-";
		if(gsnowItem.getComparedValue() != null){
			pValue = Double.toString(gsnowItem.getComparedValue());
		}
		if(gsnowItem.getRawValue() != null){
			rawValue = Double.toString(gsnowItem.getRawValue());
		}
		if(gsnowItem.getScore() != null){
			score = Double.toString(gsnowItem.getScore());
		}
		String separator = ",";
		StringBuilder sb = new StringBuilder();
		sb.append(gsnowItem.getNodes().size()).append(tab);
		sb.append(gsnowItem.getNodesIds(separator)).append(tab);
		sb.append(pValue).append(tab);
		sb.append(rawValue).append(tab);
		sb.append(score);
		return sb.toString();
	}
	public String gsnowItemToJson(Map<Integer, GSnowItem> gsnowItems, int numberOfStartingNodes){
		Integer[] array = gsnowItems.keySet().toArray(new Integer[0]);
		List<Integer> orderedList = ArrayUtils.toList(array);
		Collections.sort(orderedList);
		List<GSnowItem> list = new ArrayList<GSnowItem>();
		for(int i : orderedList){
			if(i >= numberOfStartingNodes){
				list.add(gsnowItems.get(i));
			}
		}
		Gson gson = new Gson();
		return gson.toJson(list);
	}
	public String printNodesList(List<Node> nodes){
		StringBuilder sb = new StringBuilder();
		sb.append("#id\tinput_id\tvalue").append(lineSeparator);
		
		for(Node node : nodes){
			sb.append(node.toString()).append(lineSeparator);
		}
		sb = deleteLastCh(sb);
		return sb.toString();
	}
	public String printNodes(Set<String> nodes){
		StringBuilder sb = new StringBuilder();
		sb.append("#id").append(lineSeparator);
		
		for(String node : nodes){
			sb.append(node.toString()).append(lineSeparator);
		}
		sb = deleteLastCh(sb);
		return sb.toString();
	}
}


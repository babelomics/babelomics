package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.GeneDBManager;
import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.feature.Gene;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.infrared.core.feature.XRef.XRefItem;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;

public class InfraredUtils {

	// get annotations from a filter
	public static FeatureList<AnnotationItem> getAnnotations(DBConnector dbConnector, List<String> ids, FunctionalFilter filter) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException{
		AnnotationDBManager annotationMng = new AnnotationDBManager(dbConnector);		
		FeatureList<AnnotationItem> annots = annotationMng.getAnnotation(ids, filter);
		return annots;
		//		if(filter instanceof GOFilter){
		//			return annotationMng.getGOAnnotation(ids, (GOFilter) filter);
		//		} else if(filter instanceof KeggFilter){			
		//			return annotationMng.getKeggAnnotation(ids, (KeggFilter) filter);		
		//		} else if(filter instanceof ){
		//			return annotationMng.getKeggAnnotation(ids, (KeggFilter) filter);
		//		}
		//		System.err.println("is unknown filter");
		//		return null;
	}	

	// get genome (ensembl gene format)
	public static List<String> getGenome(DBConnector dbConnector){
		try {
			//			return new GeneDBManager(dbConnector).getAllEnsemblIds();
			System.err.println("all genes: " +  new GeneDBManager(dbConnector).getAllEnsemblIds().size());
			System.err.println("protein_coding: " +  new GeneDBManager(dbConnector).getAllByBiotype("protein_coding").getFeaturesIds().size());
			return new GeneDBManager(dbConnector).getAllByBiotype("protein_coding").getFeaturesIds();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	// get chromosome region genes (ensembl gene format)
	public static List<String> getChromosomeRegionGenes(DBConnector dbConnector, String chromosome, int start, int end) {
		try {
			FeatureList<Gene> genes = new GeneDBManager(dbConnector).getAllByRegion(chromosome, start, end);
			List<String> list = new ArrayList<String> (genes.size());
			for (Gene gene: genes) {
				list.add(gene.getId());
			}
			return list;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	// translate a list of ids to ensembl format
	public static List<String> toEnsemblId(DBConnector dbConnector, List<String> ids) {
		try {	
			//			List<List<String>> list = new XRefDBManager(dbConnector).getIdsByDBName(ids, "ensembl_gene");
			List<XRef> list = new XRefDBManager(dbConnector).getByDBName(ids, "ensembl_gene");
			List<String> ensemblIds = null;
			if(list != null) {				
				ensemblIds = new ArrayList<String>();
				//				for(List<String> stringList: list) {					
				//					if(stringList!=null){
				//						for(String ensemblId: stringList) {						
				//							ensemblIds.add(ensemblId);
				//						}
				//					}
				//				}
				for(XRef xref: list) {					
					if(xref != null){
						for(XRefItem xrefItem: xref.getXrefItems().get("ensembl_gene")) {						
							ensemblIds.add(xrefItem.getDisplayName());
						}
					}
				}
				ensemblIds = ListUtils.unique(ensemblIds);
			}
			return ensemblIds;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	// translate a list of ids to ensembl format
	public static Map<String, List<String>> getEnsemblMap(DBConnector dbConnector, List<String> ids) {
		try {
			//			List<List<String>> lists = new XRefDBManager(dbConnector).getIdsByDBName(ids, "ensembl_gene");
			List<XRef> list = new XRefDBManager(dbConnector).getByDBName(ids, "ensembl_gene");
			Map<String, List<String>> map = null;			
			if ( list != null ) {
				map = new HashMap<String, List<String>>();
				//				for(int i=0 ; i<list.size() ; i++) {
				//					if ( list.get(i) == null || list.get(i).size() <= 0) {
				//						map.put(ids.get(i), null);
				//					} else {
				//						if ( !map.containsKey(ids.get(i)) ) {
				//							map.put(ids.get(i), new ArrayList<String>());
				//						}
				//						map.get(ids.get(i)).addAll(list.get(i));
				//					}
				//				}
				for(int i=0 ; i<list.size() ; i++) {
					if(list.get(i) == null || list.get(i).getXrefItems().get("ensembl_gene").size() <= 0) {
						map.put(ids.get(i), null);
					} else {
						if(!map.containsKey(ids.get(i))) {
							map.put(ids.get(i), new ArrayList<String>());
						}
						for(XRefItem xrefItem: list.get(i).getXrefItems().get("ensembl_gene")) {
							map.get(ids.get(i)).add(xrefItem.getDisplayName());
						}
					}
				}
			}
			return map;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	// translate a list of ids to unigene format
	public static Map<String, List<String>> getUnigeneMap(DBConnector dbConnector, List<String> ids) {
		try {
			//			List<List<String>> lists = new XRefDBManager(dbConnector).getIdsByDBName(ids, "unigene");
			List<XRef> list  = new XRefDBManager(dbConnector).getByDBName(ids, "unigene");
			Map<String, List<String>> map = null;			
			if(list != null) {
				map = new HashMap<String, List<String>>();
				//				for(int i=0 ; i<list.size() ; i++) {
				//					if ( list.get(i) == null || list.get(i).size() <= 0) {
				//						map.put(ids.get(i), null);
				//					} else {
				//						if ( !map.containsKey(ids.get(i)) ) {
				//							map.put(ids.get(i), new ArrayList<String>());
				//						}
				//						map.get(ids.get(i)).addAll(list.get(i));
				//					}
				//				}
				for(int i=0 ; i<list.size() ; i++) {
					if(list.get(i) == null || list.get(i).getXrefItems().get("unigene").size() <= 0) {
						map.put(ids.get(i), null);
					} else {
						if(!map.containsKey(ids.get(i))) {
							map.put(ids.get(i), new ArrayList<String>());
						}
						for(XRefItem xrefItem: list.get(i).getXrefItems().get("unigene")) {
							map.get(ids.get(i)).add(xrefItem.getDisplayName());
						}
					}
				}
			}
			return map;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

}

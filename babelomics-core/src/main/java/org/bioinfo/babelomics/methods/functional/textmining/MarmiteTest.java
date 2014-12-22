package org.bioinfo.babelomics.methods.functional.textmining;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.Config;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.db.DBConnection;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.KolmogorovSmirnovTest;

public class MarmiteTest {
	

	DBConnection dbConn = null; 
	Config config = new Config();

	List<String> geneList1;
	List<String> geneList2;
	
	Map<String, List<Score>> entityMap1;
	Map<String, List<Score>> entityMap2;
	
	List<String> excludedEntities1 = new ArrayList<String>();
	List<String> excludedEntities2 = new ArrayList<String>();
	
	List<String> validEntities;

	public MarmiteTest(String propertiesFile) throws FileNotFoundException, IOException {
		config.append(new File(propertiesFile));
		dbConn = new DBConnection("mysql", config.getProperty("INFRARED.HOST"), config.getProperty("INFRARED.PORT"), config.getProperty("INFRARED.HSA.DATABASE"), config.getProperty("INFRARED.USER"), config.getProperty("INFRARED.PASSWORD"));
	}
	
	public TestResultList<KolmogorovSmirnovTestResult> run (List<String> list1, List<String> list2, String bioentity, int scoreFilter) throws InvalidParameterException, IOException {
		
		
		geneList1 = list1;
		geneList2 = list2;
		
		// accessing to the db to get entities
		//
		entityMap1 = MarmiteUtils.getEntityMap(geneList1, bioentity, dbConn);
		entityMap2 = MarmiteUtils.getEntityMap(geneList2, bioentity, dbConn);	


		if ( entityMap1.size() == 0 ) {
			throw new InvalidParameterException("No entities found for list #1");
		}
		
		if ( entityMap2.size() == 0 ) {
			throw new InvalidParameterException("No entities found for list #2");
		}
		
		List<String> entities = MapUtils.getKeys(entityMap1);
		entities.addAll(MapUtils.getKeys(entityMap2));
		System.out.println("------->>> sizes (map1, map2, all) = (" + entityMap1.size() + ", " + entityMap2.size() + ", " + entities.size() + ")");
		entities = ListUtils.unique(entities);
		System.out.println("------->>> sizes (unique all) = (" + entities.size() + ")");

		// filtering entities by score filter
		//
		String entity;
		List<Score> scoreList1, scoreList2;
		validEntities = new ArrayList<String>();
		for(int i=0 ; i<entities.size() ; i++) {
			entity = entities.get(i);
			scoreList1 = entityMap1.get(entity); 
			scoreList2 = entityMap2.get(entity); 
			if ( scoreList1 != null && scoreList1.size() >= scoreFilter && scoreList2 != null && scoreList2.size() >= scoreFilter ) {
				validEntities.add(entity);
			}
		}
		entities.clear();
	
		// saving excluded entities from both lists
		//
		for(String key: entityMap1.keySet()) {
			if ( ! validEntities.contains(key) ) excludedEntities1.add(key);
		}
		for(String key: entityMap2.keySet()) {
			if ( ! validEntities.contains(key) ) excludedEntities2.add(key);
		}

		System.out.println("------->>> sizes (filter all, excluded map1, exclude map2) = (" + validEntities.size() + ", " + excludedEntities1.size() + ", " + excludedEntities2.size() + ")");
		
		// creating score matrixes
		//
		double[][] scoreMatrix1 = new double[validEntities.size()][];
		double[][] scoreMatrix2 = new double[validEntities.size()][];
		for(int i=0 ; i<validEntities.size() ; i++) {
			entity = validEntities.get(i);
			scoreList1 = entityMap1.get(entity); 
			scoreList2 = entityMap2.get(entity); 

			scoreMatrix1[i] = new double[scoreList1.size()];
			for (int j=0 ; j<scoreList1.size(); j++) {
				scoreMatrix1[i][j] = scoreList1.get(j).getValue();
			}
			scoreMatrix2[i] = new double[scoreList2.size()];
			for (int j=0 ; j<scoreList2.size(); j++) {
				scoreMatrix2[i][j] = scoreList2.get(j).getValue();
			}	
		}
		
		// executing kolmogorov-smirnov test  
		//
		return new KolmogorovSmirnovTest().compute(scoreMatrix1, scoreMatrix2);
	}

	
	public List<String> getGeneList1() {
		return geneList1;
	}

	public List<String> getGeneList2() {
		return geneList2;
	}

	public Map<String, List<Score>> getEntityMap1() {
		return entityMap1;
	}

	public Map<String, List<Score>> getEntityMap2() {
		return entityMap2;
	}

	public List<String> getValidEntities() {
		return validEntities;
	}


	public List<String> getExcludedEntities1() {
		return excludedEntities1;
	}


	public List<String> getExcludedEntities2() {
		return excludedEntities2;
	}

}

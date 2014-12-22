package org.bioinfo.babelomics.methods.expression.clustering;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.data.tree.multiway.MultipleTree;
import org.bioinfo.math.data.DoubleMatrix;

public class Som extends Cluster {

	public Som(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance, String home) {
		super(matrix, rowNames, colNames, distance, home);
	}

	@Override
	public String run(boolean createClusterFiles) throws Exception {
		MultipleTree nw = null;

		File tmpDir = File.createTempFile("input", ".dir");
		tmpDir.delete();
		tmpDir.mkdir();

		int dimension = heuristica(rowNames.size()); 

		File inputFile = new File(tmpDir.getAbsoluteFile() + "/in.txt");
		File outputFile = new File(tmpDir.getAbsoluteFile() + "/in_SOM_G" + dimension + "-" + dimension + ".gnf");

		//		System.out.println("(infile, outfile) = (" + inputFile.getAbsolutePath() + ", " + outputFile.getAbsolutePath() + ")");

		List<String> lines = new ArrayList<String>(rowNames.size() + 1);
		lines.add("NAMES\t" + ListUtils.toString(colNames, "\t"));
		for(int i=0 ; i<rowNames.size() ; i++) {
			lines.add(rowNames.get(i) + "\t" + ArrayUtils.toString(matrix.getRow(i), "\t"));
		}
		IOUtils.write(inputFile, lines);


		String cmdStr = home + "/bin/clustering/cluster3 -f " + inputFile.getAbsolutePath() + " -g " + getDistance(distance) + " -s -x " + dimension + " -y " + dimension;
		System.out.println("clustering command SOM: " + cmdStr);		

		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();

		if ( outputFile.exists() && outputFile.getTotalSpace() > 0 ) {
			int index;
			double value;
			List<String> values;
			Map<String, List<String>> map = new HashMap<String, List<String>> ();
			System.out.println("colnames (" + colNames.size() + ") = " + ListUtils.toString(colNames, ","));
			for (int i=0 ; i<colNames.size() ; i++) {
				values = IOUtils.column(outputFile, i+1);
				values.remove(0);
				index = 0;
				value = Math.abs(Double.parseDouble(values.get(index)));
				for(int j=1 ; j<values.size(); j++) {
					if ( Math.abs(Double.parseDouble(values.get(j))) < value ) {
						value = Math.abs(Double.parseDouble(values.get(j)));
						index = j;
					}
				}
				if ( !map.containsKey(index) ) {
					map.put(Integer.toString(index), new ArrayList<String>());
				}
				map.get(Integer.toString(index)).add(colNames.get(i));
			}
			StringBuilder sb = new StringBuilder();
			sb.append("(");
			List<String> keys = MapUtils.getKeys(map);
			int i=0;
			for(i=0 ; i<keys.size()-1 ; i++) {
				sb.append("(").append(ListUtils.toString(map.get(keys.get(i)), ",")).append("),");
			}
			sb.append("(").append(ListUtils.toString(map.get(keys.get(i)), ",")).append("));");
			//			System.out.println("newick " + sb.toString());
			//nw = new NewickParser().parse(sb.toString());
			return sb.toString();
		} else {
			return null;
			//throw new Exception("Impossible to generate newick");
		}

		//return nw;
	}

	private int getDistance(String distance) {
		if ( "none".equalsIgnoreCase(distance) ) {
			return 1;
		} else if ( "uncentered".equalsIgnoreCase(distance) ) {
			return 1;
		} else if ( "pearson".equalsIgnoreCase(distance) ) {
			return 2;
		} else if ( "spearman".equalsIgnoreCase(distance) ) {
			return 5;
		} else if ( "kendall".equalsIgnoreCase(distance) ) {
			return 6;
		} else if ( "euclidean".equalsIgnoreCase(distance) ) {
			return 7;
		}
		return 7;
	}

	private int heuristica(int row){
		int N,x;        
		N = (int)Math.round(Math.sqrt((double)row/3));
		if(N<=2){                           
			x=2;                        
		}else{                              
			if(N <= 5){                 
				x = 3;              
			}else{                      
				if(N <= 12){        
					x=4;        
				}else{              
					if(N <= 20){
						x = (int)Math.round(Math.sqrt((double)row/15));
					}else{                                       
						x = 10;                              
					}                                            
				}                                                    
			}                                                            
		}                                                                    
		return x;       
	}                    

}

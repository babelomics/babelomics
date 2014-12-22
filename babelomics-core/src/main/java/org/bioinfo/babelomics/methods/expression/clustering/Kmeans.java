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
import org.bioinfo.math.data.DoubleMatrix;

public class Kmeans extends Cluster {
	private int kvalue;
	
	public Kmeans(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance, int kvalue, String outdir, String home) {
		super(matrix, rowNames, colNames, distance, outdir, home);
		this.kvalue = kvalue;
	}

	@Override
	public String run(boolean createClusterFiles) throws Exception {
		//MultipleTree nw = null;
		
		File tmpDir = null; 
		if ( outdir == null) {
			tmpDir = File.createTempFile("input", ".dir");
			tmpDir.delete();
			tmpDir.mkdir();
		} else {
			tmpDir = new File(outdir);
		}
		
		System.out.println("---> create cluster files = " + createClusterFiles + " at " + tmpDir.getAbsolutePath());
		
		File inputFile = new File(tmpDir.getAbsoluteFile() + "/in.txt");
		File outputFile = new File(tmpDir.getAbsoluteFile() + "/in_K_G" + kvalue + ".kgg");

//		System.out.println("(infile, outfile) = (" + inputFile.getAbsolutePath() + ", " + outputFile.getAbsolutePath() + ")");
		
		List<String> lines = new ArrayList<String>(rowNames.size() + 1);
		lines.add("NAMES\t" + ListUtils.toString(colNames, "\t"));
		for(int i=0 ; i<rowNames.size() ; i++) {
			lines.add(rowNames.get(i) + "\t" + ArrayUtils.toString(matrix.getRow(i), "\t"));
		}
		IOUtils.write(inputFile, lines);
		
		String cmdStr = home + "/bin/clustering/cluster3 -f " + inputFile.getAbsolutePath() + " -g " + getDistance(distance) + " -k " + kvalue;
		System.out.println("clustering command KMEANS: " + cmdStr);		
		
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();

		if ( outputFile.exists() && outputFile.getTotalSpace() > 0 ) {
			lines = IOUtils.readLines(outputFile);
			lines.remove(0);
			Map<String, List<String>> map = new HashMap<String, List<String>> ();
			String[] values = null;
			for(String line: lines) {
				values = line.split("\t");
				if ( !map.containsKey(values[1]) ) {
					map.put(values[1], new ArrayList<String>());
				}
				map.get(values[1]).add(values[0]);
			}
			StringBuilder sb = new StringBuilder();
			sb.append("(");
			List<String> keys = MapUtils.getKeys(map);
			int i=0, nbCluster=0;
			for(i=0 ; i<keys.size()-1 ; i++) {
				if ( createClusterFiles ) {
					nbCluster = Integer.parseInt(keys.get(i)) + 1;
					IOUtils.write(tmpDir + "/cluster_" + nbCluster + ".txt", map.get(keys.get(i)));
				}
				
				sb.append("(").append(ListUtils.toString(map.get(keys.get(i)), ",")).append("),");
			}
			sb.append("(").append(ListUtils.toString(map.get(keys.get(i)), ",")).append("));");
			if ( createClusterFiles ) {
				IOUtils.write(tmpDir + "/rownames.txt", rowNames);
			}
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

}

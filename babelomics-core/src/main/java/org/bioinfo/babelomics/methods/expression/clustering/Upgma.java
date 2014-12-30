package org.bioinfo.babelomics.methods.expression.clustering;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.tree.multiway.MultipleTree;
import org.bioinfo.math.data.DoubleMatrix;

public class Upgma extends Cluster {

	public Upgma(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance, String home) {
		super(matrix, rowNames, colNames, distance, home);
	}
	
	@Override
	public String run(boolean createClusterFiles) throws Exception {
		MultipleTree nw = null;

		File inputFile = File.createTempFile("input", null);
		File outputFile = File.createTempFile("output", null);
		
		List<String> lines = new ArrayList<String>(rowNames.size() + 1);
		lines.add("#NAMES\t" + ListUtils.toString(colNames, "\t"));
		for(int i=0 ; i<rowNames.size() ; i++) {
			lines.add(rowNames.get(i) + "\t" + ArrayUtils.toString(matrix.getRow(i), "\t"));
		}
		IOUtils.write(inputFile, lines);
	
		String cmdStr = home + "/bin/clustering/cluster " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath() + " UPGMA " + (distance.equalsIgnoreCase("pearson") ? "correlation" : distance ) ;
		System.out.println("clustering command UPGMS: " + cmdStr);
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();

		System.out.println(sp.getRunnableProcess().getError());
		System.out.println(sp.getRunnableProcess().getException());
		System.out.println(sp.getRunnableProcess().getOutput());
		
		if ( outputFile.exists() && outputFile.getTotalSpace() > 0 ) {
			//nw = new NewickParser().parse(IOUtils.toString(outputFile));
			return IOUtils.toString(outputFile);
		} else {
			return null;
			//throw new Exception("Impossible to generate newick");
		}
		
		//return nw;
	}	
}

package org.bioinfo.babelomics.methods.expression.clustering;

import java.util.List;

import org.bioinfo.math.data.DoubleMatrix;

public abstract class Cluster {
	
	protected DoubleMatrix matrix;
	protected List<String> rowNames;
	protected List<String> colNames;
	protected String distance;
	protected String outdir;
	protected String home;

	public Cluster(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance, String home) {
		this(matrix, rowNames, colNames, distance, null, home);
	}

	public Cluster(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance, String outdir, String home) {
		this.matrix = matrix;
		this.rowNames = rowNames;
		this.colNames = colNames;
		this.distance = distance;
		this.outdir = outdir;
		this.home = home;
	}

	//public abstract MultipleTree run(boolean createClusterFiles) throws Exception;
	public abstract String run(boolean createClusterFiles) throws Exception;

}

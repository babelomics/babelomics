package org.bioinfo.babelomics.tools.expression;


import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class ClusteringTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void notest() {		
	}

	@Test
	public void Test() {
		//String dataset = "/mnt/commons/test/biodata/example/cyano.txt";
		String dataset = "/mnt/commons/test/biodata/example/preprocessed.txt";
		//String dataset = "/mnt/commons/test/biodata/example/merged_normalized_data.txt";
		String outdir = "/tmp/ClusteringTest";
		new File(outdir).mkdir();

		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "upgma", "--sample-clustering", "true", "--gene-clustering", "true", "--distance", "pearson", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----- clustering : UPGMA ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test1() {
		String dataset = "/mnt/commons/test/biodata/example/cyano.txt";
		String outdir = "/tmp/ClusteringTest1";
		new File(outdir).mkdir();

		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "sota", "--distance", "euclidean", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----- clustering : SOTA ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test2() {
		String dataset = "/mnt/commons/test/biodata/example/cyano.txt";
		String outdir = "/tmp/ClusteringTest2";
		new File(outdir).mkdir();

		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "som", "--distance", "euclidean", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----- clustering : SOM ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test3() {
		String dataset = "/mnt/commons/test/biodata/example/cyano.txt";
		String outdir = "/tmp/ClusteringTest3";
		new File(outdir).mkdir();

		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "kmeans", "--distance", "euclidean", "--kvalue", "10", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----- clustering : KMEANS ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
	

	public void Test4() {
		String dataset = "//mnt/commons/babelomics/tests/clustering/hypoxia_norm.txt";
		//String dataset = "/mnt/commons/test/biodata/example/merged_normalized_data.txt";
		String outdir = "/tmp/ClusteringTest";
		new File(outdir).mkdir();

		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "upgma", "--sample-clustering", "true", "--distance", "pearson", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----- clustering : UPGMA ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
	
}

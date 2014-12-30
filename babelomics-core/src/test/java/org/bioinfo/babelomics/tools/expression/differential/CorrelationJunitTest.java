package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class CorrelationJunitTest {

	@Test
	public void Test() {   
		String dataset = "/mnt/commons/test/biodata/example/correlation.txt";
		String outdir = "/tmp/CorrelationTest";
		new File(outdir).mkdir();
		
		System.out.println("-----   correlation : pearson ------");
		String []args = {"--tool", "correlation", "--dataset", dataset, "-o", outdir, "--test", "pearson", "--class-name", "indep", "--correction", "fdr", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	@Test
	public void Test1() {
		String dataset = "/mnt/commons/test/biodata/example/correlation.txt";
		String outdir = "/tmp/CorrelationTest1";
		new File(outdir).mkdir();

		System.out.println("-----  correlation : spearman ------");
		String []args = {"--tool", "correlation", "--dataset", dataset, "-o", outdir, "--test", "spearman", "--class-name", "indep", "--correction", "holm", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	@Test
	public void Test2() {
		String dataset = "/mnt/commons/test/biodata/example/correlation.txt";
		String outdir = "/tmp/CorrelationTest2";
		new File(outdir).mkdir();

		System.out.println("-----  correlation : regression  ------");
		String []args = {"--tool", "correlation", "--dataset", dataset, "-o", outdir, "--test", "regression", "--class-name", "indep", "--correction", "bonferroni", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}

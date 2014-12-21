package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class SurvivalTest {
	
	@Test
	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/survival.txt";
		String outdir = "/tmp/SurvivalTest";
		new File(outdir).mkdir();
		
		System.out.println("----- survival : cox ------");
		String []args = {"--tool", "survival", "--dataset", dataset, "-o", outdir, "--test", "cox", "--time-class", "time", "--censored-class", "censored", "--correction", "fdr", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}

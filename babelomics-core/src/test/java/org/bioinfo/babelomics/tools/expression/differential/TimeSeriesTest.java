package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class TimeSeriesTest {
	
	@Test
	public void Test() {
		
		String dataset = "/mnt/commons/test/biodata/example/masigpro.dataset";
		String outdir = "/tmp/TimeSeriesTest";
		new File(outdir).mkdir();
		
		System.out.println("----- time series test : masigpro ------");
		String []args = {"--tool", "time-dosage-series", "--dataset", dataset, "-o", outdir, "--test", "masigpro", "--contin-class", "contin", "--test", "masigpro", "--series-class", "series", "--degree", "2", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}

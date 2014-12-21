package org.bioinfo.babelomics.tools.genomic.copynumber;


import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class AgilentCGH2CNormalizationTest {

	@Test
	public void test0() {
	}


	public void test() {
		String outDirName = "/tmp/AgilentCGH2CNormalizationTest";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/cgh/agilent/normalization/dataset1/agilent_cgh_2c_dataset1.zip";
		
		System.out.println("----- agilent cgh normalization from ZIP file ------");
		String []args = { "--tool", "agilent-cgh-two-colors-normalization","--log-level", "2", "--compressed-file", dataset, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "quantile", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void test1() {
		String outDirName = "/tmp/AgilentCGH2CNormalizationTest1";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/cgh/agilent/normalization/dataset2/agilent_cgh_2c_dataset2.zip";
		
		System.out.println("----- agilent cgh normalization from ZIP file ------");
		String []args = { "--tool", "agilent-cgh-two-colors-normalization","--log-level", "2", "--compressed-file", dataset, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "quantile", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}

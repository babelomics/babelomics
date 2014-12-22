package org.bioinfo.babelomics.tools.genomic.copynumber;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class CopyNumberAnalysisTest {
	
	@Test
	public void test0() {
	}

	public void test() {
		String outDirName = "/tmp/CopyNumberAnalysisTest";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/cgh/agilent/cgh.normalized.txt";
		
		// the feature data file name is = dataset + ".featdata"
		
		System.out.println("----- copynumberanalysis dnacopy ------");
		String []args = {"--tool", "copy-number", "--normalized-file", dataset, "-o", outDirName, "--segmentation-method", "dnacopy", "--cgh-mcr", "--gap-allowed", "400", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}


	public void test1() {	
		String outDirName = "/tmp/CopyNumberAnalysisTest1";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/cgh/agilent/cgh.normalized.txt";
		
		// the feature data file name is = dataset + ".featdata"
		
		System.out.println("----- copynumberanalysis glad ------");
		String []args = {"--tool", "copy-number", "--normalized-file", dataset, "-o", outDirName, "--segmentation-method", "glad", "--cgh-mcr", "--gap-allowed", "400", "--home", System.getenv("BABELOMICS_HOME")};
		//String []args = {"-tool", "copy-number", "-normalized-file", "/mnt/commons/test/biodata/example/cgh/agilent/segmentation/data1000.txt", "-o", "/tmp/copynumber-glad", "-segmentation-method", "glad", "-cgh-mcr", "-gap-allowed", "400"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}

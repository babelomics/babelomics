package org.bioinfo.babelomics.tools.functional.tissues;


import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class AffyTmtTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	public void notest() throws Exception {
	}

	@Test
	public void Test1() throws Exception {
//		String list1 = "/mnt/commons/test/tools/tmt/list1_brain.txt"; // /ensembl1.txt
//		String list2 = "/mnt/commons/test/tools/tmt/list2_liver.txt"; // /ensembl2.txt
		String list1 = "/mnt/commons/test/tools/tmt/ensembl1.txt";
		String list2 = "/mnt/commons/test/tools/tmt/ensembl2.txt";
		String outdir = "/tmp/AffyTmtTest";
		new File(outdir).mkdir();

		System.out.println("----- affy tmt test  ------");
//		String []args = {"--tool", "tmt-affy", "--list1", list1,  "-o", outdir, "--species", "human", "--tissues", "all tissues"};
//		String []args = {"--tool", "tmt-affy", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "alltissues", "--home", System.getenv("BABELOMICS_HOME")};
		String []args = {"--tool", "tmt-affy", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "721 B lymphoblasts,Adipocyte,Adrenal Cortex", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

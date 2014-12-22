package org.bioinfo.babelomics.tools.functional.textmining;


import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class MarmiteScanTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void test0() {
		
	}
	
	@Test
	public void Test() {
		String list = "/mnt/commons/test/tools/marmitescan/marmitescan100.txt";
		//String list = "/mnt/commons/test/tools/marmitescan/marmitescan.txt";
		String outdir = "/tmp/MarmiteScanTest";
		new File(outdir).mkdir();

		System.out.println("----- marmite scan test  ------");
//		String []args = {"-tool", "marmite", "-list1", "/mnt/commons/test/tools/marmite/marmite1.txt",  "-list2", "/mnt/commons/test/tools/marmite/marmite2.txt", "-o", "/mnt/commons/test/tools/marmite/out/", "-bioentity-name", "diseases", "-bioentity-score-filter", "5"};
		String []args = {"--tool", "marmitescan", "--list", list,  "-o", outdir, "-bioentity-name", "disease", "-bioentity-score-filter", "5", "-bioentity-number-filter", "50", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}

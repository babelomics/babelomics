package org.bioinfo.babelomics.tools.functional.textmining;


import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class MarmiteTest {

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
		String list1 = "/mnt/commons/test/tools/marmite/marmite1a.txt"; ///marmite1a.txt /marmite1_40.txt
		String list2 = "/mnt/commons/test/tools/marmite/marmite2a.txt"; ///marmite2a.txt /marmite2_40.txt
		String outdir = "/tmp/MarmiteTest";
		new File(outdir).mkdir();

		System.out.println("----- marmite test  ------");
//		String []args = {"-tool", "marmite", "-list1", "/mnt/commons/test/tools/marmite/marmite1.txt",  "-list2", "/mnt/commons/test/tools/marmite/marmite2.txt", "-o", "/mnt/commons/test/tools/marmite/out/", "-bioentity-name", "diseases", "-bioentity-score-filter", "5"};
		String []args = {"--tool", "marmite", "--list1", list1,  "--list2", list2, "-o", outdir, "--bioentity-name", "disease", "--bioentity-score-filter", "5", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}

package org.bioinfo.babelomics.tools.functional.tissues;


import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class TMTTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void test0() throws Exception {
	}
		
	public void test1() throws Exception {
	System.out.println("-----     ------");
//	String []args = {"-tool", "tmt", "-list1", "/mnt/commons/test/tools/tmt/list1_brain.txt",  "-o", "/mnt/commons/test/tools/tmt/out/", "-organism", "human", "-platform", "affy", "-tissues", "all tissues"};
	String []args = {"-tool", "tmt", "-list1", "/mnt/commons/test/tools/tmt/list1_brain.txt",  "-list2", "/mnt/commons/test/tools/tmt/list2_liver.txt", "-o", "/mnt/commons/test/tools/tmt/out/", "-organism", "human", "-platform", "affy", "-tissues", "all tissues", "--home", System.getenv("BABELOMICS_HOME")};
//	String []args = {"-tool", "tmt", "-list1", "/mnt/commons/test/tools/tmt/ensembl1.txt",  "-list2", "/mnt/commons/test/tools/tmt/ensembl2.txt", "-o", "/mnt/commons/test/tools/tmt/out/", "-organism", "human", "-platform", "affy", "-tissues", "721 B lymphoblasts,Adipocyte,Adrenal Cortex"};
	
	try {
		BabelomicsMain.main(args); 
	} catch (Exception e) {
		e.printStackTrace();
	}		
}


}

package org.bioinfo.babelomics.tools.functional;


import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class FatiScanTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void Test0() {
		
	}

	public void Test1() {
		////
		String []args = {"-list", "/mnt/commons/test/biodata/example/fatiscan_input.txt", "-o", "/tmp", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--nb-partitions", "4", "--home", System.getenv("BABELOMICS_HOME")};

//		try {
//			FatiScanOld fatiscan = new FatiScanOld();
//			fatiscan.execute();
//		} catch (Exception e) {
//			e.printStackTrace();
//			//System.out.println(e.toString());
//		}
	}	

}

package org.bioinfo.babelomics.tools.preprocessing;


import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class CreateAnnotationTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	public void Test1() {
		String outdir = "/tmp/CreateAnnotationTest1";
		new File(outdir).mkdir();

		String []args = { "--tool", "create-annotation","--log-level", "2", "--species", "hsa", "--all-genome", "true", "--go-cc", "--go-bp", "--go-mf", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("CreateAnnotationTest1, args : " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test2() {
		String outdir = "/tmp/CreateAnnotationTest2";
		new File(outdir).mkdir();

		String []args = { "--tool", "create-annotation","--log-level", "2", "--species", "hsa", "--list", "AATK,BIRC7,PARM1", "--interpro", "interpro", "--kegg", "kegg", "--jaspar", "jaspar", "--output-format", "extended", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("CreateAnnotationTest2, args : " + Arrays.toString(args));

		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test3() {
		String outdir = "/tmp/CreateAnnotationTest3";
		new File(outdir).mkdir();

		String []args = { "--tool", "create-annotation","--log-level", "2", "--species", "hsa", "--list", "AATK,BIRC7,PARM1,kkkkkk", "--go-bp", "--kegg", "--go-cc", "--output-format", "compact", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("CreateAnnotationTest3, args : " + Arrays.toString(args));

		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	@Test
	public void Test4() {
		String outdir = "/tmp/CreateAnnotationTest4";
		new File(outdir).mkdir();

		String []args = { "--tool", "create-annotation","--log-level", "2", "--species", "sce", "--list", "IMP2',AATK,BIRC7,AAP1'-YHR048W", "--kegg", "--go-cc", "--output-format", "compact", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("CreateAnnotationTest3, args : " + Arrays.toString(args));

		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
}

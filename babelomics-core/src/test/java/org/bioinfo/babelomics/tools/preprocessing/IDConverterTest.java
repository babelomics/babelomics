package org.bioinfo.babelomics.tools.preprocessing;


import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class IDConverterTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	public void Test() {
		String outdir = "/tmp/IDConverterTest";
		new File(outdir).mkdir();

		String []args = { "--tool", "id-converter","--log-level", "2", "--species", "hsa", "--list", "AATK,BIRC7,PARM1,kkkkkk", "--db-names", "go,entrezgene,interpro", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("ID Converter Test, args : " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test1() {
		//String dataset = "/mnt/commons/test/biodata/example/genes.txt";
		String dataset = "/mnt/commons/test/tools/tmt/list2_liver.txt";
		String outdir = "/tmp/IDConverterTest1";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "id-converter","--log-level", "2", "--species", "hsa", "--listfile", dataset, "--go", "--kegg", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};
		//String []args = { "--tool", "id-converter","--log-level", "2", "--listfile", dataset, "--go", "--entrezgene", "--interpro", "--havana_gene", "--ensembl_gene", "--codelink", "--kegg", "-o", outdir};

		System.out.println("ID Converter Test1, args : " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test2() {
		//String dataset = "/mnt/commons/test/tools/tmt/list2_liver.txt";
		String dataset = "/mnt/commons/babelomics/tests/idconverter/genes_sce_notilde.txt";
		String outdir = "/tmp/IDConverterTest2";
		new File(outdir).mkdir();

		String []args = { "--tool", "id-converter","--log-level", "2", "--species", "sce", "--id-file", dataset, "--db-names", "go,entrezgene,interpro", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("ID Converter Test, args : " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
	
	@Test
	public void Test3() {
		//String dataset = "/mnt/commons/test/tools/tmt/list2_liver.txt";
		String dataset = "/mnt/commons/babelomics/tests/preprocessing/ids_sce_con_comillas.txt";
		String outdir = "/tmp/IDConverterTest3";
		new File(outdir).mkdir();

		String []args = { "--tool", "id-converter","--log-level", "2", "--species", "sce", "--id-file", dataset, "--db-names", "ensembl_gene", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("ID Converter Test, args : " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
	
	
}

package org.bioinfo.babelomics.tools.functional;


import static org.junit.Assert.fail;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class GenecodisTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void Test0() {
		String outdir = "/tmp/genecodis";
		new File(outdir).mkdir();
		
//		options.addOption(OptionFactory.createOption("datalist", "the ranked list"));
//		options.addOption(OptionFactory.createOption("support", "Minimum number of genes",true));
//		options.addOption(OptionFactory.createOption("support-for-random", "Minimum number of genes for correcting p-values",true));
//		options.addOption(OptionFactory.createOption("analysis", "singular_analysis or concurrence analysis",true));
//		options.addOption(OptionFactory.createOption("hypergeometric", "",false));
//		options.addOption(OptionFactory.createOption("chi-square", "",false));
//		options.addOption(OptionFactory.createOption("correction", "default correction",true));
//		
		
		String list1 = "/mnt/commons/test/example.motor";
		String list2 = "/mnt/commons/test/example.apoptosis";
//		String []args = {"--tool", "fatigo" ,"--list1", list1, "--list2", list2, "-o", "/tmp/fatigo"};
		
		
		String []args = {"--tool", "genecodis" ,"--datalist", list1, "--list2", list2, "--kegg","kegg","--go-bp", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "5", "--species", "hsa","--support", "3","--support-for-random", "3","--chi-square","true","--hypergeometric","true","--correction","fdr","--duplicates","each", "-o", outdir, "--home", System.getenv("BABELOMICS_HOME")};
		
		//String []args = {"--tool", "genecodis" ,"--datalist", list1,"--genome","true","--kegg","kegg","--go-bp", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "5", "--species", "hsa","--support", "3","--support-for-random", "3","--analysis", "concurrence","--hypergeometric","true","--correction","fdr", "-o", outdir};
		
		try {
			BabelomicsMain.main(args); 
//			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
//			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/histogram.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}

//	public void Test1() {
//		////
//		String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-rest-of-genome", "1", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref"};
//		//String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref"};
//		//String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "5"};
////		try {
////			FatiGOTool fatigo = new FatiGOTool();
////			fatigo.execute();			
////		} catch (Exception e) {
////			e.printStackTrace();
////			//System.out.println(e.toString());
////		}
//	}

}

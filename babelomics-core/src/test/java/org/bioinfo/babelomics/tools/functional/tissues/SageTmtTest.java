package org.bioinfo.babelomics.tools.functional.tissues;


import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class SageTmtTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void Test1() throws Exception {
//		String list1 = "/mnt/commons/test/tools/tmt/ensembl1.txt";
//		String list2 = "/mnt/commons/test/tools/tmt/ensembl2.txt";
		
		String list1 = "/mnt/commons/test/tools/tmt/list1_brain.txt";
		String list2 = "/mnt/commons/test/tools/tmt/list2_liver.txt";
		
		String outdir = "/tmp/SageTmtTest";
		new File(outdir).mkdir();

		System.out.println("----- sage tmt test  ------");
//		String []args = {"--tool", "tmt-sage", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "bone,brain,colon", "--histologies", "normal,cancer", "--perc-null-libraries", "0"};
//		String []args = {"--tool", "tmt-sage", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "all tissues", "--histologies", "normal,cancer", "--perc-null-libraries", "0"};
		//String []args = {"--tool", "tmt-sage", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "brain", "--histologies", "normal,cancer", "--tag-type", "long", "--min-tags", "5000", "--perc-null-libraries", "80", "perc-null-genes", "80"};
		//String []args = {"--tool", "tmt-sage", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "brain", "--histologies", "cancer", "--tag-type", "long", "--min-tags", "5000", "--perc-null-libraries", "80", "perc-null-genes", "80"};
		//String []args = {"--tool", "tmt-sage", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "alltissues", "--histologies", "cancer", "--tag-type", "long", "--min-tags", "5000", "--perc-null-libraries", "80", "perc-null-genes", "80"};
		String []args = {"--tool", "tmt-sage", "--list1", list1,  "--list2", list2, "-o", outdir, "--species", "human", "--tissues", "brain,liver", "--histologies", "cancer", "--tag-type", "long", "--min-tags", "5000", "--perc-null-libraries", "80", "--perc-null-genes", "90", "--home", System.getenv("BABELOMICS_HOME")};
		

		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}

package org.bioinfo.babelomics.tools.expression;


import static org.junit.Assert.fail;

import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PredictorTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void notest() {		
	}

	public void Test1() {
		String dataset = "/mnt/commons/test/biodata/example/predictor_input.txt";
		String outdir = "/tmp/";
		
		String []args = { "--tool", "predictor","--log-level", "2", "--dataset", dataset, "-o", outdir, "--knn", "--knn-tune", "--svm", "--svm-tune", "--random-forest", "--random-forest-tune", "--class", "class", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
//			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
//			System.out.println("cluster of genes:\n" + IOUtils.toString(new File(outdir + "/genes.nw")));
//			System.out.println("cluster of samples:\n" + IOUtils.toString(new File(outdir + "/samples.nw")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

}

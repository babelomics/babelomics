package org.bioinfo.babelomics.tools.preprocessing;


import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.IOUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PreprocessingTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	public void notest(){
		System.out.println("----------------> ");
	}
	

	public void Test0() {
		String dataset = "/mnt/commons/babelomics/tests/preprocessing/maria.jaime1000.txt";
		String outdir = "/tmp/PreprocessingTest";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--merge-replicates", "mean", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			//System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			//System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String outdir = "/tmp/PreprocessingTest";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--merge-replicates", "mean", "--logarithm-base", "2", "--impute-missing", "zero", "--filter-missing", "90", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	public void Test1() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String outdir = "/tmp/PreprocessingTest1";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir, "--impute-missing", "zero", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args);
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	


	public void Test2() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String outdir = "/tmp/PreprocessingTest2";
		new File(outdir).mkdir();
		String filename = "/mnt/commons/test/biodata/example/known_genes.txt";
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir, "--gene-file-filter", filename, "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			BabelomicsMain.main(args);
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	public void Test3() {
		//String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String dataset = "/mnt/commons/test/biodata/example/datamatrix_with_nan.txt";
		String outdir = "/tmp/PreprocessingTest3";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir, "--impute-missing", "knn", "--kvalue", "3", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}
	
	@Test
	public void Test4() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		//String dataset = "/mnt/commons/babelomics/tests/preprocessing/paco_preprocessing.txt";
		String outdir = "/tmp/PreprocessingTest";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--filter-missing", "15", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	public void Test5() {
		//String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String dataset = "/mnt/commons/babelomics/tests/preprocessing/paco_preprocessing_30.txt";
		String outdir = "/tmp/PreprocessingTest5";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--species", "hsa", "--convert-ids", "ensembl_gene", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
//			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
//			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}

	public void Test6() {
		String dataset = "/mnt/commons/babelomics/tests/preprocessing/normdata.txt";
		String outdir = "/tmp/PreprocessingTest6";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--merge-replicates", "mean", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	
	
}	
	
	
	
	
	
	
	
	
	

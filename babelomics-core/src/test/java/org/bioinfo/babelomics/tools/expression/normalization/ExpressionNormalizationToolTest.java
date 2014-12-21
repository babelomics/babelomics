package org.bioinfo.babelomics.tools.expression.normalization;


import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class ExpressionNormalizationToolTest {

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
		String outDirName = "/tmp/GenePixExpression1CNormalizationTest";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/onecolor_genepix_expression.zip";
		String tags = "genepix,one-channel";
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "normexp", "--ba-normalization", "quantiles", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}		
	
	public void Test15() {
		String outDirName = "/tmp/GenepixExpression1CNormalizationOneSingleFile";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/babelomics/tests/normalization/genepix/one-color/array1.gpr";
		String tags = "genepix,one-channel";
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "none", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test2() {
		String outDirName = "/tmp/GenePixExpression2CNormalizationTest";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/twocolor_genepix_expression.zip";
		String tags = "genepix,two-channels";
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "normexp", "--ba-normalization", "scale", "--wa-normalization", "loess", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}		

	public void Test25() {
		String outDirName = "/tmp/GenepixExpression2CNormalizationOneSingleFile";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/babelomics/tests/normalization/genepix/two-colors/array1.gpr";
		String tags = "genepix,two-channels";
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "none", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test3() {
		String outDirName = "/tmp/AgilentExpression1CNormalizationTest";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/agilent_one_color_expression.zip";
		String tags = "agilent,one-channel";
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "none", "--ba-normalization", "quantiles", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	
	
	public void Test35() {
		String outDirName = "/tmp/AgilentExpression1CNormalizationOneSingleFile";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/babelomics/tests/normalization/agilent/one-color/array_1.txt";
		String tags = "agilent,one-channel";
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "none", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	
	
	
	public void Test4() {
		String outDirName = "/tmp/AgilentExpression2CNormalizationTest2";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/GSE11968_RAW/agilent.zip";
		String tags = "agilent,two-channels";
		
		String []args = { "--tool", "agilent-expression-two-colors-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "quantiles", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	
	
	public void Test45() {
		String outDirName = "/tmp/AgilentExpression2CNormalizationOneSingleFile";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/babelomics/tests/normalization/agilent/two-colors/GSM302995.txt";
		String tags = "agilent,two-channels";
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", tags, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "none", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test5() {
		String dataset = "/mnt/commons/test/biodata/example/cel.zip";
		String outdir = "/tmp/AffyExpressionNormalizationTest";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", "affymetrix,one-channel", "-o", outdir, "--rma", "--plier", "--calls", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("executing ----------------> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	@Test
	public void Test55() {
		String dataset = "/mnt/commons/babelomics/tests/normalization/affy/Heart.CEL";
		String outdir = "/tmp/AffyExpressionNormalizationTestOneSingleFile";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", "affymetrix,one-channel", "-o", outdir, "--rma", "--plier", "--calls", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("executing ----------------> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test6() {
		String dataset = "/mnt/commons/test/biodata/example/CEL.tar.gz";
		String outdir = "/tmp/AffyExpressionNormalizationTest1";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", "affymetrix,one-channel", "-o", outdir, "--rma", "--plier", "--calls", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command ling parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test7() {
		String dataset = "/mnt/commons/test/biodata/example/cels_in_subdir.tar.gz";
		String outdir = "/tmp/AffyExpressionNormalizationTest7";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", "affymetrix,one-channel", "-o", outdir, "--rma", "--plier", "--calls", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command ling parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test8() {
		String dataset = "/mnt/commons/test/biodata/example/cels_in_subdir.zip";
		String outdir = "/tmp/AffyExpressionNormalizationTest8";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", "affymetrix,one-channel", "-o", outdir, "--rma", "--plier", "--calls", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command ling parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test9() {
		String dataset = "/mnt/commons/test/affy/GSE10245_RAW.tar";
		String outdir = "/tmp/AffyExpressionNormalizationTest9";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "expression-normalization","--log-level", "2", "--compressed-file", dataset, "--compressed-file-tags", "affymetrix,one-channel", "-o", outdir, "--rma", "--plier", "--calls", "--home", System.getenv("BABELOMICS_HOME")};

		System.out.println("command ling parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	
	
}

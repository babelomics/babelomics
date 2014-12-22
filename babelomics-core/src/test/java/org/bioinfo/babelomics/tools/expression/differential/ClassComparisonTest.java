package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.util.List;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.junit.Test;


public class ClassComparisonTest {
	
	public void notest() {
	}

	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/twoclasses.txt";
		String outdir = "/tmp/ClassComparisonTest";
		new File(outdir).mkdir();

		System.out.println("----- one class - limma ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "limma", "--class-name", "class", "--class-values", "luminal", "--correction", "fdr", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	@Test
	public void Test1() {
	    
		String dataset = "/mnt/commons/test/biodata/example/twoclasses.txt";
		String outdir = "/tmp/ClassComparisonTest1";
		new File(outdir).mkdir();

		System.out.println("----- two classes - limma ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "limma", "--class-name", "class", "--class-values", "basal,luminal", "--correction", "bonferroni", "--home", System.getenv("BABELOMICS_HOME")};
		System.out.println(ArrayUtils.toString(args, " "));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void Test2() {
	    
		String dataset = "/mnt/commons/test/biodata/example/twoclasses.txt";
		String outdir = "/tmp/ClassComparisonTest2";
		new File(outdir).mkdir();

		System.out.println("----- two classes - ttest ------");
		String test;
		//test = "t";
		test = "limma";
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", test, "--class-name", "class", "--class-values", "basal,luminal,", "--p-value", "0.04", "--correction", "hochberg", "--home", System.getenv("BABELOMICS_HOME")};
		
		System.out.println(ArrayUtils.toString(args, " "));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void Test3() {
	    
		String dataset = "/mnt/commons/test/biodata/example/twoclasses100.txt";
		String outdir = "/tmp/ClassComparisonTest3";
		new File(outdir).mkdir();

		System.out.println("----- two classes - fold change ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "fold-change", "--fold-change-value", "0.05", "--class-name", "class", "--class-values", "basal,luminal", "--correction", "holm", "--home", System.getenv("BABELOMICS_HOME")};
		System.out.println(ArrayUtils.toString(args, " "));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
		
	public void Test4() {
	    
		String dataset = "/mnt/commons/test/biodata/example/multiclasses.txt";
		String outdir = "/tmp/ClassComparisonTest4";
		new File(outdir).mkdir();

		System.out.println("----- multi classes - limma ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "limma", "--class-name", "indep", "--class-values", "1,3,5,7,9", "--correction", "by", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void Test5() {
	    
		String dataset = "/mnt/commons/test/biodata/example/multiclasses.txt";
		String outdir = "/tmp/ClassComparisonTest5";
		new File(outdir).mkdir();

		System.out.println("----- multi classes - anova ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "anova", "--class-name", "indep", "--class-values", "1,3,5,7,9", "--correction", "by", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	public void Test6() {
	    
		String dataset = "/mnt/commons/test/biodata/example/multiclasses.txt";
		String outdir = "/tmp/ClassComparisonTest6";
		new File(outdir).mkdir();

		System.out.println("----- two classes - ttest ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "t", "--class-name", "indep", "--class-values", "3,9", "--correction", "by", "--home", System.getenv("BABELOMICS_HOME")};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void Test7() {
		List<String> adjPvaluesList = StringUtils.toList("4,8,5,3,6");
		List<String> statList = StringUtils.toList("8,7,2,3,1");
		
		int[] sigOrder = ListUtils.order(adjPvaluesList);
		
		System.out.println("adj pvalues list = " + ListUtils.toString(adjPvaluesList));
		//System.out.println("stat list = " + ListUtils.toString(statList));
		System.out.println("order = " + ArrayUtils.toString(sigOrder, ","));
		System.out.println("adj pvalues ordered = " + ListUtils.toString(ListUtils.ordered(adjPvaluesList, sigOrder)));
		//System.out.println("stat ordered = " + ListUtils.toString(ListUtils.ordered(statList, sigOrder)));
		
	}
}

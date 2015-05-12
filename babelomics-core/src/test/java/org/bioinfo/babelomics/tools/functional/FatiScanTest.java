package org.bioinfo.babelomics.tools.functional;


import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsFactory;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

public class FatiScanTest {

    @Before
    public void setUp() throws Exception {
    }

    @After
    public void tearDown() throws Exception {
    }

    @Test
    public void Test0() {
        String outdir = "/tmp/fatiscan";
        new File(outdir).mkdir();

//        String list1 = "/home/ralonso/babelomics5/data/fatiscanmini.txt";
        String list1 = "/opt/babelomics/example/fatiscan_diabetesAffy.txt";

//		String []args = {"--tool", "fatigo" ,"--list1", list1, "--list2", list2, "-o", "/tmp/fatigo"};

        String cli = "--tool fatiscan --outdir " + outdir + " --species hsa --ranked-list " + list1 + " --method logistic --go-bp go-bp --go-bp-min-num-genes 5 --go-bp-max-num-genes 1000 --home " + System.getenv("BABELOMICS_HOME");
//        String[] args = {"--tool", "fatigo", "--list1", list1, "--list2", list2, "--go-bp", "--kegg", "-o", outdir, "--species", "hsa", "--home", System.getenv("BABELOMICS_HOME")};
        String[] args = cli.split(" ");
        try {
            FatiScanTool tool = (FatiScanTool) BabelomicsFactory.createTool("fatiscan");
            tool.parse(args);
            tool.run();
        } catch (ParseException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void Test1() {
        ////
        String[] args = {"-list", "/mnt/commons/test/biodata/example/fatiscan_input.txt", "-o", "/tmp", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--nb-partitions", "4", "--home", System.getenv("BABELOMICS_HOME")};

//		try {
//			FatiScanOld fatiscan = new FatiScanOld();
//			fatiscan.execute();
//		} catch (Exception e) {
//			e.printStackTrace();
//			//System.out.println(e.toString());
//		}
    }

}

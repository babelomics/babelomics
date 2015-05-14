package org.bioinfo.babelomics.tools.functional;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.InfraredUtils;
import org.bioinfo.babelomics.tools.BabelomicsFactory;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class FatiGOTest {

//    @Before
//    public void setUp() throws Exception {
//    }
//
//    @After
//    public void tearDown() throws Exception {
//    }

//    @Test
//    public void getKegg() {
//        try {
//            List<String> lines = IOUtils.readLines("/tmp/ensg_kegg.tsv");
//            for (String line : lines) {
//                String []fields = line.split("\t");
//                String ensg = fields[0];
//                String kegg = fields[1];
//                org.bioinfo.babelomics.utils.XrefManager xref = new org.bioinfo.babelomics.utils.XrefManager(ensg,"hsa");
//                Map<String, List<String>> ensenbl_transcript = xref.getXrefs("ensembl_transcript");
//                System.out.println("ensenbl_transcript = " + ensenbl_transcript);
//            }
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//    }

    @Test
    public void TestList1List2Hsa() {
        String outdir = "/tmp/fatigo";
        new File(outdir).mkdir();


        String list1 = "/home/ralonso/opt/babelomics/example/example.motor";
        String list2 = "/home/ralonso/opt/babelomics/example/example.apoptosis";

        String[] args = {"--tool", "fatigo", "--list1", list1, "--list2", list2, "--go-bp", "--go-bp-propagation","propagate","--go-bp-min-num-genes", "5", "--go-bp-max-num-genes", "1000", "-o", outdir, "--species", "hsa", "--home", System.getenv("BABELOMICS_HOME")};
        try {
            org.bioinfo.babelomics.tools.functional.FatiGOTool fatigo = (org.bioinfo.babelomics.tools.functional.FatiGOTool) BabelomicsFactory.createTool("fatigo");
            fatigo.parse(args);
            fatigo.run();
        } catch (ParseException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // @Test
    public void TestList1List2Mmu() {
        String outdir = "/tmp/fatigo";
        new File(outdir).mkdir();


//        String list1 = "/mnt/commons/test/example.motor";
        String list1 = "/home/ralonso/appl/babelomics-old/babelomics-old/example/mmu_list1.txt";
//		String list2 = "/mnt/commons/test/example.apoptosis";
        String list2 = "/home/ralonso/appl/babelomics-old/babelomics-old/example/mmu_list2.txt";

        String[] args = {"--tool", "fatigo", "--list1", list1, "--list2", list2, "--go-bp", "--go-mf", "--go-bp-min-num-genes", "5", "--go-bp-max-num-genes", "1000", "-o", outdir, "--species", "mmu", "--home", System.getenv("BABELOMICS_HOME")};
        try {
            FatiGOTool fatigo = (FatiGOTool) BabelomicsFactory.createTool("fatigo");
            fatigo.parse(args);
            fatigo.run();
        } catch (ParseException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //    @Test
    public void TestList1GenomeHsa() {
        String outdir = "/tmp/fatigo";
        new File(outdir).mkdir();


        String list1 = "/home/ralonso/appl/babelomics-old/babelomics-old/example/example.motor";

        String[] args = {"--tool", "fatigo", "--list1", list1, "--genome", "--go-bp", "--go-bp-min-num-genes", "5", "--go-bp-max-num-genes", "1000", "-o", outdir, "--species", "hsa", "--home", System.getenv("BABELOMICS_HOME")};
        try {
            FatiGOTool fatigo = (FatiGOTool) BabelomicsFactory.createTool("fatigo");
            fatigo.parse(args);
            fatigo.run();
        } catch (ParseException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    //	@Test
    public void Test0() {
        String outdir = "/tmp/fatigo";
        new File(outdir).mkdir();

//		String list1 = "/mnt/commons/test/example.motor";
        String list1 = "/home/ralonso/appl/babelomics-old/babelomics-old/example/example.motor";
//		String list2 = "/mnt/commons/test/example.apoptosis";
        String list2 = "/home/ralonso/appl/babelomics-old/babelomics-old/example/example.apoptosis";

//		String []args = {"--tool", "fatigo" ,"--list1", list1, "--list2", list2, "-o", "/tmp/fatigo"};
        String[] args = {"--tool", "fatigo", "--list1", list1, "--list2", list2, "--go-bp", "--kegg", "-o", outdir, "--species", "hsa", "--home", System.getenv("BABELOMICS_HOME")};
        try {
            FatiGOTool fatigo = (FatiGOTool) BabelomicsFactory.createTool("fatigo");
            fatigo.parse(args);
            fatigo.run();
        } catch (ParseException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void Test1() {
        ////
        String[] args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-rest-of-genome", "1", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref", "--home", System.getenv("BABELOMICS_HOME")};
        //String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref"};
        //String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "5"};
//		try {
//			FatiGOTool fatigo = new FatiGOTool();
//			fatigo.execute();			
//		} catch (Exception e) {
//			e.printStackTrace();
//			//System.out.println(e.toString());
//		}
    }

    //	@Test
    public void Test() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
        DBConnector dbConnector = new DBConnector("hsa", new File("/opt/babelomics/conf/infrared.properties"));
        System.err.println("hssss");
        GOFilter goc = new GOFilter("cellular_component");
        goc.setPropagated(false);
        System.out.print(InfraredUtils.getAnnotations(dbConnector, Arrays.asList("ACTA1"), goc).toString());
        GOFilter gob = new GOFilter("biological_process");
        gob.setPropagated(false);
        System.out.print(InfraredUtils.getAnnotations(dbConnector, Arrays.asList("ACTA1"), gob).toString());
        GOFilter gom = new GOFilter("molecular_function");
        gom.setPropagated(false);
        System.out.print(InfraredUtils.getAnnotations(dbConnector, Arrays.asList("ACTA1"), gom).toString());
    }

    //    //@Test
    public void miTest() throws FileNotFoundException {
//        DBConnector dbConnector = new DBConnector("hsa", new File("/opt/babelomics/conf/infrared.properties"));
//        List<String> list2 = InfraredUtils.getGenome(dbConnector);
//        System.out.println("list2.size() = " + list2.size());
//        String content = ListUtils.toString(list2);
//
//        PrintWriter out = new PrintWriter("/tmp/filename.txt");
//        out.println(content);
//        out.close();
//        FatiGOTool fatigoTool = (FatiGOTool)BabelomicsFactory.createTool("fatigo");
//        fatigoTool.parse(args);
//        fatigoTool.run();

    }
}

package org.bioinfo.babelomics.tools.interactome;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class SnowTest {
	String BABELOMICS_HOME = System.getenv("BABELOMICS_HOME");
	@Test
	public void testTest(){
	 System.out.println("I am a test!!");	
	}
	@Test
	public void testExample1(){
		//opt/babelomics/babelomics.sh --tool snow --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5974 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5974/job.log --list2 none --list1 /opt/babelomics/example/brca1_overexp_up.txt --side less --images  --randoms 500 --components true --intermediate 1 --interactome hsa --group curated --type genes --o-name result
		
		String outdir = "/tmp/snow/test1";
		//String outdir = "/tmp/snow/testIntermediate";
		new File(outdir).mkdirs();
		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--randoms", "500",
				"--o-name","result",
				"--interactome","hsa",
				"--group","curated",
				"--type", "genes",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/genes/brca1_overexp_up.txt",
//				"--list1",BABELOMICS_HOME+"/example/K562_symbol.txt",
				"--side", "less",
				"--intermediate","1",
				"--components","true",
				//"--xml",
				"--images",
//				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	
	//@Test
	public void ownMMUTest(){
		//opt/babelomics/babelomics.sh --tool snow --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5959 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5959/job.log --bicomponents true --list2 none --sif-file /httpd/bioinfo/wum_sessions_v0.7/4164/data/29908/mmu_alldb_proteins_interactome_nr.sif --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/29907/clean_prots_up_ratio_MMP14.txt --side less --images  --randoms 500 --components true --intermediate 1 --interactome own --type proteins --o-name result
		String outdir = "/tmp/snow/test1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir,
				"--sif-file", "/httpd/bioinfo/wum_sessions_v0.7/4164/data/29908/mmu_alldb_proteins_interactome_nr.sif", 
				"--randoms", "5",
				"--o-name","result",
				"--interactome","own",
				"--type", "proteins",
				"--list1","/httpd/bioinfo/wum_sessions_v0.7/4164/data/29907/clean_prots_up_ratio_MMP14.txt",
				"--side", "less",
				"--intermediate","1",
				"--components","1",
				//"--xml",
				"--images",
				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	//@Test
	public void btaTest(){
		//opt/babelomics/babelomics.sh --tool snow --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5944 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5944/job.log --list2 none --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/29892/bta_genes.txt --side less --images  --randoms 500 --components true --intermediate 0 --interactome bta --group curated --type genes --o-name result
		String outdir = "/tmp/snow/test1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--randoms", "100",
				"--o-name","result",
				"--interactome","bta",
				"--group","curated",
				"--type", "genes",
				"--list1","/httpd/bioinfo/wum_sessions_v0.7/4164/data/29892/bta_genes.txt",
				"--side", "less",
				"--intermediate","0",
				"--components","1",
				//"--xml",
				"--images",
				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	//@Test
	public void sceTest(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow/test1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--randoms", "100",
				"--o-name","result",
				"--interactome","sce",
				"--group","all",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/sce/genes/OET_R_DNreg.txt",
				"--side", "less",
				"--intermediate","0",
				"--components","1",
				//"--xml",
				"--images",
				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	//@Test
	public void SnowExampleOneList(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow/test1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--randoms", "100",
				"--o-name","result",
				"--interactome","hsa",
				"--group","curated",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/chr_10_block225.txt",
				"--side", "less",
				"--intermediate","0",
				"--components","1",
				//"--xml",
				"--images",
				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	
	//@Test
	public void genesTest(){

		String outdir = "/tmp/snow/test2";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--randoms", "50",
				"--o-name","result",
				"--interactome","hsa",
				"--group","all",
				"--type", "genes",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/genes/list1.txt",
				"--side", "less",
				"--intermediate","1",
				"--components","1",
				"--xml",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	//@Test
	public void transcriptsTest(){

		String outdir = "/tmp/snow2/test3";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--randoms", "50",
				"--o-name","result",
				"--interactome","hsa",
				"--group","all",
				"--type", "transcripts",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/transcripts/list1.txt",
				"--side", "less",
				"--intermediate","1",
//				"--components","1",
				"--xml",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}

//	@Test
	public void SnowExampleTwoLists(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow/testTwoLists";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--o-name","result",
				"--components", "1",
				"--interactome","hsa",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/chr_10_block225.txt",
				"--list2","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/chr_11_block80.txt",
				"--side", "less",
				"--intermediate", "0",
				"--xml",
				"--images",
				"--json",
				//"--sif",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}

//	@Test
	public void createTopoFile(){
		//./babelomics.sh --tool snow2 --outdir /tmp/ --sif-file /mnt/commons/babelomics/tests/snow2/ej8/ej8.sif  --o-sif-topo-file --o-name prueba --interactome own
		String outdir = "/mnt/commons/babelomics/tests/snow/ej8/";
		new File(outdir).mkdirs();
		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--sif-file", "/mnt/commons/babelomics/tests/snow2/ej8/ej9.sif",
				"--o-sif-topo-file",
				"--interactome","own",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	//@Test
	public void ownInteractomeTest(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow/ownTest";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow", 
				"-o", outdir, 
				"--randoms", "2",
				"--o-name","result",
				"--interactome","own",
				"--group","all",
				"--sif-file","/mnt/commons/babelomics/tests/snow2/ej8/ej8.sif",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej8/list1",
				"--side", "less",
				"--intermediate","1",
				"--components","1",
				"--xml",
				"--sif",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	public void main(String []args){
	try {
		//			for(String arg : args)
		//				System.out.println(arg);
		BabelomicsMain.main(args);

	} catch (Exception e) {
		e.printStackTrace();
	}
	}
}

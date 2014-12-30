package org.bioinfo.babelomics.tools.interactome;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class RandomsSnowToolTest {
	@Test
	public void test(){}
	//@Test
	public void test1(){
		String outdir = "/tmp/snow2/test1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "randoms-snow", 
				"-o", outdir, 
				"--random-size", "200",
				"--randoms", "10",
				"--o-name","result",
				"--interactome","hsa",
				"--group","all",
				"--type", "proteins",
				"--intermediate","1",
				"--components", "1",
				"--topo-values", "1",
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


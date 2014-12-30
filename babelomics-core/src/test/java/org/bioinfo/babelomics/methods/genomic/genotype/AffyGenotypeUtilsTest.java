package org.bioinfo.babelomics.methods.genomic.genotype;

import static org.junit.Assert.fail;

import java.io.IOException;

import org.bioinfo.microarray.AffyGenotypeUtils;
import org.junit.Test;

public class AffyGenotypeUtilsTest {

	@Test
	public void testAffyGenotypeNormalization() {
//		fail("Not yet implemented");
	}

	@Test
	public void testAffyToPedAndMap() {
//		fail("Not yet implemented");
	}

	public void testAffyMapFileCreator() {
		try {
//			AffyGenotypeUtils.affyGenomeWide6ParsedAnnotationFileCreator("/home/imedina/appl/analysis/src/snp/GenomeWideSNP_6.na29.annot.csv", "/tmp/genome_parsed_GW6_file.txt");
//			AffyGenotypeUtils.affyHumanMapping500kParsedAnnotationFileCreator("/home/imedina/appl/analysis/src/snp/Mapping250K_Nsp.na29.annot.csv", "/home/imedina/appl/analysis/src/snp/Mapping250K_Sty.na29.annot.csv", "/tmp/genome_parsed_250k_file.txt");
//			
//			AffyGenotypeUtils.affyGenomeWide6MapFileCreator("/home/imedina/appl/analysis/src/snp/GenomeWideSNP_6.na29.annot.csv", "/tmp/map_GW6_file.txt");
//			AffyGenotypeUtils.affyHumanMapping500kMapFileCreator("/home/imedina/appl/analysis/src/snp/Mapping250K_Nsp.na29.annot.csv", "/home/imedina/appl/analysis/src/snp/Mapping250K_Sty.na29.annot.csv", "/tmp/map_250k_file.txt");
			
			System.err.println("********************************");
			System.err.println("********************************");
			AffyGenotypeUtils.affyGenomeWide6ChpCallsToPedAndMap("/home/imedina/appl/analysis/src/snp/GenomeWideSNP_6.na29.annot.csv", "/tmp/map_GW6_file.txt", "", "");
			System.err.println("********************************");
			System.err.println("********************************");
		} catch (IOException e) {
			e.printStackTrace();
			fail("in testAffyMapFileCreator: " + e.toString());
		}
	}

}

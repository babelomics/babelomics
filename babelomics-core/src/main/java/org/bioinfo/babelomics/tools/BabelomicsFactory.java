package org.bioinfo.babelomics.tools;

import org.bioinfo.babelomics.tools.expression.BiclusteringTool;
import org.bioinfo.babelomics.tools.expression.Clustering;
import org.bioinfo.babelomics.tools.expression.MAPlot;
import org.bioinfo.babelomics.tools.expression.OutlierLrs;
import org.bioinfo.babelomics.tools.expression.Predictor;
import org.bioinfo.babelomics.tools.expression.RawExpressionViewer;
import org.bioinfo.babelomics.tools.expression.differential.ClassComparison;
import org.bioinfo.babelomics.tools.expression.differential.Correlation;
import org.bioinfo.babelomics.tools.expression.differential.Survival;
import org.bioinfo.babelomics.tools.expression.differential.TimeSeries;
import org.bioinfo.babelomics.tools.expression.normalization.ExpressionNormalizationTool;
import org.bioinfo.babelomics.tools.functional.Blast2GoTool;
import org.bioinfo.babelomics.tools.functional.FatiGOTool;
import org.bioinfo.babelomics.tools.functional.FatiScanTool;
import org.bioinfo.babelomics.tools.functional.GeneCodisTool;
import org.bioinfo.babelomics.tools.functional.textmining.Marmite;
import org.bioinfo.babelomics.tools.functional.textmining.MarmiteScan;
import org.bioinfo.babelomics.tools.functional.tissues.AffyTmt;
import org.bioinfo.babelomics.tools.functional.tissues.SageTmt;
import org.bioinfo.babelomics.tools.genomic.copynumber.AgilentCGH1CNormalization;
import org.bioinfo.babelomics.tools.genomic.copynumber.AgilentCGH2CNormalization;
import org.bioinfo.babelomics.tools.genomic.copynumber.CopyNumberAnalysis;
import org.bioinfo.babelomics.tools.genomic.genotype.AffyGenotypePreprocessing;
import org.bioinfo.babelomics.tools.genomic.genotype.AssociationTool;
import org.bioinfo.babelomics.tools.genomic.genotype.StratificationTool;
import org.bioinfo.babelomics.tools.graph.DescriptiveStatistics;
import org.bioinfo.babelomics.tools.graph.GoGraphViewerTool;
import org.bioinfo.babelomics.tools.interactome.RandomsSnowTool;
import org.bioinfo.babelomics.tools.interactome.Snow;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnow;
import org.bioinfo.babelomics.tools.preprocessing.CreateAnnotation;
import org.bioinfo.babelomics.tools.preprocessing.IDConverter;
import org.bioinfo.babelomics.tools.preprocessing.Preprocessing;


public class BabelomicsFactory {

	public static BabelomicsTool createTool(String toolName) {
		
		/*
		 * **************************************************************************************************
		 * *****	Preprocessing, normalization and annotation tools	*************************************
		 * **************************************************************************************************
		 */
		if(toolName.equalsIgnoreCase("expression-normalization")) {
			return new ExpressionNormalizationTool();
		}

		if(toolName.equalsIgnoreCase("affy-expression-normalization")) {
			return new ExpressionNormalizationTool();
		}
		
		if(toolName.equalsIgnoreCase("agilent-expression-one-color-normalization")) {
			return new ExpressionNormalizationTool();
		}
		
		if(toolName.equalsIgnoreCase("agilent-expression-two-colors-normalization")) {
			return new ExpressionNormalizationTool();
		}

		if(toolName.equalsIgnoreCase("genepix-expression-one-color-normalization")) {
			return new ExpressionNormalizationTool();
		}
		
		if(toolName.equalsIgnoreCase("genepix-expression-two-colors-normalization")) {
			return new ExpressionNormalizationTool();
		}

		if(toolName.equalsIgnoreCase("raw-expression-viewer")) {
			return new RawExpressionViewer();
		}
		
		if(toolName.equalsIgnoreCase("ma-plot")) {
			return new MAPlot();
		}

		if(toolName.equalsIgnoreCase("affy-genotype-preprocess")) {
			return new AffyGenotypePreprocessing();
		}
		
		if(toolName.equalsIgnoreCase("agilent-cgh-one-color-normalization")) {
			return new AgilentCGH1CNormalization();
		}

		if(toolName.equalsIgnoreCase("agilent-cgh-two-colors-normalization")) {
			return new AgilentCGH2CNormalization();
		}

		if(toolName.equalsIgnoreCase("copy-number-normalization")) {
			return new AgilentCGH2CNormalization();
		}

		if(toolName.equalsIgnoreCase("affy-snp-preprocess")) {
			return new AffyGenotypePreprocessing();
		}

		if(toolName.equalsIgnoreCase("preprocessing")) {
			return new Preprocessing();
		}
		
		if(toolName.equalsIgnoreCase("id-converter")) {
			return new IDConverter();
		}
		
		if(toolName.equalsIgnoreCase("create-annotation")) {
			return new CreateAnnotation();
		}
		
		
		/*
		 * **************************************************************************************************
		 * *****	Expression tools (differential, predictors and clustering)	*****************************
		 * **************************************************************************************************
		 */
		if(toolName.equalsIgnoreCase("class-comparison")) {
			return new ClassComparison();
		}
		
//		if(toolName.equalsIgnoreCase("differential-expression")) {
//			return new DifferentialAnalysis();
//		}

		if(toolName.equalsIgnoreCase("correlation")) {
			return new Correlation();
		}

		if(toolName.equalsIgnoreCase("survival")) {
			return new Survival();
		}

		if(toolName.equalsIgnoreCase("time-dosage-series")) {
			return new TimeSeries();
		}
		
		if(toolName.equalsIgnoreCase("class-prediction")) {
			return new Predictor();
		}

//		if(toolName.equalsIgnoreCase("regression")) {
//			return new DifferentialAnalysis();
//		}
		
		if(toolName.equalsIgnoreCase("clustering")) {
			return new Clustering();
		}

		if(toolName.equalsIgnoreCase("biclustering")) {
			return new BiclusteringTool();			
		}
		
		/*
		 * **************************************************************************************************
		 * *****	Genomic tools (Copy number, Association and Stratification)	*****************************
		 * **************************************************************************************************
		 */
		if(toolName.equalsIgnoreCase("copy-number")) {
			return new CopyNumberAnalysis();
		}
		
		if(toolName.equalsIgnoreCase("association")) {
			return new AssociationTool();
		}
		
		if(toolName.equalsIgnoreCase("stratification")) {
			return new StratificationTool();
		}
		
//		if(toolName.equalsIgnoreCase("linkage")) {
//			return new OutlierLrs();
//		}
		
		
		/*
		 * **************************************************************************************************
		 * *****	Functional profiling tools: *************************************************************
		 * **************************************************************************************************
		 */
		if(toolName.equalsIgnoreCase("fatigo")) {
			return new FatiGOTool();
		}

		if(toolName.equalsIgnoreCase("fatiscan")) {
			return new FatiScanTool();
		}

		if(toolName.equalsIgnoreCase("marmite")) {
			return new Marmite();
		}

		if(toolName.equalsIgnoreCase("marmitescan")) {
			return new MarmiteScan();
		}
		
		if(toolName.equalsIgnoreCase("tmt-affy")) {
			return new AffyTmt();
		}
		
		if(toolName.equalsIgnoreCase("tmt-sage")) {
			return new SageTmt();
		}

		if(toolName.equalsIgnoreCase("gesbap")) {
			return new OutlierLrs();
		}
		
		if(toolName.equalsIgnoreCase("snow")) {
			return new Snow();
		}
		if(toolName.equalsIgnoreCase("randoms-snow")) {
			return new RandomsSnowTool();
		}
		if(toolName.equalsIgnoreCase("network-miner")) {
			return new GSnow();
		}
		if(toolName.equalsIgnoreCase("blast2go")) {
			return new Blast2GoTool();
		}
		
		if(toolName.equalsIgnoreCase("genecodis")) {
			return new GeneCodisTool();
		}
		
		/*
		 * **************************************************************************************************
		 * *****	Graphs and displays tools	*************************************************************
		 * **************************************************************************************************
		 */
		if(toolName.equalsIgnoreCase("descriptive-statistics")) {
			return new DescriptiveStatistics();
		}
		
		if(toolName.equalsIgnoreCase("clustering-tree")) {
			return new DescriptiveStatistics();
		}
		
		if(toolName.equalsIgnoreCase("go-graph-viewer")) {
			return new GoGraphViewerTool();
		}
		
		if(toolName.equalsIgnoreCase("pca-plot")) {
			return new DescriptiveStatistics();
		}
		return null;
	}

}

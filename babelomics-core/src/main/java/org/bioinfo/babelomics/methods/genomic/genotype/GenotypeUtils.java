package org.bioinfo.babelomics.methods.genomic.genotype;

import java.io.IOException;
import java.util.List;

import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.math.util.MathUtils;

public class GenotypeUtils {

	
	public static DataFrame plinkAssocResultToDataFrame(String filename) throws InvalidIndexException, IOException {
		DataFrame dataFrame = new DataFrame();
		dataFrame.addColumn("dbsnp", IOUtils.column(filename, 1, "\\s+"));
		dataFrame.addColumn("chromosome", IOUtils.column(filename, 0, "\\s+"));
		dataFrame.addColumn("position", IOUtils.column(filename, 2, "\\s+"));
		List<String> pvalues = IOUtils.column(filename, 7, "\\s+");
		double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
		minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
		dataFrame.addColumn("p_values", pvalues);
		dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
		dataFrame.addColumn("odd_ratio", IOUtils.column(filename, 8, "\\s+"));
		dataFrame.removeRow(0);
		return dataFrame;
	}

	public static DataFrame plinkFisherResultToDataFrame(String filename) throws InvalidIndexException, IOException {
		DataFrame dataFrame = new DataFrame();
		dataFrame.addColumn("dbsnp", IOUtils.column(filename, 1, "\\s+"));
		dataFrame.addColumn("chromosome", IOUtils.column(filename, 0, "\\s+"));
		dataFrame.addColumn("position", IOUtils.column(filename, 2, "\\s+"));
		List<String> pvalues = IOUtils.column(filename, 7, "\\s+");
		double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
		minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
		dataFrame.addColumn("p_values", pvalues);
		dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
		dataFrame.addColumn("odd_ratio", IOUtils.column(filename, 8, "\\s+"));
		dataFrame.removeRow(0);
		return dataFrame;
	}

	public static DataFrame plinkLinearResultToDataFrame(String filename) throws InvalidIndexException, IOException {
		DataFrame dataFrame = new DataFrame();
		dataFrame.addColumn("dbsnp", IOUtils.column(filename, 1, "\\s+"));
		dataFrame.addColumn("chromosome", IOUtils.column(filename, 0, "\\s+"));
		dataFrame.addColumn("position", IOUtils.column(filename, 2, "\\s+"));
		List<String> pvalues = IOUtils.column(filename, 7, "\\s+");
		double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
		minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
		dataFrame.addColumn("p_values", pvalues);
		dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
		dataFrame.addColumn("odd_ratio", IOUtils.column(filename, 8, "\\s+"));
		dataFrame.removeRow(0);
		return dataFrame;
	}

	public static DataFrame plinkLogisticResultToDataFrame(String filename) throws InvalidIndexException, IOException {
		DataFrame dataFrame = new DataFrame();
		dataFrame.addColumn("dbsnp", IOUtils.column(filename, 1, "\\s+"));
		dataFrame.addColumn("chromosome", IOUtils.column(filename, 0, "\\s+"));
		dataFrame.addColumn("position", IOUtils.column(filename, 2, "\\s+"));
		List<String> pvalues = IOUtils.column(filename, 7, "\\s+");
		double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
		minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
		dataFrame.addColumn("p_values", pvalues);
		dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
		dataFrame.addColumn("odd_ratio", IOUtils.column(filename, 8, "\\s+"));
		dataFrame.removeRow(0);
		return dataFrame;
	}

	public static DataFrame plinkTdtResultToDataFrame(String filename) throws InvalidIndexException, IOException {
		DataFrame dataFrame = new DataFrame();
		dataFrame.addColumn("dbsnp", IOUtils.column(filename, 1, "\\s+"));
		dataFrame.addColumn("chromosome", IOUtils.column(filename, 0, "\\s+"));
		dataFrame.addColumn("position", IOUtils.column(filename, 2, "\\s+"));
		List<String> pvalues = IOUtils.column(filename, 7, "\\s+");
		double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
		minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
		dataFrame.addColumn("p_values", pvalues);
		dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
		dataFrame.addColumn("odd_ratio", IOUtils.column(filename, 8, "\\s+"));
		dataFrame.removeRow(0);
		return dataFrame;
	}

}

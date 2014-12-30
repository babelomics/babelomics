package org.bioinfo.babelomics.tools.functional.tissues;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.tool.OptionFactory;

public abstract class Tmt extends BabelomicsTool {

	@Override
	public void initOptions() {
//		options.addOption(OptionFactory.createOption("organism", "Organism, valid values are 'human' and 'mouse'"));
		options.addOption(OptionFactory.createOption("list1", "the feature data containig the list #1 of genes, or the feature data file"));
		options.addOption(OptionFactory.createOption("list2", "the feature data containig the list #2 of genes, or the feature data file", false));
		options.addOption(OptionFactory.createOption("tissues", "the list of tissues separated by commas. Enter 'all tissues' to take into account all available tissues"));
	}


	protected  TestResultList<TTestResult> runTtest(DoubleMatrix matrix1, DoubleMatrix matrix2) throws InvalidParameterException, MathException {		
		TTest tTest = new TTest();
		List<Double> sample1, sample2;
		TestResultList<TTestResult> res = new TestResultList<TTestResult>(matrix1.getRowDimension());
		for(int row=0 ; row < matrix1.getRowDimension() ; row++) {
				
			sample1 = new ArrayList<Double>();
			sample2 = new ArrayList<Double>();
		
			for(int col=0 ; col < matrix1.getColumnDimension() ; col++) {
				if ( ! Double.isNaN(matrix1.get(row, col)) ) {
					sample1.add(matrix1.get(row, col));
				}
			}
			for(int col=0 ; col < matrix2.getColumnDimension() ; col++) {
				if ( ! Double.isNaN(matrix2.get(row, col)) ) {
					sample2.add(matrix2.get(row, col));
				}
			}
		
			res.add(tTest.tTest(ListUtils.toDoubleArray(sample1), ListUtils.toDoubleArray(sample2)));
		}
		
		return res;
	}

	protected List<Integer> getColumns(DoubleMatrix matrix, double filter) {
		int nullValues;
		List<Integer> columnIndexes = new ArrayList<Integer>();
		double limit = matrix.getRowDimension() * filter / 100;
		//System.out.println("getColumns, limit = " + limit);
		for(int col=0 ; col<matrix.getColumnDimension() ; col++) {
			nullValues = 0;
			for(int row=0 ; row<matrix.getRowDimension() ; row++) {
				if ( Double.isNaN(matrix.get(row, col)) ) {
					nullValues++;
				}
			}
			if ( nullValues < limit ) {
				columnIndexes.add(col);
			}
		}
		return columnIndexes;
	}

	protected List<Integer> getRows(DoubleMatrix matrix, double filter) {
		int nullValues;
		List<Integer> rowIndexes = new ArrayList<Integer>();
		double limit = matrix.getColumnDimension() * filter / 100;
		System.out.println("getRow, " + filter + " % -> limit = " + limit);
		for(int row=0 ; row<matrix.getRowDimension() ; row++) {
			nullValues = 0;
			for(int col=0 ; col<matrix.getColumnDimension() ; col++) {
				if ( Double.isNaN(matrix.get(row, col)) ) {
					nullValues++;
				}
			}
			System.out.println("----> nullValues = " + nullValues);
			if ( nullValues < limit ) {
				rowIndexes.add(row);
			}
		}
		return rowIndexes;
	}

}

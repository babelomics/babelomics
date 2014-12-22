package org.bioinfo.babelomics.methods.genomic.copynumber;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.data.list.exception.InvalidIndexException;

public class CopyNumberUtils {

	public static void saveAsCGHInput(Dataset dataset, String outFileName) throws IOException, InvalidIndexException {
		dataset.load();
		
		List<String> sampleNames = dataset.getSampleNames();
		
		String header = "#ID\tChr\tStart\tEnd\t" +  ListUtils.toString(sampleNames, "\t") + "\n";
		
		StringBuilder line = new StringBuilder();;
		List<String> lines = new ArrayList<String>();

		DataFrame featureDataFrame = dataset.getFeatureData().getDataFrame();
		
		DataFrame dataFrame = new DataFrame();
		dataFrame.setRowNames(featureDataFrame.getColumn("ProbeName"));
		dataFrame.addColumn(featureDataFrame.getColumn("loc.chromosome"));
		dataFrame.addColumn(featureDataFrame.getColumn("loc.start"));
		dataFrame.addColumn(featureDataFrame.getColumn("loc.end"));
		for(int i=0 ; i<sampleNames.size() ; i++) {
			dataFrame.addColumn(ArrayUtils.toStringList(dataset.getDoubleMatrix().getColumn(i)));	
		}
		
		System.out.println("---> writting " + outFileName);
		IOUtils.write(new File(outFileName), header + dataFrame.toString(true, false));
		
//		System.out.println("feature column names = " + ListUtils.toString(dataset.getFeatureData().getDataFrame().getColumn("ProbeUID")));
		
//		line.append("#ID\tChr\tStart\tEnd\t").append(ListUtils.toString(dataset.getSampleNames(), "\t"));
//		lines.add(line.toString());
//		
//		for(int i=0 ; i<dataset.getRowDimension() ; i++) {
//			
//		}
//		
//		
//		IOUtils.write(new File(outFileName), lines);
	}
}

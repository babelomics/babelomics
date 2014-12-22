package org.bioinfo.babelomics.methods.expression.differential;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.math.result.LimmaTestResult;
import org.bioinfo.math.result.TestResultList;

public class Limma {
	
	protected String limmaBinPath;
	protected String inputFilename;
	protected List<String> classes;
	private List<String> contrast;
	protected List<String> filter;
	protected List<String> batch;
	

	private String outputFilename;

	public Limma (String limmaBinPath) {
		this.limmaBinPath = limmaBinPath;
	}
	
	public TestResultList<LimmaTestResult> compute() throws Exception {
		TestResultList<LimmaTestResult> result = null;
		
		if ( limmaBinPath == null ) {
			throw new InvalidParameterException("limma binary path missing");
		}
		
		if ( inputFilename == null ) {
			throw new InvalidParameterException("limma input filename parameter missing");
		}
		File file = FileUtils.checkFile(inputFilename);
		

		if ( classes == null ) {
			throw new InvalidParameterException("limma classes parameter missing");			
		}
		
		if ( contrast == null ) {
			throw new InvalidParameterException("limma contrast parameter missing");			
		}

		outputFilename = File.createTempFile("tmp", ".out").getAbsolutePath();
		
		List<String> env = new ArrayList<String>();
		env.add("datafile=" + inputFilename);		
		env.add("outfile=" + outputFilename);
		env.add("class=" + ListUtils.toString(classes, ","));
		env.add("contrast=" + ListUtils.toString(contrast, ","));
		if ( batch != null ) env.add("batch=" + ListUtils.toString(batch, ","));
		if ( filter != null ) env.add("filter=" + ListUtils.toString(filter, ","));
	
		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + limmaBinPath + " " + outputFilename + ".log", env);
		
		System.out.println("cmd = " + cmd.getCommandLine());
		System.out.println("env = " + ListUtils.toString(env, " "));
		
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
		
		FileUtils.checkFile(outputFilename);
		file = new File(outputFilename);
//		if ( ! file.exists() ) {
//			throw new Exception("error executing limma test (output file does not created)");
//		}

		Dataset ds = new Dataset(file, true);
		
//		new File(file.getAbsoluteFile() + ".log").delete();
//		file.delete();
				
		if ( ds == null || ds.getDoubleMatrix() == null || ds.getDoubleMatrix().getColumnDimension() != 2) {
			throw new Exception("error running limma test, output dataset can not be created");
		}
		
		result = new TestResultList<LimmaTestResult>(ds.getRowDimension());
		for(int i=0 ; i < ds.getRowDimension() ; i++) {
			result.add(new LimmaTestResult(ds.getDoubleMatrix().getRow(i)[0], ds.getDoubleMatrix().getRow(i)[1], 0.0));			
		}
		
		return result;
	}

	public String getInputFilename() {
		return inputFilename;
	}

	public void setInputFilename(String inputFilename) {
		this.inputFilename = inputFilename;
	}

	public String getLimmaBinPath() {
		return limmaBinPath;
	}

	public void setLimmaBinPath(String limmaBinPath) {
		this.limmaBinPath = limmaBinPath;
	}

	public List<String> getClasses() {
		return classes;
	}

	public void setClasses(List<String> classes) {
		this.classes = classes;
	}

	public List<String> getFilter() {
		return filter;
	}

	public void setFilter(List<String> filter) {
		this.filter = filter;
	}

	public List<String> getBatch() {
		return batch;
	}

	public void setBatch(List<String> batch) {
		this.batch = batch;
	}

	public String getOutputFilename() {
		return outputFilename;
	}

	public void setOutputFilename(String outputFilename) {
		this.outputFilename = outputFilename;
	}

	public void setContrast(List<String> contrast) {
		this.contrast = contrast;
	}

	public List<String> getContrast() {
		return contrast;
	}
}

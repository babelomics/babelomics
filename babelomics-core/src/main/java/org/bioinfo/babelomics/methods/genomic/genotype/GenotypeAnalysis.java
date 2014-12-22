package org.bioinfo.babelomics.methods.genomic.genotype;

import java.io.File;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.math.util.MathUtils;

public class GenotypeAnalysis {

	private File pedFile;
	private File mapFile;
	private String plinkPath;
	private String outdir;

	public GenotypeAnalysis() {
		this.pedFile = null;
		this.mapFile = null;
	}

	public GenotypeAnalysis(String pedFilePath, String mapFilePath) {
		this.pedFile = new File(pedFilePath);
		this.mapFile = new File(mapFilePath);
	}

	public GenotypeAnalysis(File pedFile, File mapFile) {
		this.pedFile = pedFile;
		this.mapFile = mapFile;
	}

	/*
	 * 
	 * PLINK TESTSlist
	 * 
	 */
	public void association(String assocTest) throws IOException {
		association(assocTest, 0.02);
	}

	public void association(String assocTest, double maf) throws IOException {
		checkParameters();
		if(assocTest == null) {
			throw new InvalidParameterException("association test is null");
		}
		if(!assocTest.equalsIgnoreCase("assoc") && !assocTest.equalsIgnoreCase("fisher") && !assocTest.equalsIgnoreCase("linear") && !assocTest.equalsIgnoreCase("logistic") && !assocTest.equalsIgnoreCase("tdt")) {
			throw new InvalidParameterException("association test is not valid, valid options are: 'assoc', 'fisher', 'linear', 'logistic' or 'tdt', parameter: " + assocTest);
		}
		// executing plink binary
		StringBuilder plinkCommandLine = createBasicPlinkCommand();
		plinkCommandLine.append(" --maf " + maf + " --" + assocTest.toLowerCase() + " ");
		executePlinkCommand(plinkCommandLine.toString());

		// saving the data
		FeatureData featureData;
		DataFrame dataFrame = new DataFrame();
		if(assocTest.equalsIgnoreCase("assoc")) {
			// columns:  CHR, SNP, BP, A1, F_A, F_U, A2, CHISQ, P, OR
			FileUtils.checkFile(outdir+"/plink.assoc");

		}
		if(assocTest.equalsIgnoreCase("fisher")) {
			// columns:  CHR, SNP, BP, A1, F_A, F_U, A2, P, OR
			try {
				FileUtils.checkFile(outdir+"/plink.assoc.fisher");
				dataFrame.addColumn("dbsnp", IOUtils.column(outdir+"/plink.assoc.fisher", 1, "\\s+"));
				dataFrame.addColumn("chromosome", IOUtils.column(outdir+"/plink.assoc.fisher", 0, "\\s+"));
				dataFrame.addColumn("position", IOUtils.column(outdir+"/plink.assoc.fisher", 2, "\\s+"));
				
				List<String> pvalues = IOUtils.column(outdir+"/plink.assoc.fisher", 7, "\\s+");
				double[] minusPvalueLog = MathUtils.log(ListUtils.toDoubleArray(pvalues), 2);
				minusPvalueLog = MathUtils.scalarMultiply(minusPvalueLog, -1);
				dataFrame.addColumn("p_values", pvalues);
				dataFrame.addColumn("log_p_values", ArrayUtils.toStringList(minusPvalueLog));
				dataFrame.addColumn("odd_ratio", IOUtils.column(outdir+"/plink.assoc.fisher", 8, "\\s+"));
				dataFrame.removeRow(0);
			} catch (InvalidIndexException e) {
				e.printStackTrace();
			}
		}
		if(assocTest.equalsIgnoreCase("linear")) {

		}
		if(assocTest.equalsIgnoreCase("logistic")) {
			// columns: CHR, SNP, BP, A1, TEST, NMISS, OR, STAT, P
			
		}
		if(assocTest.equalsIgnoreCase("tdt")) {

		}
		
		featureData = new FeatureData(dataFrame);
		featureData.save(new File(outdir+"/plink.featdata"));
	}

	public void stratification(String test, boolean ibsTest, double ppc, int mc, double maf) throws IOException {
		checkParameters();
		if(test == null) {
			throw new InvalidParameterException("association test is null");
		}
		// executing plink binary
		StringBuilder plinkCommandLine = createBasicPlinkCommand();
		if(test.equalsIgnoreCase("ibs-cluster")) {
			if(maf > 0) {
				plinkCommandLine.append(" --maf "+maf);
			}
			plinkCommandLine.append(" --cluster ");
			if(ibsTest) {
				plinkCommandLine.append(" --ibs-test ");
			}
			if(ppc > 0) {
				plinkCommandLine.append(" --ppc "+ppc);
			}
			if(mc > 0) {
				plinkCommandLine.append(" --mc "+mc);
			}
			executePlinkCommand(plinkCommandLine.toString());
		}
		
	}


	/*
	 * 
	 * PRIVATE METHODS
	 * 
	 */

	private StringBuilder createBasicPlinkCommand() {
		StringBuilder plinkCommandLine = new StringBuilder(plinkPath);
		plinkCommandLine.append(" --ped " + pedFile.getAbsolutePath());
		plinkCommandLine.append(" --map " + mapFile.getAbsolutePath());
		plinkCommandLine.append(" --out " + outdir + "/plink ");
		return plinkCommandLine;
	}

	private void executePlinkCommand(String command) throws IOException {
		System.out.println("executing: "+command);
//		List<String> params= new ArrayList<String>();
//		params.add(plinkPath);
//		params.add(" --ped "+pedFile.getAbsolutePath());
//		params.add(" --map "+mapFile.getAbsolutePath());
//		params.add(" --out "+outdir+"/plink");
//		params.add(" --assoc");
//		Process p = Runtime.getRuntime().exec(ListUtils.toArray(params));
		
//		Process p = Runtime.getRuntime().exec(command);
//		System.out.println(IOUtils.toString(p.getInputStream()));
//		System.out.println(IOUtils.toString(p.getErrorStream()));
//		System.out.println("exitValue: "+p.exitValue());
		
		Command plinkCommand = new Command(command);
		SingleProcess sp = new SingleProcess(plinkCommand);
		sp.runSync();
		
//		System.out.println(sp.getRunnableProcess().getError());
//		System.out.println(sp.getRunnableProcess().getException());
//		System.out.println(sp.getRunnableProcess().getExitValue());
//		System.out.println(sp.getRunnableProcess().getOutput());
	}

	private void checkParameters() throws IOException {
		FileUtils.checkFile(pedFile);
		FileUtils.checkFile(mapFile);
		FileUtils.checkFile(plinkPath);
		FileUtils.checkDirectory(outdir);
	}


	/**
	 * @param pedFile the pedFile to set
	 */
	public void setPedFile(File pedFile) {
		this.pedFile = pedFile;
	}


	/**
	 * @return the pedFile
	 */
	public File getPedFile() {
		return pedFile;
	}


	/**
	 * @param mapFile the mapFile to set
	 */
	public void setMapFile(File mapFile) {
		this.mapFile = mapFile;
	}


	/**
	 * @return the mapFile
	 */
	public File getMapFile() {
		return mapFile;
	}

	/**
	 * @param plinkPath the plinkPath to set
	 */
	public void setPlinkPath(String plinkPath) {
		this.plinkPath = plinkPath;
	}

	/**
	 * @return the plinkPath
	 */
	public String getPlinkPath() {
		return plinkPath;
	}

	/**
	 * @param outdir the outdir to set
	 */
	public void setOutdir(String outdir) {
		this.outdir = outdir;
	}

	/**
	 * @return the outdir
	 */
	public String getOutdir() {
		return outdir;
	}

}

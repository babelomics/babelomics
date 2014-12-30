package org.bioinfo.babelomics.methods.expression.differential;

import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.utils.ListUtils;

public class maSigPro {

	protected String maSigProBinPath;
	protected String inputFilename;
	protected String outdir;

	public maSigPro(String maSigProBinPath) {
		this.maSigProBinPath = maSigProBinPath;
	}

	public void compute() {
		compute(2, 0.05, "BH", 0.05, "hclust", 9);
	}

	public void compute(int degree, double q, String correction, double alfa, String clustering, int kvalue) {
		List<String> env = new ArrayList<String>();
		env.add("data=" + inputFilename);
		env.add("degree=" + degree);
		env.add("Q=" + q);
		env.add("adjust=" + correction);
		env.add("alfa=" + alfa);
		env.add("clustermethod=" + clustering);
		env.add("k=" + kvalue);
		env.add("main=");
		env.add("outdir=" + outdir);

		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + maSigProBinPath + " " + outdir + "/r.log", env);
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();

		System.out.println("cmd line = " + cmd.getCommandLine());
		System.out.println("cmd env = " + ListUtils.toString(env, " "));
	}

	public String getMaSigProBinPath() {
		return maSigProBinPath;
	}

	public void setMaSigProBinPath(String maSigProBinPath) {
		this.maSigProBinPath = maSigProBinPath;
	}

	public String getInputFilename() {
		return inputFilename;
	}

	public void setInputFilename(String inputFilename) {
		this.inputFilename = inputFilename;
	}

	public String getOutdir() {
		return outdir;
	}

	public void setOutdir(String outdir) {
		this.outdir = outdir;
	}

}

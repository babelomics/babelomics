package org.bioinfo.babelomics.utils;

import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.utils.ListUtils;

public class RCommand {

	private String scriptPath;
	private String outdir;
	private List<String> params;
	private String logFilename;
	
	public RCommand(String scriptPath, String outdir) {
		this(scriptPath, outdir, "r.log");
	}
	
	public RCommand(String scriptPath, String outdir, String logFilename) {
		this.scriptPath = scriptPath;
		this.outdir = outdir;
		this.logFilename = logFilename;
		this.params = new ArrayList<String>();
	}

	public void addParam(String key, String value){
		params.add(key + "=" + value);
	}
	public void addParam(String key, double value){
		params.add(key + "=" + value);
	}
	public void addParam(String key, int value){
		params.add(key + "=" + value);
	}
	public void addParam(String key, boolean value){
		params.add(key + "=" + value);
	}
	
	public void exec(){
		addParam("outdir", outdir);
		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + scriptPath + " " + outdir + "/" + logFilename, params);
		System.out.println("=====cmd line = " + cmd.getCommandLine());
		System.out.println("=====cmd env = " + ListUtils.toString(params, " "));
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();	
	}
	
}

package org.bioinfo.babelomics.methods.interactome;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;

public class SnowTest {

	String mode = "unknown";
	String interactome = "unknown";
	String list1FileName = null;
	String list2FileName = null;
	String interactionsFileName = null;
	boolean checkInteractions = true;
	String idNature = "protein";
	int interactionsNumber = 1;
	
	String snowBinPath = null;
	String outDirName = null;
	
	public SnowTest(String snowBinPath, String list1FileName, String interactome, String idNature, int interactionsNumber, String outDirName) {
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.interactome = interactome;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;

		this.mode = "one-list";
	}
	
	public SnowTest(String snowBinPath, String list1FileName, String interactionsFileName, boolean checkInteractions, String idNature, int interactionsNumber, String outDirName) {
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.interactome = "own";
		this.interactionsFileName = interactionsFileName;
		this.checkInteractions = checkInteractions;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;

		this.mode = "one-list";
	}

	public SnowTest(String snowBinPath, String list1FileName, String list2FileName, String interactome, String idNature, int interactionsNumber, String outDirName) {
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.list2FileName = list2FileName;
		this.interactome = interactome;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;
		
		this.mode = "two-lists";
	}

	public SnowTest(String snowBinPath, String list1FileName, String list2FileName, String interactionsFileName, boolean checkInteractions, String idNature, int interactionsNumber, String outDirName) {		
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.list2FileName = list2FileName;
		this.interactome = "own";
		this.interactionsFileName = interactionsFileName;
		this.checkInteractions = checkInteractions;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;
		
		this.mode = "two-lists";
	}


	public void run() throws InvalidParameterException {
		if ( snowBinPath == null ) {
			throw new InvalidParameterException("Missing binary path to snow script");
		}
		
		String cmd = createSnowCommand();
		
		System.out.println("SnowTest.java:run, cmd = " + cmd);

		executeSnowCommand(cmd);		
	}
	
	private String createSnowCommand() {
		StringBuilder cmd = new StringBuilder(snowBinPath);

		//die "usage: ./snow.pl <interactome> <specie> <nature> <mode> <limit> <jobname> <list-file-1> <list-file-2> <interactome-file> <same-id> <out-dir>\n";
		cmd.append(" ").append(this.interactome);
		cmd.append(" ").append("hsa");
		cmd.append(" ").append(this.idNature);
		cmd.append(" ").append(this.mode);
		cmd.append(" ").append(this.interactionsNumber);
		cmd.append(" ").append("snow");
		cmd.append(" ").append(list1FileName);
		cmd.append(" ").append("two-lists".equalsIgnoreCase(mode) ? this.list2FileName : "unknown");
		cmd.append(" ").append("own".equalsIgnoreCase(interactome) ? this.interactionsFileName : "unknown");
		cmd.append(" ").append("own".equalsIgnoreCase(interactome) ? (this.checkInteractions ? "on" : "off") : "unknown");
		cmd.append(" ").append(this.outDirName);
			
		return cmd.toString();
	}

	private void executeSnowCommand(String cmd) {
		Command snowCommand = new Command(cmd);
		SingleProcess sp = new SingleProcess(snowCommand);
		sp.runSync();
	}
}

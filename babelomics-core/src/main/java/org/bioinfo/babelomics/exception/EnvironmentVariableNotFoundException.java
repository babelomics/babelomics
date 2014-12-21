package org.bioinfo.babelomics.exception;

@SuppressWarnings("serial")
public class EnvironmentVariableNotFoundException extends Exception {

	public EnvironmentVariableNotFoundException(String environmentVariable) {
		super("Environment variable " + environmentVariable + " is not defined");
	}
	
}

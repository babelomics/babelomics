package org.bioinfo.babelomics.exception;

public class EmptyAnnotationException extends Exception{
	public EmptyAnnotationException() {		
		this("Annotations is empty for current ids");
	}
	public EmptyAnnotationException(String msg) {		
		super(msg);
	}
}

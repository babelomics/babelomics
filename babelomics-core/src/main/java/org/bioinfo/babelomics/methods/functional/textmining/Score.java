package org.bioinfo.babelomics.methods.functional.textmining;

public class Score {
	private String entity;
	private String gene;
	private double value;

	public Score(String entity, String gene, Double value) {
		this.entity = entity;
		this.gene = gene;
		this.value = value;
	}

	public String getEntity() {
		return entity;
	}

	public void setEntity(String entity) {
		this.entity = entity;
	}

	public String getGene() {
		return gene;
	}

	public void setGene(String gene) {
		this.gene = gene;
	}

	public void setValue(double value) {
		this.value = value;
	}

	public double getValue() {
		return value;
	}
	
	public String toString() {
		return "[" + entity + ", " + gene + ", " + value + "]";
	}
	
}

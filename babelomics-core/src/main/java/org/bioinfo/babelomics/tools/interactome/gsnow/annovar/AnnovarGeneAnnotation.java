package org.bioinfo.babelomics.tools.interactome.gsnow.annovar;

public class AnnovarGeneAnnotation {

	private String consequenceType;
	private String gene;
	private String chromosome;
	private int start;
	private int end;
	private String ref;
	private String alt;
	private String strand;
	
	/**
	 * @param consequenceType
	 * @param gene
	 * @param chromosome
	 * @param start
	 * @param end
	 * @param ref
	 * @param alt
	 * @param strand
	 */
	public AnnovarGeneAnnotation(String consequenceType,String gene,String chromosome,int start,int end,String ref,String alt,String strand) {
		this.consequenceType = consequenceType;
		this.gene = gene;
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
		this.ref = ref;
		this.alt = alt;
		this.strand = strand;
	}

	
	
	/**
	 * @return the consequenceType
	 */
	public String getConsequenceType() {
		return consequenceType;
	}

	/**
	 * @param consequenceType the consequenceType to set
	 */
	public void setConsequenceType(String consequenceType) {
		this.consequenceType = consequenceType;
	}

	/**
	 * @return the gene
	 */
	public String getGene() {
		return gene;
	}

	/**
	 * @param gene the gene to set
	 */
	public void setGene(String gene) {
		this.gene = gene;
	}

	/**
	 * @return the chromosome
	 */
	public String getChromosome() {
		return chromosome;
	}

	/**
	 * @param chromosome the chromosome to set
	 */
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}

	/**
	 * @param start the start to set
	 */
	public void setStart(int start) {
		this.start = start;
	}

	/**
	 * @return the end
	 */
	public int getEnd() {
		return end;
	}

	/**
	 * @param end the end to set
	 */
	public void setEnd(int end) {
		this.end = end;
	}

	/**
	 * @return the ref
	 */
	public String getRef() {
		return ref;
	}

	/**
	 * @param ref the ref to set
	 */
	public void setRef(String ref) {
		this.ref = ref;
	}

	/**
	 * @return the alt
	 */
	public String getAlt() {
		return alt;
	}

	/**
	 * @param alt the alt to set
	 */
	public void setAlt(String alt) {
		this.alt = alt;
	}

	/**
	 * @return the strand
	 */
	public String getStrand() {
		return strand;
	}

	/**
	 * @param strand the strand to set
	 */
	public void setStrand(String strand) {
		this.strand = strand;
	}
	
	
	
		
	
}

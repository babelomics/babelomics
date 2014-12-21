package org.bioinfo.babelomics.tools.interactome.gsnow.annovar;

import java.util.HashMap;

public class Variant{
	
	private String chromosome;
	private String position;
	private String id;
	private String ref;
	private String alt;
	private String qual;
	private String filter;
	private HashMap<String,String> info;
	private HashMap<String,String> genotype;
	private boolean missingValue;
	private boolean filtered;
	
	public Variant(String record){
		
		// raw read
		String[] fields = record.split("\t");

		// org.bioinfo.ngs.genomizer.main fields
		this.chromosome = fields[0];
		this.position = fields[1];
		this.id = fields[2];
		this.ref = fields[3];
		this.alt = fields[4];
		this.qual = fields[5];
		this.filter = fields[6];			

		// info
		info = loadInfo(fields[7]);
		
//		// genotype
//		if(fields.length==10){
//			genotype = loadFormat(fields[8],fields[9]);
//		} else {			
//			genotype = loadFormat(fields[8],fields[10]);
//		}
	
		if(alt.equals(".")){
			missingValue = true; 
		}		
		if(!filter.equals("PASS") && !filter.equals(".") ){
			filtered = true;
		}
		
	}
	
	public Variant(String chromosome, String position,String id, String ref, String alt, String qual,String filter, String canonicalAllele,int coverage,int refCount, int altCount){
		
		this.chromosome = chromosome;
		this.position = position;
		this.id = id;
		this.ref = ref;
		this.alt = alt;
		this.qual = qual;
		this.filter = filter;
		
		this.genotype = new HashMap<String,String>();
		genotype.put("GT", ""+canonicalAllele);
		genotype.put("DP", ""+coverage);
		genotype.put("AD", refCount + "," + altCount);
	
		if(genotype.get("GT").equals("./.")){
			missingValue = true; 
		}		
		if(!filter.equals("PASS") && !filter.equals(".") ){
			filtered = true;
		}
		
	}
	
	public Variant(String chromosome, String position, String ref){
		this(chromosome,position,".", ref, ".", "0",".", "./.",0,0,0);		
	}
	
			
	public HashMap<String,String> loadInfo(String raw){
		HashMap<String,String> info = new HashMap<String, String>();
		String[] fields;
		if(raw.contains(";")){
			fields = raw.split(";");
		} else {
			fields = new String[]{raw};	
		}
		String field;
		String[] tag;
		for(int i=0; i<fields.length; i++){
			field = fields[i];
			if(field.contains("=")){
				tag = field.split("=");
				info.put(tag[0],tag[1]);
			} else {
				info.put(field,"");
			}
		}
		return info;
	}
	
	public HashMap<String,String> loadFormat(String rawFormat,String rawGenotype){
		HashMap<String,String> format = new HashMap<String, String>();		
		if(rawFormat.contains(":") && rawGenotype.contains(":")){
			String[] keys;
			String[] values;
			keys = rawFormat.split(":");
			values = rawGenotype.split(":");
			for(int i=0; i<keys.length; i++){				
				format.put(keys[i].trim(),values[i].trim());
			}
		} else {
			format.put(rawFormat,rawGenotype);
		}
		return format;
	}

	/*
	 * 
	 * extended fields
	 * 
	 */
	public String getCanonicalAllele(){
		return genotype.get("GT");
	}
	
	public String getAllele(){
		String canonical = genotype.get("GT");		
		return canonical.replaceAll("/", "").replaceAll("0",ref).replaceAll("1",alt);
	}
	
	public boolean isHomozygous(){
		String canonical = genotype.get("GT");
		return (canonical.equalsIgnoreCase("0/0") || canonical.equalsIgnoreCase("1/1")); 
	}
	
	public int getCoverage(){
		if(info.containsKey("DP")){
			return Integer.parseInt(info.get("DP"));	
		} else {
			return Integer.parseInt(genotype.get("DP"));
		}
		
	}
	
	public void setCoverage(int coverage){
		info.put("DP", ""+coverage);		
	}
	
	public int getReferenceAlleleCount(){		
		if(genotype.containsKey("AD")){
			String[]fields = genotype.get("AD").split(",");
			return Integer.parseInt(fields[0]);
		} else {
			return -1;
		}
	}
	
	public int getAlternativeAlleleCount(){
		if(genotype.containsKey("AD")){
			String[]fields = genotype.get("AD").split(",");
			return Integer.parseInt(fields[1]);
		} else {
			return -1;
		}
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
	 * @return the position
	 */
	public String getPosition() {
		return position;
	}

	/**
	 * @param position the position to set
	 */
	public void setPosition(String position) {
		this.position = position;
	}

	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}

	/**
	 * @param id the id to set
	 */
	public void setId(String id) {
		this.id = id;
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
	 * @return the qual
	 */
	public String getQual() {
		return qual;
	}

	/**
	 * @param qual the qual to set
	 */
	public void setQual(String qual) {
		this.qual = qual;
	}

	/**
	 * @return the filter
	 */
	public String getFilter() {
		return filter;
	}

	/**
	 * @param filter the filter to set
	 */
	public void setFilter(String filter) {
		this.filter = filter;
	}

	/**
	 * @return the info
	 */
	public HashMap<String, String> getInfo() {
		return info;
	}

	/**
	 * @param info the info to set
	 */
	public void setInfo(HashMap<String, String> info) {
		this.info = info;
	}

	/**
	 * @return the genotype
	 */
	public HashMap<String, String> getGenotype() {
		return genotype;
	}

	/**
	 * @param genotype the genotype to set
	 */
	public void setGenotype(HashMap<String, String> genotype) {
		this.genotype = genotype;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "Vcf [alt=" + alt + ", chromosome=" + chromosome + ", filter="
				+ filter + ", genotype=" + genotype + ", id=" + id + ", info="
				+ info + ", position=" + position + ", qual=" + qual + ", ref="
				+ ref + "]";
	}

	/**
	 * @return the missingValue
	 */
	public boolean isMissingValue() {
		return missingValue;
	}

	/**
	 * @param missingValue the missingValue to set
	 */
	public void setMissingValue(boolean missingValue) {
		this.missingValue = missingValue;
	}

	/**
	 * @return the filtered
	 */
	public boolean isFiltered() {
		return filtered;
	}

	/**
	 * @param filtered the filtered to set
	 */
	public void setFiltered(boolean filtered) {
		this.filtered = filtered;
	}
	
}

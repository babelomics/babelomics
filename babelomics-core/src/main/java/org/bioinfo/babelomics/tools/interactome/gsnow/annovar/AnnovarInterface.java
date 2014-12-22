package org.bioinfo.babelomics.tools.interactome.gsnow.annovar;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;

public class AnnovarInterface {

	private String annovarEnginePath;
	private String humanDbPath;
	private String annovarFile;
		
	public AnnovarInterface(String annovarEnginePath, String humanDbPath){
		this.annovarEnginePath = annovarEnginePath;
		this.humanDbPath = humanDbPath;
	}

	public List<String> annotate(Vcf vcf, String outputFolder) throws SQLException, IOException{		
		annovarFile = prepareInput(vcf, outputFolder);
		annotate(annovarFile,outputFolder);
		HashMap<String,AnnovarExonAnnotation> exonAnnotations = loadExonicAnnotations(outputFolder  + "/family_variant_set_annovar.temp.exonic_variant_function");
		return getGenes(exonAnnotations, "nonsynonymous SNV");
		//insertAnnotations(vcf ,outputFolder);
	}
	
	public void annotate(String annovarInputFile, String outputFolder){
//		String commandLine = "perl " + annovarEnginePath + "/annotate_variation.pl --buildver hg19 " + annovarInputFile + " " + humanDbPath + "/ > /dev/null 2> /dev/null";
		String commandLine = "perl " + annovarEnginePath + "/annotate_variation.pl --buildver hg19 " + annovarInputFile + " " + humanDbPath + "/";

		Command command = new Command(commandLine);
		SingleProcess process = new SingleProcess(command);	
		process.getRunnableProcess().run();
	}
	
	private String prepareInput(Vcf vcf, String outputFolder) throws FileNotFoundException, SQLException{
		
		// preparing annovar input		
		String annovarFile = outputFolder  + "/family_variant_set_annovar.temp";
		PrintWriter annovarWriter = new PrintWriter(annovarFile);
		
//		List<Integer> variantIds = familyVariantDb.getVariantIds();
//		FamilyBasedVariantRecord variant;
		for (Variant variant : vcf.getVariantList()) {
			annovarWriter.print(variant.getChromosome() + "\t");
			annovarWriter.print(variant.getPosition() + "\t");
			annovarWriter.print(variant.getPosition() + "\t");
			annovarWriter.print(variant.getRef() + "\t");
			annovarWriter.print(variant.getAlt() + "\t");				
			annovarWriter.println("+");
		}
		annovarWriter.close();
		
		return annovarFile;		
	}
	
//	private void insertAnnotations(Vcf vcf, String outputFolder) throws IOException, SQLException{
//		
//		// load exonic annotations		
//		HashMap<String,AnnovarExonAnnotation> exonicAnnotations =  loadExonicAnnotations(outputFolder  + "/family_variant_set_annovar.temp.exonic_variant_function");
//		
//		// load gene annotations
//		String geneAnnotationFile = outputFolder  + "/family_variant_set_annovar.temp.variant_function";
//		BufferedReader geneReader = new BufferedReader(new FileReader(geneAnnotationFile));
//		String line,key;
//		String[] fields;
//		while((line=geneReader.readLine())!=null){
//			if(!line.startsWith("#") && line.contains("\t")){
//				fields = line.split("\t");
//				
//				// load geneAnnotation
//				AnnovarGeneAnnotation geneAnnotation = new AnnovarGeneAnnotation(fields[0],fields[1],fields[2],Integer.parseInt(fields[3]),Integer.parseInt(fields[4]),fields[5],fields[6],fields[7]);
//				
//				// get variant
//				int variantId = familyVariantDb.getVariantId(geneAnnotation.getChromosome(),geneAnnotation.getStart(),geneAnnotation.getEnd(),geneAnnotation.getRef(),geneAnnotation.getAlt());				
//				FamilyBasedVariantRecord record = familyVariantDb.getVariant(variantId);				
//
//				// add annotations
//				record.getFunctionalProfilingDescriptor().setGene(geneAnnotation.getGene());
//				key = geneAnnotation.getChromosome() + ":" + geneAnnotation.getStart() + "_" + geneAnnotation.getEnd() + "__" + geneAnnotation.getRef() + "_" + geneAnnotation.getAlt();
//				if(exonicAnnotations.containsKey(key)){
//					AnnovarExonAnnotation exonAnnotation = exonicAnnotations.get(key);
//					record.getFunctionalProfilingDescriptor().setConsequenceType(exonAnnotation.getConsequenceType());
//				} else {
//					record.getFunctionalProfilingDescriptor().setConsequenceType(geneAnnotation.getConsequenceType());
//				}		
//				
//				// update variant
//				familyVariantDb.saveVariant(record);
//				
//			}
//		}
//			
//		
//	}
	
	private HashMap<String,AnnovarExonAnnotation> loadExonicAnnotations(String exonAnnotationFile) throws IOException{
		HashMap<String,AnnovarExonAnnotation> exonAnnotations = new HashMap<String, AnnovarExonAnnotation>();
		BufferedReader exonReader = new BufferedReader(new FileReader(exonAnnotationFile));
		String line,key;
		String[] fields;
		while((line=exonReader.readLine())!=null){
			if(!line.startsWith("#") && line.contains("\t")){
				fields = line.split("\t");
				AnnovarExonAnnotation exonAnnotation = new AnnovarExonAnnotation(fields[1],fields[2],fields[3],Integer.parseInt(fields[4]),Integer.parseInt(fields[5]),fields[6],fields[7],fields[8]);
				key = exonAnnotation.getChromosome() + ":" + exonAnnotation.getStart() + "_" + exonAnnotation.getEnd() + "__" + exonAnnotation.getRef() + "_" + exonAnnotation.getAlt();
				exonAnnotations.put(key,exonAnnotation);
			}
		}
		return exonAnnotations;
	}
	private List<String> getGenes(HashMap<String,AnnovarExonAnnotation> exonAnnotations, String type){
		List<String> genes = new ArrayList<String>();
		Set<String> setGenes = new HashSet<String>();
		
		for(String key : exonAnnotations.keySet()) {
			AnnovarExonAnnotation exonAnnotation = exonAnnotations.get(key);
			if(exonAnnotation.getConsequenceType().contains(type)){
				String fields[] = exonAnnotation.getExons().split(":");
				if(fields.length > 1){
					setGenes.add(fields[0]);
				}
					
			}
		}
		for (String gene : setGenes) {
			genes.add(gene);
		}
		return genes;
	}


	/**
	 * @return the annovarFile
	 */
	public String getAnnovarFile() {
		return annovarFile;
	}

	/**
	 * @param annovarFile the annovarFile to set
	 */
	public void setAnnovarFile(String annovarFile) {
		this.annovarFile = annovarFile;
	}
	
	
}

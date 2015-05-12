package org.bioinfo.babelomics.tools.functional;

import org.bioinfo.infrared.funcannot.filter.BiocartaFilter;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.GOSlimFilter;
import org.bioinfo.infrared.funcannot.filter.InterproFilter;
import org.bioinfo.infrared.funcannot.filter.JasparFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
import org.bioinfo.infrared.funcannot.filter.MiRnaTargetFilter;
import org.bioinfo.infrared.funcannot.filter.OregannoFilter;
import org.bioinfo.infrared.funcannot.filter.ReactomeFilter;
import org.bioinfo.babelomics.utils.filters.ReconFilter;


public class FunctionalDbDescriptor {

	private String name;
	private String title;
	private String prefix;
	private String description;
	
	public FunctionalDbDescriptor(FunctionalFilter filter){
		initFromFunctionalFilter(filter);
	}
	
	public FunctionalDbDescriptor(String name, String title, String prefix, String description){
		this.name = name;
		this.title = title;
		this.prefix = prefix;
		this.description = description;		
	}
	
	protected void initFromFunctionalFilter(final FunctionalFilter filter){
		//String levels;
		// Untitled
		title = "Untitled";
		// GO
		if(filter instanceof GOFilter) {						
			GOFilter goFilter = (GOFilter) filter;
//			if(goFilter.getMinLevel()==goFilter.getMaxLevel()) {
//				levels = "(level " + goFilter.getMinLevel() + ")";
//			} else{
//				levels = "(levels from " + goFilter.getMinLevel() + " to " + goFilter.getMaxLevel() + ")";
//			}
			title = "GO " + goFilter.getNamespace().replace("_", " ");// + " " + levels;
//			name = "go_" + goFilter.getNamespace() + "_" + goFilter.getMinLevel() + "_" + goFilter.getMaxLevel();
			name = "go_" + goFilter.getNamespace();
			prefix = "go";
			description = getDefaultDescription(filter, title);			
			if(goFilter.isPropagated()) {
				description+= ", each term parent within levels has been included";
			}
			if(goFilter.getKeywords()!=null && goFilter.getKeywords().size()>0) {
				description+= " and terms have been filtered by using the keywords " + goFilter.getKeywords().toString();
			}
		}
		// GOSlim GOA
		else if(filter instanceof GOSlimFilter) {
			title = "GOSlim GOA";
			name = "go-slim";
			prefix = "go-slim";
			description = getDefaultDescription(filter, title);
		}
		// Kegg
		else if(filter instanceof KeggFilter) {
			KeggFilter keggFilter = (KeggFilter) filter;
			title = "KEGG";
			name = "kegg";
			prefix = "kegg";
			description = getDefaultDescription(filter, title);
			if(keggFilter.getKeywords()!=null && keggFilter.getKeywords().size()>0) {
				description+= " and terms have been filtered by using the keywords " + keggFilter.getKeywords().toString();
			}
		} 
		// InterPro
		else if(filter instanceof InterproFilter) {
			title = "InterPro";
			name = "interpro";
			prefix = "interpro";
			description = getDefaultDescription(filter, title);
		}
		// Reactome
		else if(filter instanceof ReactomeFilter) {
			title = "Reactome";
			name = "reactome";
			prefix = "reactome";
			description = getDefaultDescription(filter, title);
		}
		// BioCarta
		else if(filter instanceof BiocartaFilter) {
			title = "Biocarta";
			name = "biocarta";
			prefix = "biocarta";
			description = getDefaultDescription(filter, title);
		}
		// miRNA target
		else if(filter instanceof MiRnaTargetFilter) {
			title = "miRNA target";
			name = "mirna";
			prefix = "mirna";
			description = getDefaultDescription(filter, title);
		}
		// Jaspar
		else if(filter instanceof JasparFilter) {
			title = "Jaspar TFBS";
			name = "jaspar";
			prefix = "jaspar";
			description = getDefaultDescription(filter, title);
		}
		// Oreganno
		else if(filter instanceof OregannoFilter) {
			title = "ORegAnno";
			name = "oreganno";
			prefix = "oreganno";
			description = getDefaultDescription(filter, title);
		}
		else if(filter instanceof ReconFilter) {
			title = "Recon";
			name = "recon";
			prefix = "recon";
			description = getDefaultDescription(filter, title);
		}
	
	}
	
	private String getDefaultDescription(FunctionalFilter filter, String title){
		StringBuilder description = new StringBuilder();
		description.append(title);
		description.append(", allowed range of term annotations among " + filter.getMinNumberGenes() + " to " + filter.getMaxNumberGenes());
		if(filter.isGenomicNumberOfGenes()) description.append(" (from genome)");
		else description.append(" (from your input data)");
		return description.toString();
	}


	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}


	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}


	/**
	 * @return the title
	 */
	public String getTitle() {
		return title;
	}


	/**
	 * @param title the title to set
	 */
	public void setTitle(String title) {
		this.title = title;
	}


	/**
	 * @return the prefix
	 */
	public String getPrefix() {
		return prefix;
	}


	/**
	 * @param prefix the prefix to set
	 */
	public void setPrefix(String prefix) {
		this.prefix = prefix;
	}


	/**
	 * @return the description
	 */
	public String getDescription() {
		return description;
	}


	/**
	 * @param description the description to set
	 */
	public void setDescription(String description) {
		this.description = description;
	}
}


package org.bioinfo.babelomics.tools.interactome.gsnow.annovar;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


public class Vcf {
	
	private List<Variant> variantList;

	public Vcf(String vcfPath) {
		variantList = new ArrayList<Variant>();
		try {
			BufferedReader br = getBufferedReader(vcfPath);
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#"))
					continue;
				Variant variant = new Variant(line);
				variantList.add(variant);
			}
		} catch (Exception e) {
			System.err.println(e.getMessage());
		}
	}

	public List<Variant> getVariantList() {
		return variantList;
	}

	public void setVariantList(List<Variant> variantList) {
		this.variantList = variantList;
	}
	
	private BufferedReader getBufferedReader(String file){
		FileInputStream fstream;
		try {
			fstream = new FileInputStream(file);
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			return br;
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		
	}
}

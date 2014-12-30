package org.bioinfo.babelomics.methods.functional;

import java.text.DecimalFormat;
import java.util.List;

import org.bioinfo.commons.utils.ListUtils;

public class TwoListFisherTestResult {
	
	private DecimalFormat percentageFormatter = new DecimalFormat("##.##");
	private DecimalFormat pvalueFormatter = new DecimalFormat("#.#####E0");
	private DecimalFormat oddsFormatter = new DecimalFormat("######.#####E0");
	
	private String term;
	private int termSize;
	private int termSizeInGenome;
	private int list1Positives;
	private int list1Negatives;
	private double list1Percentage;
	private int list2Positives;
	private int list2Negatives;
	private double list2Percentage;
	private List<String> list1Ids;
	private List<String> list2Ids;
	private double oddsRatio;
	private double pValue;
	private double adjPValue;
		
	
	/**
	 * @param term
	 * @param list1Positives
	 * @param list1Negatives
	 * @param list2Positives
	 * @param list2Negatives
	 * @param list1Ids
	 * @param list2Ids
	 * @param value
	 * @param adjPValue
	 */
	public TwoListFisherTestResult(String term, int termSize, int termSizeInGenome, int list1Positives, int list1Negatives, int list2Positives, int list2Negatives, List<String> list1Ids, List<String> list2Ids, double oddsRatio, double pValue, double adjPValue) {
		this.term = term;
		this.termSize = termSize;
		this.termSizeInGenome = termSizeInGenome;
		this.list1Positives = list1Positives;
		this.list1Negatives = list1Negatives;		
		this.list1Percentage = 100.0*((double)list1Positives/(double)(list1Positives+list1Negatives));		
		this.list2Positives = list2Positives;
		this.list2Negatives = list2Negatives;
		this.list2Percentage = 100.0*((double)list2Positives/(double)(list2Positives+list2Negatives));
		this.list1Ids = list1Ids;
		this.list2Ids = list2Ids;
		this.oddsRatio = oddsRatio;
		this.pValue = pValue;
		this.adjPValue = adjPValue;
	}

	public String toString(){
		StringBuilder out = new StringBuilder();
		out.append(this.term).append("\t");
		out.append(this.termSize).append("\t");
		out.append(this.termSizeInGenome).append("\t");
		out.append(this.list1Positives).append("\t");
		out.append(this.list1Negatives).append("\t");
		out.append(percentageFormatter.format(this.list1Percentage)).append("\t");
		out.append(this.list2Positives).append("\t");
		out.append(this.list2Negatives).append("\t");
		out.append(percentageFormatter.format(this.list2Percentage)).append("\t");
		out.append(ListUtils.toString(this.list1Ids,",")).append("\t");
		out.append(ListUtils.toString(this.list2Ids,",")).append("\t");
		double oddsRatioLog = Math.log(oddsRatio);
		double odds;
		if(oddsRatioLog == Double.POSITIVE_INFINITY){
			odds=Double.MAX_VALUE/10000;
		} else if(oddsRatioLog == Double.NEGATIVE_INFINITY){
			odds=-Double.MAX_VALUE/10000;		
		} else odds=oddsRatioLog;
		out.append(oddsFormatter.format(odds)).append("\t");
		out.append(pvalueFormatter.format(this.pValue)).append("\t");
		out.append(pvalueFormatter.format(this.adjPValue));
		return out.toString();
	}

	public static String header(){
		StringBuilder out = new StringBuilder();
		out.append("#term").append("\t");
		out.append("term_size").append("\t");
		out.append("term_size_in_genome").append("\t");
		out.append("list1_positives").append("\t");
		out.append("list1_negatives").append("\t");
		out.append("list1_percentage").append("\t");
		out.append("list2_positives").append("\t");
		out.append("list2_negatives").append("\t");
		out.append("list2_percentage").append("\t");
		out.append("list1_positive_ids").append("\t");
		out.append("list2_positive_ids").append("\t");
		out.append("odds_ratio_log").append("\t");
		out.append("pvalue").append("\t");
		out.append("adj_pvalue");
		return out.toString();
	}
	
	/**
	 * @return the term
	 */
	public String getTerm() {
		return term;
	}


	/**
	 * @param term the term to set
	 */
	public void setTerm(String term) {
		this.term = term;
	}


	/**
	 * @return the list1Positives
	 */
	public int getList1Positives() {
		return list1Positives;
	}


	/**
	 * @param list1Positives the list1Positives to set
	 */
	public void setList1Positives(int list1Positives) {
		this.list1Positives = list1Positives;
	}


	/**
	 * @return the list1Negatives
	 */
	public int getList1Negatives() {
		return list1Negatives;
	}


	/**
	 * @param list1Negatives the list1Negatives to set
	 */
	public void setList1Negatives(int list1Negatives) {
		this.list1Negatives = list1Negatives;
	}


	/**
	 * @return the list2Positives
	 */
	public int getList2Positives() {
		return list2Positives;
	}


	/**
	 * @param list2Positives the list2Positives to set
	 */
	public void setList2Positives(int list2Positives) {
		this.list2Positives = list2Positives;
	}


	/**
	 * @return the list2Negatives
	 */
	public int getList2Negatives() {
		return list2Negatives;
	}


	/**
	 * @param list2Negatives the list2Negatives to set
	 */
	public void setList2Negatives(int list2Negatives) {
		this.list2Negatives = list2Negatives;
	}


	/**
	 * @return the list1Ids
	 */
	public List<String> getList1Ids() {
		return list1Ids;
	}


	/**
	 * @param list1Ids the list1Ids to set
	 */
	public void setList1Ids(List<String> list1Ids) {
		this.list1Ids = list1Ids;
	}


	/**
	 * @return the list2Ids
	 */
	public List<String> getList2Ids() {
		return list2Ids;
	}


	/**
	 * @param list2Ids the list2Ids to set
	 */
	public void setList2Ids(List<String> list2Ids) {
		this.list2Ids = list2Ids;
	}


	/**
	 * @return the pValue
	 */
	public double getPValue() {
		return pValue;
	}


	/**
	 * @param value the pValue to set
	 */
	public void setPValue(double value) {
		pValue = value;
	}


	/**
	 * @return the adjPValue
	 */
	public double getAdjPValue() {
		return adjPValue;
	}


	/**
	 * @param adjPValue the adjPValue to set
	 */
	public void setAdjPValue(double adjPValue) {
		this.adjPValue = adjPValue;
	}

	/**
	 * @return the termSize
	 */
	public int getTermSize() {
		return termSize;
	}

	/**
	 * @param termSize the termSize to set
	 */
	public void setTermSize(int termSize) {
		this.termSize = termSize;
	}

	/**
	 * @return the termSizeInGenome
	 */
	public int getTermSizeInGenome() {
		return termSizeInGenome;
	}

	/**
	 * @param termSizeInGenome the termSizeInGenome to set
	 */
	public void setTermSizeInGenome(int termSizeInGenome) {
		this.termSizeInGenome = termSizeInGenome;
	}

	/**
	 * @return the oddsRatio
	 */
	public double getOddsRatio() {
		return oddsRatio;
	}

	/**
	 * @param oddsRatio the oddsRatio to set
	 */
	public void setOddsRatio(double oddRatio) {
		this.oddsRatio = oddRatio;
	}

	/**
	 * @return the list1Percentage
	 */
	public double getList1Percentage() {
		return list1Percentage;
	}

	/**
	 * @param list1Percentage the list1Percentage to set
	 */
	public void setList1Percentage(double list1Percentage) {
		this.list1Percentage = list1Percentage;
	}

	/**
	 * @return the list2Percentage
	 */
	public double getList2Percentage() {
		return list2Percentage;
	}

	/**
	 * @param list2Percentage the list2Percentage to set
	 */
	public void setList2Percentage(double list2Percentage) {
		this.list2Percentage = list2Percentage;
	}
	
}

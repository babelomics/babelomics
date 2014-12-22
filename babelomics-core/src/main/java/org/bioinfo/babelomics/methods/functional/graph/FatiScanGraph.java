package org.bioinfo.babelomics.methods.functional.graph;

import java.awt.Color;
import java.awt.Font;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.methods.functional.GeneSetAnalysisTestResult;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.AnnotationFeature;
import org.bioinfo.graphics.canvas.panel.AnnotationPanel;
import org.bioinfo.graphics.canvas.track.AnnotationTrack;

public class FatiScanGraph {

	public static void fatiScanGraph(List<GeneSetAnalysisTestResult> significants, String title, String fileName){
		
		Canvas canvas = new Canvas(title);
		canvas.setBorderWidth(0);
		canvas.setTitleFont(new Font(Font.SANS_SERIF,Font.BOLD,12));
		
		AnnotationPanel ap = new AnnotationPanel("", 0, 0, 400, 100);
		
		ap.setPanelBorderWidth(0);
		
		// compute percentages
		List<Double> percentages = new ArrayList<Double>(significants.size());
		
		// get term order
		for(GeneSetAnalysisTestResult testResult: significants){
			percentages.add(getPercentage(testResult.getList1Positives(),testResult.getList1Negatives(),testResult.getList2Positives(), testResult.getList2Negatives()));
		}
		int[] termOrder = ListUtils.order(percentages);
						
		Color leftColor,rightColor;
		
		double percentage;
		GeneSetAnalysisTestResult testResult;
		int i;
		
		for(int si=0; si<termOrder.length; si++){
			i = termOrder[termOrder.length - si - 1];
			//i = termOrder.length - termOrder[si] - 1;	
			testResult = significants.get(i);
			percentage = percentages.get(i);
			
			AnnotationTrack track1 = new AnnotationTrack(testResult.getTerm(), "");
			
			// compute colors
			if(percentages.get(i)>=0.5) {
				leftColor = Color.red;
				rightColor = Color.gray;
			} else {
				leftColor = Color.gray;
				rightColor = Color.blue;
			}
			
			// left
			AnnotationFeature leftAnnotationFeature = new AnnotationFeature("feat1", "desc1", 0, (int)Math.round(400*percentage));
			leftAnnotationFeature.setBackgroundColor(leftColor);
			track1.add(leftAnnotationFeature);
			
			//right
			AnnotationFeature rightAnnotationFeature = new AnnotationFeature("feat1", "desc1", (int)Math.round(400*percentage), 400);
			rightAnnotationFeature.setBackgroundColor(rightColor);
			track1.add(rightAnnotationFeature);
			
			ap.addAnnotationTrack(track1);
		}
				
		ap.setHeight(significants.size()*30);
		
		canvas.addPanel(ap);
		canvas.render();
		
		try {			
			canvas.save(fileName);
		} catch (IOException e) {
			e.printStackTrace();	
		}
		
	}
	
	private static double getPercentage(int p1, int n1, int p2, int n2){		
		double per1 = (double)p1/(p1+n1);
		double per2 = (double)p2/(p2+n2);
		return per1/(per1+per2);
	}
}

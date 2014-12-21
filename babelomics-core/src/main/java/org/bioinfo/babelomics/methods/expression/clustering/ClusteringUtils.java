package org.bioinfo.babelomics.methods.expression.clustering;

import java.awt.Color;
import java.io.IOException;
import java.util.List;

import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.tree.multiway.MultipleTree;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.panel.NewickPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.math.data.DoubleMatrix;

public class ClusteringUtils {


	public static void saveImageTree(MultipleTree tree, String title, String imgFilename, boolean vertical, boolean createMap) throws IOException {

		int width = vertical ? (tree.getNumberOfLevels() * 10 + 30) : (tree.getNumberOfLeaves() * 10);
		int height = vertical ? (tree.getNumberOfLeaves() * 10) : (tree.getNumberOfLevels() * 10 + 30);

		// inferior bounds, not too small
		//
		if ( width < 256 ) width = 256;
		if (height < 256) height = 256;

		// superior bounds, not too big
		//
		if (width > 8000) width = 8000;
		if (height > 8000) height = 8000;

		NewickPanel newickPanel = new NewickPanel(title, 0, 0, width, height);
		newickPanel.setNewick(tree);
		newickPanel.setShowLabels(true);
		newickPanel.setVertical(vertical);

		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderPadding(0);
		canvas.setSpaceSeparator(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);
		canvas.setHeight(newickPanel.getHeight());
		canvas.setWidth(newickPanel.getWidth());

		canvas.addPanel(newickPanel);

		canvas.render(false);

		try {
			System.out.println("saving tree in : " + imgFilename);
			System.out.println("tree labels : " + ListUtils.toString(tree.getLabels(), ","));
			canvas.save(imgFilename, createMap);
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}



	public static void saveImageTree(DoubleMatrix matrix, MultipleTree vTree, MultipleTree hTree, List<String> columnLabels, List<String> rowLabels, String imgFilename, boolean createMap) throws IOException {

		int cellSide = 20;
		int rowLabelsWidth = getMaxStringLengh(rowLabels) * 9;
		int colLabelsWidth = getMaxStringLengh(columnLabels) * 9;
		int infoWidth = 0;

		int rowDimension = 0, columnDimension = 0;

		rowDimension = matrix.getRowDimension();
		columnDimension = matrix.getColumnDimension();

		//		System.out.println("sizes from trees: rows = " + rowDimension + ", cols = " + columnDimension);
		//		System.out.println("sizes from matrix: rows = " + matrix.getRowDimension() + ", cols = " + matrix.getColumnDimension());
		//
		//		System.out.println("width grid = " + (columnDimension * cellSide));
		//		System.out.println("height grid = " + (rowDimension * cellSide));
		//
		//		System.out.println("division (width / hTree.getNumberOfLeaves) = " + (1.0 * columnDimension * cellSide) / (1.0 * hTree.getNumberOfLeaves()));

		//		int hPanelX = rowLabelsWidth + (vTree.getNumberOfLevels() * cellSide);
		//		int hPanelY = 0;
		//		int hPanelWidth = columnDimension * cellSide + 10; 
		//		int hPanelHeight = hTree.getNumberOfLevels() * cellSide;
		//
		//		System.out.println("hpanel (x, y) = (" + hPanelX + ", " + hPanelY + ")");
		//		System.out.println("hpanel (width, height) = (" + hPanelWidth + ", " + hPanelHeight + ")");
		//		System.out.println("hpanel (levels, leaves) = (" + (hTree.getNumberOfLevels()) + ", " + (hTree.getNumberOfLeaves()) + ")");
		//		
		//
		//		NewickPanel newickHPanel = new NewickPanel("", hPanelX, hPanelY, hPanelWidth, hPanelHeight);

		int x, y, width, height;

		NewickPanel newickHPanel = null;
		if ( hTree != null ) {
			x = (vTree == null) ? rowLabelsWidth : rowLabelsWidth + (vTree.getNumberOfLevels() * cellSide);
			y = 0;
			width = columnDimension * cellSide + 10;
			height = hTree.getNumberOfLevels() * cellSide;

			newickHPanel = new NewickPanel("", x, y, width, height);

			newickHPanel.setNewick(hTree);
			newickHPanel.setLevelSeparation(cellSide);
			newickHPanel.setLeafSeparation(cellSide);
			newickHPanel.setShowLabels(false);
			newickHPanel.setVertical(false);
		}


		NewickPanel newickVPanel = null;
		if ( vTree != null ) {
			x = 0;
			y = (hTree == null) ? colLabelsWidth : colLabelsWidth + (hTree.getNumberOfLevels() * cellSide);
			width = vTree.getNumberOfLevels() * cellSide;
			height = rowDimension * cellSide; // before there was: (rowDimension * cellSide) + colLabelsWidth

			newickVPanel = new NewickPanel("", x, y, width, height);
			
			newickVPanel.setNewick(vTree);
			newickVPanel.setLevelSeparation(cellSide);
			newickVPanel.setLeafSeparation(cellSide);
			newickVPanel.setShowLabels(false);
			newickVPanel.setVertical(true);
		}

		GridPanel gridPanel = new GridPanel("", 
											(vTree == null) ? 0 : vTree.getNumberOfLevels() * cellSide, 
											(hTree == null) ? 0 : hTree.getNumberOfLevels() * cellSide, 
											(columnDimension * cellSide) + rowLabelsWidth + infoWidth, 
											(rowDimension * cellSide) + colLabelsWidth);

		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
		gridTrack.setName("");
		gridTrack.setColumnLabels(columnLabels);
		gridTrack.setRowLabels(rowLabels);
		gridTrack.setTopRegion(colLabelsWidth);
		gridTrack.setLeftRegion(rowLabelsWidth);
		gridTrack.setRightRegion(infoWidth);
		ScoreFeature feature;

		double mean, deviation, min, max, offset, standard;
		double[] values = new double[gridTrack.getColumnDimension()];
		for(int row=0 ; row<gridTrack.getRowDimension() ; row++) {
			mean = matrix.getRowMean(row);
			deviation = matrix.getRowStdDeviation(row);
			min = Double.MAX_VALUE;
			max = Double.MIN_NORMAL;
			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
				values[column] = (deviation == 0) ? Double.NaN : (matrix.get(row, column)-mean)/(deviation);
				if ( min > values[column] ) min = values[column];
				if ( max < values[column] ) max = values[column];
			}
			offset = ( min <= 0 ) ? Math.abs(min) : (-1 * min);
			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
				standard = (values[column] + offset) / ( max + offset);

				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, standard);
				gridTrack.setFeature(row, column, feature);
			}
		}
		gridPanel.add(gridTrack);

		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderPadding(4);
		canvas.setSpaceSeparator(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);

		int canvasHeight = gridPanel.getHeight() + (newickHPanel == null ? 0 : newickHPanel.getHeight()) + (cellSide * 2);
		int canvasWidth = gridPanel.getWidth() + (newickVPanel == null ? 0 : newickVPanel.getWidth()) + (cellSide * 2);
		System.out.println("canvas height = " + canvasHeight);
		System.out.println("canvas width  = " + canvasWidth);

		canvas.setHeight(canvasHeight);
		canvas.setWidth(canvasWidth);

		canvas.addPanel(gridPanel);
		if ( newickHPanel != null ) { canvas.addPanel(newickHPanel); }
		if ( newickVPanel != null ) { canvas.addPanel(newickVPanel); }

		canvas.render();
		canvas.save(imgFilename, createMap);		
	}


//	public static void saveImageTreeOK(DoubleMatrix matrix, MultipleTree vTree, MultipleTree hTree, String imgFilename, boolean createMap) throws IOException {
//
//		int cellSide = 20;
//		int rowLabelsWidth = getMaxStringLengh(vTree.getLabels()) * 9;
//		int colLabelsWidth = getMaxStringLengh(hTree.getLabels()) * 9;
//		int infoWidth = 0;
//
//
//		int rowDimension = vTree.getNumberOfLeaves();
//		int columnDimension = hTree.getNumberOfLeaves();
//
//		System.out.println("sizes from trees: rows = " + rowDimension + ", cols = " + columnDimension);
//		System.out.println("sizes from matrix: rows = " + matrix.getRowDimension() + ", cols = " + matrix.getColumnDimension());
//
//		System.out.println("width grid = " + (columnDimension * cellSide));
//		System.out.println("height grid = " + (rowDimension * cellSide));
//
//		System.out.println("division (width / hTree.getNumberOfLeaves) = " + (1.0 * columnDimension * cellSide) / (1.0 * hTree.getNumberOfLeaves()));
//
//		//		NewickPanel newickHPanel = new NewickPanel("", hTree.getNumberOfLeaves() * cellSide, 
//		//													   hTree.getNumberOfLevels() * cellSide, 
//		//													   rowLabelsWidth + (vTree.getNumberOfLevels() * cellSide), 
//		//													   0);
//
//		int hPanelX = rowLabelsWidth + (vTree.getNumberOfLevels() * cellSide);
//		int hPanelY = 0;
//		int hPanelWidth = columnDimension * cellSide + 10; 
//		int hPanelHeight = hTree.getNumberOfLevels() * cellSide;
//
//		System.out.println("hpanel (x, y) = (" + hPanelX + ", " + hPanelY + ")");
//		System.out.println("hpanel (width, height) = (" + hPanelWidth + ", " + hPanelHeight + ")");
//		System.out.println("hpanel (levels, leaves) = (" + (hTree.getNumberOfLevels()) + ", " + (hTree.getNumberOfLeaves()) + ")");
//
//
//		NewickPanel newickHPanel = new NewickPanel("", hPanelX, hPanelY, hPanelWidth, hPanelHeight);
//
//		newickHPanel.setNewick(hTree);
////		newickHPanel.setLevelSeparation(cellSide);
////		newickHPanel.setLeafSeparation(cellSide);
//		newickHPanel.setShowLabels(false);
//		newickHPanel.setVertical(false);
//
//
//		NewickPanel newickVPanel = new NewickPanel("", 0, colLabelsWidth + (hTree.getNumberOfLevels() * cellSide),
//				vTree.getNumberOfLevels() * cellSide, 
//				(rowDimension * cellSide) + colLabelsWidth); //vTree.getNumberOfLeaves() * cellSide);
//
//		newickVPanel.setNewick(vTree);
////		newickVPanel.setLevelSeparation(cellSide);
////		newickVPanel.setLeafSeparation(cellSide);
//		newickVPanel.setShowLabels(false);
//		newickVPanel.setVertical(true);
//
//		//		GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + infoWidth, 
//		//												(columnDimension * cellSide) + colLabelsWidth, 
//		//												vTree.getNumberOfLevels() * cellSide, 
//		//												newickHPanel.getHeight());
//
//		GridPanel gridPanel = new GridPanel("", 
//				vTree.getNumberOfLevels() * cellSide, 
//				hTree.getNumberOfLevels() * cellSide, 
//				(columnDimension * cellSide) + rowLabelsWidth + infoWidth, 
//				(rowDimension * cellSide) + colLabelsWidth);
//
//		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
//		gridTrack.setName("");
//		gridTrack.setColumnLabels(hTree.getLabels());
//		gridTrack.setRowLabels(vTree.getLabels());
//		gridTrack.setTopRegion(colLabelsWidth);
//		gridTrack.setLeftRegion(rowLabelsWidth);
//		gridTrack.setRightRegion(infoWidth);
//		ScoreFeature feature;
//
//		double mean, deviation, min, max, offset, standard;
//		double[] values = new double[gridTrack.getColumnDimension()];
//		for(int row=0 ; row<gridTrack.getRowDimension() ; row++) {
//			mean = matrix.getRowMean(row);
//			deviation = matrix.getRowStdDeviation(row);
//			min = Double.MAX_VALUE;
//			max = Double.MIN_NORMAL;
//			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
//				values[column] = (deviation == 0) ? Double.NaN : (matrix.get(row, column)-mean)/(deviation);
//				if ( min > values[column] ) min = values[column];
//				if ( max < values[column] ) max = values[column];
//			}
//			offset = ( min <= 0 ) ? Math.abs(min) : (-1 * min);
//			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
//				standard = (values[column] + offset) / ( max + offset);
//
//				//				System.out.print(matrix.get(row, column) + "\t");
//				//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, matrix.get(row, column));
//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, standard);
//				gridTrack.setFeature(row, column, feature);
//			}
//			//			System.out.println("");
//		}
//		gridPanel.add(gridTrack);
//
//		Canvas canvas = new Canvas("");
//		canvas.setBorderWidth(0);
//		canvas.setBorderPadding(4);
//		canvas.setSpaceSeparator(0);
//		canvas.setBorderColor(Color.BLACK);
//		canvas.setBackGroundColor(Color.WHITE);
//
//		int canvasHeight = gridPanel.getHeight() + newickHPanel.getHeight() + (cellSide * 2);
//		int canvasWidth = gridPanel.getWidth() + newickVPanel.getWidth() + (cellSide * 2);
//		System.out.println("canvas height = " + canvasHeight);
//		System.out.println("canvas width  = " + canvasWidth);
//
//		canvas.setHeight(canvasHeight);
//		canvas.setWidth(canvasWidth);
//
//		canvas.addPanel(gridPanel);
//		canvas.addPanel(newickHPanel);
//		canvas.addPanel(newickVPanel);
//
//		canvas.render();
//		canvas.save(imgFilename, createMap);		
//	}

	private static int getMaxStringLengh(List<String> names) {
		int max = 0;
		for(String name: names) {
			if ( name.length() > max ) {
				max = name.length();
			}
		}
		return max;
	}

}

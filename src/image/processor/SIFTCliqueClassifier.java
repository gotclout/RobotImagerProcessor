package image.processor;

//package edu.ksu.cis.robotics.vision.entry.core;

//import imagehunter.connection.HunterSender;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.imageio.ImageIO;

//import edu.ksu.cis.robotics.core.VSVectorFloat;
import edu.ksu.cis.robotics.vision.SIFTDriver;
import edu.ksu.cis.robotics.vision.SIFTDriverSpec;
import edu.ksu.cis.robotics.vision.util.ImageConversions;

import mpi.cbg.fly.Feature;
import mpi.cbg.fly.Filter;
import mpi.cbg.fly.FloatArray2D;
import mpi.cbg.fly.FloatArray2DScaleOctave;

public class SIFTCliqueClassifier {
	
	private static int steps = 3;
	// initial sigma
	private static float initial_sigma = 1.6f;
	// feature descriptor size
	private static int fdsize = 4;
	// feature descriptor orientation bins
	private static int fdbins = 8;
	// size restrictions for scale octaves, use octaves < max_size and > min_size only
	private static int min_size = 64;
	private static int max_size = 1024;
	
	
	Map<String, List<Vector<Feature>>> classFeatureMap;
	
	public SIFTCliqueClassifier(Map<String, List<Vector<Feature>>> classFeatureMap){
		this.classFeatureMap = classFeatureMap;
	}
	
	public BufferedImage annotateImage(BufferedImage img){
		Map<String, Vector<Feature>> annotMap = new HashMap<String, Vector<Feature>>();
		
		int[][] image;
		boolean upscale = img.getHeight() < 1000 || img.getWidth()< 1000;
		float[] img_data_greyscale = ImageConversions.convertRGBBuffIToFloatArr(img);
		FloatArray2D fa = new FloatArray2D(img_data_greyscale, img.getWidth(), img.getHeight());
		Filter.enhance(fa, 1.0f);
		if(upscale){
			FloatArray2D fat = new FloatArray2D( fa.width * 2 - 1, fa.height * 2 - 1 ); 
			FloatArray2DScaleOctave.upsample( fa, fat );
			fa = fat;
			fa = Filter.computeGaussianFastMirror( fa, ( float )Math.sqrt( initial_sigma * initial_sigma - 1.0 ) );
		}
		else{
			fa = Filter.computeGaussianFastMirror( fa, ( float )Math.sqrt( initial_sigma * initial_sigma - 0.25 ) );
		}
		Vector<Feature> features = SIFTDriverSpec.getFeaturesFromImage(img);
		int width, height;
		width = fa.width;
		height = fa.height;
		
		image = new int[height][width];
		
	    for (int y = 0; y < image.length; y++)
	    	for (int x = 0; x < image[0].length; x++)
	    		image[y][x] = (int)(fa.get(x, y) * 255) * (256*256 + 256 + 1);
			
		int rgbOutput[] = new int[width * height];
		
		
		for(Feature f : features){
			System.out.println("Drawing features");
			drawCircleWithFeaturePoint(f, image);
		}
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)
				rgbOutput[x+width*y] = image[y][x];
		
		BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		bi.setRGB(0, 0, width, height, rgbOutput, 0, width);
		Graphics g = bi.getGraphics();
		g.setColor(Color.GREEN);
		Set<String> keySet = classFeatureMap.keySet();
		for(String s : keySet){
			List<Vector<Feature>> exemplars = classFeatureMap.get(s);
			boolean matched = false;
			for(Vector<Feature> vect : exemplars){
				List<Feature> matchList = new ArrayList<Feature>();	
				for(Feature fExemp : vect){
					boolean subed = false;
					Feature match = checkForMatch(fExemp, features);
					if(match != null) matchList.add(match);
					//for(Feature fImg : features){
						/*if(match(fExemp, fImg) && !subed){
							subed = true;
							matchList.add(fImg);
						}*/
					//}
				}
				float numer = (float)matchList.size();
				float denom = (float)vect.size();
				float matchLev = numer / denom;
				System.out.println(matchLev);
				//if(matchLev > 0.002f && !matched) {
				System.out.println("Match list size:" + matchList.size());
				if(matchList.size() > 1){
					drawBox(bi, s, matchList);
					matched = true;
				//}
				}
			}
			
		}
		return bi;
	}
	
	public double dist(Feature f1, Feature f2){
		double preSqrt = 0.0;
		for(int i = 0; i < 128; i++){
			preSqrt += Math.pow((double)(f1.descriptor[i] - f2.descriptor[i]), 2.0);
		}
		return Math.sqrt(preSqrt);
	}
	
	public double distSqrd(Feature f1, Feature f2){
		double preSqrt = 0.0;
		for(int i = 0; i < 128; i++){
			preSqrt += Math.pow((double)(f1.descriptor[i] - f2.descriptor[i]), 2.0);
		}
		return preSqrt;
	}
	
	public boolean match(Feature f1, Feature f2){
		double dist = dist(f1, f2);
		return(dist < 0.255);
	}
	
	public Feature checkForMatch(Feature feature, Vector<Feature> features){
		double  dsq, distsq1 = 100000000d, distsq2 = 100000000d;
		Feature minKey = null;
		for(Feature f1 : features){
			dsq = distSqrd(feature, f1);
			 if (dsq < distsq1) {
					distsq2 = distsq1;
					distsq1 = dsq;
					minKey = f1;
			 }
			 else if (dsq < distsq2) {
					distsq2 = dsq;
			}
		}
		if (10.0d * 10.0d * distsq1 < 6.0d * 4.0d * distsq2)
		      return minKey;
		    else return null;
	}
	
	public void drawBox(BufferedImage img, String className, List<Feature> features){
		Graphics g = img.getGraphics();
		g.setColor(Color.GREEN);
		float minx = (float)img.getWidth();
		float maxx = 0.0f;
		float miny = (float)img.getHeight();
		float maxy = 0.0f;
		//refine the maxes and mins
		for(Feature f : features){
			if(f.location[0] < minx) minx = f.location[0];
			if(f.location[0] > maxx) maxx = f.location[0];
			if(f.location[1] < miny) miny = f.location[1];
			if(f.location[1] > maxy) maxy = f.location[1];
		}
		if(features.size() == 1){
			minx = features.get(0).location[0] - 300;
			maxx = features.get(0).location[0] + 300;
			miny = features.get(0).location[1] - 300;
			maxy = features.get(0).location[1] + 300;
		}
		g.drawRect((int)minx, (int)miny, (int)(maxx - minx), (int)(maxy - miny));
		System.out.println("X: " + minx + " Y: " + miny + " width: " + (maxx - minx) + " height: " + (maxy - miny));
		g.drawString(className, (int)(minx - 4.0f), (int)(miny - 4.0f));
	}
	
	
	
	private void drawCircleWithFeaturePoint(Feature f, int[][] image) {
		try{
			int x, y;
			
			double cos = Math.cos(f.orientation);
			double sin = Math.sin(f.orientation);
			for (int r = 0; r <= f.scale; r++)
			{
				x = (int)(cos * r) + (int)f.location[0];
				y = (int)(sin * r) + (int)f.location[1];
				
				image[y][x] = 255 * 256 * 256;
			}
			
			for (float theta = 0; theta < 2*Math.PI; theta += Math.PI / 60)
			{
				cos = Math.cos(theta);
				sin = Math.sin(theta);
				x = (int)(cos * f.scale) + (int)f.location[0];
				y = (int)(sin * f.scale) + (int)f.location[1];
				
				image[y][x] = 255 * 256 * 256;
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void main(String args[]){
		//File pe1f = new File("panda_express_logo/www_cozbaldwin_com_images_panda_express_logo_jpg.png");
		//File pe2f = new File("panda_express_logo/pecup2.jpg");
		/*HunterSender sender = new HunterSender("192.168.7.94");
		File pe1f = new File("norvig.jpg");
		List<Vector<Feature>> featureList = new ArrayList<Vector<Feature>>();
		Map<String, List<Vector<Feature>>> classes = new HashMap<String, List<Vector<Feature>>>();
		try{
			featureList.add(SIFTDriverSpec.getFeaturesFromImage(ImageIO.read(pe1f)));
			//featureList.add(SIFTDriverSpec.getFeaturesFromImage(ImageIO.read(pe2f)));
			classes.put("Artificial Intelligence a Modern Approach", featureList);
			SIFTCliqueClassifier scc = new SIFTCliqueClassifier(classes);
			BufferedImage img = ImageIO.read(new File("book_test.png"));
			BufferedImage annotated = scc.annotateImage(img);
			img = null;
			ImageIO.write(annotated, "png", new File("scc-testout.png"));
			sender.sendImageToServer(0, 0, 2.0f, "test1", annotated);
			sender.sendImageToServer(1, -2, 0.0f, "test2", annotated);
			sender.sendImageToServer(1, 2, 0.0f, "test2", annotated);
			sender.sendImageToServer(0, -3, 0.0f, "test3", annotated);
			sender.sendImageToServer(0, -6, 0.0f, "test4", annotated);
			sender.sendImageToServer(0, -8, 0.0f, "test5", annotated);
			sender.sendImageToServer(1, -12, 0.0f, "test6", annotated);
			sender.sendImageToServer(2, -15, 0.0f, "test7", annotated);
			sender.sendImageToServer(1, -20, 0.0f, "test8", annotated);
			
		}
		catch(Exception e){
			e.printStackTrace();
		}*/
	}

}

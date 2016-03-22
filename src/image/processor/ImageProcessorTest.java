package image.processor;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;

import javax.imageio.ImageIO;

import mpi.cbg.fly.Feature;
import mpi.cbg.fly.Filter;
import mpi.cbg.fly.FloatArray2D;
import mpi.cbg.fly.FloatArray2DSIFT;
import mpi.cbg.fly.FloatArray2DScaleOctave;
import edu.ksu.cis.robotics.vision.util.ImageConversions;

public class ImageProcessorTest {
//	 steps
	private static int steps = 3; //3 o
	// initial sigma
	private static float initial_sigma = 1.6f;
	// feature descriptor size
	private static int fdsize = 4;
	// feature descriptor orientation bins
	private static int fdbins = 8;
	// size restrictions for scale octaves, use octaves < max_size and > min_size only
	private static int min_size = 64; //32 o 64?
	private static int max_size = 512; //384 o 512? 1024?
	
	private static int[][] image;
	private ArrayList<Integer> kMeansValues;
	public static ImageProcessor ip = new ImageProcessor(/*"",""*/);
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//n number of files followed by n file names
		
		int n , i , j , l, k;
		n = i = j = k = l = 0;
		String fileName  = "";
		boolean useMeanPoint = false;
			
		n = Integer.valueOf( args[0] ).intValue();
		
		for(k = 0; k < 1; k++){
			
		if(k == 1)
			{
				useMeanPoint = true;
			}
			/*if(k == 0)
			{
				//siftOnly
				for(i = 0; i < n; i++)
				{
					fileName = "";
					
					l = args[i+1].length();
					
					for(j = 0; j < l; j++)
					{
						if(args[i+1].charAt(j) == '.')
							break;
						
						fileName += args[i+1].charAt(j);
		
					}
					
					fileName += "_siftOnly.jpg";
					myTest(args[i+1], fileName, false, false, true, useMeanPoint);
				}
			}*/
			
			//clustering without outlier detection
			/*for(i = 0; i < n; i++)
			{
				fileName = "";
				
				l = args[i+1].length();
				
				for(j = 0; j < l; j++)
				{	
					if(args[i+1].charAt(j) == '.')
						break;
					
					fileName += args[i+1].charAt(j);
				}

				if(useMeanPoint)
				{
					fileName += "_without_outlier.jpg";
				}
				else
				{
					fileName += "_without_outlier_actual_mean.jpg";
				}
				
				myTest(args[i+1], fileName, false, false, false, useMeanPoint);
			}*/
			
			/*//clustering with nearest neighbor outlier detection
			for(i = 0; i < n; i++)
			{
				fileName = "";
				
				l = args[i+1].length();
				
				for(j = 0; j < l; j++)
				{
					if(args[i+1].charAt(j) == '.')
						break;
					
					fileName += args[i+1].charAt(j);
				}
				
				if(useMeanPoint)
				{
					fileName += "_with_nearest_neighbor.jpg";
				}
				else
				{
					fileName += "_with_nearest_neighbor_actual_mean.jpg";
				}
				myTest(args[i+1], fileName, true, false, false, useMeanPoint);
			}
			
			*/
			//clustering with top percentage outlier detection
			for(i = 0; i < n; i++)
			{
				fileName = "";
				
				l = args[i+1].length();
				
				for(j = 0; j < l; j++)
				{
					if(args[i+1].charAt(j) == '.')
						break;
					
					fileName += args[i+1].charAt(j);
				}

				if(useMeanPoint)
				{
					fileName += "_with_top_percentage.jpg";
				}
				else
				{
					fileName += "_with_top_percentage_actual_mean.jpg";
				}
				
				myTest(args[i+1], fileName, false, true, false, useMeanPoint);
			}
			
			//Tracking without outlier detection
			//trackObjects(args, false, false, useMeanPoint, false);
			//Tracking with nearest neighbor detection
		//	trackObjects(args, true, false, useMeanPoint, false);
			//Tracking with top % outlier detection
			//trackObjects(args, false, true, useMeanPoint, false);
			
			
			//use these to test the same image against it'self
			//it's more of a sanity check
			//Tracking without outlier detection
			//trackObjects(args, false, false, useMeanPoint, true);
			//Tracking with nearest neighbor detection
			//trackObjects(args, true, false, useMeanPoint, true);
			//Tracking with top % outlier detection
			//trackObjects(args, false, true, useMeanPoint, true);
			
		}

	}
	
	public static void trackObjects(String[] files, boolean useNearestNeighbor, boolean useTopPercentage, 
			boolean useMeanPoint, boolean testSame)
	{
		FloatArray2D fa;
		if(useNearestNeighbor)
		{
			System.out.println("Tracking images with nearest neighbor outlier detection\n");
		}
		else if(useTopPercentage)
		{
			System.out.println("Tracking images with top % outlier detection\n");
		}
		else
		{
			System.out.println("Tracking images without outlier detection\n");
		}
		
		Vector< Feature > feats = null; 
		int numFiles = Integer.valueOf( files[0] ).intValue();
		
		//ImageProcessor ip = new ImageProcessor(/*"",""*/);
		for(int i = 1; i < numFiles + 1; i++)
		{
			fa = getFloatArray2D(files[i]);
			//feats = getFeaturesFromFile(files[i]);
			feats = getFeaturesFromFile(fa);
			if(useNearestNeighbor)
			{
				feats = ip.detectOutliers(feats, 0);
			}
			else if(useTopPercentage)
			{
				feats = ip.detectOutliers(feats, 1);
			}
			
			BufferedImage img = getImageFromFile(files[i]);
			System.out.println("Beginning Image Processor Tracking Test For Image File " + 
						files[i] + "\n");
			ip.processPhoto(img, feats, useMeanPoint); 
			if(testSame)
			{
				ip.processPhoto(img, feats, useMeanPoint, testSame);
			}
			String outfile = files[i].substring(0, files[i].length() - 4) +
			"_tracking_" + Integer.toString(i) + ".jpg";
			
			Vector<Feature> feats2;
			feats2 = new Vector<Feature>(); 
			int k = ip.getKMeans().getK();
			for(int j = 0; j < k; j++)
			{
				//Feature f1 = ip.getKMeans().getClusterMeanPoint(j);
				//Feature f2 = ip.getKMeans().getClusterMeanFeat(j);
				feats2.add(ip.getKMeans().getClusterMeanFeat(j));
			}
			displayFeatures(fa, feats2, outfile);
		}
	}
	
	public static BufferedImage getImageFromFile(String infile)
	{
		BufferedImage img = null;
		try {
			File picture = new File(infile);
			img = ImageIO.read(picture);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return img;
	}
	
	public static FloatArray2D getFloatArray2D(String infile)
	{
		BufferedImage img = getImageFromFile(infile);
		
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
		
		return fa;
	}
	
	public static Vector< Feature > getFeaturesFromFile(FloatArray2D fa)
	{
		//FloatArray2D fa = getFloatArray2D(infile);
		Vector< Feature > fs1;
		FloatArray2DSIFT sift = new FloatArray2DSIFT( fdsize, fdbins );
		
		System.out.print( "processing SIFT ..." );
		sift.init( fa, steps, initial_sigma, min_size, max_size );
		fs1 = sift.run( max_size );
		Collections.sort(fs1);
		return fs1;
	}
	
	public static Vector< Feature > getFeaturesFromFile(String infile)
	{
		FloatArray2D fa = getFloatArray2D(infile);
		Vector< Feature > fs1;
		FloatArray2DSIFT sift = new FloatArray2DSIFT( fdsize, fdbins );
		
		System.out.print( "processing SIFT ..." );
		sift.init( fa, steps, initial_sigma, min_size, max_size );
		fs1 = sift.run( max_size );
		Collections.sort(fs1);
		return fs1;
	}
	
	public static void myTest(String infile, String outfile, 
			boolean useNearestNeighbor, boolean useTopPercentage, 
			boolean siftOnly, boolean useMeanPoint)
	{
		BufferedImage img = getImageFromFile(infile);
		
		boolean upscale = img.getHeight() < 1000 || img.getWidth()< 1000;
		float[] img_data_greyscale = ImageConversions.convertRGBBuffIToFloatArr(img);
		FloatArray2D fa = new FloatArray2D(img_data_greyscale, img.getWidth(), img.getHeight());
		Filter.enhance(fa, 1.0f);
		Vector< Feature > fs1;
		FloatArray2DSIFT sift = new FloatArray2DSIFT( fdsize, fdbins );
		if(upscale){
			FloatArray2D fat = new FloatArray2D( fa.width * 2 - 1, fa.height * 2 - 1 ); 
			FloatArray2DScaleOctave.upsample( fa, fat );
			fa = fat;
			fa = Filter.computeGaussianFastMirror( fa, ( float )Math.sqrt( initial_sigma * initial_sigma - 1.0 ) );
		}
		else{
			fa = Filter.computeGaussianFastMirror( fa, ( float )Math.sqrt( initial_sigma * initial_sigma - 0.25 ) );
		}
		System.out.print( "processing SIFT ..." );
		sift.init( fa, steps, initial_sigma, min_size, max_size );
		fs1 = sift.run( max_size );
		Collections.sort(fs1);
		Vector<Feature> feats2;
		
		if(!siftOnly)
		{	
			if(useNearestNeighbor)
			{
				fs1 = ip.detectOutliers(fs1, 0);
			}
			else if(useTopPercentage)
			{
				fs1 = ip.detectOutliers(fs1, 1);
			}
			
			if(useNearestNeighbor)
			{
				System.out.println("Tracking images with nearest neighbor outlier detection\n");
			}
			else if(useTopPercentage)
			{
				System.out.println("Tracking images with top % outlier detection\n");
			}
			else
			{
				System.out.println("Tracking images without outlier detection\n");
			}
			System.out.println("The current file is " + infile + "\n");
			ip.processPhoto(img, fs1, useMeanPoint); 
			feats2 = new Vector<Feature>(); 
			int k = ip.getKMeans().getK();
			for(int i = 0; i < k; i++)
			{
				feats2.add(ip.getKMeans().getClusterMeanFeat(i));
			}
		}
		else
			feats2 = fs1;
		
		Collections.sort( feats2 );
		for(Feature f : feats2){
			System.out.println("Feature at: X: " + f.location[0] + " Y: " + f.location[1] + ", Scale: " + f.scale);
		}
		System.out.println("Num Features: " + feats2.size());
		displayFeatures(fa, feats2, outfile);
	}
	
	private static void displayFeatures(FloatArray2D fa, Vector<Feature> fs1, String outputFile) {
		int width, height;
		width = fa.width;
		height = fa.height;
		
		image = new int[height][width];
	    for (int y = 0; y < image.length; y++)
	    	for (int x = 0; x < image[0].length; x++)
	    		image[y][x] = (int)(fa.get(x, y) * 255) * (256*256 + 256 + 1);
			
		int rgbOutput[] = new int[width * height];
		
		for(Feature f : fs1){
			drawCircleWithFeaturePoint(f);
		}
		
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)
				rgbOutput[x+width*y] = image[y][x];
		
		BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		bi.setRGB(0, 0, width, height, rgbOutput, 0, width);
		
		File f = new File (outputFile);
		String parentDir = f.getParent();
		
		if (parentDir != null && new File(parentDir).mkdirs())
			System.out.println("Created directory " + parentDir);
		
		try
		{
			ImageIO.write (bi, "bmp", f);
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	private static void drawCircleWithFeaturePoint(Feature f) {
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
			
			if(x < image[0].length && y < image.length && y >= 0 && x >= 0)
				image[y][x] = 255 * 256 * 256;
		}
	}
}

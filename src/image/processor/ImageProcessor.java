/******************************
 * Author: Robert L. Foster Jr.
 * Processes objects in images from the robot camera into
 * clusters.
 ******************************/
package image.processor;

//import edu.ksu.cis.robotics.vision.*;
//import edu.ksu.cis.robotics.vision.test.SIFTTest;
import robotserver.camera.CameraDriver;
import imagehunter.connection.HunterSender;
import javaclient.PlannerInterface;
import javaclient.PlayerClient;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Vector;
import java.util.List;
import java.util.ArrayList;
import javax.imageio.ImageIO;
import com.sun.tools.javac.util.Pair;
import edu.ksu.cis.robotics.vision.SIFTDriver;
import mpi.cbg.fly.Feature;

public class ImageProcessor{
	
	public int imageCount = 12;
	CameraDriver	camera;
	kMeans			km = new kMeans();
	ArrayList<ClusterData> cDataList;
	final long processTime = 6000;
	SIFTCliqueClassifier classifier;
	RobotConnection con;
	BufferedImage	image = null;
	final String changeToRobotDir     = "cd /Users/rfoster/";
	final String deleteAllImgsFromCam = "gphoto2 -f /store_00010001/DCIM/100CANON -D";
	final String getImageFromCam      = "gphoto2 -f /store_00010001/DCIM/100CANON --get-file 1";
	final String takePicture          = "gphoto2 --capture-image";
	final String picOnDisk			  = "/Users/rfoster/IMG_0001.JPG";
	//HunterSender sender; currently using player joy to send remove to change
	//should probably use this to get the map and update it
	
	public BufferedImage gphoto2TakePhoto()
	{
		BufferedImage img = null;
		File picture = new File("/Users/rfoster/RobotPhotos/IMG_0001.JPG");
		//working directory
		File wd = new File("/Users/rfoster");
		
		//delete old from disk
		picture.delete();
		
		Process proc;
		try{
			//delete the old one from cam
			proc = Runtime.getRuntime().exec(deleteAllImgsFromCam, null, wd);
			Thread.sleep(300);
			//take a new one
			proc = Runtime.getRuntime().exec(takePicture, null, wd);
			Thread.sleep(300);
			//save to disk
			proc = Runtime.getRuntime().exec(getImageFromCam, null, wd);
			Thread.sleep(1000);
		}
		catch (Exception e){
			e.printStackTrace();
		}
		
		picture = new File(picOnDisk);
		
		try {
			img = ImageIO.read(picture);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return img;
	}
	
	//indicates the direction of travel in a moving object
	enum Direction{
		NORTH, 
		SOUTH, 
		EAST, 
		WEST, 
		NORTHEAST, 
		SOUTHEAST, 
		NORTHWEST, 
		SOUTHWEST, 
		STATIONARY
	}
	
	public class ClusterData
	{
		cluster aCluster;
		boolean inMotion;
		//number of points matching in this cluster
		public int matches;
		int id;
		Direction direction;
		//list of possible matches for caching and removing
		ArrayList<cluster> matchList = new ArrayList<cluster>();
		ClusterData()
		{
			;
		}
		ClusterData(cluster C, boolean motion, int s)
		{
			this.aCluster = C;
			this.inMotion = motion;
			id = s;
			matches = 0;
		}
		
		Direction getDirection()
		{
			return direction;
		}
		
		void setDirection(Direction d)
		{
			this.direction = d;
		}
		
		cluster getCluster()
		{
			return aCluster;
		}
		
		boolean getMotion()
		{
			return inMotion;
		}
		
		void setMotion(boolean motion)
		{
			this.inMotion = motion;
		}
	}

	public ImageProcessor(/*String conAdd, String sendAdd*/)
	{
		cDataList = new ArrayList<ClusterData>();
		classifier = new SIFTCliqueClassifier(null);
		//con = new RobotConnection(conAdd);
		//sender = new HunterSender(sendAdd);
	}
	
	public kMeans getKMeans()
	{
		return this.km;
	}
	
	public boolean updateMap()
	{
		boolean retVal = true;
		
		return retVal;
	}
	public boolean takePhotoAndProcess()
	{
		boolean			retVal = true;
		Vector<Feature>	features;
		int 			numClusters = 0;
		cluster 		currentCluster = null;
		ClusterData		cData;
		
		//image = camera.getPicture(); //obsolete?
		imageCount++;
		image = con.getImg();
		//save it to disk for future analysis
		try {
			ImageIO.write(image, "jpeg", new File("/Users/rfoster/RobotPhotos/IMG_" 
					+ imageCount + ".jpg" ));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//image = gphoto2TakePhoto();
		if(image != null)
		{
			features = SIFTDriver.getFeaturesFromImage(image);
			features = detectOutliers(features, 0);
			List<kMeansPoint> points = this.generatePointsList(features);
			numClusters = points.size();
			km = new kMeans(((int)Math.sqrt(features.size()/2)), points);
			km.runKMeans();
			boolean firstSet = cDataList.isEmpty();
			ArrayList<ClusterData> newGuys = new ArrayList<ClusterData>();
			for(int i = 0; i < km.getK(); i++)
			{
				currentCluster = km.getCluster(i);
				if(currentCluster != null)
				{
					if(!firstSet)
					{
						ClusterData match = trackClusters(currentCluster, false);
						if(match == null)
						{
							cData = new ClusterData(currentCluster, false, cDataList.size());
							newGuys.add(cData); 
						}
						else
						{
							cData = new ClusterData(currentCluster, false, cDataList.size());
							//if new cluster contains more points than old cluster 
							//add it and remove the new one? maybe cache it
							if(cData.getCluster().numPoints > 
							   match.getCluster().numPoints)
							{
								newGuys.add(cData);
								cDataList.remove(match);
							}
							setDirection(match, cData, newGuys.size());
						}
					}
					else //add them all the first time, there is nothing to map to
					{
						cData = new ClusterData(currentCluster, false, cDataList.size());
						cDataList.add(cData);
					}
					//System.out.println("Cluster " + i + " mean is " + currentCluster.getMean() + '\n' );
				}
			}
			if(!newGuys.isEmpty())
				cDataList.addAll(newGuys);		
		}
		else
		{
			retVal = false;
		}
		
		return retVal;
			
	} 
	
	public ClusterData mapToExistingCluster2(cluster current, boolean useMeanPoint, boolean printTables)
	{
		Vector<Feature> feats = new Vector<Feature>();
		int i;
		for(i = 0; i < cDataList.size(); i++)
		{
			feats.add(cDataList.get(i).getCluster().getActualMean().getFeature());
		}
		Feature match = classifier.checkForMatch(current.getActualMean().getFeature(), feats);
		if(match != null)
		{
			Feature f = null;
			for(i = 0; i < cDataList.size(); i++)
			{
				f = cDataList.get(i).getCluster().getActualMean().getFeature();
			
				if(classifier.match(match, f))
				{
					return cDataList.get(i);
				}
			}
		}
		return null;
	}
	public ClusterData mapToExistingCluster(cluster current, boolean useMeanPoint, boolean printTables)
	{
		ClusterData theMatch = null;
		int matches = 0;
		int max = 0;
		kMeansPoint currentMean = null;
		int numClusters;
		ArrayList < Pair<ClusterData, Integer> > clusterMap = null;
		if(current != null)
		{
			numClusters = cDataList.size();
			clusterMap = new ArrayList< Pair<ClusterData, Integer> >();
			for(int i = 0; i < numClusters; i++)
			{
				//max = 0;
				if(useMeanPoint)
				{
					currentMean = current.getMean();
					if(kMeansPointEqual(currentMean, cDataList.get(i).getCluster().getMean()))
					{
						theMatch = cDataList.get(i);
						theMatch.matches = max = matches;
					}
				}
				else
				{
					Vector<kMeansPoint> oldPoints = cDataList.get(i).getCluster().getClusterPoint();
					Vector<kMeansPoint> newPoints = current.getClusterPoint();
					kMeansPoint last = null;
					kMeansPoint prev = null;
					int mp = 0;
					matches = 0;
					for(int j = 0; j < oldPoints.size(); j++)
					{
						for(int k = 0; k < newPoints.size(); k++)
						{
							if(kMeansPointEqual(oldPoints.get(j), newPoints.get(k)))
							{
								mp++;
								prev = last;
								last = newPoints.get(k);
								matches++;
							}
						}
						mp = 0;
					}
					if(matches > max)
					{
						if(theMatch != null)
						{
							float oldRatio, newRatio;
							oldRatio = theMatch.getCluster().numPoints/(float)theMatch.matches;
							newRatio = cDataList.get(i).getCluster().numPoints/(float)matches;
							if((theMatch.getCluster().numPoints == cDataList.get(i).getCluster().numPoints))
							{
								//if they have the same num points we want to let it be the match 
								//regardless of the num matches
							}
							else if(oldRatio < newRatio)
							{
								//more points doesn't make a better match if the ratio higher 
								//keep the one with the lowest ratio
								boolean b = false;
							}
							else
							{
								theMatch.matches = 0;
								theMatch = cDataList.get(i);
								theMatch.matches = max = matches;
							}
						}
						else
						{
							if(theMatch != null)
							{
								theMatch.matches = 0;
							}
							theMatch = cDataList.get(i);
							theMatch.matches = max = matches;
						}
					}
					else if(matches == max && theMatch != null)
					{
						boolean b = false;
						//map to the one with fewest points
						if(theMatch.getCluster().numPoints < cDataList.get(i).getCluster().numPoints)
						{
							//this one has fewer points so keep the smaller ratio
						}
						else if(theMatch.getCluster().numPoints == cDataList.get(i).getCluster().numPoints)
						{
							//map to the smallest area
							if(theMatch.getCluster().getClusterArea() > 
								cDataList.get(i).getCluster().getClusterArea())
							{
								theMatch.matches = 0;
								theMatch = cDataList.get(i);
								theMatch.matches = max = matches;
							}
						}
						else
						{
							theMatch.matches = 0;
							theMatch = cDataList.get(i);
							theMatch.matches = max = matches;
						}
					}
					Pair<ClusterData, Integer> p = new Pair<ClusterData, Integer>(cDataList.get(i), matches);
					clusterMap.add(p);
				}
			}
			if(printTables)
			{
				System.out.println("\nPrinting match table for new cluster " + current.getClusterNumber() + ":\n" + 
						"New Cluster\tOld Cluster\tNew Cluster Points\tOld Cluster Points\tNew Cluster Area\t" +
						"Old Cluster Area\tNum Matches\n");
				for(int idx = 0; idx < clusterMap.size(); idx++)
				{
					System.out.println(current.getClusterNumber() + "\t\t"+ clusterMap.get(idx).fst.id + "\t\t" +
							current.numPoints+ "\t\t\t" + clusterMap.get(idx).fst.getCluster().numPoints
							 + "\t\t\t" + current.getClusterArea() + "\t" + clusterMap.get(idx).fst.getCluster().getClusterArea() +
							 "\t" + clusterMap.get(idx).snd + "\n");	
				}
			}
		}
		
		if(theMatch != null && theMatch.matches == 0)
		{
			theMatch = null;
		}
		
		return theMatch;
	}
	
	public boolean identifyMatchAsObstruction(ClusterData oldObj, ClusterData newObj)
	{
		boolean retVal = true;
		cluster oldCluster = oldObj.getCluster();
		cluster newCluster = newObj.getCluster();
		
		int matchPoints = oldObj.matches;
		int oldClusterPoints = oldCluster.numPoints;
		int newClusterPoints = newCluster.numPoints;
		
		double oldArea = oldCluster.getClusterArea();
		double newArea = newCluster.getClusterArea();
		
		System.out.println("Identified old cluster " + oldObj.id + " as a match to new cluster " +
				newObj.id + ", with " + oldObj.matches + " match points.\n" + "The old cluster area is "
				+ oldArea + " the new area is " + newArea + ".\n");
		
		if(oldArea > newArea)
		{
			if(oldClusterPoints > newClusterPoints)
			{
				//the new one may have been blocked by an object in the old one
				//keep both
			}
			else
			{
				//the old one may have had outliers that cause the area to be larger
				//it's likely that the new one is a better representation of the same object
				//keep the new one ditch the old one
				cDataList.remove(oldObj);
			}
		}
		else // oldarea <= newArea
		{
			if(oldClusterPoints > newClusterPoints)
			{
				//the old one is stronger probably the same object 
				//keep the old one
				retVal = false;
			}
			else
			{
				//the new one may be a combination of new objects that contain the old one 
				//keep both
			}
		}
		return retVal;
	}
	
	public boolean kMeansPointEqual(kMeansPoint left, kMeansPoint right)
	{
		boolean orientationEqual, locationEqual = false, scaleEqual, retVal = true;
		if(left != null && right != null)
		{
			Feature lf = left.getFeature();
			Feature rf =  right.getFeature();
			//arc distance
			retVal = classifier.match(lf, rf);
			/*this allows for various types of comparisons, but arc distance seems
			 * to work best 
			 * orientationEqual = lf.orientation == rf.orientation;
			scaleEqual = lf.scale == rf.scale;
			locationEqual = (lf.location[0] == rf.location[0] && 
					lf.location[1] == rf.location[1]);
			
			if(scaleEqual && orientationEqual && locationEqual)
			{*/
			/*	float [] lDims = null;
				float [] rDims = null;
				
				lDims = lf.descriptor;
				rDims = rf.descriptor;
				
				for(int i = 0; i < 128; i++)
				{
					if(lDims[i] != rDims[i])
					{
						retVal = false;
						break;
					}
				}
			/*}
			else
			{
				retVal = false;
			}*/
		}

		/*boolean retVal = true;
		float [] lDims = null;
		float [] rDims = null;
		if(left != null && right != null)
		{
			lDims = left.getDims();
			rDims = right.getDims();
		
			for(int i = 0; i < 128; i++)
			{
				if(lDims[i] != rDims[i])
				{
					retVal = false;
					break;
				}
			}
		}
		else
		{
			retVal = false;
		}
		return retVal;*/
		//return (lf == rf);
		//retVal = locationEqual;
		return retVal;
	}
	public boolean processPhoto(BufferedImage	image, Vector<Feature> feat, 
			boolean useMeanPoint/*, int numClusters */)
	{
		return processPhoto(image, feat, 
				useMeanPoint/*, int numClusters */, false);
		
	}
	public boolean processPhoto(BufferedImage	image, Vector<Feature> feat, 
			boolean useMeanPoint/*, int numClusters */, boolean testSame)
	{
		boolean			retVal = true;
		Vector<Feature>	features;
		cluster 		currentCluster = null;
		ClusterData		cData;
		int 			numClusters = 0;
		
		if(image != null)
		{
			//features = SIFTDriver.getFeaturesFromImage(image);
			//features = detectOutliers(features, 0);
			List<kMeansPoint> points = this.generatePointsList(feat);
			numClusters = ((int)Math.sqrt(feat.size()/2));
			//km = new kMeans(numClusters, points);
			if(testSame)
			{
				km.InitializekMeans(numClusters, points);
				km.runKMeansWithPriors(km.randomMeanValues);
			}
			else
			{
				km.InitializekMeans(numClusters, points);
				km.runKMeans();
			}
			boolean firstSet = cDataList.isEmpty();
			ArrayList<ClusterData> newGuys = new ArrayList<ClusterData>();
			
			for(int i = 0; i < km.getK(); i++)
			{
				currentCluster = km.getCluster(i);
				if(currentCluster != null)
				{
					if(!firstSet)
					{
						ClusterData match = trackClusters(currentCluster, useMeanPoint);
						if(match == null)
						{
							cData = new ClusterData(currentCluster, false, cDataList.size());
							newGuys.add(cData); 
						}
						else
						{
							cData = new ClusterData(currentCluster, false, cDataList.size());
							//if new cluster contains more points than old cluster 
							//add it and remove the new one? maybe cache it
							if(identifyMatchAsObstruction(match, cData))
							{
								newGuys.add(cData);
								//cDataList.remove(match);
							}
							setDirection(match, cData, newGuys.size());
						}
					}
					else //add them all the first time, there is nothing to map to
					{
						cData = new ClusterData(currentCluster, false, cDataList.size());
						cDataList.add(cData);
					}
					//System.out.println("Cluster " + i + " mean is " + currentCluster.getMean() + '\n' );
				}
			}
			if(!newGuys.isEmpty())
				cDataList.addAll(newGuys);		
		}
		else
		{
			retVal = false;
		}
		
		return retVal;
	}
	
	public List<kMeansPoint> generatePointsList(Vector<Feature> features)
	{
		int i;
		List<kMeansPoint> kMeansPoints;
		kMeansPoint	newPoint;
		kMeansPoints = new ArrayList<kMeansPoint>();
		Feature feat;
		
 		for(i = 0; i < features.size(); i++)
		{
			feat = features.get(i);
			newPoint = new kMeansPoint(feat, feat.descriptor, 0, 0);
			newPoint.numNeighbors = calcNumNeighbors(feat, features);
			kMeansPoints.add(newPoint);
		}
		return kMeansPoints;
		
	}
	
	public void setDirection(ClusterData oldClusterData, ClusterData newClusterData, int s)
	{
		double slope;
		double rise;
		double run;
		
		run = oldClusterData.getCluster().getActualMean().getFeature().location[0] -
			  newClusterData.getCluster().getActualMean().getFeature().location[0];
		rise = oldClusterData.getCluster().getActualMean().getFeature().location[1] -
		       newClusterData.getCluster().getActualMean().getFeature().location[1];

		//slope = rise/run;
		if(rise > 0 && run > 0)
		{
			if( rise > run) //north
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.NORTH);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.NORTH);
			}
			else //northeast
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.NORTHEAST);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.NORTHEAST);
			}
		}
		else if(rise < 0 && run < 0)
		{
			if(rise < run) //south
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.SOUTH);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.SOUTH);
			}
			else //southwest
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.SOUTHWEST);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.SOUTHWEST);
			}
		}
		else if(rise > 0 && run < 0)
		{
			if(Math.abs(rise) > Math.abs(run)) //north
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.NORTH);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.NORTH);
			}
			else //northwest
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.NORTHWEST);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.NORTHWEST);
			}
		}
		else if(rise < 0 && run > 0)
		{
			if(Math.abs(rise) > Math.abs(run)) //south
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.SOUTH);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.SOUTH);
			}
			else //southeast
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.SOUTHEAST);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.SOUTHEAST);
			}
		}
		else if(rise == 0 && run != 0)
		{
			if(run > 0) //east
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.EAST);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.EAST);
			}
			else // run < 0 west
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.WEST);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.WEST);
			}
		}
		else if(run == 0 && rise != 0)
		{
			if(rise > 0) //north
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.NORTH);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.NORTH);
			}
			else // rise < 0 south
			{
				oldClusterData.setMotion(true);
				oldClusterData.setDirection(Direction.SOUTH);
				newClusterData.setMotion(true);
				newClusterData.setDirection(Direction.SOUTH);
			}
		}
		else if(run == 0 && rise == 0) //no motion probably shouldn't happen
		{
			oldClusterData.setMotion(false);
			oldClusterData.setDirection(Direction.STATIONARY);
			newClusterData.setMotion(false);
			newClusterData.setDirection(Direction.STATIONARY);
		}
		String dirString = "The object appears to be moving in a ";
		String dir ;
		if(newClusterData.getDirection() == Direction.SOUTH)
		{
			dirString += "Southbound direction";
			dir = "S";
		}
		else if(newClusterData.getDirection() == Direction.NORTH)
		{
			dirString += "Northbound direction";
			dir = "N";
		}
		else if(newClusterData.getDirection() == Direction.EAST)
		{
			dirString += "Eastbound direction";
			dir = "E";
		}
		else if(newClusterData.getDirection() == Direction.WEST)
		{
			dirString += "Westbound direction";
			dir = "W";
		}
		else if(newClusterData.getDirection() == Direction.NORTHEAST)
		{
			dirString += "Northeasterly direction";
			dir = "NE";
		}
		else if(newClusterData.getDirection() == Direction.SOUTHEAST)
		{
			dirString += "Southeasterly direction";
			dir = "SE";
		}
		else if(newClusterData.getDirection() == Direction.NORTHWEST)
		{
			dirString += "Northwesterly direction";
			dir = "NW";
		}
		else if(newClusterData.getDirection() == Direction.SOUTHWEST)
		{
			dirString += "Southwesterly direction";
			dir = "SW";
		}
		else 
		{
			dirString += "Stationary position";
			dir = "ST";
		}
		
		dirString += "\n";
		System.out.println(dirString);
		String oString = "";
		int prevMatch = 1;
		int current = cDataList.size() + s ;
		int i;
		for( i = 0; i < cDataList.size(); i++)
		{
			if(cDataList.get(i).getCluster() == oldClusterData.getCluster())
			{
				prevMatch += i;
				break;
			}
		}

		float percent = oldClusterData.matches / (float) newClusterData.getCluster().numPoints;
		if(percent > 1)
		{
			percent = 100;
		}
		else
		{
			percent *= 100;
		}
		oString += current +
			"\t" + prevMatch + "\t" + percent + "%\t" + newClusterData.getCluster().getClusterArea() +
			"\t" + oldClusterData.getCluster().getClusterArea() + "\t" + dir + "\n";
		System.out.println("***\n" + oString);
	}
	
	public ClusterData trackClusters(cluster newCluster, boolean useMean)
	{
		ClusterData retVal = null;
		double threshold = .25;
		
		ClusterData match = mapToExistingCluster(newCluster, useMean, true);
		if(match != null)
		{
			if(!(match.matches >= threshold * newCluster.numPoints))
			{
				match = null;
				//not enough matches to be the same cluster
			}
		}
		
		retVal = match;
		
		return retVal;
	}
	
	public int calcNumNeighbors(Feature f, Vector<Feature> feats)
	{
		int numNeighbors = 0;
		double distance;
		//this one is the nearest neighbor and it's distance per feature
		double distances[][] = new double[feats.size()][2];
		//this one is the distance from a given feature to each point
		float distanceVector[] = new float[feats.size() -1];
		ArrayList < Pair<Feature, Double> > featureMap = null;
		double avgDist;
		Feature current, neighbor;
		int i = 0, j, pos;
		
		for(i =0; i < feats.size(); i++)
		{
			distances[i][0] = -1;
		}
		Collections.sort(feats);

		featureMap = new ArrayList < Pair<Feature, Double> >();
		for(i = 0; i < feats.size(); i++)
		{
			current = feats.get(i);
			for(j = 0; j < feats.size(); j++)
			{
				neighbor = feats.get(j);
				distance = current.descriptorDistance(neighbor);
				if(current != neighbor)
				{
					if (distances[i][0] == -1.0)
					{
						distances[i][0] = distance;
						distances[i][1] = j;
					}
					else if(distances[i][0] > distance)
					{	
						distances[i][0] = distance;
						distances[i][1] = j;
					}
				}
			}
			Pair<Feature, Double> p = new Pair<Feature, Double>(current, distances[i][0]);
			featureMap.add(p);
		}

		avgDist = 0;
		for(i = 0; i < feats.size(); i++)
		{
			avgDist += distances[i][0];
		}
		
		avgDist /= feats.size();
		
		pos = 0;
		for(i = 0; i < feats.size(); i++)
		{
			if(feats.get(i) != f)
			{
				distanceVector[pos] = f.descriptorDistance(feats.get(i));
				pos++;
			}
		}

		for(j = 0; j < featureMap.size()-1; j++)
		{
			if(distanceVector[j] <= avgDist)
			{
				numNeighbors++;
			}
		}
		
		return numNeighbors;
	}
	
	public Vector<Feature> detectOutliers(Vector<Feature> feats, int detectionType)
	{
		//change the table to be an array of pairs node and distance
		double distance;
		double distances[][] = new double[feats.size()][2];
		ArrayList < Pair<Feature, Double> > featureMap = null;
		double avgDist;
		double max = -1;
		Feature current, neighbor;
		int maxPos;
		int i = 0, j;
		double percentOutliers = .40;
		int numOutliers = (int) (percentOutliers * feats.size());
		
		for(i =0; i < feats.size(); i++)
		{
			distances[i][0] = -1;
		}
		Collections.sort(feats);

		featureMap = new ArrayList < Pair<Feature, Double> >();
		for(i = 0; i < feats.size(); i++)
		{
			current = feats.get(i);
			for(j = 0; j < feats.size(); j++)
			{
				neighbor = feats.get(j);
				distance = current.descriptorDistance(neighbor);
				if(current != neighbor)
				{
					if (distances[i][0] == -1.0)
					{
						distances[i][0] = distance;
						distances[i][1] = j;
					}
					else if(distances[i][0] > distance)
					{	
						distances[i][0] = distance;
						distances[i][1] = j;
					}
				}
			}
			Pair<Feature, Double> p = new Pair<Feature, Double>(current, distances[i][0]);
			featureMap.add(p);
		}
		maxPos = -888;
		int fullSize = feats.size();
		int rem = 0, pos, rem2, remFinal;
		avgDist = 0;
		for(i = 0; i < feats.size(); i++)
		{
			avgDist += distances[i][0];
		}
		
		avgDist /= feats.size();
		//remove out of average
		if(detectionType == 0)
		{
			for(j = 0; j < featureMap.size(); j++)
			{
				if(featureMap.get(j).snd > avgDist)
				{
					featureMap.remove(j);
					feats.remove(j);
					if(feats.size() == 1)
					{
						//none are in the average break and leave the last one
						//think about whether or not this is possible
						break;
					}
				}
			}
		}
		//remove top %
		else if(detectionType == 1)
		{
			for(i = 0; i < numOutliers; i++)
			{
				maxPos = 0;
				for(int k = 0; k< featureMap.size(); k++)
				{
					if(featureMap.get(k).snd > max)
					{
						maxPos = k;
						max = featureMap.get(k).snd;
					}
				}
				featureMap.remove(maxPos);
				feats.remove(maxPos);
			}
		}
		// this is out of date but should perform the same as the above removing top %
		else if(detectionType == 2)
		{
			for(j = 0; j < numOutliers; j++)
			{
				remFinal = rem2 = pos = 0;
				//find the first valid distance and set it to max
				for(int firstVal = 0; firstVal < fullSize; firstVal++)
				{
					if(distances[firstVal][0] != -888)
					{
						maxPos = firstVal - rem2;
						max = distances[firstVal][0];
						remFinal = rem2;
						break;
					}
					else
						rem2++; //the position of this entry in the distances table is first val but in the 
					            //feats vector it's firstVal - the num removed
				}
	
				//find the next valid distance and compare it 
				rem = 0;
				for(int k = 0; k < fullSize; k++)
				{
					if(distances[k][0] != -888)
						pos++;
					else
					{
						rem++;
						continue;
					}
					
					if(distances[k][0] > max)
					{
						max = distances[k][0];
						maxPos = k - rem;
						remFinal = rem;
					}
				}
				try{
					feats.remove(maxPos);
				}
				catch(ArrayIndexOutOfBoundsException e){
					System.out.println("max pos = " + maxPos + " rem = " + rem + " rem2 = " + rem2 +
							" remFinal = " + remFinal + "size = " + feats.size());
				}
				distances[maxPos+ remFinal][0] = -888;	//no longer valid
			}
		}
		if(feats.size() == 0)
		{
			System.out.println("Outlier Error \n");
		}
		
		return feats;
	}
	
	public static void main(String args[]) 
	{
		boolean start = true;
		boolean ok    = true;
		//I'm storing these and checking to waste time so I can sleep less
		long pTime = 0;
		long cTime = 0;
		long dTime = 0;
		//con = new RobotConnection(conAdd); this was in the constructor
		ImageProcessor    ip    = new ImageProcessor(/*"129.130.103.713",""*/);
		ip.con = new RobotConnection("129.130.102.140");
		//ip.con.goTo(3.840, -28.8, 1.893);
		//java client used to set goal then destroy it and use only hunter
		/*PlayerClient      robot = new PlayerClient("129.130.32.111", 6665);
		PlannerInterface  plani	= robot.requestInterfacePlanner(0, 'a');
		
		if(plani == null || robot == null)
			ok = false;*/
		while(ok)
		{
			ip.takePhotoAndProcess();
		}
		//use this for receiving images wirelessly from the robot because it's slow
		/*while(ok)
		{
			if(start)
			{
				//plani.setGoal (3000, 6000, 0); //values?
				//robot.close();
				pTime = System.currentTimeMillis();
				ip.takePhotoAndProcess();
				start = false;
			}
			else
			{
				cTime = System.currentTimeMillis();
				dTime = cTime  - pTime;
				if( dTime >= ip.processTime)
				{
					pTime = System.currentTimeMillis();
					ip.takePhotoAndProcess();
				}
				else
				{
					try
					{
						Thread.sleep(ip.processTime - dTime);
					}
					catch(InterruptedException e)
					{
						e.printStackTrace();
					}
				}
			}
		}*/
	}	
}
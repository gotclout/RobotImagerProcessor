/*

 * Implements the k-means algorithm

 *

 * Robert L. Foster Jr.

 * Computer Information Sciences

 * Kansas State University

 *

 * Created: October 21, 2008

 * Last updated: October 21, 2008

 *

 */

package image.processor;

import java.io.*;

import java.util.*;

import java.lang.*; 

import mpi.cbg.fly.Feature;

/**
 * Implements the k-means algorithm
 * @author	Robert L. Foster Jr.	rfoster@ksu.edu
 */ 

class kMeansPoint implements Comparable< kMeansPoint >, Serializable{

	/** Dimensions 0-127 */
	private float [] dims;

	/** Value in dimension x */
	private float x;

	/** Value in dimension y */
	private float y;

	/** Assigned cluster */
	private int clusterNumber;
	
	/**Feature Data**/
	private Feature siftFeat;
	
	/**Outlier from other features**/
	private boolean outlier;
	
	/**Number of neighbors ing avg dist**/
	public int numNeighbors = 0;
	
	public void setOutlier(boolean value)
	{
		this.outlier = value;
	}
	
	public boolean getOutlier()
	{
		return this.outlier;
	}
	
	public int compareTo( kMeansPoint k )
	{
		//return this.getFeature().scale < k.getFeature().scale ? 1 : this.getFeature().scale == k.getFeature().scale ? 0 : -1;
		return this.numNeighbors > k.numNeighbors ? 1 : this.numNeighbors == k.numNeighbors ? 0 : -1;
	}
	
	/**
	 * Creates a new instance of data point
	 *
	 * @param	_x	value in dimension x
	 * @param	_y	value in dimension y
	 */
	public kMeansPoint(Feature feat, float [] _dims, int _x, int _y) 
	{
		this.dims = _dims;
		this.x = _x;
		this.y = _y;
		this.clusterNumber = 0;
		this.siftFeat = feat;
	} // end of kMeansPoint()

	/**
	 * Assigns the data point to a cluster
	 *
	 * @param	_clusterNumber	the cluster to which this data point is to be assigned
	 */
	public void assignToCluster(int _clusterNumber) 
	{
		this.clusterNumber = _clusterNumber;
	} // end of assignToCluster()

	/**
	 * Returns the cluster to which the data point belongs
	 *
	 * @return	the cluster number to which the data point belongs
	 */
	public int getClusterNumber() 
	{
		return this.clusterNumber;
	} // end of getClusterNumber()

	/**
	 * Returns the value of data point in x dimension
	 *
	 * @return	the value in x dimension
	 */

	public float getX() 
	{
		return this.siftFeat.location[0];
	} // end of getX()
	
	/**
	 * Returns the value of data point in x dimension
	 *
	 * @return	the value in x dimension
	 */
	public float [] getDims()
	{
		return this.dims;
	} //end of getDims()

	/**
	 * Returns the value of data point in y dimension
	 *
	 * @return	the value in y dimension
	 */
	public float getY() 
	{
		return this.siftFeat.location[1];
	} // end of getY()
	
	/**
	 * Returns the feature of the kMeanspoint
	 *
	 * @return	the feature value
	 */
	public Feature getFeature()
	{
		return this.siftFeat;
	}
	
	/**
	 * Returns the feature of the kMeanspoint
	 *
	 * @return	the feature value
	 */
	public void setFeature(Feature feat)
	{
		this.siftFeat = feat;
	}

	/**
	 * Returns the distance between two data points
	 *
	 * @param	dp1 	the first data point
	 * @param	dp2 	the second data point
	 * @return	the distance between the two data points
	 */
	public static double distance(kMeansPoint dp1, kMeansPoint dp2) 
	{
		//sum += p1[i] - p2[i] *  p1[i] - p2[i]
		double result = 0;
		double resultX = dp1.getX() - dp2.getX();
		double resultY = dp1.getY() - dp2.getY();
		double sum;
		int i;
		sum = 0;
		for(i = 0; i < 128; i++)
		{
			sum += (dp1.dims[i]  - dp2.dims[i]) * (dp1.dims[i]  - dp2.dims[i]);
		}
		result = Math.sqrt(resultX*resultX + resultY*resultY);
		result = Math.sqrt(sum);
		return result;
	} // end of distance()
	
//	public static double distance(kMeansPoint dp1, kMeansPoint dp2) {
//
//		double result = 0;
//		double resultX = dp1.getFeature().location[0] - dp2.getFeature().location[0];
//		double resultY = dp1.getFeature().location[1] - dp2.getFeature().location[1];
//		result = Math.sqrt(resultX*resultX + resultY*resultY);
//		return result;
//
//		} // end of distance()

	/**
	 * Returns a string representation of this kMeansPoint
	 *
	 * @return	a string representation of this data point
	 */
	public String toString()
	{
		int i ;
		String retVal = "";
		retVal = "[" ;
		for(i = 0; i < 128; i++)
		{
			retVal += " " + this.dims[i];
			if(i != 127)
			{
				retVal += ",";
			}
		}
		retVal += "] (" + this.x + "," + this.y + ")[" + this.clusterNumber + "]";
		return retVal;
	} // end of toString()

	/**
	 * Main method -- to test the kMeansPoint class
	 *
	 * @param	args	command line arguments
	 */
	public static void main(String[] args)
	{
		kMeansPoint dp1 = new kMeansPoint(null, null, -3, -4);
		kMeansPoint dp2 = new kMeansPoint(null, null, 0, 4);
		System.out.println(kMeansPoint.distance(dp1, dp2));
		System.out.println(dp1.getX());
		System.out.println(dp2.getY());
		dp1.assignToCluster(7);
		System.out.println(dp1.getClusterNumber());
		dp1.assignToCluster(17);
		System.out.println(dp1.getClusterNumber());
		System.out.println(dp2.getClusterNumber());
		System.out.println(dp1);

	} // end of main()
} // end of class

public class kMeans 
{
	/** Number of clusters */
	private int k;

	/** Array of clusters */
	private cluster[] clusters;

	/** Number of iterations */
	private int nIterations;

	/** Vector of data points */
	private Vector kMeansPoints;

	/** Name of the input file */
	private String inputFileName;
	
	/**Array of values chosen as mean points*/
	public ArrayList<Integer> randomMeanValues = new ArrayList<Integer>();

	/**
	 * Returns a new instance of kMeans algorithm
	 *
	 * @param	k		number of clusters
	 * @param	inputFileName	name of the file containing input data
	 */
    public kMeans() 
    {
		//this means the creator will call the initialize method before use
	} // end of kMeans()
	/**
	 * Returns a new instance of kMeans algorithm
	 *
	 * @param	k		number of clusters
	 * @param	inputFileName	name of the file containing input data
	 */
    public kMeans(int k, String inputFileName) 
    {
		this.k = k;
		this.inputFileName = inputFileName;
		this.clusters = new cluster[this.k];
		this.nIterations = 0;
		this.kMeansPoints = new Vector();
	} // end of kMeans()

	/**
	 * Returns a new instance of kMeans algorithm
	 *
	 * @param	k		number of clusters
	 * @param	kMeansPoints	List containing objects of type kMeansPoint
	 */
    public kMeans(int k, List kMeansPoints) 
    {
		this.k = k;
		this.inputFileName = inputFileName;
		this.clusters = new cluster[this.k];
		this.nIterations = 0;
		this.kMeansPoints=new Vector(kMeansPoints);
	} // end of kMeans()
    
    /**
	 * Returns a new instance of kMeans algorithm
	 *
	 * @param	k		number of clusters
	 * @param	kMeansPoints	List containing objects of type kMeansPoint
	 */
    public void InitializekMeans(int k, List kMeansPoints) 
    {
		this.k = k;
		this.inputFileName = inputFileName;
		this.clusters = new cluster[this.k];
		this.nIterations = 0;
		this.kMeansPoints=new Vector(kMeansPoints);
	} // end of kMeans()
         
	/**
	 * Reads the input data from the file and stores the data points in the vector
	 */
	public void readData() throws IOException
	{
		BufferedReader in = new BufferedReader(new FileReader(this.inputFileName));
		String line = "";
		while ((line = in.readLine()) != null )
		{
			StringTokenizer st = new StringTokenizer(line, " \t\n\r\f,");
				if (st.countTokens() == 2) 
				{
					kMeansPoint dp = new kMeansPoint(null, null, Integer.parseInt(st.nextToken()), Integer.parseInt(st.nextToken()));
					dp.assignToCluster(0);
					this.kMeansPoints.add(dp);
                }
		}
		in.close();
	} // end of readData()

	public void calculatePointsPerCluster()
	{
		Vector<kMeansPoint> points;
		for(int i = 0; i < this.k; i++)
		{
			points = new Vector<kMeansPoint>();
			Iterator it = this.kMeansPoints.iterator();
			while (it.hasNext())
			{
				kMeansPoint point = (kMeansPoint)it.next();
				if(point.getClusterNumber() == i)
				{
					clusters[i].numPoints++;
					if(clusters[i].west < point.getFeature().location[0])
						clusters[i].west = point.getFeature().location[0];
					if(clusters[i].east > point.getFeature().location[0])
						clusters[i].east = point.getFeature().location[0];
					if(clusters[i].south < point.getFeature().location[1])
						clusters[i].south = point.getFeature().location[1];
					if(clusters[i].north > point.getFeature().location[1])
						clusters[i].north = point.getFeature().location[1];
					points.add(point);
				}
			}
			clusters[i].setClusterPoints(points);
			clusters[i].setClusterArea(Math.abs((clusters[i].west - clusters[i].east) * 
					(clusters[i].north - clusters[i].south)));
			clusters[i].calculateActualMean();
		}
	}
	/**
	 * Runs the k-means algorithm over the data set
	 */
	public void runKMeans() 
	{
		// Select k points as initial means
		int r = 0;
		int s = this.kMeansPoints.size();
		int incBy = s/k;
		Collections.sort( this.kMeansPoints );
		for (int i=0; i < k; i++){
			this.clusters[i] = new cluster(i);
			r = (int)(Math.random() * this.kMeansPoints.size());
			//r = i;
			randomMeanValues.add(r);
			this.clusters[i].setMean((kMeansPoint)(this.kMeansPoints.get(r)));
			//r += incBy;
		}
		do {
			// Form k clusters
			Iterator i = this.kMeansPoints.iterator();
			while (i.hasNext())
			{
				this.assignToCluster((kMeansPoint)(i.next()));
			}
			this.nIterations++;
		}
		// Repeat while centroids do not change
		while (this.updateMeans());
		calculatePointsPerCluster();
		//generateClusterFeatures(true);
	} // end of runKMeans()
	/**
	 * Runs the k-means algorithm over the data set with the same 
	 * values for means as prior runs
	 */
	public void runKMeansWithPriors(ArrayList<Integer> priors) 
	{
		// Select k points as initial means
		Collections.sort( this.kMeansPoints );
		for (int i=0; i < k; i++){
			this.clusters[i] = new cluster(i);
			int s = this.kMeansPoints.size();
			int r = priors.get(i);
			randomMeanValues.add(r);
			this.clusters[i].setMean((kMeansPoint)(this.kMeansPoints.get(r)));
		}
		do {
			// Form k clusters
			Iterator i = this.kMeansPoints.iterator();
			while (i.hasNext())
			{
				this.assignToCluster((kMeansPoint)(i.next()));
			}
			this.nIterations++;
		}
		// Repeat while centroids do not change
		while (this.updateMeans());
		calculatePointsPerCluster();
		//generateClusterFeatures(true);
	} // end of runKMeans()

	/**
	 * Assigns a data point to one of the k clusters based on its distance from the means of the clusters
	 *
	 * @param	dp	data point to be assigned
	 */
	private void assignToCluster(kMeansPoint dp) 
	{
		int currentCluster = dp.getClusterNumber();

		double minDistance = kMeansPoint.distance(dp, this.clusters[currentCluster].getMean());

		for (int i=0; i <this.k; i++)
		{
			if (kMeansPoint.distance(dp, this.clusters[i].getMean()) < minDistance) 
			{
				minDistance = kMeansPoint.distance(dp, this.clusters[i].getMean());
				currentCluster = i;
			}
		}
		dp.assignToCluster(currentCluster);	
	} // end of assignToCluster


	/**
	 * Updates the means of all k clusters, and returns if they have changed or not
	 *
	 * @return	have the updated means of the clusters changed or not
	 */

	private boolean updateMeans() 
	{
		boolean reply = false;
		int[] x = new int[this.k];
		int[] y = new int[this.k];
		float[][] dims = new float[this.k][128];
		int[] size = new int[this.k];
		Feature feat[] = new Feature[this.k];
		kMeansPoint[] pastMeans = new kMeansPoint[this.k];

		for (int i=0; i<this.k; i++) 
		{
			x[i] = 0;
			y[i] = 0;
			size[i] = 0;
			feat[i] = new Feature();
			for(int j = 0; j < 128; j++)
			{
				dims[i][j] = 0;
				feat[i].descriptor[j] = 0;
			}
			feat[i].location[0] = 0;
			feat[i].location[1] = 0;
			feat[i].orientation = 0;
			feat[i].scale = 0;
			
			//generateClusterFeatures(true);
			pastMeans[i] = this.clusters[i].getMean();
		}
	
		int idx = 0;
		Iterator i = this.kMeansPoints.iterator();

		while (i.hasNext()) 
		{
			kMeansPoint dp = (kMeansPoint)(i.next());
			int currentCluster = dp.getClusterNumber();
			float[] _dims = dp.getDims();
			x[currentCluster] += dp.getX();
			y[currentCluster] += dp.getY();
			feat[currentCluster].location[0] += dp.getFeature().location[0];
			feat[currentCluster].location[1] += dp.getFeature().location[1];
			feat[currentCluster].orientation += dp.getFeature().orientation;
			feat[currentCluster].scale += dp.getFeature().scale;
			
			for(int j = 0; j < 128; j++)
			{
				dims[currentCluster][j] += _dims[j];
				feat[currentCluster].descriptor[j] += dp.getFeature().descriptor[j];
			}

			size[currentCluster]++;
		}

		for (int j=0; j < this.k; j++ ) 
		{
			if(size[j] != 0) 
			{
				x[j] /= size[j];
				y[j] /= size[j];
				float[] d = new float[128];
				feat[j].location[0] /= size[j];
				feat[j].location[1] /= size[j];
				feat[j].orientation /= size[j];
				feat[j].scale /= size[j];
		
				for(int k = 0; k < 128; k++)
				{
					dims[j][k] /= size[j];
					d[k] = dims[j][k];
					feat[j].descriptor[k] /= size[j];
				}
				
				kMeansPoint temp = new kMeansPoint(feat[j], d, x[j], y[j]);
				temp.assignToCluster(j);
				this.clusters[j].setMean(temp);

				if (kMeansPoint.distance(pastMeans[j], this.clusters[j].getMean()) !=0 )
				{
					reply = true;
					//break;
				}
			}
		}
		return reply;
	} // end of updateMeans()
	
	/**
	 * Returns the value of k
	 *
	 * @return	the value of k
	 */
	public int getK() 
	{
		return this.k;
	} // end of getK()
	
	/**
	 * Returns the specified cluster by index
	 *
	 * @param	index	index of the cluster to be returned
	 * @return	return the specified cluster by index
	 */
	public cluster getCluster(int index) 
	{
		return this.clusters[index];
	} // end of getCluster()

	/**
	 * Returns the specified cluster by index
	 *
	 * @param	index	index of the cluster to be returned
	 * @return	return the specified cluster by index
	 */
	public Feature getClusterFeat(int index)
	{
		return getCluster(index).getMean().getFeature();
	}
	
	/**
	 * Returns the specified cluster by index
	 *
	 * @param	index	index of the cluster to be returned
	 * @return	return the specified cluster by index
	 */
	public Feature getClusterMeanFeat(int index)
	{
		return getCluster(index).getActualMean().getFeature();
	}
	
	public Feature getClusterMeanPoint(int index)
	{
		return getCluster(index).getMean().getFeature();
	}
	
	/**
	 * Returns the string output of the data points
	 *
	 * @return  the string output of the data points
	 */
	public String toString()
	{
		return this.kMeansPoints.toString();
	} // end of toString()
	
	public boolean detectOutliers()
	{
		boolean retVal = true;
		
		return retVal;
	}

	/**
	 * Returns the data points
	 *
	 * @return  the data points
	 */
	public boolean generateClusterFeatures(boolean firstRun)
	{
		boolean retVal = true;
		int	numClusters = this.k;
		int numPoints = 0, j, distIdx = 0;
		Iterator it;
		kMeansPoint point, cMean; 
		Feature newFeat = null;
		double distance[][];
		double rDistance;
		boolean detectedOutlier = false;
		if(k > 0)
		{
			for(int i = 0; i < k; i++)
			{
				it = this.kMeansPoints.iterator();
				distance = new double[clusters[i].numPoints][2];
				while (it.hasNext())
				{
					point = (kMeansPoint)it.next();
					cMean = this.clusters[i].getMean();
					if(point.getClusterNumber() == i)
					{
							if(newFeat == null)
							{
								newFeat = new Feature();
								newFeat.scale = newFeat.orientation = 0;
								newFeat.location = new float[2];
								newFeat.descriptor = new float [128];
								for(j = 0; j < 2; j++)
									newFeat.location[j] = 0;
								for(j = 0; j < 128; j++)
									newFeat.descriptor[j] = 0;
								
							}
							
							Feature currentFeat = point.getFeature();
							newFeat.scale += currentFeat.scale;
							newFeat.orientation += currentFeat.orientation;
							for(j = 0; j < 2; j++)
								newFeat.location[j] += currentFeat.location[j];
							for(j = 0; j < 128; j++)
								newFeat.descriptor[j] += currentFeat.descriptor[j];
							numPoints++;
					}
				}
				if(newFeat != null)
				{
					newFeat.scale /= numPoints;
					newFeat.orientation /= numPoints;
					for(j = 0; j < 2; j++)
						newFeat.location[j] /= numPoints;
					for(j = 0; j < 128; j++)
						newFeat.descriptor[j] /= numPoints;
					this.clusters[i].getMean().setFeature(newFeat);
					newFeat = null;
					distIdx = 0;
				}
			}
		}
		else
			retVal = false;
		
		return retVal;
	}
	
	/**
	 * Returns the data points
	 *
	 * @return  the data points
	 */
	public Vector getDataPoints() 
	{
		return this.kMeansPoints ;
	} // end of getDataPoints()

	/**
	 * Main method -- to test the kMeans class
	 *
	 * @param   args    command line arguments
	 */
	public static void main(String[] args) 
	{
		kMeans km = new kMeans(2, "input1");

		try 
		{
			km.readData();
		} 
		catch (Exception e) 
		{
			System.err.println(e);
			System.exit(-1);
		}         

		km.runKMeans();
		System.out.println(km);          

        } // end of main()
} // end of class

/*
 * Represents an abstraction for a cluster of data points in two dimensional space
 *
 * Robert L. Foster Jr.
 * Computer Information Sciences
 * Kansas State University
 *
 * Created: October 21, 2008
 * Last updated: October 21, 2008
 *
 */

/**
 * Represents an abstraction for a cluster of data points in two dimensional space
 * @author		rfoster@ksu.edu
 */
 class cluster 
 {
	/** Cluster Number */
	private int clusterNumber;

	/** Mean data point of this cluster */
	private kMeansPoint mean;
	
	/**The actual mean of points in this cluster */
	private kMeansPoint actualMean;
	
	/** Points that belong to this cluster **/
	private Vector<kMeansPoint> clusterPoints;
	
	/**The number of points that make up the cluster**/
	public int numPoints;
	
	/** The outermost boundaries **/
	double north, south, east, west;
	
	/** The area of the cluster **/
	double clusterArea;
	
	/**
	 * Returns a new instance of cluster
	 *
	 * @param	_clusterNumber	the cluster number of this cluster
	 */
	public cluster(int _clusterNumber) 
	{
		this.clusterNumber = _clusterNumber;
		this.numPoints = 0;
		north = south = east = west = 0;
	} // end of cluster()
	
	public boolean setClusterPoints(Vector<kMeansPoint> points)
	{
		boolean retVal = true;
		/*for(int i = 0; i < clusterPoints.size(); i++)
		{
			kMeansPoint point = clusterPoints.get(i);
			if(point.getClusterNumber() == clusterNumber)
			{
				clusterPoints.add(point);
			}
		}*/
		clusterPoints = points;
		return retVal;
	}
	
	public Vector<kMeansPoint> getClusterPoint()
	{
		return clusterPoints;
	}
	/**
	 * Sets the mean data point of this cluster
	 *
	 * @param	meanDataPoint	the new mean data point for this cluster
	 */
	public void setMean(kMeansPoint meanDataPoint) 
	{
		this.mean = meanDataPoint;
	} // end of setMean()
	
	public void setActualMean(kMeansPoint theMean)
	{
		this.actualMean = theMean;
	}
	
	public kMeansPoint getActualMean()
	{
		return this.actualMean;
	}
	/**
	 * Returns the mean data point of this cluster
	 *
	 * @return	the mean data point of this cluster
	 */
	public kMeansPoint getMean()
	{
		return this.mean;
	} // end of getMean()

	/**
	 * Returns the cluster number of this cluster
	 *
	 * @return	the cluster number of this cluster
	 */
	public int getClusterNumber() 
	{
		return this.clusterNumber;
	} // end of getClusterNumber()
	
	/**
	 * Returns the cluster area of this cluster
	 *
	 * @return	the cluster area
	 */
	public double getClusterArea() 
	{
		return this.clusterArea;
	} // end of getClusterNumber()
	
	/**
	 * Sets the cluster area of this cluster
	 *
	 * @return	void
	 */
	public void setClusterArea(double value) 
	{
		this.clusterArea = value;
	} // end of getClusterNumber()
	
	public void calculateActualMean()
	{
		Feature mean = new Feature();
		int j;
		Vector<kMeansPoint> points = this.getClusterPoint();
		Feature currentFeat;
		//init
		mean.scale = mean.orientation = 0;
		mean.location = new float[2];
		mean.descriptor = new float [128];
		for(j = 0; j < 2; j++)
			mean.location[j] = 0;
		for(j = 0; j < 128; j++)
			mean.descriptor[j] = 0;
		
		for(kMeansPoint k : points)
		{
			currentFeat = k.getFeature();
			
			mean.scale += currentFeat.scale;
			mean.orientation += currentFeat.orientation;
			for(j = 0; j < 2; j++)
				mean.location[j] += currentFeat.location[j];
			for(j = 0; j < 128; j++)
				mean.descriptor[j] += currentFeat.descriptor[j];
		}
		
		mean.scale /= numPoints;
		mean.orientation /= numPoints;
		for(j = 0; j < 2; j++)
			mean.location[j] /= numPoints;
		for(j = 0; j < 128; j++)
			mean.descriptor[j] /= numPoints;
		
		kMeansPoint meanPoint = new kMeansPoint(mean, mean.descriptor, 0, 0);
		this.setActualMean(meanPoint);
	}

	/**
	 * Main method -- to test the cluster class
	 *
	 * @param	args	command line arguments
	 */
	public static void main(String[] args) 
	{
		cluster c1 = new cluster(1);

		c1.setMean(new kMeansPoint(null, null, 3,4));

		System.out.println(c1.getMean());

	} // end of main()

} // end of class


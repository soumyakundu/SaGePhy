package se.cbb.jprime.topology;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.NavigableMap;
import java.util.TreeMap;

import org.jfree.util.PublicCloneable;

import se.cbb.jprime.math.PRNG;

/**
 * Represents discretisation information on the part of a host tree spanning <b>one</b> epoch.
 * <p/>
 * An epoch is a timespan defined by two branching events on the host tree without any
 * divergence events occurring the time in-between. The definition is intuitively extended
 * to include the single-arc epoch of a top time arc, and the epoch ending with leaves. 
 * <p/>
 * Essentially, this class defines an array of discretisation points. The "lateral"
 * dimension (or columns, if you wish) of this array corresponds to the contemporary
 * arcs of the epoch. The vertical dimension (or rows) corresponds to a slicing of
 * the epoch into several equidistant time sub-intervals. A discretisation mid-point has
 * then been placed in each such interval.
 * <p/>
 * When accessing the times of discretisation points, both of the two endpoints are included.
 * Note that the distance from one these to the adjacent point is half of the normal timestep.
 * <p/>
 * See <code>EpochDiscretiser</code> for more info.
 * 
 * @author Joel Sjöstrand.
 * @author Mehmood Alam Khan
 */
public class Epoch implements PublicCloneable {
	
	/** Number ID. */
	private int no;
	
	/** Vector of intersecting arcs in the epoch (defined via head vertex). */
	private int[] m_arcs;
	
	/** List of times of epoch, both endpoint times included. */
	private double[] m_times;
	
	/** Host Tree. */
	private RBTreeEpochDiscretiser tree;
		
	/** Timestep between discretised times. Stored explicitly for speed. */
	private double m_timestep;
	
	/** Random Arc selected fro transfer*/
	private int transferedToArc = -1;
	
	/**
	 * Constructor.
	 * @param no the number identifier.
	 * @param arcs the arcs intersecting the epoch's time span.
	 * @param loTime the time of the epoch's lower divergence event.
	 * @param upTime the time of the epoch's upper divergence event.
	 * @param noOfIvs the number of sub-intervals to slice epoch into.
	 */
	public Epoch(int no, List<Integer> arcs, double loTime, double upTime, int noOfIvs, RBTreeEpochDiscretiser tree) {
		this.no = no;
		m_arcs = new int[arcs.size()];
		int i = 0;
		for (int arc : arcs) {
			m_arcs[i++] = arc;
		}
		m_times = new double[noOfIvs + 2];
		this.tree = tree;
		m_timestep = (upTime - loTime) / noOfIvs;
	
		assert(upTime > loTime);
		
		// Add times. Treat endpoints separately.
		m_times[0] = loTime;
		for (i = 0; i < noOfIvs; ++i) {
			m_times[i+1] = (loTime + m_timestep / 2.0 + i * m_timestep);
		}
		m_times[m_times.length-1] = upTime;
	}
	
	/**
	 * Copy-constructor.
	 * @param orig original.
	 */
	public Epoch(Epoch orig) {
		this.m_arcs = new int[orig.m_arcs.length];
		System.arraycopy(orig.m_arcs, 0, this.m_arcs, 0, orig.m_arcs.length);
		this.m_times = new double[orig.m_times.length];
		System.arraycopy(orig.m_times, 0, this.m_times, 0, orig.m_times.length);
		this.m_timestep = orig.m_timestep;
	}
	
	/**
	 * Returns the discretised times of the epoch including endpoints,
	 * implying that the latter coincide with end times of epochs above/below.
	 * Also please note that e.g. for a 4-interval epoch with discretization
	 * times (t0,t1,t2,t3,t4,t5) and timestep dt we have t1-t0 = t5-t4 = dt/2,
	 * while t2-t1 = t3-t2 = t4-t3 = dt.
	 * @return the times of the epoch.
	 */
	public double[] getTimes() {
		return m_times;
	}
	
	/**
	 * Returns a discretised time.
	 * @param index the time index.
	 * @return the time.
	 */
	public double getTime(int index) {
		return m_times[index];
	}
	
	/**
	 * Returns the epoch time boundary closest to the leaves.
	 * @return the lowermost time.
	 */
	public double getLowerTime() {
		return m_times[0];
	}
	
	/**
	 * Returns the epoch time boundary closest to the top.
	 * @return the uppermost time.
	 */
	public double getUpperTime() {
		return m_times[m_times.length-1];
	}
	
	/**
	 * Returns the number of discretised times of an epoch,
	 * both end points included (two greater than the
	 * number of intervals).
	 * @return the number of discretised times.
	 */
	public int getNoOfTimes() {
		return m_times.length;
	}
	
	/**
	 * Returns the number of intervals, see also getNoOfTimes().
	 * @return the number of discretisation intervals.
	 */
	public int getNoOfIntervals() {
		return (m_times.length - 2);
	}
	
	/**
	 * Returns the timestep of discretization intervals.
	 * Note: the distance from an endpoint to its neighbour
	 * is only half the timestep.
	 * @return the timestep.
	 */
	public double getTimestep() {
		return m_timestep;
	}
	
	/**
	 * Returns the entire timespan of the epoch.
	 * @return the time span.
	 */
	public double getTimespan() {
		return (m_times[m_times.length-1] - m_times[0]);
	}
	
	/**
	 * Returns the vector of host tree arcs which intersect the epoch. These appear in the order used to reference a specific arc.
	 * Arcs are defined via their lower vertex in the host tree.
	 */
	public int[] getArcs() {
		return m_arcs;
	}
	
	/**
	 * Returns the arc at a specified index from the vector
	 * of arcs which intersect the epoch.
	 * Arcs are referenced via their lower vertex.
	 * @param index the index of the arc in the epoch.
	 * @return the arc, defined by its lower vertex.
	 */
	public int getArc(int index) {
		return m_arcs[index];
	}
	
	/**
	 * Returns the number of arcs the epoch spans.
	 * @return the number of arcs.
	 */
	public int getNoOfArcs() {
		return m_arcs.length;
	}
	
	/**
	 * Returns the total number of points over all spanned
	 * arcs in the epoch (i.e. number of times multiplied with
	 * number of arcs).
	 * @return the number of points.
	 */
	public int getNoOfPoints() {
		return (m_times.length * m_arcs.length);
	}
	
	@Override
	public Object clone() throws CloneNotSupportedException {
		return new Epoch(this);
	}
	
	/**
	 * Samples an arc uniformly from the epoch.
	 * @param prng PRNG.
	 */
	public int sampleArc(PRNG prng) {
		return this.m_arcs[prng.nextInt(this.m_arcs.length)];
	}
	
	public int getTranferedToArc(){
		return this.transferedToArc;
	}
	
	public void setTranferedToArc(int arc){
		this.transferedToArc= arc;
	}
	
	/**
	 * Find index of the given arc of species tree
	 * @param arc
	 * 
	 */
	public int findIndexOfArc(int speciesArc){
		int idx=0;
		while(true){
			if (this.m_arcs[idx] == speciesArc){
				break;
			}
			++idx;
		}
		return idx;
	}
	
	/**
	 * Samples an arc from the epoch, excluding a given arc.
	 * @param prng PRNG.
	 * @param excludeArc arc to exclude.
	 */
	public int sampleArc(PRNG prng, int excludeArc, int fromArc, double eventTime, String distance_bias, ArrayList<Integer> emptyArcs) {
		int arc;
		if (this.m_arcs.length == 1 && excludeArc == m_arcs[0]) {
			throw new IllegalArgumentException("Cannot exclude arc from sampling (it's the only one!): " + excludeArc);
		}
		if (!distance_bias.equals("none")) {
			NavigableMap<BigDecimal, Integer> dist = setDistance(excludeArc, eventTime, distance_bias, emptyArcs);
			BigDecimal total = dist.lastKey();
			
			//System.out.println("Final Total: " + total);
			
			BigDecimal rand;
			
			//System.out.println("Check A");
			
			do {
				rand = total.multiply(BigDecimal.valueOf(prng.nextDouble()));
				//System.out.println("Random: " + rand);
			} while (rand == total);
			
			//System.out.println("Check B");
			
			arc = dist.higherEntry(rand).getValue();
		} else {
			
			//System.out.println(this.m_arcs);
			//System.out.println("No Distance Bias");
			
			int to = prng.nextInt(this.m_arcs.length);
			arc = this.m_arcs[to];
			while (true) {
				if (to != fromArc && arc != excludeArc) {
					break;
				}
				to = prng.nextInt(this.m_arcs.length);
				arc = this.m_arcs[to];
			}
		}
		this.setTranferedToArc(arc);
		return arc;
	}
	
	/** Returns HashMap of phylogenetic distances to all lineages in the epoch eligible to receive transfers from given recipient. */	
	
	public NavigableMap<BigDecimal, Integer> setDistance(int excludeArc, double eventTime, String distance_bias, ArrayList<Integer> emptyArcs) {
		//System.out.println("Distance Bias");
		BigDecimal total = new BigDecimal(0);
		BigDecimal value;
		NavigableMap<BigDecimal, Integer> phy_dist = new TreeMap<BigDecimal, Integer>();
		for (int x : m_arcs) {
			if (x != excludeArc) {
				if (emptyArcs == null || !emptyArcs.contains(x)) {
					if (distance_bias.equals("exponential")) {
						//System.out.println("Exponential");
						double doubleValue = (Math.pow(10.0, (-1.0 * tree.getDistance(excludeArc, x, eventTime))));
						if (Double.isInfinite(doubleValue)) {
							doubleValue = Double.MAX_VALUE;
						}
						value = BigDecimal.valueOf(doubleValue);
					} else {
						//System.out.println("Simple");
						double doubleValue = (1.0 / tree.getDistance(excludeArc, x, eventTime));
						if (Double.isInfinite(doubleValue)) {
							doubleValue = Double.MAX_VALUE;
						}
						value = BigDecimal.valueOf(doubleValue);
					}
					total = total.add(value);
					phy_dist.put(total, x);
					
					//System.out.println("Distance: " + tree.getDistance(excludeArc, x, eventTime));
					//System.out.println("Value: " + value);
					//System.out.println("Total: " + total);
					
				}
			}
		}
		
		//System.out.println(phy_dist);
		
		return phy_dist;
	}
	
	/**
	 * Returns the epoch ID number.
	 * @return the number.
	 */
	public int getNo() {
		return this.no;
	}
}



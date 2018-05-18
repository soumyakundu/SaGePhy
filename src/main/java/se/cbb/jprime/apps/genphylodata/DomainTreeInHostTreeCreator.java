package se.cbb.jprime.apps.genphylodata;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import se.cbb.jprime.apps.genphylodata.GuestVertex.Event;
import se.cbb.jprime.io.NewickIOException;
import se.cbb.jprime.io.NewickVertex;
import se.cbb.jprime.io.PrIMENewickTree;
import se.cbb.jprime.math.ExponentialDistribution;
import se.cbb.jprime.math.NumberManipulation;
import se.cbb.jprime.math.PRNG;
import se.cbb.jprime.topology.Epoch;
import se.cbb.jprime.topology.NamesMap;
import se.cbb.jprime.topology.RBTree;
import se.cbb.jprime.topology.RBTreeEpochDiscretiser;
import se.cbb.jprime.topology.TimesMap;
import se.cbb.jprime.topology.TopologyException;

/**
 * Creates unpruned trees evolving over a host tree.
 *
 * @author Joel Sjöstrand.
 * @author Mehmood Alam Khan
 */

public class DomainTreeInHostTreeCreator implements UnprunedGuestTreeCreator {

	/** Host tree. */
	private RBTreeEpochDiscretiser hostTree;

	/** Host names. */
	private NamesMap hostNames;

	/** Duplication rate. */
	private double lambda;

	/** Loss rate. */
	private double mu;

	/** Transfer rate. */
	private double tau;

	/** Replacing Transfer rate. */
	private double theta;

	/** Transfer distance bias. */
	private String distance_bias;
	
	/** Gene birth sampling bias. */
	private double gene_birth;

	/** Host Tree or Guest Tree. */
	private boolean ishost;
	
	/** Intra-gene transfer rate. */
	private double intra_gene;
	
	/** Intra-species transfer rate. */
	private double intra_species;

	/** Sampling probability. */
	private double rho;

	/**
	 * Constructor.
	 * @param host host tree.
	 * @param lambda duplication rate.
	 * @param mu loss rate.
	 * @param tau transfer rate.
	 * @param rho probability of sampling leaf.
	 * @throws TopologyException.
	 * @throws NewickIOException.
	 */

	public DomainTreeInHostTreeCreator(PrIMENewickTree host, ArrayList<PrIMENewickTree> guests, double lambda, double mu, double tau, double theta, String distance_bias, double gene_birth, boolean ishost, double intra_gene, double intra_species, double rho, Double stem) throws TopologyException, NewickIOException {

		// Host tree.
		RBTree S = new RBTree(host, "HostTree");
		TimesMap hostTimes = host.getTimesMap("HostTimes");
		this.hostNames = host.getVertexNamesMap(true, "HostNames");
		if (stem != null) {
			hostTimes.getArcTimes()[S.getRoot()] = stem;
		}
		// Hack: Set 0 stem to eps.
		if (hostTimes.getArcTime(S.getRoot()) <= 0.0) {
			hostTimes.getArcTimes()[S.getRoot()] = 1.0e-64;
		}
		this.hostTree = new RBTreeEpochDiscretiser(S, hostNames, hostTimes);

		// generate epoch and arc id information for each edge of the species tree. this is useful when doing sampling realization from DLTRS model

		// Rates.
		this.lambda = lambda;
		this.mu = mu;
		this.tau = tau;
		this.theta = theta;
		this.distance_bias = distance_bias;
		this.gene_birth = gene_birth;
		this.ishost = ishost;
		this.intra_gene = intra_gene;
		this.intra_species = intra_species;
		this.rho = rho;
		if (lambda < 0 || mu < 0 || tau < 0) {
			throw new IllegalArgumentException("Cannot have rate less than 0.");
		}
		if (rho < 0 || rho > 1) {
			throw new IllegalArgumentException("Cannot have leaf sampling probability outside [0,1].");
		}
	}

	@Override
	public GuestVertex createUnprunedTree(PRNG prng) {

		// Currently processed lineages.
		LinkedList<GuestVertex> alive = new LinkedList<GuestVertex>();

		// Used to sort replacing transfer lineages by decreasing event times.
		NavigableMap<Double, GuestVertex> sortedlist = new TreeMap<Double, GuestVertex>();

		GuestVertex root;
		int myRoot;

		// Single lineage at tip.
		if (this.ishost == true) {
			myRoot = hostTree.getRoot();
			root = this.createGuestVertex(myRoot, hostTree.getTipToLeafTime(), prng);
		} else if (this.gene_birth != 0) {
			myRoot = hostTree.sampleRoot(prng, this.gene_birth);
			root = this.createGuestVertex(myRoot, hostTree.getVertexUpperTime(myRoot), prng);
		} else {
			myRoot = hostTree.getRoot();
			root = this.createGuestVertex(myRoot, hostTree.getTipToLeafTime(), prng);
		}
		//myRoot = 189;
		//root = this.createGuestVertex(myRoot, hostTree.getVertexUpperTime(myRoot), prng);
		//System.out.println(hostTree.getVertexUpperTime(myRoot));
		//System.out.println(hostTree.getTipToLeafTime());
		System.out.println(myRoot);
		alive.add(root);

		// Recursively process lineages.
		while (!alive.isEmpty()) {

			// Prints out current lineages and their events.
			//this.listprint(alive);
			//System.out.println("-------------------------");

			GuestVertex lin = alive.pop();
			if (lin.event == Event.LOSS || lin.event == Event.REPLACING_LOSS || lin.event == Event.LEAF || lin.event == Event.UNSAMPLED_LEAF) {
				// Lineage ends.
				if (alive.isEmpty()) {
					if (!sortedlist.isEmpty()) {
						alive.add(sortedlist.pollLastEntry().getValue());
					}
				}
				continue;
			}

			GuestVertex lc = null;
			GuestVertex rc = null;

			if (lin.event == Event.SPECIATION) {
				lc = this.createGuestVertex(hostTree.getLeftChild(lin.sigma), lin.abstime, prng);
				rc = this.createGuestVertex(hostTree.getRightChild(lin.sigma), lin.abstime, prng);

			} else if (lin.event == Event.DUPLICATION) {
				lc = this.createGuestVertex(lin.sigma, lin.abstime, prng);
				rc = this.createGuestVertex(lin.sigma, lin.abstime, prng);

			} else if (lin.event == Event.ADDITIVE_TRANSFER) {
				if (prng.nextDouble() < 0.5) {
					lc = this.createGuestVertex(lin.sigma, lin.abstime, prng);
					int transferedToArc= lin.epoch.sampleArc(prng, lin.sigma, lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, this.distance_bias, null);
					lin.setTransferedFromArc(lin.sigma);
					rc = this.createGuestVertex(transferedToArc, lin.abstime, prng);
					lin.setTransferedToArc(lin.epoch.getTranferedToArc());

					//System.out.println("Additive Transfer to: " + transferedToArc);

				} else {
					rc = this.createGuestVertex(lin.sigma, lin.abstime, prng);
					int transferedToArc= lin.epoch.sampleArc(prng, lin.sigma, lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, this.distance_bias, null);
					lin.setTransferedFromArc(lin.sigma);
					lc = this.createGuestVertex(transferedToArc, lin.abstime, prng);
					lin.setTransferedToArc(lin.epoch.getTranferedToArc());

					//System.out.println("Additive Transfer to: " + transferedToArc);
				}

			} else if (lin.event == Event.REPLACING_TRANSFER) {
				if (prng.nextDouble() < 0.5) {
					lc = this.createGuestVertex(lin.sigma, lin.abstime, prng);
					int transferedToArc= lin.epoch.sampleArc(prng, lin.sigma, lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, this.distance_bias, null);
					lin.setTransferedFromArc(lin.sigma);
					lin.setTransferedToArc(lin.epoch.getTranferedToArc());

					GuestVertex node = findVertex(root, lin);
					ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

					while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {

						if (!emptyArcs.contains(lin.transferedToArc)) {
							emptyArcs.add(lin.transferedToArc);
						}

						if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
							lin.transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, this.distance_bias, emptyArcs);
							node = findVertex(root, lin);
						}
					}

					if (node != null) {
						rc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng);
						node.event = Event.REPLACING_LOSS;
						node.abstime = lin.abstime;
						Iterator<Map.Entry<Double, GuestVertex>> iter = sortedlist.entrySet().iterator();
						while (iter.hasNext()) {
							GuestVertex x = iter.next().getValue();
							if (isAncestor(node, x)) {
								iter.remove();
							}
						}
						node.setChildren(null);

						//System.out.println("Replacing Transfer to: " + transferedToArc);

					} else {
						lin.event = Event.ADDITIVE_TRANSFER;
						rc = this.createGuestVertex(transferedToArc, lin.abstime, prng);
						System.out.println("Could not replace from: " + lin.sigma + " at " + lin.abstime);
					}

				} else {
					rc = this.createGuestVertex(lin.sigma, lin.abstime, prng);
					int transferedToArc= lin.epoch.sampleArc(prng, lin.sigma, lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, this.distance_bias, null);
					lin.setTransferedFromArc(lin.sigma);
					lin.setTransferedToArc(lin.epoch.getTranferedToArc());

					GuestVertex node = findVertex(root, lin);
					ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

					while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {

						if (!emptyArcs.contains(lin.transferedToArc)) {
							emptyArcs.add(lin.transferedToArc);
						}

						if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
							lin.transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, this.distance_bias, emptyArcs);
							node = findVertex(root, lin);
						}
					}

					if (node != null) {
						lc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng);
						node.event = Event.REPLACING_LOSS;
						node.abstime = lin.abstime;
						Iterator<Map.Entry<Double, GuestVertex>> iter = sortedlist.entrySet().iterator();
						while (iter.hasNext()) {
							GuestVertex x = iter.next().getValue();
							if (isAncestor(node, x)) {
								iter.remove();
							}
						}
						node.setChildren(null);

						//System.out.println("Replacing Transfer to: " + transferedToArc);

					} else {
						lin.event = Event.ADDITIVE_TRANSFER;
						lc = this.createGuestVertex(transferedToArc, lin.abstime, prng);
						System.out.println("Could not replace from: " + lin.sigma + " at " + lin.abstime);
					}
				}

			} else {
				throw new UnsupportedOperationException("Unexpected event type.");
			}

			ArrayList<NewickVertex> children = new ArrayList<NewickVertex>(2);
			children.add(lc);
			children.add(rc);
			lin.setChildren(children);
			lc.setParent(lin);
			rc.setParent(lin);

			if (lc.event != Event.REPLACING_TRANSFER) {
				alive.add(lc);
			} else {
				sortedlist.put(lc.abstime, lc);
			}

			if (rc.event != Event.REPLACING_TRANSFER) {
				alive.add(rc);
			} else {
				sortedlist.put(rc.abstime, rc);
			}

			if (alive.isEmpty()) {
				if (!sortedlist.isEmpty()) {
					alive.add(sortedlist.pollLastEntry().getValue());
				}
			}
		}

		// Restore 0 length stem.
		if (root.getBranchLength() <= 1.0e-32) {
			root.setBranchLength(0.0);
		}

		//System.out.println("Still to be replaced:");
		//System.out.println(sortedlist);
		//System.out.println("-------------------------");

		return root;
	}

	//Prints all elements of a given LinkedList of GuestVertices.
	public void listprint(LinkedList<GuestVertex> lineages) {
		for(GuestVertex vertex: lineages) {
			System.out.println(vertex.sigma);
			System.out.println(vertex.event + ": " + vertex.abstime);
		}
	}

	//Finds a given vertex in a given tree.
	public GuestVertex findVertex(GuestVertex node, GuestVertex x) {
		if (node != null) {
			if (node.sigma == x.getTransferedToArc() && node != x.getRightChild() && node != x.getLeftChild() && node.abstime < x.abstime && (node.abstime + node.getBranchLength()) > x.abstime) {
				return node;
			} else {
				GuestVertex new_node = findVertex(node.getLeftChild(), x);
				if (new_node == null) {
					new_node = findVertex(node.getRightChild(), x);
				}
				return new_node;
			}
		} else {
			return null;
		}
	}

	//Determines if a given node is an ancestor of another node.
	public boolean isAncestor(GuestVertex replaced, GuestVertex node) {
		if (replaced != null) {
			if (replaced == node) {
				return true;
			} else {
				boolean new_node = isAncestor(replaced.getLeftChild(), node);
				if (new_node == false) {
					new_node = isAncestor(replaced.getRightChild(), node);
				}
				return new_node;
			}
		} else {
			return false;
		}
	}

	/**
	 * Samples a guest vertex, given a process starting in host arc X at a given time.
	 * @param X host arc.
	 * @param startTime start time of process.
	 * @param prng PRNG.
	 * @return guest vertex.
	 */
	private GuestVertex createGuestVertex(int X, double startTime, PRNG prng) {
		boolean isRoot = this.hostTree.isRoot(X);
		double sum = isRoot ? this.lambda + this.mu : this.lambda + this.mu + this.tau;
		if (sum == 0.0) { sum = 1e-48; }
		ExponentialDistribution pd = new ExponentialDistribution(sum);
		double lowerTime = this.hostTree.getVertexTime(X);
		double branchTime = pd.sampleValue(prng);
		double eventTime = startTime - branchTime;
		GuestVertex.Event event;
		Epoch epoch;
		if (eventTime <= lowerTime) {
			// LEAF OR SPECIATION.
			eventTime = lowerTime;
			branchTime = NumberManipulation.roundToSignificantFigures(startTime - eventTime, 8);
			if (this.hostTree.isLeaf(X)) {
				event = (prng.nextDouble() < this.rho ? GuestVertex.Event.LEAF : GuestVertex.Event.UNSAMPLED_LEAF);
			} else {
				event = GuestVertex.Event.SPECIATION;
			}
			epoch = this.hostTree.getEpochAbove(X);
		} else {
			// DUPLICATION, LOSS OR TRANSFER.
			double rnd = prng.nextDouble();
			if (isRoot) {
				if (rnd < this.lambda / sum) {
					event = GuestVertex.Event.DUPLICATION;
				} else {
					event = GuestVertex.Event.LOSS;
				}
			} else {
				if (rnd >= (this.mu + this.tau) / sum) {
					event = GuestVertex.Event.DUPLICATION;
				} else if (rnd < this.mu / sum) {
					event = GuestVertex.Event.LOSS;
				} else {
					double trn = prng.nextDouble();
					if (trn < this.theta) {
						event = GuestVertex.Event.REPLACING_TRANSFER;
					} else {
						event = GuestVertex.Event.ADDITIVE_TRANSFER;
					}
				}
			}
			// Find correct epoch.
			int epno = this.hostTree.getEpochNoAbove(X);
			while (this.hostTree.getEpoch(epno).getUpperTime() < eventTime) {
				epno++;
			}
			epoch = this.hostTree.getEpoch(epno);
		}
		return new GuestVertex(event, X, epoch, eventTime, branchTime);
	}

	@Override
	public List<Integer> getHostLeaves() {
		return this.hostTree.getLeaves();
	}

	@Override
	public String getInfo(GuestVertex guestRoot, boolean doML) {
		StringBuilder sb = new StringBuilder(1024);
		int noOfVertices = 0;
		int noOfLeaves = 0;
		int noOfSpecs = 0;
		int noOfDups = 0;
		int noOfLosses = 0;
		int noOfReplacing_Losses = 0;
		int noOfAdditive_Trans = 0;
		int noOfReplacing_Trans = 0;
		double totalTime = 0.0;
		double totalTimeBeneathStem = 0.0;

		LinkedList<GuestVertex> vertices = new LinkedList<GuestVertex>();
		if (guestRoot != null) {
			vertices.add(guestRoot);
		}
		double hostRootTime = hostTree.getVertexTime(hostTree.getRoot());
		while (!vertices.isEmpty()) {
			GuestVertex v = vertices.pop();
			if (!v.isLeaf()) {
				vertices.add(v.getLeftChild());
				vertices.add(v.getRightChild());
			}
			noOfVertices++;
			totalTime += v.getBranchLength();
			totalTimeBeneathStem += Math.max(Math.min(hostRootTime - v.abstime, v.getBranchLength()), 0.0);
			switch (v.event) {
			case DUPLICATION:
				noOfDups++;
				break;
			case LOSS:
				noOfLosses++;
				break;
			case REPLACING_LOSS:
				noOfReplacing_Losses++;
				break;
			case ADDITIVE_TRANSFER:
				noOfAdditive_Trans++;
				break;
			case REPLACING_TRANSFER:
				noOfReplacing_Trans++;
				break;
			case SPECIATION:
				noOfSpecs++;
				break;
			case LEAF:
			case UNSAMPLED_LEAF:
				noOfLeaves++;
				break;
			default:
				throw new UnsupportedOperationException("Unexpected event type.");
			}
		}
		totalTime = NumberManipulation.roundToSignificantFigures(totalTime, 8);
		totalTimeBeneathStem = NumberManipulation.roundToSignificantFigures(totalTimeBeneathStem, 8);

		sb.append("No. of vertices:\t").append(noOfVertices).append('\n');
		sb.append("No. of extant leaves:\t").append(noOfLeaves).append('\n');
		sb.append("No. of speciations:\t").append(noOfSpecs).append('\n');
		sb.append("No. of duplications:\t").append(noOfDups).append('\n');
		sb.append("No. of losses:\t").append(noOfLosses).append('\n');
		sb.append("No. of replacing losses:\t").append(noOfReplacing_Losses).append('\n');
		sb.append("No. of additive transfers:\t").append(noOfAdditive_Trans).append('\n');
		sb.append("No. of replacing transfers:\t").append(noOfReplacing_Trans).append('\n');
		sb.append("Total branch time:\t").append(totalTime).append('\n');
		sb.append("Total branch time beneath host stem:\t").append(totalTimeBeneathStem).append('\n');
		if (doML) {
			double dupMLEst = NumberManipulation.roundToSignificantFigures(noOfDups / totalTime, 8);
			double lossMLEst = NumberManipulation.roundToSignificantFigures(noOfLosses / totalTime, 8);
			double replacing_lossMLEst = NumberManipulation.roundToSignificantFigures(noOfReplacing_Losses / totalTime, 8);
			double additive_transMLEst = NumberManipulation.roundToSignificantFigures(noOfAdditive_Trans / totalTimeBeneathStem, 8);  // Excl. stem-spanning arcs.
			double replacing_transMLEst = NumberManipulation.roundToSignificantFigures(noOfReplacing_Trans / totalTimeBeneathStem, 8);  // Excl. stem-spanning arcs.
			sb.append("Duplication ML estimate:\t").append(dupMLEst).append('\n');
			sb.append("Loss ML estimate:\t").append(lossMLEst).append('\n');
			sb.append("Replacing Loss ML estimate:\t").append(replacing_lossMLEst).append('\n');
			sb.append("Additive Transfer ML estimate:\t").append(additive_transMLEst).append('\n');
			sb.append("Replacing Transfer ML estimate:\t").append(replacing_transMLEst).append('\n');
		}
		return sb.toString();
	}

	@Override
	public String getLeafMap(GuestVertex guestRoot) {
		StringBuilder sb = new StringBuilder(1024);
		LinkedList<GuestVertex> vertices = new LinkedList<GuestVertex>();
		if (guestRoot != null) {
			vertices.add(guestRoot);
		}
		while (!vertices.isEmpty()) {
			GuestVertex v = vertices.pop();
			if (!v.isLeaf()) {
				vertices.add(v.getLeftChild());
				vertices.add(v.getRightChild());
			} else {
				if (v.event == Event.LEAF || v.event == Event.UNSAMPLED_LEAF) {
					sb.append(v.getName()).append('\t').append(hostNames.get(v.sigma)).append('\n');
				}
			}
		}
		return sb.toString();
	}

	@Override
	public String getSigma(GuestVertex guestRoot) {
		StringBuilder sb = new StringBuilder(4096);
		sb.append("# GUEST-TO-HOST MAP\n");
		sb.append("Host tree:\t").append(hostTree.toString()).append('\n');
		sb.append("Guest vertex name:\tGuest vertex ID:\tGuest vertex type:\tGuest vertex time:\tHost vertex/arc ID:\tHost epoch ID:\n");
		LinkedList<GuestVertex> vertices = new LinkedList<GuestVertex>();
		if (guestRoot != null) {
			vertices.add(guestRoot);
		}
		while (!vertices.isEmpty()) {
			GuestVertex v = vertices.pop();
			if (!v.isLeaf()) {
				vertices.add(v.getLeftChild());
				vertices.add(v.getRightChild());
			}
			sb.append(v.getName()).append('\t');
			sb.append(v.getNumber()).append('\t');
			sb.append(v.event.toString()).append('\t');
			sb.append(v.abstime).append('\t');
			sb.append(v.sigma).append('\t');
			sb.append(v.epoch.getNo()).append('\n');

		}
		return sb.toString();
	}

	@Override
	public String getHost() {
		return this.hostTree.toString();
	}
}

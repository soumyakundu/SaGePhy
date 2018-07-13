package se.cbb.jprime.apps.genphylodata;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.NavigableMap;
import java.util.TreeMap;

import se.cbb.jprime.apps.genphylodata.GuestVertex.Event;
import se.cbb.jprime.io.NewickIOException;
import se.cbb.jprime.io.NewickVertex;
import se.cbb.jprime.io.PrIMENewickTree;
import se.cbb.jprime.math.ExponentialDistribution;
import se.cbb.jprime.math.NumberManipulation;
import se.cbb.jprime.math.PRNG;
import se.cbb.jprime.misc.Pair;
import se.cbb.jprime.topology.Epoch;
import se.cbb.jprime.topology.NamesMap;
import se.cbb.jprime.topology.RBTree;
import se.cbb.jprime.topology.RBTreeEpochDiscretiser;
import se.cbb.jprime.topology.TimesMap;
import se.cbb.jprime.topology.TopologyException;

/**
 * Creates domain tree evolving over multiple gene trees and one species tree.
 *
 * @author Soumya Kundu
 */

public class DomainTreeInHostTreeCreator implements UnprunedGuestTreeCreator {

	/** Species tree. */
	private RBTreeEpochDiscretiser speciesTree;
	
	/** Species names. */
	private NamesMap speciesNames;
	
	/** Guest trees. */
	private ArrayList<RBTreeEpochDiscretiser> guestTrees = 
			new ArrayList<RBTreeEpochDiscretiser>();

	/** Guest names. */
	private ArrayList<NamesMap> guestNames = new ArrayList<NamesMap>();

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
	
	private double distance_bias_rate;
	
	private boolean doDomainBirth;
	
	/** Gene birth sampling bias. */
	private double domain_birth;

	/** Host Tree or Guest Tree. */
	private boolean isHost;
	
	/** Intra-gene transfer rate. */
	private double inter_gene;
	
	/** Intra-species transfer rate. */
	private double inter_species;

	/** Sampling probability. */
	private double rho;
	
	/** Domain tree root. */
	private GuestVertex root;
	
	private HashMap<RBTreeEpochDiscretiser, Integer> used;
	
	/**
	 * Constructor.
	 * @param species species tree.
	 * @param guests guest trees.
	 * @param lambda duplication rate.
	 * @param mu loss rate.
	 * @param tau transfer rate.
	 * @param theta replacing transfer rate.
	 * @param distance_bias type of transfer distance bias.
	 * @param doDomainBirth sample domain birth location on gene tree.
	 * @param domain_birth coefficient for domain birth sampling.
	 * @param isHost is this the species tree.
	 * @param inter_gene inter-gene transfer probability.
	 * @param inter_species inter-species transfer probability.
	 * @param rho probability of sampling leaf.
	 * @param stem 
	 * @throws TopologyException.
	 * @throws NewickIOException.
	 */

	public DomainTreeInHostTreeCreator(PrIMENewickTree species, 
			ArrayList<PrIMENewickTree> guests, double lambda, double mu, 
			double tau, double theta, String distance_bias, 
			double distance_bias_rate, boolean doDomainBirth, 
			double domain_birth, boolean isHost, double inter_gene, 
			double inter_species, double rho, Double stem) 
				throws TopologyException, NewickIOException {

		// Species tree.
		RBTree H = new RBTree(species, species.getTreeName());
		TimesMap speciesTimes = species.getTimesMap("SpeciesTimes");
		this.speciesNames = species.getVertexNamesMap(true, "SpeciesNames");
		if (stem != null) {
			speciesTimes.getArcTimes()[H.getRoot()] = stem;
		}
		if (speciesTimes.getArcTime(H.getRoot()) <= 0.0) {
			speciesTimes.getArcTimes()[H.getRoot()] = 1.0e-64;
		}
		this.speciesTree = new RBTreeEpochDiscretiser(H, speciesNames, speciesTimes);
		
		// Gene trees.
		for (int i = 0; i < guests.size(); i++) {
			RBTree S = new RBTree(guests.get(i), guests.get(i).getTreeName());
			TimesMap guestTimes = guests.get(i).getTimesMap("GuestTimes" + i);
			this.guestNames.add(guests.get(i).getVertexNamesMap(true, "GuestNames" + i));
			if (guestTimes.getArcTime(S.getRoot()) <= 0.0) {
				guestTimes.getArcTimes()[S.getRoot()] = 1.0e-64;
			}
			this.guestTrees.add(new RBTreeEpochDiscretiser(S, guestNames.get(i), guestTimes));
		}
		
		// generate epoch and arc id information for each edge of the species 
		// tree. this is useful when doing sampling realization from DLTRS model.

		// Parameters.
		this.lambda = lambda;
		this.mu = mu;
		this.tau = tau;
		this.theta = theta;
		this.distance_bias = distance_bias;
		this.distance_bias_rate = distance_bias_rate;
		this.doDomainBirth = doDomainBirth;
		this.domain_birth = domain_birth;
		this.isHost = isHost;
		this.inter_gene = inter_gene;
		this.inter_species = inter_species;
		this.rho = rho;
		if (lambda < 0 || mu < 0 || tau < 0 || domain_birth < 0) {
			throw new IllegalArgumentException("Cannot have rate less than 0.");
		}
		if (theta < 0 || theta > 1 || inter_gene < 0 || inter_gene > 1 ||
			inter_species < 0 || inter_species > 1 || rho < 0 || rho > 1) {
			throw new IllegalArgumentException("Cannot have probability outside [0,1].");
		}
	}

	@Override
	public GuestVertex createUnprunedTree(PRNG prng) throws NewickIOException {

		// Currently processed lineages.
		LinkedList<GuestVertex> alive = new LinkedList<GuestVertex>();

		// Used to sort replacing transfer lineages by decreasing event times.
		NavigableMap<Double, GuestVertex> sortedlist = new TreeMap<Double, GuestVertex>();

		GuestVertex root;
		int myRoot;
		
		int first_guest = prng.nextInt(guestTrees.size());

		this.used = new HashMap<RBTreeEpochDiscretiser, Integer>();
		for (RBTreeEpochDiscretiser k : guestTrees) {
			used.put(k, 0);
		}
		
		// Single lineage at tip.
		if (this.isHost == true) {
			myRoot = guestTrees.get(first_guest).getRoot();
			root = this.createGuestVertex(myRoot, 
					guestTrees.get(first_guest).getTipToLeafTime(), prng, 
					guestTrees.get(first_guest));
		} else if (this.doDomainBirth != false) {
			myRoot = guestTrees.get(first_guest).sampleRoot(prng, this.domain_birth);
			root = this.createGuestVertex(myRoot, 
					guestTrees.get(first_guest).getVertexUpperTime(myRoot), 
					prng, guestTrees.get(first_guest));
		} else {
			myRoot = guestTrees.get(first_guest).getRoot();
			root = this.createGuestVertex(myRoot, 
					guestTrees.get(first_guest).getTipToLeafTime(), prng, 
					guestTrees.get(first_guest));
		}
		
		//System.out.println("Root: " + myRoot);
		alive.add(root);
		this.root = root;

		// Recursively process lineages.
		while (!alive.isEmpty()) {
			GuestVertex lin = alive.pop();
			//System.out.println(lin.getHostVertex() + ": " + lin.event);
			if (lin.event == Event.LOSS || lin.event == Event.REPLACING_LOSS ||
				lin.event == Event.LEAF || lin.event == Event.UNSAMPLED_LEAF) {
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
				lc = this.createGuestVertex(lin.guestTree.getLeftChild(lin.sigma), 
						lin.abstime, prng, lin.guestTree);
				rc = this.createGuestVertex(lin.guestTree.getRightChild(lin.sigma), 
						lin.abstime, prng, lin.guestTree);

			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.DUPLICATION) {
				lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
				rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);

			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.ADDITIVE_TRANSFER_INTRAGENE_INTRASPECIES) {
				if (prng.nextDouble() < 0.5) {
					int transferedToArc = lin.epoch.sampleArc(prng, lin.sigma,
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime,
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							true);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(rc.getGuestTree());
					}
						
				} else {
					int transferedToArc = lin.epoch.sampleArc(prng, lin.sigma,
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							true);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(lc.getGuestTree());
					}
				}
				
			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.ADDITIVE_TRANSFER_INTRAGENE_INTERSPECIES) {
				if (prng.nextDouble() < 0.5) {
					int transferedToArc= lin.epoch.sampleArc(prng, lin.sigma, 
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							false);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(rc.getGuestTree());
					}

				} else {
					int transferedToArc= lin.epoch.sampleArc(prng, lin.sigma, 
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							false);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(lc.getGuestTree());
					}
				}
				
			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.ADDITIVE_TRANSFER_INTERGENE_INTRASPECIES) {
				if (prng.nextDouble() < 0.5) {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, true, null, prng);
					int transferedToArc = igtransfer.first;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, igtransfer.second);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(rc.getGuestTree());
					}
						
				} else {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, true, null, prng);
					int transferedToArc = igtransfer.first;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, igtransfer.second);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(lc.getGuestTree());
					}
				}
				
			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.ADDITIVE_TRANSFER_INTERGENE_INTERSPECIES) {
				if (prng.nextDouble() < 0.5) {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, false, null, prng);
					int transferedToArc = igtransfer.first;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, igtransfer.second);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(rc.getGuestTree());
					}

				} else {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, false, null, prng);
					int transferedToArc = igtransfer.first;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, igtransfer.second);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());
						lin.setTransferedFromGuest(lin.getGuestTree());
						lin.setTransferedToGuest(lc.getGuestTree());
					}
				}

			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.REPLACING_TRANSFER_INTRAGENE_INTRASPECIES) {
				if (prng.nextDouble() < 0.5) {
					int transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							true);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, lin.guestTree);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) 
								&& lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								lin.transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
										lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
										this.distance_bias, this.distance_bias_rate, emptyArcs, 
										lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
										true);
								node = findVertex(this.root, lin, lin.guestTree);
							}
						}

						if (node != null) {
							rc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, lin.guestTree);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());

						} else {
							rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
							lin.event = Event.ADDITIVE_TRANSFER_INTRAGENE_INTRASPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}

				} else {
					int transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							true);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, lin.guestTree);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) && 
								lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								lin.transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
										lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
										this.distance_bias, this.distance_bias_rate, emptyArcs, 
										lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
										true);
								node = findVertex(this.root, lin, lin.guestTree);
							}
						}

						if (node != null) {
							lc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, lin.guestTree);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());

						} else {
							lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
							lin.event = Event.ADDITIVE_TRANSFER_INTRAGENE_INTRASPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}
				}
				
			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.REPLACING_TRANSFER_INTRAGENE_INTERSPECIES) {
				if (prng.nextDouble() < 0.5) {
					int transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							false);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, lin.guestTree);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) 
								&& lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								lin.transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
										lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
										this.distance_bias, this.distance_bias_rate, emptyArcs, 
										lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
										false);
								node = findVertex(this.root, lin, lin.guestTree);
							}
						}

						if (node != null) {
							rc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, lin.guestTree);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());

						} else {
							rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
							lin.event = Event.ADDITIVE_TRANSFER_INTRAGENE_INTERSPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}

				} else {
					int transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
							lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
							this.distance_bias, this.distance_bias_rate, null, 
							lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
							false);
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, lin.guestTree);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) && 
								lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								lin.transferedToArc = lin.epoch.sampleArc(prng, lin.sigma, 
										lin.epoch.findIndexOfArc(lin.sigma), lin.abstime, 
										this.distance_bias, this.distance_bias_rate, emptyArcs, 
										lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex(), 
										false);
								node = findVertex(this.root, lin, lin.guestTree);
							}
						}

						if (node != null) {
							lc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, lin.guestTree);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());

						} else {
							lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, lin.guestTree);
							lin.event = Event.ADDITIVE_TRANSFER_INTRAGENE_INTERSPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}
				}
				
			//--------------------------------------------------------------//
								
			} else if (lin.event == Event.REPLACING_TRANSFER_INTERGENE_INTRASPECIES) {
				if (prng.nextDouble() < 0.5) {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, true, null, prng);
					int transferedToArc = igtransfer.first;
					RBTreeEpochDiscretiser orig_gene = igtransfer.second;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, igtransfer.second);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) && 
								lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								igtransfer = this.sampleInterGeneArc(lin, true, emptyArcs, prng);
								lin.transferedToArc = igtransfer.first;
								node = findVertex(this.root, lin, igtransfer.second);
							}
						}

						if (node != null) {
							rc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, igtransfer.second);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());
							

						} else {
							rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, orig_gene);
							lin.event = Event.ADDITIVE_TRANSFER_INTERGENE_INTRASPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}

				} else {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, true, null, prng);
					int transferedToArc = igtransfer.first;
					RBTreeEpochDiscretiser orig_gene = igtransfer.second;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, igtransfer.second);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) && 
								lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								igtransfer = this.sampleInterGeneArc(lin, true, emptyArcs, prng);
								lin.transferedToArc = igtransfer.first;
								node = findVertex(this.root, lin, igtransfer.second);
							}
						}

						if (node != null) {
							lc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, igtransfer.second);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());

						} else {
							lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, orig_gene);
							lin.event = Event.ADDITIVE_TRANSFER_INTERGENE_INTRASPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}
				}
				
			//--------------------------------------------------------------//
				
			} else if (lin.event == Event.REPLACING_TRANSFER_INTERGENE_INTERSPECIES) {
				if (prng.nextDouble() < 0.5) {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, false, null, prng);
					int transferedToArc = igtransfer.first;
					RBTreeEpochDiscretiser orig_gene = igtransfer.second;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, lc, prng);
						lin = swapped.first;
						lc = swapped.second;
					} else {
						lc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, igtransfer.second);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) && 
								lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								igtransfer = this.sampleInterGeneArc(lin, false, emptyArcs, prng);
								lin.transferedToArc = igtransfer.first;
								node = findVertex(this.root, lin, igtransfer.second);
							}
						}

						if (node != null) {
							rc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, igtransfer.second);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());
							

						} else {
							rc = this.createGuestVertex(transferedToArc, lin.abstime, prng, orig_gene);
							lin.event = Event.ADDITIVE_TRANSFER_INTERGENE_INTERSPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(rc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}

				} else {
					Pair<Integer, RBTreeEpochDiscretiser> igtransfer = this.sampleInterGeneArc(lin, false, null, prng);
					int transferedToArc = igtransfer.first;
					RBTreeEpochDiscretiser orig_gene = igtransfer.second;
					if (transferedToArc == -5) {
						Pair<GuestVertex, GuestVertex> swapped = swap(lin, rc, prng);
						lin = swapped.first;
						rc = swapped.second;
					} else {
						rc = this.createGuestVertex(lin.sigma, lin.abstime, prng, lin.guestTree);
						lin.setTransferedFromArc(lin.sigma);
						lin.setTransferedToArc(lin.epoch.getTransferedToArc());

						GuestVertex node = findVertex(this.root, lin, igtransfer.second);
						ArrayList<Integer> emptyArcs = new ArrayList<Integer>();

						while (node == null && emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1) && 
								lin.transferedToArc != -5) {

							if (!emptyArcs.contains(lin.transferedToArc)) {
								emptyArcs.add(lin.transferedToArc);
							}

							if (emptyArcs.size() != (lin.epoch.getNoOfArcs() - 1)) {
								igtransfer = this.sampleInterGeneArc(lin, false, emptyArcs, prng);
								lin.transferedToArc = igtransfer.first;
								node = findVertex(this.root, lin, igtransfer.second);
							}
						}

						if (node != null) {
							lc = this.createGuestVertex(lin.transferedToArc, lin.abstime, prng, igtransfer.second);
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
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());

						} else {
							lc = this.createGuestVertex(transferedToArc, lin.abstime, prng, orig_gene);
							lin.event = Event.ADDITIVE_TRANSFER_INTERGENE_INTERSPECIES;
							lin.setTransferedToArc(transferedToArc);
							lin.epoch.setTransferedToArc(transferedToArc);
							lin.setTransferedFromGuest(lin.getGuestTree());
							lin.setTransferedToGuest(lc.getGuestTree());
							//System.out.println("Could not replace from - Tree: " + 
							//		lin.guestTreeIndex + ", Lineage: " + lin.sigma + 
							//		" at " + lin.abstime);
						}
					}
				}
				
			//--------------------------------------------------------------//

			} else {
				throw new UnsupportedOperationException("Unexpected event type.");
			}

			if (lin != null) {
				ArrayList<NewickVertex> children = new ArrayList<NewickVertex>(2);
				children.add(lc);
				children.add(rc);
				lin.setChildren(children);
				lc.setParent(lin);
				rc.setParent(lin);
			} else {
				//System.out.println("Null Lineage -------- Null Lineage");
			}
				
			if (lc != null) {
				if (lc.event != Event.REPLACING_TRANSFER_INTRAGENE_INTRASPECIES && 
					lc.event != Event.REPLACING_TRANSFER_INTRAGENE_INTERSPECIES && 
					lc.event != Event.REPLACING_TRANSFER_INTERGENE_INTRASPECIES && 
					lc.event != Event.REPLACING_TRANSFER_INTERGENE_INTERSPECIES) {
					alive.add(lc);
				} else {
					sortedlist.put(lc.abstime, lc);
				}
			}

			if (rc != null) {
				if (rc.event != Event.REPLACING_TRANSFER_INTRAGENE_INTRASPECIES && 
					rc.event != Event.REPLACING_TRANSFER_INTRAGENE_INTERSPECIES &&
					rc.event != Event.REPLACING_TRANSFER_INTERGENE_INTRASPECIES && 
					rc.event != Event.REPLACING_TRANSFER_INTERGENE_INTERSPECIES) {
					alive.add(rc);
				} else {
					sortedlist.put(rc.abstime, rc);
				}
			}

			if (alive.isEmpty()) {
				if (!sortedlist.isEmpty()) {
					alive.add(sortedlist.pollLastEntry().getValue());
				}
			}
		}

		// Restore 0 length stem.
		if (this.root.getBranchLength() <= 1.0e-32) {
			this.root.setBranchLength(0.0);
		}

		return this.root;
	}

	//Prints all elements of a given LinkedList of GuestVertices.
	public void listprint(LinkedList<GuestVertex> lineages) {
		for(GuestVertex vertex: lineages) {
			System.out.println(vertex.sigma);
			System.out.println(vertex.event + ": " + vertex.abstime);
		}
	}

	//Finds a given vertex in a given tree.
	public GuestVertex findVertex(GuestVertex node, GuestVertex x, RBTreeEpochDiscretiser guestTree) {
		if (node != null) {
			if (node.sigma == x.getTransferedToArc() && node.guestTree == guestTree && 
					node != x.getRightChild() && node != x.getLeftChild() && 
					node.abstime < x.abstime && (node.abstime + node.getBranchLength()) > x.abstime) {
				return node;
			} else {
				GuestVertex new_node = findVertex(node.getLeftChild(), x, guestTree);
				if (new_node == null) {
					new_node = findVertex(node.getRightChild(), x, guestTree);
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
	
	//Replaces parent vertex with child vertex.
	public Pair<GuestVertex, GuestVertex> swap(GuestVertex lin, GuestVertex child, PRNG prng) {
		child = this.createGuestVertex(lin.sigma, lin.abstime + lin.branchtime, prng, lin.guestTree);
		if (lin != this.root) {
			if (lin.getParent().getChildren().get(0) == lin) {
				lin.getParent().getChildren().set(0, child);
			} else {
				lin.getParent().getChildren().set(1, child);
			}
		} else {
			this.root = child;
		}
		if (child != this.root) {
			child.setParent(lin.getParent());
		}
		lin.setParent(null);
		lin = null;
		Pair<GuestVertex, GuestVertex> swapped = new Pair<GuestVertex, GuestVertex>(lin, child);
		return swapped;
	}
	
	public Pair<Integer, RBTreeEpochDiscretiser> sampleInterGeneArc(GuestVertex lin, 
			Boolean intraSpecies, ArrayList<Integer> emptyArcs, PRNG prng) throws NewickIOException{
		ArrayList<RBTreeEpochDiscretiser> otherTrees = new ArrayList<RBTreeEpochDiscretiser>();
		ArrayList<Pair<Integer, RBTreeEpochDiscretiser>> candidates = 
				new ArrayList<Pair<Integer, RBTreeEpochDiscretiser>>();
		for (RBTreeEpochDiscretiser tree : this.guestTrees) {
			if (tree != lin.guestTree) {
				otherTrees.add(tree);
			}
		}
		for (RBTreeEpochDiscretiser tree : otherTrees) {
			int k = tree.getRBTree().getNoOfVertices();
			for (int i = 0; i < k; i++) {
				int arcHost = tree.getRBTree().getNewickTree().getVertex(i).getHostVertex();
				int hostArc = lin.guestTree.getRBTree().getNewickTree().getVertex(lin.sigma).getHostVertex();
				if (i != lin.sigma && ((emptyArcs == null) ? true : !emptyArcs.contains(i)) && 
					tree.getVertexTime(i) < lin.abstime && tree.getVertexUpperTime(i) > lin.abstime && 
					((intraSpecies) ? (arcHost == hostArc) : (arcHost != hostArc))) {
					candidates.add(new Pair<Integer, RBTreeEpochDiscretiser>(i, tree));
				}
			}
		}
		if (candidates.size() < 1) {
			return new Pair<Integer, RBTreeEpochDiscretiser>(-5, null);
		} else {
			int candidate = prng.nextInt(candidates.size());
			lin.epoch.setTransferedToArc(candidates.get(candidate).first);
			return candidates.get(candidate);
		}
	}

	/**
	 * Samples a guest vertex, given a process starting in host arc X at a given time.
	 * @param X host arc.
	 * @param startTime start time of process.
	 * @param prng PRNG.
	 * @return guest vertex.
	 */
	private GuestVertex createGuestVertex(int X, double startTime, PRNG prng, RBTreeEpochDiscretiser guestTree) {
		boolean isRoot = guestTree.isRoot(X);
		double sum = isRoot ? this.lambda + this.mu : this.lambda + this.mu + this.tau;
		if (sum == 0.0) { sum = 1e-48; }
		ExponentialDistribution pd = new ExponentialDistribution(sum);
		double lowerTime = guestTree.getVertexTime(X);
		double branchTime = pd.sampleValue(prng);
		double eventTime = startTime - branchTime;
		GuestVertex.Event event;
		Epoch epoch;
		if (eventTime <= lowerTime) {
			// LEAF OR SPECIATION.
			eventTime = lowerTime;
			branchTime = NumberManipulation.roundToSignificantFigures(startTime - eventTime, 8);
			if (guestTree.isLeaf(X)) {
				event = (prng.nextDouble() < this.rho ? GuestVertex.Event.LEAF : GuestVertex.Event.UNSAMPLED_LEAF);
			} else {
				event = GuestVertex.Event.SPECIATION;
			}
			epoch = guestTree.getEpochAbove(X);
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
					double trn_gen = prng.nextDouble();
					double trn_spe = prng.nextDouble();
					if (trn < this.theta) {
						if (trn_gen < this.inter_gene) {
							if (trn_spe < this.inter_species) {
								event = GuestVertex.Event.REPLACING_TRANSFER_INTERGENE_INTERSPECIES;
							} else {
								event = GuestVertex.Event.REPLACING_TRANSFER_INTERGENE_INTRASPECIES;
							}
						} else {
							if (trn_spe < this.inter_species) {
								event = GuestVertex.Event.REPLACING_TRANSFER_INTRAGENE_INTERSPECIES;
							} else {
								event = GuestVertex.Event.REPLACING_TRANSFER_INTRAGENE_INTRASPECIES;
							}
						}
					} else {
						if (trn_gen < this.inter_gene) {
							if (trn_spe < this.inter_species) {
								event = GuestVertex.Event.ADDITIVE_TRANSFER_INTERGENE_INTERSPECIES;
							} else {
								event = GuestVertex.Event.ADDITIVE_TRANSFER_INTERGENE_INTRASPECIES;
							}
						} else {
							if (trn_spe < this.inter_species) {
								event = GuestVertex.Event.ADDITIVE_TRANSFER_INTRAGENE_INTERSPECIES;
							} else {
								event = GuestVertex.Event.ADDITIVE_TRANSFER_INTRAGENE_INTRASPECIES;
							}
						}
					}
				}
			}
			// Find correct epoch.
			int epno = guestTree.getEpochNoAbove(X);
			while (guestTree.getEpoch(epno).getUpperTime() < eventTime) {
				epno++;
			}
			epoch = guestTree.getEpoch(epno);
		}
		this.used.put(guestTree, used.get(guestTree) + 1);
		return new GuestVertex(event, X, epoch, eventTime, branchTime, guestTree);
	}

	@Override
	public List<Integer> getHostLeaves() {
		return this.speciesTree.getLeaves();
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
		int noOfAdditive_Trans_Intragene_Intraspecies = 0;
		int noOfAdditive_Trans_Intragene_Interspecies = 0;
		int noOfAdditive_Trans_Intergene_Intraspecies = 0;
		int noOfAdditive_Trans_Intergene_Interspecies = 0;
		int noOfReplacing_Trans_Intragene_Intraspecies = 0;
		int noOfReplacing_Trans_Intragene_Interspecies = 0;
		int noOfReplacing_Trans_Intergene_Intraspecies = 0;
		int noOfReplacing_Trans_Intergene_Interspecies = 0;
		double totalTime = 0.0;
		double totalTimeBeneathStem = 0.0;

		LinkedList<GuestVertex> vertices = new LinkedList<GuestVertex>();
		if (guestRoot != null) {
			vertices.add(guestRoot);
		}
		//System.out.println(guestRoot.guestTree.getNoOfArcs(1));
		double hostRootTime = 0.0;
		hostRootTime = this.speciesTree.getVertexTime(this.speciesTree.getRoot()); //guestRoot.guestTree.getVertexTime(guestRoot.guestTree.getRoot());
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
			case ADDITIVE_TRANSFER_INTRAGENE_INTRASPECIES:
				noOfAdditive_Trans_Intragene_Intraspecies++;
				break;
			case ADDITIVE_TRANSFER_INTRAGENE_INTERSPECIES:
				noOfAdditive_Trans_Intragene_Interspecies++;
				break;
			case ADDITIVE_TRANSFER_INTERGENE_INTRASPECIES:
				noOfAdditive_Trans_Intergene_Intraspecies++;
				break;
			case ADDITIVE_TRANSFER_INTERGENE_INTERSPECIES:
				noOfAdditive_Trans_Intergene_Interspecies++;
				break;
			case REPLACING_TRANSFER_INTRAGENE_INTRASPECIES:
				noOfReplacing_Trans_Intragene_Intraspecies++;
				break;
			case REPLACING_TRANSFER_INTRAGENE_INTERSPECIES:
				noOfReplacing_Trans_Intragene_Interspecies++;
				break;
			case REPLACING_TRANSFER_INTERGENE_INTRASPECIES:
				noOfReplacing_Trans_Intergene_Intraspecies++;
				break;
			case REPLACING_TRANSFER_INTERGENE_INTERSPECIES:
				noOfReplacing_Trans_Intergene_Interspecies++;
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
		sb.append("No. of additive transfers intragene intraspecies:\t").append(noOfAdditive_Trans_Intragene_Intraspecies).append('\n');
		sb.append("No. of additive transfers intragene interspecies:\t").append(noOfAdditive_Trans_Intragene_Interspecies).append('\n');
		sb.append("No. of additive transfers intergene intraspecies:\t").append(noOfAdditive_Trans_Intergene_Intraspecies).append('\n');
		sb.append("No. of additive transfers intergene interspecies:\t").append(noOfAdditive_Trans_Intergene_Interspecies).append('\n');
		sb.append("No. of replacing transfers intragene intraspecies:\t").append(noOfReplacing_Trans_Intragene_Intraspecies).append('\n');
		sb.append("No. of replacing transfers intragene interspecies:\t").append(noOfReplacing_Trans_Intragene_Interspecies).append('\n');
		sb.append("No. of replacing transfers intergene intraspecies:\t").append(noOfReplacing_Trans_Intergene_Intraspecies).append('\n');
		sb.append("No. of replacing transfers intergene interspecies:\t").append(noOfReplacing_Trans_Intergene_Interspecies).append('\n');
		sb.append("Total branch time:\t").append(totalTime).append('\n');
		sb.append("Total branch time beneath host stem:\t").append(totalTimeBeneathStem).append('\n');
		if (doML) {
			double dupMLEst = NumberManipulation.roundToSignificantFigures(noOfDups / totalTime, 8);
			double lossMLEst = NumberManipulation.roundToSignificantFigures(noOfLosses / totalTime, 8);
			double replacing_lossMLEst = NumberManipulation.roundToSignificantFigures(noOfReplacing_Losses / totalTime, 8);
			double additive_trans_intragene_intraspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfAdditive_Trans_Intragene_Intraspecies / totalTime, 8);
			double additive_trans_intragene_interspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfAdditive_Trans_Intragene_Interspecies / totalTime, 8);
			double additive_trans_intergene_intraspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfAdditive_Trans_Intergene_Intraspecies / totalTime, 8);
			double additive_trans_intergene_interspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfAdditive_Trans_Intergene_Interspecies / totalTime, 8);
			double replacing_trans_intragene_intraspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfReplacing_Trans_Intragene_Intraspecies / totalTime, 8);
			double replacing_trans_intragene_interspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfReplacing_Trans_Intragene_Interspecies / totalTime, 8);
			double replacing_trans_intergene_intraspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfReplacing_Trans_Intergene_Intraspecies / totalTime, 8);
			double replacing_trans_intergene_interspecies_MLEst = NumberManipulation.roundToSignificantFigures(noOfReplacing_Trans_Intergene_Interspecies / totalTime, 8);
			sb.append("Duplication ML estimate:\t").append(dupMLEst).append('\n');
			sb.append("Loss ML estimate:\t").append(lossMLEst).append('\n');
			sb.append("Replacing Loss ML estimate:\t").append(replacing_lossMLEst).append('\n');
			sb.append("Additive Transfer Intragene Intraspecies ML estimate:\t").append(additive_trans_intragene_intraspecies_MLEst).append('\n');
			sb.append("Additive Transfer Intragene Interspecies ML estimate:\t").append(additive_trans_intragene_interspecies_MLEst).append('\n');
			sb.append("Additive Transfer Intergene Intraspecies ML estimate:\t").append(additive_trans_intergene_intraspecies_MLEst).append('\n');
			sb.append("Additive Transfer Intergene Interspecies ML estimate:\t").append(additive_trans_intergene_interspecies_MLEst).append('\n');
			sb.append("Replacing Transfer Intragene Intraspecies ML estimate:\t").append(replacing_trans_intragene_intraspecies_MLEst).append('\n');
			sb.append("Replacing Transfer Intragene Interspecies ML estimate:\t").append(replacing_trans_intragene_interspecies_MLEst).append('\n');
			sb.append("Replacing Transfer Intergene Intraspecies ML estimate:\t").append(replacing_trans_intergene_intraspecies_MLEst).append('\n');
			sb.append("Replacing Transfer Intergene Interspecies ML estimate:\t").append(replacing_trans_intergene_interspecies_MLEst).append('\n');
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
					sb.append(v.getName()).append('\t').append(guestNames.get(guestTrees.indexOf(v.guestTree)).get(v.sigma)).append('\n');
				}
			}
		}
		return sb.toString();
	}

	@Override
	public String getSigma(GuestVertex guestRoot) {
		StringBuilder sb = new StringBuilder(4096);
		sb.append("# GUEST-TO-HOST MAP\n");
		sb.append("Host tree:\t").append(guestRoot.guestTree.toString()).append('\n');
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
		return this.speciesTree.toString();
	}
	
	@Override
	public HashMap<RBTreeEpochDiscretiser, Integer> getUsed() {
		return this.used;
	}
}

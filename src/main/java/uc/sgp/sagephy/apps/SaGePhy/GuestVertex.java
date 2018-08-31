package uc.sgp.sagephy.apps.SaGePhy;

import java.util.LinkedList;
import uc.sgp.sagephy.io.NewickVertex;
import uc.sgp.sagephy.topology.Epoch;
import uc.sgp.sagephy.topology.RBTreeEpochDiscretiser;

public class GuestVertex extends NewickVertex {

	/** Event types. */
	public enum Event {
		SPECIATION,
		LEAF,             // Sampled leaf.
		UNSAMPLED_LEAF,   // Unsampled leaf.
		DUPLICATION,
		LOSS,
		REPLACING_LOSS,
		ADDITIVE_TRANSFER,
		REPLACING_TRANSFER,
		ADDITIVE_TRANSFER_INTRAGENE_INTRASPECIES,
		ADDITIVE_TRANSFER_INTRAGENE_INTERSPECIES,
		ADDITIVE_TRANSFER_INTERGENE_INTRASPECIES,
		ADDITIVE_TRANSFER_INTERGENE_INTERSPECIES,
		REPLACING_TRANSFER_INTRAGENE_INTRASPECIES,
		REPLACING_TRANSFER_INTRAGENE_INTERSPECIES,
		REPLACING_TRANSFER_INTERGENE_INTRASPECIES,
		REPLACING_TRANSFER_INTERGENE_INTERSPECIES,
		HYBRID_DONATION,
		HYBRID_DONATION_FROM_EXTINCT_DONOR,
		ALLOPLOIDIC_HYBRID_RECEPTION,  // This is only one of the lineages of the polyploidisation.
		AUTOPLOIDIC_HYBRID_RECEPTION,  // Obligate duplication due to polyploidisation.
	}

	/** Prunability status. */
	public enum Prunability {
		UNPRUNABLE,
		PRUNABLE,
		COLLAPSABLE,
		UNKNOWN
	}

	/** Type of event of vertex. */
	Event event;

	/** Absolute time. */
	double abstime;

	/** Host vertex/arc. */
	int sigma;

	/** Tranfered from arc */
	int transferedFromArc = -1;

	/** Tranfered to arc */
	int transferedToArc = -1;

	String transferedFromGuest = null;

	String transferedToGuest = null;

	/** Epoch. Not always applicable. */
	public Epoch epoch = null;

	/** Prunability status. */
	Prunability prunability = Prunability.UNKNOWN;

	/** Guest tree. */
	public RBTreeEpochDiscretiser guestTree;

	public String treeName;

	double branchtime;

	/**
	 * Host arcs that the arc (where this vertex is head) passes by.
	 * Only applicable for pruned trees.
	 */
	//List<Integer> hostArcs = null;

	/**
	 * Constructor.
	 * @param event type of event of this vertex.
	 * @param sigma enclosing host arc/vertex ID.
	 * @param epoch enclosing epoch.
	 * @param abstime absolute time of the this vertex.
	 * @param branchtime arc time of this vertex.
	 */
	public GuestVertex(Event event, int sigma, Epoch epoch, double abstime, double branchtime, RBTreeEpochDiscretiser guestTree) {
		super(-1, "", branchtime, "");
		this.event = event;
		this.sigma = sigma;
		this.host_vertex = sigma;
		this.epoch = epoch;
		this.abstime = abstime;
		this.guestTree = guestTree;
		this.treeName = guestTree.getRBTree().getName();
		this.guest_tree = guestTree.getRBTree().getName();
		//System.out.println(guestTree.getRBTree().getName());
		this.branchtime = branchtime;
	}

	/**
	 * Shallow copy constructor. References to parent or children are not included.
	 * @param orig the original.
	 */
	public GuestVertex(GuestVertex orig) {
		super(orig);
		this.event = orig.event;
		this.abstime = orig.abstime;
		this.sigma = orig.sigma;
		this.host_vertex = orig.sigma;
		this.epoch = orig.epoch;
		this.prunability = orig.prunability;
		this.guestTree = orig.guestTree;
		this.treeName = orig.treeName;
		this.guest_tree = orig.treeName;
		this.branchtime = orig.branchtime;
	}

	public void setTransferedToArc(int x){
		this.transferedToArc= x;
	}

	public int getTransferedToArc(){
		return this.transferedToArc;
	}

	public void setTransferedFromArc(int x){
		this.transferedFromArc= x;
	}

	public int getTransferedFromArc(){
		return this.transferedFromArc;
	}

	public void setTransferedToGuest(String x){
		this.transferedToGuest = x;
	}

	public String getTransferedToGuest(){
		return this.transferedToGuest;
	}

	public void setTransferedFromGuest(String x){
		this.transferedFromGuest= x;
	}

	public String getTransferedFromGuest(){
		return this.transferedFromGuest;
	}

	/**
	 * Helper. Returns the left child. If there is a single child, returns that one.
	 * @return left child.
	 */
	public GuestVertex getLeftChild() {
		if (this.getChildren() != null) {
		if (this.isLeaf()) { return null; }
			return (GuestVertex) this.getChildren().get(0);
		} else {
			return null;
		}
	}

	/**
	 * Helper. Returns the right child. If there is a single child, null is returned.
	 * @return right child.
	 */
	public GuestVertex getRightChild() {
		if (this.getChildren() != null) {
		if (this.getChildren().size() == 1) { return null; }
			return (GuestVertex) this.getChildren().get(1);
		} else {
			return null;
		}
	}

	/**
	 * Helper.
	 * @param root root.
	 */
	public static void setMeta(GuestVertex root, Boolean isSpecies) {
		LinkedList<NewickVertex> vertices = new LinkedList<NewickVertex>();
		vertices.add(root);
		while (!vertices.isEmpty()) {
			GuestVertex v = (GuestVertex) vertices.pop();
			StringBuilder sb = new StringBuilder(1024);
			//sb.append("[&&PRIME");
			sb.append("[ID=").append(v.getNumber());
			if (!isSpecies) {
				sb.append(" HOST=").append(v.getHostVertex());
				sb.append(" GUEST=").append(v.getGuestTree());
			}
			switch (v.event) {
			case DUPLICATION:
				double [] disTimes= v.epoch.getTimes();
				int j=0;
				while (true && j < disTimes.length){
					if (disTimes[j] >= v.abstime){
						break;
					}
					++j;
				}
				String dispt= "DISCPT=(" + v.epoch.getNo() + "," + j +")";
				if (isSpecies) {
					sb.append(" VERTEXTYPE=Speciation" + " "+ dispt);
				} else {
					sb.append(" VERTEXTYPE=Duplication" + " "+ dispt);
				}
				break;

			case LOSS:
				sb.append(" VERTEXTYPE=Loss");
				break;
			case REPLACING_LOSS:
				sb.append(" VERTEXTYPE=Replacing Loss");
				break;
			case ADDITIVE_TRANSFER:
				sb.append(" VERTEXTYPE=Additive Transfer");
				String fromToArc= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+")";

				double [] discTimes= v.epoch.getTimes();
				int i=0;
				while (true && i < discTimes.length){
					if (discTimes[i] >= v.abstime){
						break;
					}
					++i;
				}
				String discpt= "DISCPT=(" + v.epoch.getNo() + "," + i +")";

				//String speciesEdge= "SPECIES_EDGE=("+ v.getTransferedFromArc() +","+ v.epoch.getNoOfArcs() +")";
				//sb.append(" FROMTOLINEAGE="+ fromToArc +" "+ speciesEdge + " "+ discpt);
				sb.append(" FROMTOLINEAGE="+ fromToArc +" "+ discpt);
				break;
			case REPLACING_TRANSFER:
				sb.append(" VERTEXTYPE=Replacing Transfer");
				String fromToArc2= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+")";

				double [] discTimes2= v.epoch.getTimes();
				int r=0;
				while (true && r < discTimes2.length){
					if (discTimes2[r] >= v.abstime){
						break;
					}
					++r;
				}
				String discpt2= "DISCPT=(" + v.epoch.getNo() + "," + r +")";

				//String speciesEdge= "SPECIES_EDGE=("+ v.getTransferedFromArc() +","+ v.epoch.getNoOfArcs() +")";
				//sb.append(" FROMTOLINEAGE="+ fromToArc2 +" "+ speciesEdge + " "+ discpt2);
				sb.append(" FROMTOLINEAGE="+ fromToArc2 +" "+ discpt2);
				break;
			case ADDITIVE_TRANSFER_INTRAGENE_INTRASPECIES:
				sb.append(" VERTEXTYPE=Additive Transfer Intragene Intraspecies");
				String fromToArc1= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes1= v.epoch.getTimes();
				int i1=0;
				while (true && i1 < discTimes1.length){
					if (discTimes1[i1] >= v.abstime){
						break;
					}
					++i1;
				}
				String discpt1= "DISCPT=(" + v.epoch.getNo() + "," + i1 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc1 +" "+ discpt1);
				break;
			case ADDITIVE_TRANSFER_INTRAGENE_INTERSPECIES:
				sb.append(" VERTEXTYPE=Additive Transfer Intragene Interspecies");
				String fromToArc11= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes11= v.epoch.getTimes();
				int i11=0;
				while (true && i11 < discTimes11.length){
					if (discTimes11[i11] >= v.abstime){
						break;
					}
					++i11;
				}
				String discpt11= "DISCPT=(" + v.epoch.getNo() + "," + i11 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc11 +" "+ discpt11);
				break;
			case ADDITIVE_TRANSFER_INTERGENE_INTRASPECIES:
				sb.append(" VERTEXTYPE=Additive Transfer Intergene Intraspecies");
				String fromToArc111= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes111= v.epoch.getTimes();
				int i111=0;
				while (true && i111 < discTimes111.length){
					if (discTimes111[i111] >= v.abstime){
						break;
					}
					++i111;
				}
				String discpt111= "DISCPT=(" + v.epoch.getNo() + "," + i111 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc111 +" "+ discpt111);
				break;
			case ADDITIVE_TRANSFER_INTERGENE_INTERSPECIES:
				sb.append(" VERTEXTYPE=Additive Transfer Intergene Interspecies");
				String fromToArc1111= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes1111= v.epoch.getTimes();
				int i1111=0;
				while (true && i1111 < discTimes1111.length){
					if (discTimes1111[i1111] >= v.abstime){
						break;
					}
					++i1111;
				}
				String discpt1111= "DISCPT=(" + v.epoch.getNo() + "," + i1111 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc1111 +" "+ discpt1111);
				break;
			case REPLACING_TRANSFER_INTRAGENE_INTRASPECIES:
				sb.append(" VERTEXTYPE=Replacing Transfer Intragene Intraspecies");
				String fromToArc21= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes21= v.epoch.getTimes();
				int r1=0;
				while (true && r1 < discTimes21.length){
					if (discTimes21[r1] >= v.abstime){
						break;
					}
					++r1;
				}
				String discpt21= "DISCPT=(" + v.epoch.getNo() + "," + r1 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc21 +" "+ discpt21);
				break;
			case REPLACING_TRANSFER_INTRAGENE_INTERSPECIES:
				sb.append(" VERTEXTYPE=Replacing Transfer Intragene Interspecies");
				String fromToArc211= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes211= v.epoch.getTimes();
				int r11=0;
				while (true && r11 < discTimes211.length){
					if (discTimes211[r11] >= v.abstime){
						break;
					}
					++r11;
				}
				String discpt211= "DISCPT=(" + v.epoch.getNo() + "," + r11 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc211 +" "+ discpt211);
				break;
			case REPLACING_TRANSFER_INTERGENE_INTRASPECIES:
				sb.append(" VERTEXTYPE=Replacing Transfer Intergene Intraspecies");
				String fromToArc2111= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes2111= v.epoch.getTimes();
				int r111=0;
				while (true && r111 < discTimes2111.length){
					if (discTimes2111[r111] >= v.abstime){
						break;
					}
					++r111;
				}
				String discpt2111= "DISCPT=(" + v.epoch.getNo() + "," + r111 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc2111 +" "+ discpt2111);
				break;
			case REPLACING_TRANSFER_INTERGENE_INTERSPECIES:
				sb.append(" VERTEXTYPE=Replacing Transfer Intergene Interspecies");
				String fromToArc21111= "("+v.getTransferedFromArc()+","+ v.getTransferedToArc()+";"+v.getTransferedFromGuest()+","+v.getTransferedToGuest()+")";

				double [] discTimes21111= v.epoch.getTimes();
				int r1111=0;
				while (true && r1111 < discTimes21111.length){
					if (discTimes21111[r1111] >= v.abstime){
						break;
					}
					++r1111;
				}
				String discpt21111= "DISCPT=(" + v.epoch.getNo() + "," + r1111 +")";
				sb.append(" FROMTOLINEAGE="+ fromToArc21111 +" "+ discpt21111);
				break;
			case SPECIATION:
				double [] disctTimes= v.epoch.getTimes();
				int k=0;
				while (true && k < disctTimes.length){
					if (disctTimes[k] >= v.abstime){
						break;
					}
					++k;
				}
				String disctpt= "DISCPT=(" + v.epoch.getNo() + "," + k +")";
				sb.append(" VERTEXTYPE=Speciation"+ " "+ disctpt);
				break;
			case LEAF:
				sb.append(" VERTEXTYPE=Leaf");
				break;
			case UNSAMPLED_LEAF:
				sb.append(" VERTEXTYPE=UnsampledLeaf");
				break;
			case HYBRID_DONATION:
				sb.append(" VERTEXTYPE=HybridDonation");
				break;
			case HYBRID_DONATION_FROM_EXTINCT_DONOR:
				sb.append(" VERTEXTYPE=HybridDonationFromExtinctDonor");
				break;
			case ALLOPLOIDIC_HYBRID_RECEPTION:
				sb.append(" VERTEXTYPE=AlloploidicHybridReception");
				break;
			case AUTOPLOIDIC_HYBRID_RECEPTION:
				sb.append(" VERTEXTYPE=AutoploidicHybridReception");
				break;
			default:
				throw new UnsupportedOperationException("Invalid vertex event type.");
			}
			sb.append("]");
			v.setMeta(sb.toString());
			if (!v.isLeaf()) {
				vertices.addAll(v.getChildren());
			}
		}
	}
}

package uc.sgp.sagephy.apps.SaGePhy;

import java.util.HashMap;
import java.util.List;

import uc.sgp.sagephy.io.NewickIOException;
import uc.sgp.sagephy.math.PRNG;
import uc.sgp.sagephy.topology.RBTreeEpochDiscretiser;

/**
 * Interface of the all-mighty unpruned tree creators.
 * 
 * @author Joel Sjöstrand.
 */
public interface UnprunedGuestTreeCreator {
	
	public HashMap<RBTreeEpochDiscretiser, Integer> getUsed();

	/**
	 * Creates an unpruned tree. Appropriate labelling of vertices may not be in place in the returned tree.
	 * @param prng PRNG.
	 * @return the root.
	 * @throws NewickIOException 
	 */
	public GuestVertex createUnprunedTree(PRNG prng) throws NewickIOException;

	/**
	 * Returns the host leaves.
	 * @return the host leaves.
	 */
	public List<Integer> getHostLeaves();
	
	/**
	 * Returns the host as a string.
	 * @return the host.
	 */
	public String getHost();

	/**
	 * Creates auxiliary info on a tree (pruned or unpruned) created by the same instance.
	 * @param guestRoot guest tree root.
	 * @param doML include ML estimates.
	 * @return the info.
	 */
	public String getInfo(GuestVertex guestRoot, boolean doML);

	/**
	 * Creates a leaf map on a tree (pruned or unpruned) created by the same instance.
	 * @param guestRoot the guest root. May be null.
	 * @return the leaf map.
	 */
	public String getLeafMap(GuestVertex guestRoot);

	/**
	 * Creates a guest-to-host mapping info on a tree (pruned or unpruned) created by the same instance.
	 * @param guestRoot the guest root. May be null.
	 * @return the info.
	 */
	public String getSigma(GuestVertex guestRoot);
	
}

package uc.sgp.sagephy.apps.SaGePhy;

import uc.sgp.sagephy.mcmc.GenerativeModel;
import uc.sgp.sagephy.topology.DoubleMap;
import uc.sgp.sagephy.topology.NamesMap;
import uc.sgp.sagephy.topology.RootedTree;

/**
 * Interface for branch length relaxation models.
 * 
 * @author Joel Sjöstrand.
 */
public interface RateModel extends GenerativeModel {

	/**
	 * Returns rates.
	 * @param t the tree.
	 * @param names leaf/vertex names.
	 * @param origLengths the original lengths
	 * @return rates.
	 */
	public DoubleMap getRates(RootedTree t, NamesMap names, DoubleMap origLengths);
	
	/**
	 * Ultrametricity requirement on original lengths.
	 * @return true if required; false if arbitrary.
	 */
	public boolean lengthsMustBeUltrametric();
}

package uc.sgp.sagephy.apps.SaGePhy;

import java.util.HashMap;
import java.util.Map;

import uc.sgp.sagephy.math.ExponentialDistribution;
import uc.sgp.sagephy.math.PRNG;
import uc.sgp.sagephy.topology.DoubleMap;
import uc.sgp.sagephy.topology.NamesMap;
import uc.sgp.sagephy.topology.RootedTree;

/**
 * IID exponential rate model.
 * 
 * @author Joel Sjöstrand.
 */
public class IIDExponentialRateModel implements RateModel {

	/** Probability distribution. */
	private ExponentialDistribution pd;

	/** PRNG. */
	private PRNG prng;
	
	/**
	 * Constructor.
	 * @param lambda rate.
	 * @param prng PRNG.
	 */
	public IIDExponentialRateModel(double lambda, PRNG prng) {
		pd = new ExponentialDistribution(lambda);
		this.prng = prng;
	}
	
	@Override
	public Map<String, String> getModelParameters() {
		HashMap<String, String> kv =  new HashMap<String, String>(1);
		kv.put("lambda", ""+this.pd.getRate());
		return kv;
	}

	@Override
	public String getModelName() {
		return "IIDExponentialRates";
	}

	@Override
	public DoubleMap getRates(RootedTree t, NamesMap names, DoubleMap origLengths) {
		int n = t.getNoOfVertices();
		DoubleMap rates = new DoubleMap("Rates", n);
		for (int x = 0; x < n; ++x) {
			rates.set(x, pd.sampleValue(this.prng));
		}
		return rates;
	}

	@Override
	public boolean lengthsMustBeUltrametric() {
		return false;
	}

}

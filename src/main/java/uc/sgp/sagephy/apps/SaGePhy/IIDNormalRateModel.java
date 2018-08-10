package uc.sgp.sagephy.apps.SaGePhy;

import java.util.HashMap;
import java.util.Map;

import uc.sgp.sagephy.math.NormalDistribution;
import uc.sgp.sagephy.math.PRNG;
import uc.sgp.sagephy.topology.DoubleMap;
import uc.sgp.sagephy.topology.NamesMap;
import uc.sgp.sagephy.topology.RootedTree;

/**
 * IID normal rate model.
 * 
 * @author Joel Sjöstrand.
 */
public class IIDNormalRateModel implements RateModel {

	/** Probability distribution. */
	private NormalDistribution pd;

	/** PRNG. */
	private PRNG prng;
	
	/**
	 * Constructor.
	 * @param mu mu.
	 * @param sigma2 sigma2.
	 * @param prng PRNG.
	 */
	public IIDNormalRateModel(double mu, double sigma2, PRNG prng) {
		pd = new NormalDistribution(mu, sigma2);
		this.prng = prng;
	}
	
	@Override
	public Map<String, String> getModelParameters() {
		HashMap<String, String> kv =  new HashMap<String, String>(2);
		kv.put("mu", ""+this.pd.getMean());
		kv.put("sigma2", ""+this.pd.getVariance());
		return kv;
	}

	@Override
	public String getModelName() {
		return "IIDNormalRates";
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

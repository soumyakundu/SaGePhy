package uc.sgp.sagephy.apps.SaGePhy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import uc.sgp.sagephy.io.GMLFactory;
import uc.sgp.sagephy.io.GMLGraph;
import uc.sgp.sagephy.io.GMLIOException;
import uc.sgp.sagephy.io.GMLReader;
import uc.sgp.sagephy.io.NewickIOException;
import uc.sgp.sagephy.io.PrIMENewickTree;
import uc.sgp.sagephy.io.PrIMENewickTreeReader;
import uc.sgp.sagephy.topology.TopologyException;

import com.beust.jcommander.Parameter;


/**
 * Contains user settings.
 *
 * @author Soumya Kundu.
 */
public class DomainTreeGenParameters {
	
	/** Required parameters. */
	@Parameter(description = "<host tree file or string> <guest tree directory> <dup rate> <loss rate> <trans rate> <out prefix>")
	public List<String> args = new ArrayList<String>();

	/** Help. */
	@Parameter(names = {"-h", "--help"}, description = "Display help.")
	public Boolean help = false;

	/** Output. */
	@Parameter(names = {"-q", "--quiet"}, description = "Suppress creation of auxiliary files and write only pruned tree directly to stdout.")
	public Boolean doQuiet = false;

	/** PRNG seed. */
	@Parameter(names = {"-s", "--seed"}, description = "PRNG seed. Default: Random seed.")
	public String seed = null;

	/** Min leaves. */
	@Parameter(names = {"-min", "--min-leaves"}, description = "Minimum number of extant guest leaves required.")
	public Integer min = 2;

	/** Max leaves. */
	@Parameter(names = {"-max", "--max-leaves"}, description = "Maximum number of extant guest leaves allowed.")
	public Integer max = 64;

	/** Min leaves per host tree leaf. */
	@Parameter(names = {"-minper", "--min-leaves-per-host-leaf"}, description = "Minimum number of extant guest leaves per host leaf required.")
	public Integer minper = 0;

	/** Max leaves per host tree leaf. */
	@Parameter(names = {"-maxper", "--max-leaves-per-host-leaf"}, description = "Maximum number of extant guest leaves per host leaf allowed.")
	public Integer maxper = 10;

	/** Leaf size samples. */
	@Parameter(names = {"-sizes", "--leaf-sizes-file"}, description = "Samples the desired number of extant leaves uniformly from a single-column file. This is suitable for mimicking a known leaf size distribution. Default: No sampling.")
	public String leafSizes = null;

	/** Do meta. */
	@Parameter(names = {"-nox", "--no-auxiliary-tags"}, description = "Exclude auxiliary PrIME tags in output trees and suppress creation of files for the mappings.")
	public Boolean excludeMeta = false;

	/** Leaf sampling. */
	@Parameter(names = {"-p", "--leaf-sampling-probability"}, description = "Governs the probability of observing a guest tree leaf. Lineages that fail to be observed will be pruned away similarly to lineages lost during the evolutionary process.")
	public String leafSamplingProb = "1.0";

	/** Enforce parsimony. */
	@Parameter(names = {"-mpr", "--enforce-most-parsimonious-reconciliation"}, description = "Requires the pruned guest tree to not have any duplication vertex that occurs on a host edge other than" +
			" what a most-parsimonious reconciliation of the pruned guest tree topology would suggest, and that cannot be mistaken for a speciation.")
	public Boolean enforceParsimony = false;

	/** Attempts. */
	@Parameter(names = {"-a", "--max-attempts"}, description = "Maximum number of attempts at creating random tree that meets requirements. If not met, no tree is output.")
	public Integer maxAttempts = 10000;

	/** Stem edge. */
	@Parameter(names = {"-stem", "--override-host-stem"}, description = "If set, overrides the stem edge of the host tree by the specified value. If no stem edge is desired, this can be set to 0. Default: Value in host tree.")
	public String stem = null;

	/** Vertex prefix. */
	@Parameter(names = {"-vp", "--vertex-prefix"}, description = "Vertex prefix.")
	public String vertexPrefix = "D";

	/** Sigma. */
	@Parameter(names = {"-vph", "--vertex-prefix-host-map"}, description = "Append host vertex/edge belonging to vertex prefix.")
	public Boolean appendSigma = true;

	/** Hybrid network. */
	@Parameter(names = {"-hybrid", "--hybrid-host-graph"}, arity = 3, description = "<post-hyb timespan> <post-hyb dup fact> <post-hyb loss fact>. Assumes that the input host tree is in fact a hybrid DAG." +
			" The transfer rate parameter will be ignored (i.e., set to 0). The additional parameters refer to a change-factor applied to the parameters for a limited time following a hybrid speciation.")
	public List<String> hybrid = null;

	/** Rate of replacing transfers. */
	@Parameter(names = {"-rt", "--replacing-transfers"}, description = "Probability of replacing horizontal gene transfers.")
	public String replacing_transfers = "0.5";

	/** Type of distance bias for transfers. */
	@Parameter(names = {"-db", "--distance-bias"}, description = "Type of distance-bias for horizontal gene transfers. Options: none, simple, exponential.")
	public String distance_bias = "none";
	
	/** Coefficient for transfer distance-bias. */
	@Parameter(names = {"-dbr", "--distance-bias-rate"}, description = "Set level of bias towards recipients that are phylogenetically closer.")
	public String distance_bias_rate = "1.0";
	
	@Parameter(names = {"-dbs", "--domain-birth-sampling"}, description = "Randomly sample location of domain birth on gene tree.")
	public Boolean doDomainBirth = false;
	
	/** Coefficient for sampling domain tree birth location. */
	@Parameter(names = {"-dbc", "--domain-birth-coefficient"}, description = "Set level of bias towards root of gene tree for domain tree birth location.")
	public String domain_birth = "1.0";
	
	@Parameter(names = {"-ig", "--inter-gene-transfers"}, description = "Probability of transfers across gene trees.")
	public String inter_gene = "0.5";
	
	@Parameter(names = {"-is", "--inter-species-transfers"}, description = "Probability of transfers across species.")
	public String inter_species = "0.5";
	
	@Parameter(names = {"-all", "--all-gene-trees"}, description = "Enforce the evolution of the domain tree in all gene trees.")
	public Boolean all_genes = false;
	
	/**
	 * Returns output and info streams.
	 * @return streams.
	 */
	public BufferedWriter getOutputFile(String suffix) {
		try {
			String outprefix = args.get(5).trim();
			return new BufferedWriter(new FileWriter(outprefix + suffix));
		} catch (Exception e) {
			throw new IllegalArgumentException("Invalid output file prefix.", e);
		}
	}

	public UnprunedGuestTreeCreator getHostTreeCreator() throws TopologyException, NewickIOException, IOException {
		PrIMENewickTree host;
		File f = new File(this.args.get(0));
		if (f.exists()) {
			host = PrIMENewickTreeReader.readTree(f, false, true);
		} else {
			host = PrIMENewickTreeReader.readTree(args.get(0), false, true, null);
		}
		ArrayList<PrIMENewickTree> guests = new ArrayList<PrIMENewickTree>();
		File dir = new File(this.args.get(1));
		File[] genes = dir.listFiles();	
		for (File j : genes) {
			guests.add(PrIMENewickTreeReader.readTree(j, false, true));
		}
		return new DomainTreeInHostTreeCreator(host, guests, this.getDuplicationRate(), this.getLossRate(), this.getTransferRate(), this.getReplacingTransferRate(), this.distance_bias, this.getDistanceBiasRate(), this.doDomainBirth, this.getDomainBirthSampling(), false, this.getInterGeneTransfers(), this.getInterSpeciesTransfers(), this.getLeafSamplingProb(), this.getStem());
	}

	public double getDuplicationRate() {
		return Double.parseDouble(this.args.get(2));
	}

	public double getLossRate() {
		return Double.parseDouble(this.args.get(3));
	}

	public double getTransferRate() {
		return Double.parseDouble(this.args.get(4));
	}

	public ArrayList<Integer> getLeafSizes() throws FileNotFoundException {
		// Read file.
		if (this.leafSizes == null) { return null; }
		ArrayList<Integer> samples = new ArrayList<Integer>(8192);
		Scanner sc = new Scanner(new File(this.leafSizes));
		while (sc.hasNextLine()) {
			samples.add(Integer.parseInt(sc.nextLine()));
		}
		sc.close();
		return samples;
	}

	public Double getStem() {
		return (this.stem == null ? 0.0 : Double.parseDouble(this.stem));
	}

	public Double getLeafSamplingProb() {
		return Double.parseDouble(this.leafSamplingProb);
	}

	public Double getReplacingTransferRate() {
		return Double.parseDouble(this.replacing_transfers);
	}
	
	public Double getDistanceBiasRate() {
		return Double.parseDouble(this.distance_bias_rate);
	}
	
	public Double getDomainBirthSampling() {
		return Double.parseDouble(this.domain_birth);
	}
	
	public Double getInterGeneTransfers() {
		return Double.parseDouble(this.inter_gene);
	}
	
	public Double getInterSpeciesTransfers() {
		return Double.parseDouble(this.inter_species);
	}
	
	public GuestTreeInHybridGraphCreator getHostHybridGraphCreator() throws GMLIOException, IOException {
		String str = this.args.get(0);
		File f = new File(str);
		double dup = this.getDuplicationRate();
		double loss = this.getLossRate();
		double timespan = Double.parseDouble(this.hybrid.get(0));
		double dupfact = Double.parseDouble(this.hybrid.get(1));
		double lossfact = Double.parseDouble(this.hybrid.get(2));

		GMLGraph host;
		if (f.exists()) {
			host = GMLFactory.getGraphs(GMLReader.readGML(f)).get(0);
		} else {
			host = GMLFactory.getGraphs(GMLReader.readGML(str)).get(0);
		}
		return new GuestTreeInHybridGraphCreator(host, dup, loss, dupfact, lossfact, timespan, this.getLeafSamplingProb(), this.getStem());
	}

}



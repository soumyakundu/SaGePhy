package uc.sgp.sagephy.apps.SaGePhy;

import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.LinkedList;

import com.beust.jcommander.JCommander;

import uc.sgp.sagephy.apps.SaGePhyApp;
import uc.sgp.sagephy.io.JCommanderUsageWrapper;
import uc.sgp.sagephy.io.NewickTreeWriter;
import uc.sgp.sagephy.io.PrIMENewickTree;
import uc.sgp.sagephy.io.PrIMENewickTreeReader;
import uc.sgp.sagephy.math.NumberManipulation;
import uc.sgp.sagephy.misc.Pair;

public class HostTreeGen implements SaGePhyApp {

	@Override
	public String getAppName() {
		return "HostTreeGen";
	}

	@Override
	public void main(String[] args) throws Exception {

		try {
			// ================ PARSE USER OPTIONS AND ARGUMENTS ================
			HostTreeGenParameters params = new HostTreeGenParameters();
			JCommander jc = new JCommander(params, args);
			int noargs = params.doQuiet ? 3 : 4;
			if (args.length == 0 || params.help || params.args.size() != noargs) {
				StringBuilder sb = new StringBuilder(65536);
				sb.append("Usage:\n" +
						"    java -jar sagephy-X.Y.Z.jar HostTreeGen [options] <time interval> <birth rate> <death rate> <out prefix>\n");
				JCommanderUsageWrapper.getUnsortedUsage(jc, params, sb);
				System.out.println(sb.toString());
				return;
			}

			// Create hostTree.
			double T = params.getT();
			String S;
			if (params.bifurcationStart) {
				// Bifurcating host tree.
				S = "(H0:" + T + ",H1:" + T + "):0.0;";
			} else {
				S = "H0:" + T + ";";
			}
			PrIMENewickTree host = PrIMENewickTreeReader.readTree(S, false, true, null);

			// Machine.
			int minper = (params.bifurcationStart ? 1 : 0);
			GuestTreeMachina machine = new GuestTreeMachina(params.seed, params.min, params.max, minper, Integer.MAX_VALUE, params.getLeafSizes(), params.maxAttempts,
						params.vertexPrefix, params.excludeMeta, false, true);

			// Primus motor.
			GuestTreeInHostTreeCreator motor = new GuestTreeInHostTreeCreator(host, params.getBirthRate(), params.getDeathRate(), 0.0, 0.0, "", 1.0, false, 0.0, true, params.getLeafSamplingProb(), null);

			Pair<PrIMENewickTree, PrIMENewickTree> guestTree = null;
			try {
				guestTree = machine.sampleGuestTree(motor, null);
			} catch (MaxAttemptsException ex) {
				if (!params.doQuiet) {
					BufferedWriter outinfo = params.getOutputFile(".info");
					outinfo.write("Failed to produce valid pruned tree within max allowed attempts.\n");
					outinfo.close();
				}
				System.err.println("Failed to produce valid pruned tree within max allowed attempts.");
				System.exit(0);
			}

			// Fix stem override.
			if (params.stem != null) {
				guestTree.second.getRoot().setBranchLength(Double.parseDouble(params.stem));
				if (guestTree.first != null) {
					guestTree.first.getRoot().setBranchLength(Double.parseDouble(params.stem));
				}
			}

			// Print output.
			if (params.doQuiet) {
				if (params.excludeMeta) {
					System.out.println(guestTree.first == null ? ";" : NewickTreeWriter.write(guestTree.first));
				} else {
					System.out.println(guestTree.first == null ? "[NAME=PrunedTree];" : NewickTreeWriter.write(guestTree.first));
				}
			} else {
				BufferedWriter out = params.getOutputFile(".unpruned.tree");
				out.write(NewickTreeWriter.write(guestTree.second) + '\n');
				out.close();
				out = params.getOutputFile(".pruned.tree");
				if (params.excludeMeta) {
					out.write(guestTree.first == null ? ";\n" : NewickTreeWriter.write(guestTree.first) + '\n');
				} else {
					out.write(guestTree.first == null ? "[NAME=PrunedTree];" : NewickTreeWriter.write(guestTree.first) + '\n');
				}
				out.close();
				out = params.getOutputFile(".unpruned.info");
				out.write("# HOSTTREEGEN\n");
				out.write("Arguments:\t" +  Arrays.toString(args) + '\n');
				out.write("Attempts:\t" + machine.getAttempts() + '\n');
				out.write(this.getInfo(guestTree.second, true));
				out.close();
				out = params.getOutputFile(".pruned.info");
				out.write("# HOSTTREEGEN\n");
				out.write("Arguments:\t" +  Arrays.toString(args) + '\n');
				out.write("Attempts:\t" + machine.getAttempts() + '\n');
				out.write(this.getInfo(guestTree.first, false));
				out.close();
			}
		} catch (Exception e) {
			e.printStackTrace(System.err);
			System.err.print("\nUse option -h or --help to show usage.\n");
		}
	}

	/**
	 * Creates info.
	 * @param tree tree.
	 * @param doML include ML estimates.
	 * @return the info.
	 */
	private String getInfo(PrIMENewickTree tree, boolean doML) {
		StringBuilder sb = new StringBuilder(1024);
		int noOfVertices = 0;
		int noOfBirths = 0;
		int noOfDeaths = 0;
		int noOfLeaves = 0;
		double totalTime = 0.0;

		LinkedList<GuestVertex> vertices = new LinkedList<GuestVertex>();
		if (tree != null) {
			vertices.add((GuestVertex) tree.getRoot());
		}
		while (!vertices.isEmpty()) {
			GuestVertex v = vertices.pop();
			if (!v.isLeaf()) {
				vertices.add(v.getLeftChild());
				vertices.add(v.getRightChild());
			}
			noOfVertices++;
			totalTime += v.getBranchLength();
			switch (v.event) {
			case DUPLICATION:
			case SPECIATION:   // For the forced bifurcation case.
				noOfBirths++;
				break;
			case LOSS:
				noOfDeaths++;
				break;
			case LEAF:
			case UNSAMPLED_LEAF:
				noOfLeaves++;
				break;
			default:
				break;
			}
		}
		totalTime = NumberManipulation.roundToSignificantFigures(totalTime, 8);

		sb.append("No. of vertices:\t").append(noOfVertices).append('\n');
		sb.append("No. of extant leaves:\t").append(noOfLeaves).append('\n');
		sb.append("No. of births:\t").append(noOfBirths).append('\n');
		sb.append("No. of deaths:\t").append(noOfDeaths).append('\n');
		sb.append("Total branch time:\t").append(totalTime).append('\n');
		if (doML) {
			double birthMLEst = NumberManipulation.roundToSignificantFigures(noOfBirths / totalTime, 8);
			double deathMLEst = NumberManipulation.roundToSignificantFigures(noOfDeaths / totalTime, 8);
			sb.append("Birth ML estimate:\t").append(birthMLEst).append('\n');
			sb.append("Death ML estimate:\t").append(deathMLEst).append('\n');
		}
		return sb.toString();
	}
}

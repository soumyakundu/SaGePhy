package uc.sgp.sagephy.apps.SaGePhy;

import java.io.BufferedWriter;
import java.util.Arrays;
import com.beust.jcommander.JCommander;

import uc.sgp.sagephy.apps.SaGePhyApp;
import uc.sgp.sagephy.io.JCommanderUsageWrapper;
import uc.sgp.sagephy.io.NewickTreeWriter;
import uc.sgp.sagephy.io.PrIMENewickTree;
import uc.sgp.sagephy.misc.Pair;

public class GuestTreeGen implements SaGePhyApp {

	@Override
	public String getAppName() {
		return "GuestTreeGen";
	}

	@Override
	public void main(String[] args) throws Exception {

		try {
			// ================ PARSE USER OPTIONS AND ARGUMENTS ================
			GuestTreeGenParameters params = new GuestTreeGenParameters();
			JCommander jc = new JCommander(params, args);
			int noargs = params.doQuiet ? 4 : 5;
			if (args.length == 0 || params.help || params.args.size() != noargs) {
				StringBuilder sb = new StringBuilder(65536);
				sb.append("Usage:\n" +
						"    java -jar sagephy-X.Y.Z.jar GuestTreeGen [options] <host tree> <dup rate> <loss rate> <trans rate> <out prefix>\n");
				JCommanderUsageWrapper.getUnsortedUsage(jc, params, sb);
				System.out.println(sb.toString());
				return;
			}

			// Machine.
			GuestTreeMachina machina = new GuestTreeMachina(params.seed, params.min, params.max, params.minper, params.maxper, params.getLeafSizes(), params.maxAttempts,
					params.vertexPrefix, params.excludeMeta, params.appendSigma, false);

			// Machine motor.
			UnprunedGuestTreeCreator motor = (params.hybrid == null || params.hybrid.isEmpty()) ? params.getHostTreeCreator() : params.getHostHybridGraphCreator();

			// Create guest tree.
			Pair<PrIMENewickTree, PrIMENewickTree> guestTree = null;
			try {
				guestTree  = machina.sampleGuestTree(motor, null);
			} catch (MaxAttemptsException ex) {
				if (!params.doQuiet) {
					BufferedWriter outinfo = params.getOutputFile(".pruned.info");
					outinfo.write("Failed to produce valid pruned tree within max allowed attempts.\n");
					outinfo.close();
				}
				System.err.println("Failed to produce valid pruned tree within max allowed attempts.");
				System.exit(0);
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
					out.write(guestTree.first == null ? "[NAME=PrunedTree];\n" : NewickTreeWriter.write(guestTree.first) + '\n');
				}
				out.close();
				out = params.getOutputFile(".unpruned.info");
				out.write("# GUESTTREEGEN\n");
				out.write("Arguments:\t" +  Arrays.toString(args) + '\n');
				out.write("Attempts:\t" + machina.getAttempts() + '\n');
				out.write(motor.getInfo((GuestVertex) guestTree.second.getRoot(), true));
				out.close();
				out = params.getOutputFile(".pruned.info");
				out.write("# GUESTTREEGEN\n");
				out.write("Arguments:\t" +  Arrays.toString(args) + '\n');
				out.write("Attempts:\t" + machina.getAttempts() + '\n');
				out.write(motor.getInfo((GuestVertex) guestTree.first.getRoot(), true));
				out.close();
				if (!params.excludeMeta) {
					out = params.getOutputFile(".unpruned.guest2host");
					out.write(motor.getSigma((GuestVertex) guestTree.second.getRoot()));
					out.close();
					out = params.getOutputFile(".pruned.guest2host");
					out.write(motor.getSigma((GuestVertex) guestTree.first.getRoot()));
					out.close();
					out = params.getOutputFile(".unpruned.leafmap");
					out.write(motor.getLeafMap((GuestVertex) guestTree.second.getRoot()));
					out.close();
					out = params.getOutputFile(".pruned.leafmap");
					out.write(motor.getLeafMap((GuestVertex) guestTree.first.getRoot()));
					out.close();
				}
			}
		} catch (Exception e) {
			e.printStackTrace(System.err);
			System.err.print("\nUse option -h or --help to show usage.\n");
		}
	}
}

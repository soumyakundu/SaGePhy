package uc.sgp.sagephy.apps;

import java.util.ArrayList;
import java.util.Set;
import java.util.TreeMap;

import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;
import org.reflections.util.FilterBuilder;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;

/**
 * Starter class for all apps. Uses reflection to search for main methods in the
 * <code>uc.sgp.sagephy.apps</code> package, and uses the corresponding class name
 * as a program name identifier. Apps are thus started as an ordinary executable:
 * <code>./jprime-X.Y.Z.jar &lt;application&gt; [options] &lt;arguments&gt;</code> etc.,
 * or alternatively: <code>java -jar jprime-X.Y.Z.jar &lt;application&gt; [options] &lt;arguments&gt;</code> etc.
 * 
 * @author Joel Sjöstrand.
 */
public class SaGePhyStarter {

	/**
	 * Starts a SaGePhy application located in the <code>uc.sgp.sagephy.apps</code> folder
	 * (or sub-folder).
	 * @param args program name, followed by arguments (and necessarily in that order).
	 * @throws Exception.
	 */
	public static void main(String[] args) throws Exception {
		
		// Find all starter methods in apps package.
		Logger logger = (Logger) Reflections.log;
		Level oldLvl = logger.getLevel();
		logger.setLevel(Level.OFF);
		Reflections reflections = new Reflections(
			    new ConfigurationBuilder()
			        .filterInputsBy(new FilterBuilder().include(FilterBuilder.prefix("uc.sgp.sagephy.apps")))
			        .setUrls(ClasspathHelper.forPackage("uc.sgp.sagephy.apps"))
			        .setScanners(new SubTypesScanner())
			        
			);
		Set<Class<? extends SaGePhyApp>> apps = reflections.getSubTypesOf(SaGePhyApp.class);
		logger.setLevel(oldLvl);
		
		// Create map linking name and app used for launching apps (we assume uniqueness).
		TreeMap<String, SaGePhyApp> map = new TreeMap<String, SaGePhyApp>();
		// A list of application names used to list available apps
		ArrayList<String> appNames = new ArrayList<>();
		for (Class<? extends SaGePhyApp> c : apps) {
			try {
				// NOTE: Empty constructor assumed!!!!
				SaGePhyApp app = c.newInstance();
				String name = app.getAppName();
				// app can have multiple aliases
				String[] names = name.split("\\s+");
				for (String n : names) {
					map.put(n, app);
				}
				String listedAppName = names[0];
				for (int i = 1; i < names.length; i++) {
					listedAppName += " (" + names[i] + ")";
				}
				appNames.add(listedAppName);
			} catch (InstantiationException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		
		if (args.length < 1 || args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("--help")) {
			// No app specified.
			System.out.println(
					"================================================================================\n" +
					"SaGePhy is a Java phylogenetic simulation framework for gene and sub-gene evolution\n" +
					"developed at the University of Connecticut Computational Biology Laboratory\n" +
					"through modification of the GenPhyloData simulation framework.\n\n" +
					"Releases, source code and tutorial: https://github.com/soumyakundu/simulation_framework\n\n" +
					"License: SaGePhy is available under the New BSD License.\n"
					);
			System.out.println("Usage: sagephy-x.y.z.jar <application> [options] <arguments>\n");
			System.out.println("List of available applications:");
			for (String k : appNames) {
				System.out.println("    " + k);
			}
			System.out.println("You can usually obtain app-specific help thus:\n" +
					"    java -jar sagephy-x.y.z.jar <application> -h");
			System.out.println("================================================================================\n");
		} else if (!map.containsKey(args[0])) {
			System.out.println("Unknown application. Use -h to show help.");
			System.out.println("Usage: java -jar sagephy-x.y.z.jar <application> [options] <arguments>");
		} else {
			// Start app. Remove app name first, though.
			String[] appArgs = new String[args.length - 1];
			System.arraycopy(args, 1, appArgs, 0, appArgs.length);
			SaGePhyApp app = map.get(args[0]);
			app.main(appArgs);
		}
	}

}

package uc.sgp.sagephy.apps;

/**
 * Interface for SaGePhy apps. Considered using annotations too, but that seemed less practical.
 * We assume empty constructors for apps at the moment.
 * 
 * @author Joel Sjöstrand.
 */
public interface SaGePhyApp {
	
	/**
	 * Returns the app's name. This can be a string of several names separated by whitespace.
	 * These names will then be the accepted command line parameters for launching the App.
	 * @return the name(s).
	 */
	public String getAppName();
	
	/**
	 * Starts the application.
	 * @param args the application arguments.
	 */
	public void main(String[] args) throws Exception;
	
}

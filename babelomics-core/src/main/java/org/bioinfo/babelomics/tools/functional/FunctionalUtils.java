package org.bioinfo.babelomics.tools.functional;

public class FunctionalUtils {

	public static boolean isEnsemblID(String id) {
		if ( id != null && id.startsWith("ENSG") && id.length() == 15 ) {
			return true;
		}
		return false;
	}
}

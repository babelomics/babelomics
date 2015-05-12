package org.bioinfo.babelomics.utils.filters;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;

/**
 * Created by ralonso on 1/2/15.
 */
public class ReconFilter extends FunctionalFilter {

    public ReconFilter() {
        this(2, 500);
    }

    public ReconFilter(int minNumberGenes, int maxNumberGenes) {
        super(minNumberGenes, maxNumberGenes);
    }

    public String getSQLWhereClause(String prefixSqlField) {
        return "";
    }
}


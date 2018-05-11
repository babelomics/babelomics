# Babelomics

## Babelomics is an integrative platform for the analysis of Transcriptomics, Proteomics and Genomics data with advanced functional profiling

This version of Babelomics integrates primary (normalization, calls, etc.) and secondary (signatures, predictors, associations, TDTs, clustering, etc.) analysis tools within an environment that allows relating genomic data and/or interpreting them by means of different functional enrichment or gene set methods. Such interpretation is made using functional definitions, protein-protein interactions...

### How to install

        git clone https://github.com/babelomics/babelomics.git
        cd babelomics/babelomics-web/
        bower install
        git submodule update --init
        cd lib/jsorolla/
        bower install

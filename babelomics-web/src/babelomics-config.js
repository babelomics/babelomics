CELLBASE_HOST = "http://www.ebi.ac.uk/cellbase/webservices/rest";
CELLBASE_VERSION = "v3";
// OPENCGA_HOST = "http://ws-beta.bioinfo.cipf.es/opencga-staging/rest";
OPENCGA_HOST = "http://test.babelomics.org/opencga/rest";
OPENCGA_VERSION = "v1";

if (
    window.location.host.indexOf("localhost") != -1 ||
    window.location.host.indexOf("fsalavert") != -1 ||
    window.location.host.indexOf("rsanchez") != -1 ||
    window.location.host.indexOf("imedina") != -1 ||
    window.location.host.indexOf("aaleman") != -1 ||
    window.location.protocol === "file:"
    ) {

    CELLBASE_HOST = "http://www.ebi.ac.uk/cellbase/webservices/rest";
    CELLBASE_VERSION = "v3";
    OPENCGA_HOST = "http://test.babelomics.org/opencga/rest";
    //OPENCGA_HOST = "http://mem18:8080/opencga/rest";
    OPENCGA_VERSION = "v1";
//    OPENCGA_HOST = "http://localhost:8080/opencga/rest";
}

var AVAILABLE_SPECIES = {
    "text": "Species",
    "items": [
        {
            "text": "Vertebrates",
            "items": [
                {
                    "text": "Homo sapiens",
                    "assembly": "GRCh37.p10",
                    "region": {"chromosome": "13", "start": 32889611, "end": 32889611},
                    "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"],
                    "url": "ftp://ftp.ensembl.org/pub/release-75/"
                },
            ]
        }
    ]
};

/** Reference to a species from the list to be shown at start **/
var DEFAULT_SPECIES = AVAILABLE_SPECIES.items[0].items[0];
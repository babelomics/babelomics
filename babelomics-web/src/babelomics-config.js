CELLBASE_HOST = "http://wwwdev.ebi.ac.uk/cellbase/webservices/rest";
CELLBASE_VERSION = "v3";
// OPENCGA_HOST = "http://ws-beta.bioinfo.cipf.es/opencga-staging/rest";
OPENCGA_HOST = "http://ws.babelomics.org/opencga/rest";
//OPENCGA_HOST = "http://web02:8080/opencga/rest";
OPENCGA_VERSION = "v1";

//if (
//    window.location.host.indexOf("localhost") != -1 ||
//    window.location.host.indexOf("fsalavert") != -1 ||
//    window.location.host.indexOf("rsanchez") != -1 ||
//    window.location.host.indexOf("imedina") != -1 ||
//    window.location.host.indexOf("aaleman") != -1 ||
//    window.location.protocol === "file:"
//    ) {
//
//    CELLBASE_HOST = "http://wwwdev.ebi.ac.uk/cellbase/webservices/rest";
//    CELLBASE_VERSION = "v3";
//    OPENCGA_HOST = "http://test.babelomics.org/opencga/rest";
//    //OPENCGA_HOST = "http://mem18:8080/opencga/rest";
//    OPENCGA_VERSION = "v1";
////    OPENCGA_HOST = "http://localhost:8080/opencga/rest";
//}


/** Config file for genome maps**/
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

/** available species for protein viewer **/
var SPECIES = [
    {name: "Homo sapiens", value: "hsapiens", dbname: "hgnc_symbol"},
    {name: "Mus musculus", value: "mmusculus", dbname: "mgi_symbol"},
    {name: "Rattus norvegicus", value: "rnorvegicus", dbname: "EntrezGene"},
    {name: "Drosophila melanogaster", value: "dmelanogaster", dbname: "EntrezGene"},
    {name: "Caenorhabditis elegans", value: "celegans", dbname: "sgd"},
    {name: "Saccharomyces cerevisiae", value: "scerevisiae", dbname: "sgd"},
    {name: "Danio rerio", value: "drerio", dbname: "EntrezGene"},
    {name: "Arabidopsis thaliana", value: "athaliana", dbname: "EntrezGene"}
];

STUDY_NAME = "WorkSpace";

var TOOLS = ["affy-expression-normalization", "agilent-expression-one-color-normalization", "agilent-expression-two-colors-normalization", "association", "burden",
    "class-comparison", "class-prediction", "clustering", "communities-structure-detection", "correlation", "fatigo", "fatiscan", "genepix-expression-one-color-normalization",
    "genepix-expression-two-colors-normalization", "network-miner", "oncodriveclust", "oncodrivefm", "preprocessing", "rnaseq-diffexpr", "rnaseq-norm", "snow", "stratification",
    "survival"];

var MAINTENANCE_MODE = false;

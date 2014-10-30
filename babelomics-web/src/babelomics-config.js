CELLBASE_HOST = "http://www.ebi.ac.uk/cellbase/webservices/rest";
CELLBASE_VERSION = "v3";
OPENCGA_HOST = "http://ws-beta.bioinfo.cipf.es/opencga-staging/rest";
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
    OPENCGA_VERSION = "v1";
    OPENCGA_HOST = "http://localhost:8080/opencga/rest";
}
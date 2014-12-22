/**
 * Created by ralonso on 11/28/14.
 */

function ParseMatrix(contentFile) {
    this.contentFile = contentFile;
}

ParseMatrix.prototype.getDataset = function () {
    var datasetHeader = {
        "#NUMBER_FEATURES": "",
        "#NUMBER_SAMPLES": "",
        "#VARIABLE": new Array(),
        "#NAMES": new Array()
    };
    var lines = this.contentFile.split("\n");
    var featNumber = 0;
    var firstLine = 0;
    for (var i = 0; i < lines.length; i++) {
        var line = lines[i];
        if (line.trim() == "") {
            firstLine++;
            continue
        }
        if (line.indexOf("#") == 0) {
            var fields = line.split("\t");
            if (fields[0] == "#VARIABLE") {
                var arr = datasetHeader[fields[0]]
                arr.push(fields);
                datasetHeader[fields[0]] = arr;
            }
        }
        else {
            featNumber++;
            /** Samples **/
            var fields = new Array();
            if (datasetHeader["#NAMES"].length == 0) {
                if (i > firstLine) {
                    fields = lines[i - 1].split("\t");
                    datasetHeader["#NAMES"] = fields.slice(1);
                }
                else {
                    fields = lines[i].split("\t");
                    var names = new Array();
                    for (var z = 1; z < fields.length; z++) {
                        names.push("Sample" + z);
                    }
                    datasetHeader["#NAMES"] = names;
                }
            }
        }
    }
    datasetHeader["#NUMBER_FEATURES"] = featNumber;
    datasetHeader["#NUMBER_SAMPLES"] = datasetHeader["#NAMES"].length;
    return datasetHeader;
};
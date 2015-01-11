/**
 * Created by ralonso on 11/28/14.
 */

function FormUtils(scope) {
    this.scope = scope;
}

FormUtils.prototype.basicValidationForm = function () {
    var validated = true;
    var msg = "";
    if( this.scope.$.outdir.selectedFile === undefined ||  this.scope.$.outdir.selectedFile.type != "FOLDER"){
        msg += "Error: Please select an output folder.\n";
        validated=false;
    }
    if( this.scope.$.inputFile.selectedFile === undefined ||  this.scope.$.inputFile.selectedFile.type != "FILE"){
        msg += "Error: Please select an input file.\n";
        validated=false;
    }
     if( this.scope.$.jobName.value == "" ){
            msg += "Error: Please add a job name.\n";
            validated=false;
        }
    if(!validated){
        alert(msg)
    }
    return validated;
}
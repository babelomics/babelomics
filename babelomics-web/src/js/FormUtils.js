/**
 * Created by ralonso on 11/28/14.
 */

function FormUtils(scope) {
    this.scope = scope;
}

FormUtils.prototype.basicValidationForm = function () {
    var validated = true;
    var msg = "";
    if( this.scope.$.outdir.selectedFile === undefined){
        msg += "Missing output folder. Please select one.\n";
        validated=false;
    }
    if( this.scope.$.inputFile.selectedFile === undefined ||  this.scope.$.inputFile.selectedFile.type != "FILE"){
        msg += "Missing input file. Please select one.\n";
        validated=false;
    }
     if( this.scope.$.jobName.value == "" ){
            msg += "Missing job name. Please select one.\n";
            validated=false;
        }
    if(!validated){
        alert(msg)
    }
    return validated;
}
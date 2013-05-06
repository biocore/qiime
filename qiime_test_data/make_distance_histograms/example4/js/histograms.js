function visibilityAndOpacity(chkbox, group) {
    toggleCheckbox(chkbox, group);
    setAllOpacity();
}

function setAllOpacity() {
    var visible = document.getElementsByName("visible");
    var opacity_value = (1.0/visible.length)*100.0
    for(i=0; i<visible.length; i++) {
        SetOpacity(visible[i],opacity_value)
    }
}

function toggleCheckbox(chkbox, group) { 
    var visSetting = (chkbox.checked) ? "visible" : "hidden"; 
    document.getElementById(group).style.visibility = visSetting;
    document.getElementById(group).name = visSetting;
}

function makeVisible(group) { 
    document.getElementById(group).style.visibility = "visible";
    document.getElementById(group).name = "visible";
}

function makeHidden(group) { 
    document.getElementById(group).style.visibility = "hidden";
    document.getElementById(group).name = "hidden";
}

function mouseoverVisible(group) {
    var checkname = "check_" + group
    if (! document.getElementById(checkname).checked){
        makeVisible(group);
        setAllOpacity();
    }
}

function mouseoverHidden(group) {
    var checkname = "check_" + group
    if (! document.getElementById(checkname).checked){
        makeHidden(group);
        setAllOpacity();
    }
}

function SetOpacity(elem, opacityAsInt)
{
	var opacityAsDecimal = opacityAsInt;
	
	if (opacityAsInt > 100)
		opacityAsInt = opacityAsDecimal = 100; 
	else if (opacityAsInt < 0)
		opacityAsInt = opacityAsDecimal = 0; 
	
	opacityAsDecimal /= 100;
	if (opacityAsInt < 1)
		opacityAsInt = 1; // IE7 bug, text smoothing cuts out if 0
	
	elem.style.opacity = opacityAsDecimal;
	elem.style.filter  = "alpha(opacity=" + opacityAsInt + ")";
}

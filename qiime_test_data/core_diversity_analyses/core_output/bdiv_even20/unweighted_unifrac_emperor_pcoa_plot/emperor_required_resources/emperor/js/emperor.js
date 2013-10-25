/*
 * __author__ = "Meg Pirrung"
 * __copyright__ = "Copyright 2013, Emperor"
 * __credits__ = ["Meg Pirrung","Antonio Gonzalez Pena","Yoshiki Vazquez Baeza","Jackson Chen"]
 * __license__ = "BSD"
 * __version__ = "0.9.1"
 * __maintainer__ = "Meg Pirrung"
 * __email__ = "meganap@gmail.com"
 * __status__ = "Release"
 */

// spheres and ellipses that are being displayed on screen
var g_plotSpheres = {};
var g_plotEllipses = {};
var g_plotTaxa = {};
var g_plotVectors = {};
var g_plotEdges = {};
var g_parallelPlots = []

// sample identifiers of all items that are plotted
var g_plotIds = [];

// line objects used to represent the axes
var g_xAxisLine;
var g_yAxisLine;
var g_zAxisLine;
var g_viewingAxes = [0, 1, 2];

// scene elements for the webgl plot
var g_mainScene;
var g_sceneCamera;
var g_sceneLight;
var g_mainRenderer;
var g_sceneControl;

// general multipurpose variables
var g_elementsGroup; // group that holds the plotted shapes
var g_categoryIndex = 0; // current coloring category index
var g_genericSphere; // generic sphere used for plots
var g_categoryName = ""; // current coloring category
var g_foundId = ""; // id of currently located point
var g_time;
var g_visiblePoints = 0;
var g_sphereScaler = 1.0;
var g_keyBuilt = false;
var g_useDiscreteColors = true;
var g_screenshotBind;
var g_separator_left;
var g_separator_history;

// valid ascii codes for filename
var g_validAsciiCodes = new Array();
// adding: .-_
g_validAsciiCodes = g_validAsciiCodes.concat([45,46,95]);
// adding: 0->9
g_validAsciiCodes = g_validAsciiCodes.concat([48,49,50,51,52,53,54,55,56,57]);
// adding: A->Z
g_validAsciiCodes = g_validAsciiCodes.concat([65,66,67,68,68,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90]);
// adding: a->z
g_validAsciiCodes = g_validAsciiCodes.concat([97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122]);

// taken from the qiime/colors.py module; a total of 29 colors
var k_QIIME_COLORS = [
"0xFF0000", // red1
"0x0000FF", // blue1
"0xF27304", // orange1
"0x008000", // green
"0x91278D", // purple1
"0xFFFF00", // yellow1
"0x7CECF4", // cyan1
"0xF49AC2", // pink1
"0x5DA09E", // teal1
"0x6B440B", // brown1
"0x808080", // gray1
"0xF79679", // red2
"0x7DA9D8", // blue2
"0xFCC688", // orange2
"0x80C99B", // green2
"0xA287BF", // purple2
"0xFFF899", // yellow2
"0xC49C6B", // brown2
"0xC0C0C0", // gray2
"0xED008A", // red3
"0x00B6FF", // blue3
"0xA54700", // orange3
"0x808000", // green3
"0x008080"] // teal3

// Taken from http://stackoverflow.com/questions/18082/validate-numbers-in-javascript-isnumeric
function isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

/*This function recenter the camera to the initial position it had*/
function resetCamera() {
	// We need to reset the camera controls first before modifying the values of the camera (this is the reset view!)
	g_sceneControl.reset();
	
	g_sceneCamera.aspect = document.getElementById('pcoaPlotWrapper').offsetWidth/document.getElementById('pcoaPlotWrapper').offsetHeight;
	g_sceneCamera.rotation.set( 0, 0, 0);
	g_sceneCamera.updateProjectionMatrix();
	
	if ($('#scale_checkbox').is(':checked'))
		g_sceneCamera.position.set(0 , 0, (g_maximum*4.2) + g_radius);
	else
		g_sceneCamera.position.set(0 , 0, (g_maximum*4.8) + g_radius);	
}

/*Removes duplicates from a list of samples*/
function dedupe(list) {
	var set = {};
	for (var i = 0; i < list.length; i++){
		set[list[i]] = true;
	}
	list = [];
	for (var obj in set){
		list.push(obj);
	}
	return list;
}

/*Toggle between the scaled and unscaled version of the plot

  This function will change multiple elements in the plot, as described by the
  percentage explained in each of the PC axes.

  Note that this function will also change the position of the camera, light and
  adds a scaling value to the sphere size slider to make the size consistent
  between scaled and unscaled versions of the plot.
*/
function toggleScaleCoordinates(element) {
	var axesLen;
	var operation;

	// used only for vector and edges re-drawing
	var currentPosition = [], currentColor = 0x000000;

	if (!isNumeric(g_fractionExplained[g_viewingAxes[0]])) {
		alert("PC" + (g_viewingAxes[0]+1) + " is too small for this feature, change your selection.");
		return;
	} else if (!isNumeric(g_fractionExplained[g_viewingAxes[1]])) {
		alert("PC" + (g_viewingAxes[1]+1) + " is too small for this feature, change your selection.");
		return;
	} else if (!isNumeric(g_fractionExplained[g_viewingAxes[2]])) {
		alert("PC" + (g_viewingAxes[2]+1) + " is too small for this feature, change your selection.");
		return;
	}

	// XOR operation for the checkbox widget, this will select an operation
	// to perform over various properties, either a multiplication or a division
	if(element.checked == true){
		operation = function(a, b){return a*b};
		g_sphereScaler = g_fractionExplained[g_viewingAxes[0]];
	}
	else{
		operation = function(a, b){return a/b};
		g_sphereScaler = 1;
	}

	// force an update of the size of the spheres
	$("#sradiusslider").slider("value",$("#sradiusslider").slider("value"));

	// scale other properties
	g_xMaximumValue = operation(g_xMaximumValue,g_fractionExplained[g_viewingAxes[0]]);
	g_yMaximumValue = operation(g_yMaximumValue,g_fractionExplained[g_viewingAxes[1]]);
	g_zMaximumValue = operation(g_zMaximumValue,g_fractionExplained[g_viewingAxes[2]]);
	g_xMinimumValue = operation(g_xMinimumValue,g_fractionExplained[g_viewingAxes[0]]);
	g_yMinimumValue = operation(g_yMinimumValue,g_fractionExplained[g_viewingAxes[1]]);
	g_zMinimumValue = operation(g_zMinimumValue,g_fractionExplained[g_viewingAxes[2]]);
	g_maximum = operation(g_maximum, g_fractionExplained[g_viewingAxes[0]])

	// scale the position of the light
	g_sceneLight.position.set(
		operation(g_sceneLight.position.x, g_fractionExplained[g_viewingAxes[0]]),
		operation(g_sceneLight.position.y, g_fractionExplained[g_viewingAxes[0]]),
		operation(g_sceneLight.position.z, g_fractionExplained[g_viewingAxes[0]]));

	// scale the position of the camera according to pc1
	g_sceneCamera.position.set(
		operation(g_sceneCamera.position.x, g_fractionExplained[g_viewingAxes[0]]),
		operation(g_sceneCamera.position.y, g_fractionExplained[g_viewingAxes[0]]),
		operation(g_sceneCamera.position.z, g_fractionExplained[g_viewingAxes[0]]))
	// scale the axis lines
	drawAxisLines();

	// set the new position of each of the sphere objects
	for (sample_id in g_plotSpheres){
		// scale the position of the spheres
		g_plotSpheres[sample_id].position.set(
			operation(g_plotSpheres[sample_id].position.x,g_fractionExplained[g_viewingAxes[0]]),
			operation(g_plotSpheres[sample_id].position.y,g_fractionExplained[g_viewingAxes[1]]),
			operation(g_plotSpheres[sample_id].position.z,g_fractionExplained[g_viewingAxes[2]]));
	}

	// ellipses won't always be available hence the two separate loops
	for (sample_id in g_plotEllipses){
		// scale the dimensions of the positions of each ellipse
		g_plotEllipses[sample_id].position.set(
			operation(g_plotEllipses[sample_id].position.x, g_fractionExplained[g_viewingAxes[0]]),
			operation(g_plotEllipses[sample_id].position.y, g_fractionExplained[g_viewingAxes[1]]),
			operation(g_plotEllipses[sample_id].position.z, g_fractionExplained[g_viewingAxes[2]]));

		// scale the dimensions of the ellipse
		g_plotEllipses[sample_id].scale.set(
			operation(g_plotEllipses[sample_id].scale.x, g_fractionExplained[g_viewingAxes[0]]),
			operation(g_plotEllipses[sample_id].scale.y, g_fractionExplained[g_viewingAxes[1]]),
			operation(g_plotEllipses[sample_id].scale.z, g_fractionExplained[g_viewingAxes[2]]));
	}

	for (index in g_plotTaxa){
		//scale the dimensions of the positions of each taxa-sphere
		g_plotTaxa[index].position.set(
			operation(g_plotTaxa[index].position.x, g_fractionExplained[g_viewingAxes[0]]),
			operation(g_plotTaxa[index].position.y, g_fractionExplained[g_viewingAxes[1]]),
			operation(g_plotTaxa[index].position.z, g_fractionExplained[g_viewingAxes[2]]));

		//scale the dimensions of each taxa-sphere
		g_plotTaxa[index].scale.set(
			operation(g_plotTaxa[index].scale.x, g_fractionExplained[g_viewingAxes[0]]),
			operation(g_plotTaxa[index].scale.y, g_fractionExplained[g_viewingAxes[0]]),
			operation(g_plotTaxa[index].scale.z, g_fractionExplained[g_viewingAxes[0]]));
	}

	// each line is indexed by a sample, creating in turn TOTAL_SAMPLES-1 lines
	for (sample_id in g_plotVectors){

		// the color has to be formatted as an hex number for makeLine to work
		currentColor = g_plotVectors[sample_id].material.color.getHex();

		// updating the position of a vertex in a line is a really expensive
		// operation, hence we just remove it from the group and create it again
		g_elementsGroup.remove(g_plotVectors[sample_id]);

		for (vertex in g_plotVectors[sample_id].geometry.vertices){
			currentPosition[vertex] = g_plotVectors[sample_id].geometry.vertices[vertex];

			// scale the position of each of the vertices
			currentPosition[vertex].x = operation(currentPosition[vertex].x,
				g_fractionExplained[g_viewingAxes[0]])
			currentPosition[vertex].y = operation(currentPosition[vertex].y,
				g_fractionExplained[g_viewingAxes[1]])
			currentPosition[vertex].z = operation(currentPosition[vertex].z,
				g_fractionExplained[g_viewingAxes[2]])

			// create an array we can pass to makeLine
			currentPosition[vertex] = [currentPosition[vertex].x,
				currentPosition[vertex].y, currentPosition[vertex].z]
		}

		// add the element to the main vector array and to the group
		g_plotVectors[sample_id] = makeLine(currentPosition[0],
			currentPosition[1], currentColor, 2);
		g_elementsGroup.add(g_plotVectors[sample_id]);
	}

	// support scaling of edges in plot comparisons
	for(var sample_id in g_plotEdges){
		// each edge is composed of two separate lines
		for(var section in g_plotEdges[sample_id]){

			// the color has to be formatted as an hex number for makeLine to work
			currentColor = g_plotEdges[sample_id][section].material.color.getHex();

			// remove them completely from the group and scene we no longer need
			// these objects as re-creating them is as expensive as modifying
			// most of their features
			g_elementsGroup.remove(g_plotEdges[sample_id][section])
			g_mainScene.remove(g_plotEdges[sample_id][section])

			for (vertex in g_plotEdges[sample_id][section].geometry.vertices){
				currentPosition[vertex] = g_plotEdges[sample_id][section].geometry.vertices[vertex];

				// scale the position of each of the vertices
				currentPosition[vertex].x = operation(currentPosition[vertex].x,
					g_fractionExplained[g_viewingAxes[0]])
				currentPosition[vertex].y = operation(currentPosition[vertex].y,
					g_fractionExplained[g_viewingAxes[1]])
				currentPosition[vertex].z = operation(currentPosition[vertex].z,
					g_fractionExplained[g_viewingAxes[2]])

				// create an array we can pass to makeLine
				currentPosition[vertex] = [currentPosition[vertex].x,
					currentPosition[vertex].y, currentPosition[vertex].z]
			}

			// add the element to the main vector array and to the group
			g_plotEdges[sample_id][section] = makeLine(currentPosition[0],
				currentPosition[1], currentColor, 2);
			g_elementsGroup.add(g_plotEdges[sample_id][section]);
			g_mainScene.add(g_plotEdges[sample_id][section]);

		}
	}
}

/* Toggle between discrete and continuous coloring for samples and labels */
function toggleContinuousAndDiscreteColors(element){
	g_useDiscreteColors = !element.checked;

	// re-coloring the samples and labels now will use the appropriate coloring
	colorByMenuChanged();
	labelMenuChanged();
}

/*Generate a list of colors that corresponds to all the samples in the plot

  This function will generate a list of coloring values depending on the
  coloring scheme that the system is currently using (discrete or continuous).
*/
function getColorList(vals) {
	var colors = {};

	// cases with one or two categories are basically the same no matter if the
	// coloring scheme is continuous or discrete; choose red or red and blue
	if(vals.length == 1){
		colors[vals[0]] = new THREE.Color();
		colors[vals[0]].setHex("0xff0000");
	}
	else if (vals.length == 2) {
		colors[vals[0]] = new THREE.Color();
		colors[vals[0]].setHex("0xff0000");
		colors[vals[1]] = new THREE.Color();
		colors[vals[1]].setHex("0x0000ff");
	}
	else {
		for(var index in vals){
			colors[vals[index]] = new THREE.Color();
			if(g_useDiscreteColors){
				// get the next available color
				colors[vals[index]].setHex(getDiscreteColor(index)*1);
			}
			else{
				// multiplying the value by 0.66 makes the colormap go R->G->B
				THREE.ColorConverter.setHSV(colors[vals[index]], index*.66/vals.length, 1, 1)
			}
		}
	}
	return colors;
}

/* Retrieve one of the discrete colors from the list

  This function will return the color at the requested index, if this value
  value is greater than the number of colors available, the function will just
  rollover and retrieve the next available color.
*/
function getDiscreteColor(index){
	var size = k_QIIME_COLORS.length;
	if(index >= size){
		index = index - (Math.floor(index/size)*size)
	}

	return k_QIIME_COLORS[index]
}

/*Start timer (for debugging)*/
function startTimer() {
	var d=new Date()
	g_time = d.getTime();
}

/*End timer (for debugging)*/
function stopTimer(info) {
	var d=new Date()
	g_time = d.getTime() - g_time;
	console.log("time to " +info +":"+g_time+"ms")
}

/* Sorting function that deals with alpha and numeric elements

  This function takes a list of strings, divides it into two new lists, one
  that's alpha-only and one that's numeric only. The resulting list will have
  sorted all alpha elements at the beginning & all numeric elements at the end.
 */
function _splitAndSortNumericAndAlpha(list){
	var numericPart = [], alphaPart = [], result = [];

	// separate the numeric and the alpha elements of the array
	for(var index = 0; index < list.length; index++){
		if(isNaN(parseFloat(list[index]))){
			alphaPart.push(list[index])
		}
		else{
			numericPart.push(list[index])
		}
	}

	// sort each of the two parts, numeric part is ascending order
	alphaPart.sort();
	numericPart.sort(function(a,b){return parseFloat(a)-parseFloat(b)})

	return result.concat(alphaPart, numericPart);
}

/*This function is called when a new value is selected in the colorBy menu */
function colorByMenuChanged() {
	// set the new current category and index
	g_categoryName = document.getElementById('colorbycombo')[document.getElementById('colorbycombo').selectedIndex].value;
	g_categoryIndex = g_mappingFileHeaders.indexOf(g_categoryName);

	// get all values of this category from the mapping
	var vals = [];
	for(var i in g_plotIds){
		vals.push(g_mappingFileData[g_plotIds[i]][g_categoryIndex]);
	}

	vals = _splitAndSortNumericAndAlpha(dedupe(vals));
	colors = getColorList(vals);
	
	// build the colorby table in HTML
	var lines = "<table id='colorbylist_table'>";
	for(var i in vals){
		// each field is identified by the value it has in the deduplicated
		// list of values and by the number of the column in the mapping file
		// if this is done otherwise, weird characters have to be extemped etc.
		var idString = "r"+i+"c"+g_categoryIndex;

		// set the div id so that we can reference this div later
		lines += "<tr><td><div id=\""+idString+"\"class=\"colorbox\" name=\""+vals[i]+"\"></div></td><td title=\""+vals[i]+"\">";

		if(vals[i].length > 25){
			lines+= vals[i].substring(0,25) + "..."
		}
		else{
			lines += vals[i];
		}

		lines+= "</td></tr>";
	}
	lines += "</table>";
	document.getElementById("colorbylist").innerHTML = lines;

	for(var i in vals){
		// each field is identified by the value it has in the deduplicated
		// list of values and by the number of the column in the mapping
		var idString = "r"+i+"c"+g_categoryIndex;

		// get the div built earlier and turn it into a color picker
		$('#'+idString).css('backgroundColor',"#"+colors[vals[i]].getHexString());
		$("#"+idString).spectrum({
			localStorageKey: 'key',
			color: colors[vals[i]].getHexString(),
			showInitial: true,
			showInput: true,
			preferredFormat: "hex6",
			change:
				function(color) {
					$(this).css('backgroundColor', color.toHexString());
					var c = color.toHexString();
					if(c.length == 4){
						c = "#"+c.charAt(1)+c.charAt(1)+c.charAt(2)+c.charAt(2)+c.charAt(3)+c.charAt(3);
					}
					colorChanged($(this).attr('name'), c);
					colors[$(this).attr('name')] = c;
					colorParallelPlots(vals, colors);
				}
		});
	}

	colorParallelPlots(vals, colors);
	setKey(vals, colors);
}

function colorParallelPlots(vals,colors) 
{
	pwidth = document.getElementById('parallelPlotWrapper').offsetWidth
	pheight = document.getElementById('parallelPlotWrapper').offsetHeight
	
	document.getElementById('parallelPlotWrapper').innerHTML = '<div id="parallelPlot" class="parcoords" style="width:'+pwidth+'px;height:'+pheight+'px"></div>'
	
	var color = function(d) {
		var sid = d[0];
		var divid = sid.replace(/\./g,'')+"_key";
		var catValue = g_mappingFileData[sid][g_categoryIndex];
		var catColor = colors[catValue];

		try {
			var hex = '#'+catColor.getHexString();
		}catch(TypeError) {
			var hex = catColor;
		}
		return hex;
	}

	var pc = d3.parcoords()("#parallelPlot")
	  .data(g_parallelPlots)
	  .color(color)
	  .margin({ top: 40, left: 50, bottom: 40, right: 0 })
	  .mode("queue")
	  .render()
	  .brushable();

	$('.parcoords text').css('stroke', $('#axeslabelscolor').css('backgroundColor'));
	$('.parcoords .axis line, .parcoords .axis path').css('stroke', $('#axescolor').css('backgroundColor'));
}

/*Callback when the scaling by drop-down menu changes

  This function will create one slider and a label for each category.
*/
function scalingByMenuChanged(){
	var scalingByCategoryName = document.getElementById('scalingbycombo')[document.getElementById('scalingbycombo').selectedIndex].value;
	var scalingByCategoryIndex = g_mappingFileHeaders.indexOf(scalingByCategoryName);
	var values = [], lines, idString;

	// get all values of this category from the mapping
	for(var i in g_plotIds){
		values.push(g_mappingFileData[g_plotIds[i]][scalingByCategoryIndex]);
	}
	values = _splitAndSortNumericAndAlpha(dedupe(values));

	lines = "<table width=\"90%\">"
	for(var i in values){
		// construct a sanitized category name for all HTML elements
		idString = "r"+i+"c"+scalingByCategoryIndex;

		lines += "<tr><td>";
		// add a label with the name of the category
		lines +=" <label for=\""+values[i]+"\" class=\"text\">"
		// do not make the category name too long
		if(values[i].length > 25){
			lines+= values[i].substring(0,25) + "..."
		}
		else{
			lines += values[i];
		}
		lines +="</label>"

		// // add a slider and a current-value-label to the table
		lines += "<label id=\""+idString+"scalingvalue\" name=\""+values[i]+"\" class=\"slidervalue\"></label>"
		lines += "<div id=\""+idString+"scalingslider\" name=\""+values[i]+"\" class=\"slider-range-max\"></div>"
		lines +="</td></tr>";
	}
	lines += "</table>";
	document.getElementById("scalingbylist").innerHTML = lines;

	// set all sliders to the default value of 5, that's reflected as no scaling
	for(var i in values){
		var idString = "r"+i+"c"+scalingByCategoryIndex;
		$("#"+idString+"scalingslider").slider({
			range: "max",
			min: 1,
			max: 20,
			value: 5,
			slide: function( event, ui ) {
				sphereRadiusChange(ui, $(this).attr('name'));
			},
			change: function( event, ui ) {
				sphereRadiusChange(ui, $(this).attr('name'));
			}
		});
		document.getElementById(idString+"scalingvalue").innerHTML = $("#"+idString+"scalingslider").slider("value")/5;
	}
}


/*This function is called when a new value is selected in the showBy menu*/
function showByMenuChanged() {
	g_categoryName = document.getElementById('showbycombo')[document.getElementById('showbycombo').selectedIndex].value;
	var showByMenuIndex = g_mappingFileHeaders.indexOf(g_categoryName);
	var vals = [];

	for(var i in g_plotIds){
		var sid = g_plotIds[i];
		var divid = sid.replace(/\./g,'');
		// get all of the values for the selected category
		vals.push(g_mappingFileData[sid][showByMenuIndex]);
		// set everything to visible
		try {
			g_elementsGroup.add(g_plotEllipses[sid])
		}
		catch(TypeError){}
		try {
			g_elementsGroup.add(g_plotSpheres[sid])
		}
		catch(TypeError){}
		try {
			g_elementsGroup.add(g_plotVectors[sid])
		}
		catch(TypeError){}
		$('#'+divid+"_label").css('display','block');
	}

	g_visiblePoints = g_plotIds.length;
	changePointCount();

	vals = _splitAndSortNumericAndAlpha(dedupe(vals));

	// build the showby checkbox table in HTML
	var lines = "<form name=\"showbyform\"><table width=\"90%\">"
	for(var i in vals){
		// tag each slider & percent label with the idString to avoid conflicts
		var idString = "r"+i+"c"+showByMenuIndex;

		// make the size of the checkmark fixed so proportions don't look weird
		lines += "<tr><td width=\"10px\">";
		lines +="<input name=\""+vals[i]+"_show\" value=\""+vals[i]+"\" type=\"checkbox\" checked=\"yes\" onClick=\"toggleVisible(\'"+vals[i]+"\')\">";
		lines +="</input></td><td title=\""+vals[i]+"\">";
		if(vals[i].length > 25){
			lines+= vals[i].substring(0,25) + "..."
		}
		else{
			lines += vals[i];
		}

		// add a slider and a current-value-label to the table
		lines +="</td></tr>";
		lines += "<tr><td></td><td>";
		lines += "<label id=\""+idString+"opacityvalue\" name=\""+vals[i]+"\" class=\"slidervalue\"></label>"
		lines += "<div id=\""+idString+"opacityslider\" name=\""+vals[i]+"\" class=\"slider-range-max\"></div>"
		lines +="</td></tr>";
	}
	lines += "</table></form>";
	document.getElementById("showbylist").innerHTML = lines;

	// set all the sliders to 100 % and to respond to the sphereOpacityChange
	// with the name of the category they have in the mapping file
	for(var i in vals){
		var idString = "r"+i+"c"+showByMenuIndex;
		$("#"+idString+"opacityslider").slider({
			range: "max",
			min: 0,
			max: 100,
			value: 100,
			slide: function( event, ui ) {
				sphereOpacityChange(ui, $(this).attr('name'));
			},
			change: function( event, ui ) {
				sphereOpacityChange(ui, $(this).attr('name'));
			}
		});
		document.getElementById(idString+"opacityvalue").innerHTML = $("#"+idString+"opacityslider").slider("value")+"%";
	}

	// change the value of the general opacity for all the spheres, this action
	// has to take place after the creation of the sliders for all categories;
	// that is the for loop right befor this statement, as this will produce
	// a callback that will require to change all the sliders in available
	$("#sopacityslider").slider("value", 100)
}

/*Toggle plot items by category selected in showby menu*/
function toggleVisible(value) {

	var hidden = !document.showbyform.elements[value+'_show'].checked;
	g_categoryName = document.getElementById('showbycombo')[document.getElementById('showbycombo').selectedIndex].value;

	//change visibility of points depending on metadata category
	for(var i in g_plotIds){
	var sid = g_plotIds[i];
	var divid = sid.replace(/\./g,'');
	var mappingVal = g_mappingFileData[sid][g_mappingFileHeaders.indexOf(g_categoryName)]
		if(mappingVal == value && hidden){
			try{
				g_elementsGroup.remove(g_plotEllipses[sid]);
			}
			catch(TypeError){}
			try{
				g_elementsGroup.remove(g_plotSpheres[sid]);
				g_visiblePoints--
			}
			catch(TypeError){}
			try{
				g_elementsGroup.remove(g_plotVectors[sid]);
			}
			catch(TypeError){}
			$('#'+divid+"_label").css('display','none');
		}
		else if(mappingVal == value && !hidden)
		{
			try {
				g_elementsGroup.add(g_plotEllipses[sid]);
			}
			catch(TypeError){}
			try {
				g_elementsGroup.add(g_plotSpheres[sid]);
				g_visiblePoints++;
			}
			catch(TypeError){}
			try{
				g_elementsGroup.add(g_plotVectors[sid]);
			}
			catch(TypeError){}
			$('#'+divid+"_label").css('display','block');
		}
	}
	changePointCount()

}

/*Build the plot legend in HTML*/
function setKey(values, colors) {
	if(g_keyBuilt){
		for(var i = 0; i < values.length; i++){
			colorChanged(values[i], '#'+colors[values[i]].getHexString());
		}
	}
	else {
		var keyHTML = "<table class=\"key\">";
		for(var i in g_plotIds){
			var sid = g_plotIds[i];
			var divid = sid.replace(/\./g,'')+"_key";
			var catValue = g_mappingFileData[sid][g_categoryIndex];
			var catColor = colors[catValue];
			keyHTML += "<tr id=\""+divid+"row\"><td><div id=\""+divid+"\" name=\""+sid+"\" class=\"colorbox\" style=\"background-color:#";
			keyHTML += catColor.getHexString();
			keyHTML += ";\"></div>";
			keyHTML +="</td><td>";
			keyHTML += sid;
			keyHTML += "</td></tr>";

			try {
				g_plotEllipses[g_plotIds[i]].material.color.setHex("0x"+catColor.getHexString());
			}
			catch(TypeError){}
			try {
				g_plotSpheres[g_plotIds[i]].material.color.setHex("0x"+catColor.getHexString());
			}
			catch(TypeError){}
			try {
				g_plotVectors[g_plotIds[i]].material.color.setHex("0x"+catColor.getHexString());
			}
			catch(TypeError){}
		}
		keyHTML += "</table>";
		document.getElementById("key").innerHTML = keyHTML;

		for(var i in g_plotIds){
			var sid = g_plotIds[i];
			var divid = sid.replace(/\./g,'')+"_key";
			$('#'+divid).attr('name',sid);
			$('#'+divid).dblclick(function () {
			toggleFinder($(this), $(this).attr('name'));
			});
		}
		g_keyBuilt = true;
	}
}

/*Toggle an arrow to locate a sample by double clicking the box @ the key menu*/
function toggleFinder(div, divName) {
	if(g_foundId != divName) {
		$('.colorbox').css('border','1px solid black');
		div.css('border','1px solid white');
		$('#finder').css('opacity',1);
		var coords = toScreenXY(g_plotSpheres[divName].position, g_sceneCamera, $('#main_plot'));
		$('#finder').css('left',coords['x']-15);
		$('#finder').css('top',coords['y']-5);
		g_foundId = divName;
	}
	else {
		if($('#finder').css('opacity') == 1) {
			$('#finder').css('opacity',0);
			div.css('border','1px solid black');
			g_foundId = null
		}
		else {
			$('#finder').css('opacity',1);
			div.css('border','1px solid white');
		}
	}
}

/*Callback for the colorChanged event as triggered by the color picker*/
function colorChanged(catValue,color) {
	for(var i in g_plotIds)
	{
		var sid = g_plotIds[i]
		if(g_mappingFileData[g_plotIds[i]][g_categoryIndex] == catValue)
		{
			// get the valid divId for the key and set its color
			$("#"+sid.replace(/\./g,'')+"_key").css('backgroundColor',color);
			// set the color of the corresponding sphere and ellipse 
			try {
				g_plotEllipses[sid].material.color.setHex(color.replace('#','0x'));
			}
			catch(TypeError){}
			try {
				g_plotSpheres[sid].material.color.setHex(color.replace('#','0x'));
			}
			catch(TypeError){}
			try{
				g_plotVectors[sid].material.color.setHex(color.replace('#','0x'));
			}
			catch(TypeError){}
		}
	}
}

/* This function is called when q new color is selected for #taxaspherescolor */
function colorChangedForTaxaSpheres(color){
	for (index in g_plotTaxa){
		g_plotTaxa[index].material.color.setHex(color)
	}
}

/* This function is called when a new color is selected for the edges

 The two input parameters are a hexadecimal formatted color and an index, the
 index indicates which side of the edges are going to be re-colored.
*/
function colorChangedForEdges(color, index){
	for(var sample_id in g_plotEdges){
		currentColor = g_plotEdges[sample_id][index].material.color.setHex(color);
	}
}

/*This function is called when a new value is selected in the label menu*/
function labelMenuChanged() {
	if(document.getElementById('labelcombo').selectedIndex == 0){
		document.getElementById("labellist").innerHTML = "";
		return;
	}

	// set the new current category and index
	var labelCategory = document.getElementById('labelcombo')[document.getElementById('labelcombo').selectedIndex].value;
	var labelCatIndex = g_mappingFileHeaders.indexOf(labelCategory);

	// get all values of this category from the mapping
	var vals = [];
	for(var i in g_plotIds){
		vals.push(g_mappingFileData[g_plotIds[i]][labelCatIndex]);
	}

	vals = _splitAndSortNumericAndAlpha(dedupe(vals));
	colors = getColorList(vals);

	// build the label table in HTML
	var lines = "<form name=\"labels\" id=\"labelForm\"><table>";
	for(var i in vals){
		// each field is identified by the value it has in the deduplicated
		// list of values and by the number of the column in the mapping file
		// if this is done otherwise, weird characters have to be extemped etc.
		var idString = "r"+i+"c"+g_categoryIndex;

		// set the div id, checkbox name so that we can reference this later
		lines += "<tr><td><input name=\""+vals[i]+"\" type=\"checkbox\" checked=\"true\" onClick=\"toggleLabels()\" ></input></td><td><div id=\""+idString+"Label\" class=\"colorbox\" name=\""+vals[i]+"\"></div></td><td title=\""+vals[i]+"\">";

		if(vals[i].length > 25){
			lines+= vals[i].substring(0,25) + "..."
		}
		else{
			lines += vals[i];
		}

		lines+= "</td></tr>";
	}

	lines += "</table></form>";
	document.getElementById("labellist").innerHTML = lines;

	for(var i in vals){
		// each field is identified by the value it has in the deduplicated
		// list of values and by the number of the column in the mapping file
		// if this is done otherwise, weird characters have to be extemped etc.
		var idString = "r"+i+"c"+g_categoryIndex;

		// get the div built earlier and turn it into a color picker
		$('#'+idString+'Label').css('backgroundColor',"#"+colors[vals[i]].getHexString());
		labelColorChanged(vals[i], "#"+colors[vals[i]].getHexString());

		$("#"+idString+'Label').spectrum({
			color: colors[vals[i]].getHexString(),
			showInitial: true,
			showPalette: true,
			preferredFormat: "hex6",
			palette: [['red', 'green', 'blue']],
			change:
				function(color) {
					$(this).css('backgroundColor', color.toHexString());
					labelColorChanged($(this).attr('name'), color.toHexString());
				}
		});
	}
}

/*This function is called when a label color is changed*/
function labelColorChanged(value, color) {
	g_categoryName = document.getElementById('labelcombo')[document.getElementById('labelcombo').selectedIndex].value;
	value = value.replace('_','');

	for(var i in g_plotIds){
		var sid = g_plotIds[i];
		var divid = sid.replace(/\./g,'');
		if(g_mappingFileData[sid][g_mappingFileHeaders.indexOf(g_categoryName)] == value){
			$('#'+divid+"_label").css('color', color);
		}
	}
}

/*This function turns the labels on and off*/
function toggleLabels() {
	if(document.plotoptions.elements[0].checked){
		$('#labelForm').css('display','block');
		$('#labels').css('display','block');
		$("#lopacityslider").slider('enable');
		$("#labelColor").spectrum('enable');
		document.getElementById('labelcombo').disabled = false;

		if(document.labels == null){
			return;
		}

		// get the current category name to show the labels
		g_categoryName = document.getElementById('labelcombo')[document.getElementById('labelcombo').selectedIndex].value;

		// for each of the labels check if they are enabled or not
		for(var i = 0; i < document.labels.elements.length; i++){
			var hidden = !document.labels.elements[i].checked;
			var value = document.labels.elements[i].name;

			for(var j in g_plotIds){
				var sid = g_plotIds[j];
				var divid = sid.replace(/\./g,'');

				if(g_mappingFileData[sid][g_mappingFileHeaders.indexOf(g_categoryName)] == value && hidden){
					$('#'+divid+"_label").css('display', 'none');
				}
				else if(g_mappingFileData[sid][g_mappingFileHeaders.indexOf(g_categoryName)] == value && !hidden){
					$('#'+divid+"_label").css('display', 'block');
				}
			}
		}
	}
	else{
		$('#labels').css('display','none');
	}
}

/*This function turns the labels with the lineages on and off*/
function toggleTaxaLabels(){
	// present labels if the visibility checkbox is marked
	if(document.biplotoptions.elements[0].checked){
		$('#taxalabels').css('display','block');

		for(var key in g_taxaPositions){
			var taxa_label = g_taxaPositions[key]['lineage'];
			var divid = taxa_label.replace(/\./g,'');
			$('#'+key+"_taxalabel").css('display', 'block');
		}
	}
	else{
		$('#taxalabels').css('display','none');
	}
}

/* Turn on and off the spheres representing the biplots on screen */
function toggleBiplotVisibility(){
	// reduce the opacity to zero if the element should be off or to 0.5
	// if the element is supposed to be present; 0.5 is the default value
	if(!document.biplotsvisibility.elements[0].checked){
		for (index in g_plotTaxa){
			g_mainScene.remove(g_plotTaxa[index]);
		}
	}
	else{
		for (index in g_plotTaxa){
			g_mainScene.add(g_plotTaxa[index])
		}
	}
}

/* Turn on and off the lines connecting the samples being compared */
function toggleEdgesVisibility(){

	// each edge is really composed of two lines and those elements are stored
	// in each of the keys that are stored for each sample comparison
	if(!document.edgesvisibility.elements[0].checked){
		for (index in g_plotEdges){
			g_mainScene.remove(g_plotEdges[index][0]);
			g_mainScene.remove(g_plotEdges[index][1]);
		}
	}
	else{
		for (index in g_plotEdges){
			g_mainScene.add(g_plotEdges[index][0]);
			g_mainScene.add(g_plotEdges[index][1]);
		}
	}
}

/*This function finds the screen coordinates of any position in the current plot.

  The main purpose of this function is to be used for calculating the placement
  of the labels.
*/
function toScreenXY( position, camera, jqdiv ) {

	var screenPosition = position.clone();
	var screenProjectionMatrix = new THREE.Matrix4();

	// multiply the matrices and aply the vector to the projection matrix
	screenProjectionMatrix.multiplyMatrices( camera.projectionMatrix,
		camera.matrixWorldInverse);
	screenPosition.applyProjection(screenProjectionMatrix);

	return { x: (screenPosition.x + 1)*jqdiv.width()/2 + jqdiv.offset().left,
		y: (-screenPosition.y+1)*jqdiv.height()/2 + jqdiv.offset().top};
}

/*This function is used to filter the key to a user's provided search string*/
function filterKey() {
	var searchVal = document.keyFilter.filterBox.value.toLowerCase();

	for(var i in g_plotIds){
		var sid = g_plotIds[i];
		var divid = sid.replace(/\./g,'')+"_keyrow";

		if(sid.toLowerCase().indexOf(searchVal) != -1){
			$('#'+divid).css('display','block');
		}
		else{
			$('#'+divid).css('display','none');
		}
	}
}

/*This function handles events from the ellipse opacity slider*/
function ellipseOpacityChange(ui) {
	document.getElementById('ellipseopacity').innerHTML = ui.value + "%";
	ellipseOpacity = ui.value/100;

	for(var sid in g_plotEllipses){
		g_plotEllipses[sid].material.opacity = ellipseOpacity;
	}
}

/*This function handles events from the sphere opacity sliders

  Note that there are two type of sliders, the master opacity slider in the
  options tab and the 'per-category' slider in the visibility tab. The later
  one controls the opacity only for the spheres belonging to a specific category
  in the mapping file.

  As for the parameters 'ui' is the jQuery slider element and category the
  string with the value of that category of the mapping file or null when the
  callback is originated from the master opacity slider.
*/
function sphereOpacityChange(ui, category) {
	var sphereOpacity = ui.value/100;
	var showByCategoryName = document.getElementById('showbycombo')[document.getElementById('showbycombo').selectedIndex].value;
	var showByCategoryIndex = g_mappingFileHeaders.indexOf(showByCategoryName);
	var vals = [], idString, newValue;

	// get all values of this category from the mapping
	for(var i in g_plotIds){
		vals.push(g_mappingFileData[g_plotIds[i]][showByCategoryIndex]);
	}
	vals = _splitAndSortNumericAndAlpha(dedupe(vals));

	// category as null means that it's the general opacity slider (the on in the options tab)
	if (category == null) {
		for (index in vals){
			idString = "r"+index+"c"+showByCategoryIndex;
			$("#"+idString+"opacityslider").slider("value", ui.value);
			document.getElementById(idString+"opacityvalue").innerHTML = $("#"+idString+"opacityslider").slider("value")+"%";
		}
		document.getElementById('sphereopacity').innerHTML = ui.value + "%";
	}
	else{
		// each field is identified by the value it has in the deduplicated
		// list of values and by the number of the column in the mapping file
		// if this is done otherwise, weird characters have to be extemped etc.
		idString = "r"+vals.indexOf(category)+"c"+showByCategoryIndex;

		for(var i in g_plotIds){
			if(g_mappingFileData[g_plotIds[i]][showByCategoryIndex] == category){
				g_plotSpheres[g_plotIds[i]].material.opacity = sphereOpacity;
			}
			document.getElementById(idString+"opacityvalue").innerHTML = $("#"+idString+"opacityslider").slider("value")+"%";
		}
	}
}

/*This function handles events from the vectors opacity slider*/
function vectorsOpacityChange(ui) {
	document.getElementById('vectorsopacity').innerHTML = ui.value + "%";
	var vectorsOpacity = ui.value/100;

	for(var sample_id in g_plotVectors){
		g_plotVectors[sample_id].material.opacity = vectorsOpacity;
	}
}

/*This function handles events from the label opacity slider*/
function labelOpacityChange(ui) {
	document.getElementById('labelopacity').innerHTML = ui.value + "%";
	labelOpacity = ui.value/100;

	$('#labels').css('opacity', labelOpacity);
}

/*This function handles events from the sphere radius slider

  Note that this function will get a scaling value added depending on whether or
  not the plot being displayed is scaled by the percent explained in each axis.

  When the category argument is null, the callback comes from the master radius
  slider in any other case it refers to a specific category.
*/
function sphereRadiusChange(ui, category) {
	var scale = (ui.value/5.0)*g_sphereScaler;
	var scalingByCategoryName = document.getElementById('scalingbycombo')[document.getElementById('scalingbycombo').selectedIndex].value;
	var scalingByCategoryIndex = g_mappingFileHeaders.indexOf(scalingByCategoryName);
	var values = [], idString;

	// get all values of this category from the mapping
	for(var i in g_plotIds){
		values.push(g_mappingFileData[g_plotIds[i]][scalingByCategoryIndex]);
	}
	values = _splitAndSortNumericAndAlpha(dedupe(values));

	if (category == null){
		for (index in values){
			idString = "r"+index+"c"+scalingByCategoryIndex;
			$("#"+idString+"scalingslider").slider("value", ui.value);
			document.getElementById(idString+"scalingvalue").innerHTML = $("#"+idString+"scalingslider").slider("value")/5;
		}

		// set the value for the master scaling slider
		document.getElementById('sphereradius').innerHTML = ui.value/5;
	}
	else{
		// each field is identified by the value it has in the deduplicated
		// list of values and by the number of the column in the mapping file
		// if this is done otherwise, weird characters have to be extemped etc.
		idString = "r"+values.indexOf(category)+"c"+scalingByCategoryIndex;

		for(var i in g_plotIds){
			if(g_mappingFileData[g_plotIds[i]][scalingByCategoryIndex] == category){
				g_plotSpheres[g_plotIds[i]].scale.set(scale, scale, scale);
			}
			document.getElementById(idString+"scalingvalue").innerHTML = $("#"+idString+"scalingslider").slider("value")/5;
		}
	}
}

/*Setup the interface elements required for the sidebar of the main interface*/
function setJqueryUi() {
	$("#menutabs").tabs();
	$("#plottype").buttonset();
	$("input[name='plottype']").change(togglePlots);
	
	$("#labelColor").css('backgroundColor', '#ffffff');

	$("#labelColor").spectrum({
		color: '#fffffff',
		showInitial: true,
		showPalette: true,
		preferredFormat: "hex6",
		palette: [['black', 'red', 'green', 'blue']],
		change:
			function(color) {
				$(this).css('backgroundColor', color.toHexString());
				$('#labels').css('color', color.toHexString());
				for(var i in g_plotIds){
					var sid = g_plotIds[i];
					var divid = sid.replace(/\./g,'');
					$('#'+divid+"_label").css('color', color.toHexString());
				}
				document.getElementById('labelcombo').selectedIndex = 0;
				labelMenuChanged();
			}
	});

	// check whether or not there is an ellipse opacity slider in the plot
	if (document.getElementById('ellipseopacity')){
		$("#eopacityslider").slider({
			range: "max",
			min: 0,
			max: 100,
			value: 20,
			slide: function( event, ui ) {
				ellipseOpacityChange(ui);
			},
			change: function( event, ui ) {
				ellipseOpacityChange(ui);
			}
		});
		document.getElementById('ellipseopacity').innerHTML = $( "#eopacityslider" ).slider( "value")+"%";
	}

	// check whether or not there is a vectors opacity slider in the plot
	if (document.getElementById('vectorsopacity')){
		$("#vopacityslider").slider({
			range: "max",
			min: 0,
			max: 100,
			value: 100,
			slide: function( event, ui ) {
				vectorsOpacityChange(ui);
			},
			change: function( event, ui ) {
				vectorsOpacityChange(ui);
			}
		});
		document.getElementById('vectorsopacity').innerHTML = $( "#vopacityslider" ).slider( "value")+"%";
	}

	// check if we are presenting biplots, to decide whether or not we should
	// show the color picker for the biplot spheres, white is the default color
	if(document.getElementById('taxaspherescolor')){
		$('#taxaspherescolor').css('backgroundColor',"#FFFFFF");
		$("#taxaspherescolor").spectrum({
			localStorageKey: 'key',
			color: "#FFFFFF",
			preferredFormat: "hex6",
			showInitial: true,
			showInput: true,
			change:
				function(color) {
					// pass a boolean flag to convert to hex6 string
					var c = color.toHexString(true);
					$(this).css('backgroundColor', c);
					colorChangedForTaxaSpheres(c.replace('#', '0x'));
				}
		});
	}
	// set up the color selector for the taxa labels
	if(document.getElementById('taxalabelcolor')){
		$('#taxalabelcolor').css('backgroundColor',"#FFFFFF");
		$("#taxalabelcolor").spectrum({
			localStorageKey: 'key',
			color: "#FFFFFF",
			preferredFormat: "hex6",
			showInitial: true,
			showInput: true,
			change:
				function(color) {
				// pass a boolean flag to convert to hex6 string
				var c = color.toHexString(true);

				// set the color for the box and for the renderer
				$(this).css('backgroundColor', c);

				// get the taxonomic assignments and append '_taxalabel' to
				// retrieve all the labels belonging to a sphere in the plot
				for(var key in g_taxaPositions) {
					$('#'+key+"_taxalabel").css('color', c);
				}
			}// on-change callback
		});
	}

	// check if the plot is a comparison plot if so, setup the elements that
	// will allow the user to change the color of the two sides of the edges
	if(document.getElementById('edgecolorselector_a')){
		$('#edgecolorselector_a').css('backgroundColor',"#FFFFFF");
		$("#edgecolorselector_a").spectrum({
			localStorageKey: 'key',
			color: "#FFFFFF",
			preferredFormat: "hex6",
			showInitial: true,
			showInput: true,
			change:
				function(color) {
					// pass a boolean flag to convert to hex6 string
					var c = color.toHexString(true);
					$(this).css('backgroundColor', c);
					colorChangedForEdges(c.replace('#', '0x'), 0);
				}
		});
	}
	if(document.getElementById('edgecolorselector_b')){
		$('#edgecolorselector_b').css('backgroundColor',"#FF0000");
		$("#edgecolorselector_b").spectrum({
			localStorageKey: 'key',
			color: "#FF0000",
			preferredFormat: "hex6",
			showInitial: true,
			showInput: true,
			change:
				function(color) {
					// pass a boolean flag to convert to hex6 string
					var c = color.toHexString(true);
					$(this).css('backgroundColor', c);
					colorChangedForEdges(c.replace('#', '0x'), 1);
				}
		});
	}

	$("#sopacityslider").slider({
		range: "max",
		min: 0,
		max: 100,
		value: 100,
		slide: function( event, ui ) {
			sphereOpacityChange(ui, null);
		},
		change: function( event, ui ) {
			sphereOpacityChange(ui, null);
		}
	});
	document.getElementById('sphereopacity').innerHTML = $( "#sopacityslider" ).slider( "value")+"%";

	$("#sradiusslider" ).slider({
		range: "max",
		min: 1,
		max: 20,
		value: 5,
		slide: function( event, ui ) {
			sphereRadiusChange(ui, null);
		},
		change: function( event, ui ) {
			sphereRadiusChange(ui, null);
		}
	});
	document.getElementById('sphereradius').innerHTML = $( "#sradiusslider" ).slider( "value")/5;

	$("#lopacityslider").slider({
		range: "max",
		min: 0,
		max: 100,
		value: 100,
		slide: function( event, ui ) {
			labelOpacityChange(ui);
		},
		change: function( event, ui ) {
			labelOpacityChange(ui);
		}
	});
	document.getElementById('labelopacity').innerHTML = $( "#lopacityslider" ).slider( "value")+"%"

	//default color for axes labels is white
	$('#axeslabelscolor').css('backgroundColor',"#FFFFFF");
	$("#axeslabelscolor").spectrum({
		localStorageKey: 'key',
		color: "#FFFFFF",
		showInitial: true,
		showInput: true,
		showPalette: true,
		preferredFormat: "hex6",
		palette: [['white', 'black']],
		change:
			function(color) {
				// pass a boolean flag to convert to hex6 string
				var c = color.toHexString(true);
				// set the color for the box and for the renderer
				$(this).css('backgroundColor', c);

				//set css for the text in the parallel plot
				$('.parcoords text').css('stroke', c);

				// change the css color of the 3d plot labels
				$("#pc1_label").css('color', c);
				$("#pc2_label").css('color', c);
				$("#pc3_label").css('color', c);

			}
	});

	//default color for the axes is white
	$('#axescolor').css('backgroundColor',"#ffffff");
	$("#axescolor").spectrum({
		localStorageKey: 'key',
		color: "#ffffff",
		showInitial: true,
		showInput: true,
		showPalette: true,
		preferredFormat: "hex6",
		palette: [['white', 'black']],
		change:
			function(color) {
				// pass a boolean flag to convert to hex6 string
				var c = color.toHexString(true);
				// create a new three.color from the string
				var axesColor = new THREE.Color();
				axesColor.setHex(c.replace('#','0x'));
				g_xAxisLine.material.color = axesColor;
				g_yAxisLine.material.color = axesColor;
				g_zAxisLine.material.color = axesColor;


				// set the color for the box and for the renderer
				$(this).css('backgroundColor', c);

				//set css for the lines of the parallel cords
				$('.parcoords .axis line, .parcoords .axis path').css('stroke', c);
			}
	});

	// the default color palette for the background is black and white
	$('#rendererbackgroundcolor').css('backgroundColor',"#000000");
	$('#parallelPlotWrapper').css('backgroundColor',"#000000");
	$("#rendererbackgroundcolor").spectrum({
		localStorageKey: 'key',
		color: "#000000",
		showInitial: true,
		showInput: true,
		showPalette: true,
		preferredFormat: "hex6",
		palette: [['white', 'black']],
		change:
			function(color) {
				// pass a boolean flag to convert to hex6 string
				var c = color.toHexString(true);

				// create a new three.color from the string
				var rendererBackgroundColor = new THREE.Color();
				rendererBackgroundColor.setHex(c.replace('#','0x'));

				// set the color for the box and for the renderer
				$(this).css('backgroundColor', c);
				g_mainRenderer.setClearColor(rendererBackgroundColor, 1);
				
				$('#parallelPlotWrapper').css('backgroundColor', c);
			}
	});
}

/*Draw the ellipses in the plot as described by the g_ellipsesDimensions array

  Note that this is a function that won't always get executed since this should
  only happen when plotting a jaccknifed principal coordinates analysis
*/
function drawEllipses() {
	for(var sid in g_ellipsesDimensions) {
		//draw ellipsoid
		var emesh = new THREE.Mesh( g_genericSphere,new THREE.MeshLambertMaterial() );
		emesh.scale.x = g_ellipsesDimensions[sid]['width']/g_radius;
		emesh.scale.y = g_ellipsesDimensions[sid]['height']/g_radius;
		emesh.scale.z = g_ellipsesDimensions[sid]['length']/g_radius;
		emesh.position.set(g_ellipsesDimensions[sid]['x'],g_ellipsesDimensions[sid]['y'],g_ellipsesDimensions[sid]['z']);
		emesh.material.color = new THREE.Color()
		emesh.material.transparent = true;
		emesh.material.depthWrite = false;
		emesh.material.opacity = 0.2;
		emesh.updateMatrix();
		emesh.matrixAutoUpdate = true;
		if(g_mappingFileData[sid] != undefined){
			g_elementsGroup.add( emesh );
			g_plotEllipses[sid] = emesh;
		}
	}
}

/*Draw the spheres in the plot as described by the g_spherePositions array*/
function drawSpheres() {
	for(var sid in g_spherePositions){
		//draw ball
		var mesh = new THREE.Mesh( g_genericSphere, new THREE.MeshLambertMaterial() );
		mesh.material.color = new THREE.Color()
		mesh.material.transparent = true;
		mesh.material.depthWrite = false;
		mesh.material.opacity = 1;
		mesh.position.set(g_spherePositions[sid]['x'], g_spherePositions[sid]['y'], g_spherePositions[sid]['z']);
		mesh.updateMatrix();
		mesh.matrixAutoUpdate = true;
		if(g_mappingFileData[sid] != undefined){
			g_elementsGroup.add( mesh );
			g_plotSpheres[sid] = mesh;
			g_plotIds.push(sid);
		}
	}
}

/*Draw the taxa spheres in the plot as described by the g_taxaPositions array

  Note that this is a function that won't always have an effect because the
  g_taxaPositions array must have elements stored in it.
*/
function drawTaxa(){
	// All taxa spheres are white by default
	var whiteColor = new THREE.Color();
	whiteColor.setHex("0xFFFFFF");

	for (var key in g_taxaPositions){
		var mesh = new THREE.Mesh(g_genericSphere,
			new THREE.MeshLambertMaterial());

		// set the volume of the sphere
		mesh.scale.x = g_taxaPositions[key]['radius'];
		mesh.scale.y = g_taxaPositions[key]['radius'];
		mesh.scale.z = g_taxaPositions[key]['radius'];

		// set the position
		mesh.position.set(g_taxaPositions[key]['x'],
			g_taxaPositions[key]['y'],
			g_taxaPositions[key]['z']);

		// the legacy color of these spheres is white
		mesh.material.color = whiteColor;
		mesh.material.transparent = true;
		mesh.material.opacity = 0.5;
		mesh.updateMatrix();
		mesh.matrixAutoUpdate = true;

		// add the element to the scene and to the g_plotTaxa dictionary
		g_elementsGroup.add(mesh)
		g_mainScene.add(mesh);
		g_plotTaxa[key] = mesh;
	}
}

/*Draw the lines for the plot as described in the g_vectorPositions array

  Note that this will draw a single line between each of the samples, hence
  resulting in N-1 lines being drawn where N is the total number of samples,
  similarly to spheres or ellipsoids these elements are added to the
  g_elementsGroup variable.
*/
function drawVectors(){
	var current_vector, previous = null;

	// There are as many vectors as categories were specified
	for (var categoryKey in g_vectorPositions){
		for (var sampleKey in g_vectorPositions[categoryKey]){

			// retrieve the initial position on the first loop
			if (previous == null){
				previous = g_vectorPositions[categoryKey][sampleKey];
			}
			// if we already have the initial position, draw the line with
			// the current position and the value of previous (initial position)
			else{
				current = g_vectorPositions[categoryKey][sampleKey];
				current_vector = makeLine(previous, current, 0xFFFFFF, 2)
				previous = current;
				g_plotVectors[sampleKey] = current_vector;

				if(g_mappingFileData[sampleKey] != undefined){
					g_elementsGroup.add(current_vector);
					g_plotVectors[sampleKey] = current_vector;
				}
			}
		}

		// reset previous so the algorithms work on the next category
		previous = null;
	}
}

/*Draw the lines that connect samples being compared (see g_comparePositiosn)

  This will draw two lines between each compared sample, one with color red and
  the other one with color white, what must be noted here is that, these two
  lines visually compose a single line and are both stored in the g_plotEdges
  array in arrays of two elements where the first element is the red line and
  the second element is the white line.
  
  A dynamic value that contains the coordinates of the spheres is passed in, that way 
  to allow the negating of the values.

  In the case of a non-serial comparison plot, all edges will originate in the
  same point.
*/
function drawEdges(spherepositions){
	var previous = null, origin = null, current, middle_point, index=0, line_a, line_b;
					
	// note that this function is composed of an if-else statement with a loop
	// that's almost identical under each case. This approach was taken as
	// otherwise the comparison would need to happen N times instead of 1 time
	// (N is the number of edges*2).

	// if the comparison is serial draw one edge after the other
	if (g_isSerialComparisonPlot == true){
		for (var sampleKey in spherepositions){
			for (var edgePosition in spherepositions[sampleKey]){
			
				// if we don't have a start point store it and move along
				if (previous == null) {
					previous = spherepositions[sampleKey][edgePosition];
				}
				// if we already have a start point then draw the edge
				else{
					current = spherepositions[sampleKey][edgePosition];
					
					// the edge is composed by two lines so calculate the middle
					// point between these two lines and end the first line in this
					// point and start the second line in this point
					middle_point = [(previous[0]+current[0])/2,
						(previous[1]+current[1])/2, (previous[2]+current[2])/2];

					currentColorA_hex = $("#edgecolorselector_a").spectrum("get").toHexString(true);
					currentColorB_hex = $("#edgecolorselector_b").spectrum("get").toHexString(true);
					line_a = makeLine(previous, middle_point, currentColorA_hex, 2);
					line_b = makeLine(middle_point, current, currentColorB_hex, 2);
					line_a.transparent = false;
					line_b.transparent = false;

					// index the two lines by the name of the sample plus a suffix
					g_plotEdges[sampleKey+'_'+index.toString()] = [line_a, line_b];

					g_elementsGroup.add(line_a);
					g_elementsGroup.add(line_b);
					g_mainScene.add(line_a);
					g_mainScene.add(line_b);

					// the current line becomes the previous line for the next
					// iteration as all samples must be connected
					previous = spherepositions[sampleKey][edgePosition];
				}
				index = index+1;
			}

			// if we've finished with the connecting lines let a new line start
			previous = null;
		}
	}
	// if the comparison is not serial, originate all edges in the same coords
	else{
		for (var sampleKey in spherepositions){
			for (var edgePosition in spherepositions[sampleKey]){
				if (origin == null) {
					origin = spherepositions[sampleKey][edgePosition];
				}
				else{
					current = spherepositions[sampleKey][edgePosition];

					// edges are composed of two lines so use the start and
					// the end point to calculate the position of the vertices
					middle_point = [(origin[0]+current[0])/2,
						(origin[1]+current[1])/2, (origin[2]+current[2])/2];

					// in the case of centered comparisons the origins are
					// painted in color white one one side and red on the other
					currentColorA_hex = $("#edgecolorselector_a").spectrum("get").toHexString(true);
					currentColorB_hex = $("#edgecolorselector_b").spectrum("get").toHexString(true);
					line_a = makeLine(origin, middle_point, currentColorA_hex, 2);
					line_b = makeLine(middle_point, current, currentColorB_hex, 2);
					line_a.transparent = false;
					line_b.transparent = false;

					// given that these are just sample repetitions just
					// just add a suffix at the end of the sample id
					g_plotEdges[sampleKey+'_'+index.toString()] = [line_a, line_b];

					g_elementsGroup.add(line_a);
					g_elementsGroup.add(line_b);
					g_mainScene.add(line_a);
					g_mainScene.add(line_b);

				}
				index = index + 1;
			}
			origin = null;
		}
	}
}

/*Save the current view to SVG
  This will take the current webGL renderer, convert it to SVG and then generate 
  a file to download. Additionally it will create the labels if this option is selected.
*/
function saveSVG(button){
    // add a name subfix for the filenames
    if ((g_segments<=8 && g_visiblePoints>=10000) || (g_segments>8 && g_visiblePoints>=5000)) {
        var res = confirm("The number of segments (" + g_segments + ") combined with the number " +
            "of samples could take a long time and in some computers the browser will crash. " +
            "If this happens we suggest to lower the number of segments or use the png " +
            "implementation. Do you want to continue?");
        if (res==false) return;
    }
 
    $('body').css('cursor','progress');
    
    var width = document.getElementById('pcoaPlotWrapper').offsetWidth;
    var height = document.getElementById('pcoaPlotWrapper').offsetHeight;
    
    var color = $("#rendererbackgroundcolor").spectrum("get").toHexString(true);
    var rendererBackgroundColor = new THREE.Color();
    rendererBackgroundColor.setHex(color.replace('#','0x'));

	var svgRenderer = new THREE.SVGRenderer({ antialias: true, preserveDrawingBuffer: true }); 
	// this is the proper way to set the color of the background but it doesn't work   
    svgRenderer.setClearColor(rendererBackgroundColor, 1);
    svgRenderer.setSize(width, height);
    svgRenderer.render(g_mainScene, g_sceneCamera);
    svgRenderer.sortObjects = true;
        
    // converting svgRenderer to string: http://stackoverflow.com/questions/17398134/three-svgrenderer-save-text-of-image
    var XMLS = new XMLSerializer(); 
    var svgfile = XMLS.serializeToString(svgRenderer.domElement);
    
    // hacking the color to the svg
    var index = svgfile.indexOf('viewBox="')+9;
    var viewBox = svgfile.substring(index, svgfile.indexOf('"',index))
    viewBox = viewBox.split(" ");
    var background = '<rect id="background" height="' + viewBox[3] + '" width="' + viewBox[2] + '" y="' + 
        viewBox[1] + '" x="' + viewBox[0] + '" stroke-width="0" stroke="#000000" fill="' + color + '"/>'
    index = svgfile.indexOf('>',index)+1;
    svgfile = svgfile.substr(0, index) + background + svgfile.substr(index);
    
    // adding xmlns header to open in the browser 
    svgfile = svgfile.replace('viewBox=', 'xmlns="http://www.w3.org/2000/svg" viewBox=')
    saveAs(new Blob([svgfile], {type: "text/plain;charset=utf-8"}), 
         $('#saveas_name').val() + ".svg");
    
    if ($('#saveas_legends').is(':checked')) {
        var labels_text = '', pos_y = 1, increment = 40, max_len = 0, font_size = 12;
        $('#colorbylist_table tr div').each(function() {
            if ($(this).attr('name').length > max_len) max_len = $(this).attr('name').length 
            
            // adding rectangle
            labels_text += '<rect height="27" width="27" y="' + pos_y + 
                '" x="1" stroke-width="1" ' + 'stroke="#FFFFFF" fill="' + 
                $("#" + $(this).attr('id')).spectrum("get").toHexString(true) + '"/>';
            // adding text
            labels_text += '<text xml:space="preserve" y="' + (pos_y+20) + '" x="35" ' + 
                'font-family="Monospace" font-size="' + font_size + '" stroke-width="0" ' +
                'stroke="#000000" fill="#000000">' + $(this).attr('name') + '</text>';
            pos_y += increment;
        });
        labels_text = '<svg width="' + ((font_size*max_len) + 10) + '" height="' + 
            (pos_y-10) + '" xmlns="http://www.w3.org/2000/svg"><g>' + labels_text + 
            '</g></svg>';
        
        saveAs(new Blob([labels_text], {type: "text/plain;charset=utf-8"}), 
            $('#saveas_name').val() + "_labels.svg");
    }
    
    $('body').css('cursor','default');
}

/*Utility function to draw two-vertices lines at a time

  This function allows you to create a line with only two vertices i. e. the
  start point and the end point, plus the color and width of the line. The
  start and end point must be 3 elements array. The color must be a hex-string
  or a hex number.
*/
function makeLine(coords_a, coords_b, color, width){
	// based on the example described in:
	// https://github.com/mrdoob/three.js/wiki/Drawing-lines
	var material, geometry, line;

	// make the material transparent and with full opacity
	material = new THREE.LineBasicMaterial({color:color, linewidth:width});
	material.matrixAutoUpdate = true;
	material.transparent = true;
	material.opacity = 1.0;

	// add the two vertices to the geometry
	geometry = new THREE.Geometry();
	geometry.vertices.push(new THREE.Vector3(coords_a[0], coords_a[1], coords_a[2]));
	geometry.vertices.push(new THREE.Vector3(coords_b[0], coords_b[1], coords_b[2]));

	// the line will contain the two vertices and the described material
	line = new THREE.Line(geometry, material);

	return line;
}

/*Draw each of the lines that represent the X, Y and Z axes in the plot

  The length of each of these axes depend on the ranges that the data being
  displayed uses.
*/
function drawAxisLines() {
	var axesColorFromColorPicker;

	// removing axes, if they do not exist the scene doesn't complain
	g_mainScene.remove(g_xAxisLine);
	g_mainScene.remove(g_yAxisLine);
	g_mainScene.remove(g_zAxisLine);

	// value should be retrieved from the picker every time the axes are drawn
	axesColorFromColorPicker = $("#axescolor").spectrum("get").toHexString(true);
	axesColorFromColorPicker = axesColorFromColorPicker.replace('#','0x')
	axesColorFromColorPicker = parseInt(axesColorFromColorPicker, 16)

	// one line for each of the axes
	g_xAxisLine = makeLine([g_xMinimumValue, g_yMinimumValue, g_zMinimumValue],
		[g_xMaximumValue, g_yMinimumValue, g_zMinimumValue], axesColorFromColorPicker, 3);
	g_yAxisLine = makeLine([g_xMinimumValue, g_yMinimumValue, g_zMinimumValue],
		[g_xMinimumValue, g_yMaximumValue, g_zMinimumValue], axesColorFromColorPicker, 3);
	g_zAxisLine = makeLine([g_xMinimumValue, g_yMinimumValue, g_zMinimumValue],
		[g_xMinimumValue, g_yMinimumValue, g_zMaximumValue], axesColorFromColorPicker, 3);

	// axes shouldn't be transparent
	g_xAxisLine.material.transparent = false;
	g_yAxisLine.material.transparent = false;
	g_zAxisLine.material.transparent = false;

	g_mainScene.add(g_xAxisLine)
	g_mainScene.add(g_yAxisLine)
	g_mainScene.add(g_zAxisLine)
};

/* update point count label */
function changePointCount() {
	document.getElementById('pointCount').innerHTML = g_visiblePoints+'/'+g_plotIds.length+' points'
}

/* Validating and modifying the view axes */
function changeAxesDisplayed() {
	if (!jQuery.isEmptyObject(g_vectorPositions) || !jQuery.isEmptyObject(g_taxaPositions) ||
			!jQuery.isEmptyObject(g_ellipsesDimensions) || g_number_of_custom_axes!=0) {
			resetCamera();
			return;
	}
	
	// HACK: this is a work around for cases when the scale is on
	if ($('#scale_checkbox').is(':checked')) toggleScaleCoordinates({'checked': false});

	var pc1_axis = $("#pc1_axis").val(), pc2_axis = $("#pc2_axis").val(),
		pc3_axis = $("#pc3_axis").val();
	var pc1_value = parseInt(pc1_axis.substring(1))-1,
		pc2_value = parseInt(pc2_axis.substring(1))-1,
		pc3_value = parseInt(pc3_axis.substring(1))-1;

	//g_fractionExplained
	if (pc1_axis==pc2_axis || pc1_axis==pc3_axis || pc2_axis==pc3_axis) {
		$("#refresh_axes_label").html('<font color="red">Not valid values, try again.</font>');
		return;
	}
	if (pc1_value>pc2_value || pc1_value>pc3_value || pc2_value>pc3_value) {
		$("#refresh_axes_label").html('<font color="red">PC3 should be > than P2, and P2 than PC1.</font>');
		return;
	}

	for (var sid in g_spherePositions) {
		g_spherePositions[sid]['x'] = g_spherePositions[sid][pc1_axis];
		g_spherePositions[sid]['y'] = g_spherePositions[sid][pc2_axis];
		g_spherePositions[sid]['z'] = g_spherePositions[sid][pc3_axis];
	}
	
	comparisonPositionlength = Object.keys(g_comparisonPositions).length
	spherePositionslength = Object.keys(g_spherePositions).length/comparisonPositionlength
	for (var sampleKey in g_comparisonPositions) {
		for (var j=0;j<spherePositionslength;j++) {
			var sid = sampleKey + "_" + j
				g_comparisonPositions[sampleKey][j][0] = g_spherePositions[sid]['x']
				g_comparisonPositions[sampleKey][j][1] = g_spherePositions[sid]['y']
				g_comparisonPositions[sampleKey][j][2] = g_spherePositions[sid]['z']
 		}
	}
	
	checkedboxes = []
    if ($('#flip_axes_1').is(':checked')) {
		for(var sid in g_spherePositions){
			g_spherePositions[sid]['x'] = g_spherePositions[sid][pc1_axis]*(-1);
		}
 		checkedboxes.push(0);
	}
    if ($('#flip_axes_2').is(':checked')) {
		for(var sid in g_spherePositions){
			g_spherePositions[sid]['y'] = g_spherePositions[sid][pc2_axis]*(-1);
		}		
 		checkedboxes.push(1);
	}
    if ($('#flip_axes_3').is(':checked')) {
		for(var sid in g_spherePositions){
			g_spherePositions[sid]['z'] = g_spherePositions[sid][pc3_axis]*(-1);
		}		
 		checkedboxes.push(2);
	}
	flipEdges(checkedboxes);
	
	// Setting up new positions
	var max_x = Number.NEGATIVE_INFINITY, max_y = Number.NEGATIVE_INFINITY,
		max_z = Number.NEGATIVE_INFINITY, min_x = Number.POSITIVE_INFINITY,
		min_y = Number.POSITIVE_INFINITY, min_z = Number.POSITIVE_INFINITY;
	for (var sid in g_spherePositions) {
		if (g_spherePositions[sid]['x']>max_x)
			max_x=g_spherePositions[sid]['x'];
		if (g_spherePositions[sid]['y']>max_y)
			max_y=g_spherePositions[sid]['y'];
		if (g_spherePositions[sid]['z']>max_z)
			max_z=g_spherePositions[sid]['z'];
		if (g_spherePositions[sid]['x']<min_x)
			min_x=g_spherePositions[sid]['x'];
		if (g_spherePositions[sid]['y']<min_y)
			min_y=g_spherePositions[sid]['y'];
		if (g_spherePositions[sid]['z']<min_z)
			min_z=g_spherePositions[sid]['z'];
	}
		
	for (var sample_id in g_plotSpheres){
		g_plotSpheres[sample_id].position.set(g_spherePositions[sample_id]['x'],
			g_spherePositions[sample_id]['y'], g_spherePositions[sample_id]['z']);
	}

	// Setting up new axes for axes by coords explained
	g_viewingAxes = [pc1_value, pc2_value, pc3_value]
	g_pc1Label = "PC" + (g_viewingAxes[0]+1) + " (" + g_fractionExplainedRounded[g_viewingAxes[0]] + " %)";
	g_pc2Label = "PC" + (g_viewingAxes[1]+1) + " (" + g_fractionExplainedRounded[g_viewingAxes[1]] + " %)";
	g_pc3Label = "PC" + (g_viewingAxes[2]+1) + " (" + g_fractionExplainedRounded[g_viewingAxes[2]] + " %)";
			
	g_xMaximumValue = max_x + (max_x>=0 ? 6*g_radius : -6*g_radius);
	g_yMaximumValue = max_y + (max_y>=0 ? 6*g_radius : -6*g_radius);
	g_zMaximumValue = max_z + (max_z>=0 ? 6*g_radius : -6*g_radius);
	g_xMinimumValue = min_x + (min_x>=0 ? 6*g_radius : -6*g_radius);
	g_yMinimumValue = min_y + (min_y>=0 ? 6*g_radius : -6*g_radius);
	g_zMinimumValue = min_z + (min_z>=0 ? 6*g_radius : -6*g_radius);
	drawAxisLines();
	buildAxisLabels();

	// HACK: this is a work around for cases when the scale is on
	if ($('#scale_checkbox').is(':checked')) toggleScaleCoordinates({'checked': true});
	
	// Change the css color of the 3d plot labels, set colors here because buildAxesLabels reverts color to default
	axeslabelscolor = $('#axeslabelscolor').css( "background-color" );
	axeslabelscolor_hex = $("#axeslabelscolor").spectrum("get").toHexString(true);
	$("#pc1_label").css('color', axeslabelscolor);
	$("#pc2_label").css('color', axeslabelscolor);
	$("#pc3_label").css('color', axeslabelscolor);

	resetCamera();
}

/*This function flips the lines in comparison plots when the user selects the option to negate the axes.*/
function flipEdges(axis) {
		var flippedPositions2d = new Array();
		for (var sampleKey in g_comparisonPositions){
			flippedPositions1d = []
			for (var edgePosition in g_comparisonPositions[sampleKey]){
				flippedPositions = [g_comparisonPositions[sampleKey][edgePosition][0], 
									g_comparisonPositions[sampleKey][edgePosition][1],
									g_comparisonPositions[sampleKey][edgePosition][2]]
				for (var i=0;i<axis.length;i++) {
					flippedPositions[axis[i]] *= (-1);
				}
				flippedPositions1d.push(flippedPositions);
			}
			flippedPositions2d[sampleKey] = flippedPositions1d;
		}
		removeEdges();
		drawEdges(flippedPositions2d);
}

/*Removes the lines in comparison plots so the negated lines can be drawn*/
function removeEdges() {
	for(var sample_id in g_plotEdges){
		for(var section in g_plotEdges[sample_id]){
			g_mainScene.remove(g_plotEdges[sample_id][section]);
		}
	}
}

function clean_label_refresh_axes() {
	$("#refresh_axes_label").html("");
}

function togglePlots() {

	// set some interface changes for 3D visualizations
	if(document.getElementById('pcoa').checked)
	{
		document.getElementById('pcoaPlotWrapper').className = 'plotWrapper';
		document.getElementById('pcoaoptions').className = '';
		document.getElementById('pcoaviewoptions').className = '';
		document.getElementById('pcoaaxes').className = '';
		document.getElementById('parallelPlotWrapper').className += ' invisible'
		document.getElementById('paralleloptions').className += ' invisible'

		// key menu is the default
		$("#menutabs").tabs('select',0);

		// make all tabs usable
		$("#menutabs").tabs({disabled: []});
		
		// adding ctrl-p
		g_screenshotBind = THREEx.Screenshot.bindKey(g_mainRenderer, {charCode: 16});
	}
	// changes for parallel plots
	else{
		document.getElementById('parallelPlotWrapper').className = document.getElementById('parallelPlotWrapper').className.replace(/(?:^|\s)invisible(?!\S)/ , '');
		document.getElementById('paralleloptions').className = document.getElementById('paralleloptions').className.replace(/(?:^|\s)invisible(?!\S)/ , '');
		document.getElementById('pcoaPlotWrapper').className += ' invisible'
		document.getElementById('pcoaoptions').className += ' invisible'
		document.getElementById('pcoaviewoptions').className += ' invisible'
		document.getElementById('pcoaaxes').className += ' invisible'

		// switch back to the key menu
		$("#menutabs").tabs('select',0);

		// make the visibility, scaling, labels and axes tabs un-usable
		// they have no contextualized meaning in when lookin at parallel plots
		// 0 = Key, 1 = Colors, 2 = Visibility, 3 = Scaling, 4 = Labels, 5 = Axes, 6 = View, 7 = Options
		$("#menutabs").tabs({disabled: [2,3,4,5,7]});
		
		// removing the ctrl-p 
        g_screenshotBind.unbind();
		
		colorByMenuChanged();
	}
}

function setParallelPlots() {
	g_parallelPlots = []

	// get the number of axes being presented on screen but remove the ones
	// that are represented by all of the custom axes (if there are any)
	var num_axes = g_fractionExplained.length-g_number_of_custom_axes;

	for(p in g_spherePositions){
		var dataline = []
		dataline.push(g_spherePositions[p].name)
		for(var i = 1; i < num_axes+1; i++){
			dataline.push(g_spherePositions[p]['P'+i])
		}
		g_parallelPlots.push(dataline)
	}

	pwidth = document.getElementById('pcoaPlotWrapper').offsetWidth
	pheight = document.getElementById('pcoaPlotWrapper').offsetHeight

	document.getElementById('parallelPlotWrapper').innerHTML = '<div id="parallelPlot" class="parcoords" style="width:'+pwidth+'px;height:'+pheight+'px"></div>'
}

// Resets the aspect ratio after dragging and window resize
function aspectReset() {
	winWidth = Math.min(document.getElementById('pcoaPlotWrapper').offsetWidth,document.getElementById('pcoaPlotWrapper').offsetHeight);
	winAspect = document.getElementById('pcoaPlotWrapper').offsetWidth/document.getElementById('pcoaPlotWrapper').offsetHeight;                               
	resetDivSizes(g_separator_left*$(window).width());
	containmentLeft = $(window).width()*0.5;
	containmentRight = $(window).width()*0.99;
	g_sceneCamera.aspect = winAspect;
	g_sceneCamera.updateProjectionMatrix();		

}

// Makes separator draggable and implements drag function
function separator_draggable() {
	$('.separator').draggable({
		axis: 'x',
		containment: [containmentLeft, 0, containmentRight, $(window).height()],
		helper: 'clone',
		drag: function (event, ui) {
			offset = ui.offset.left;
			if (offset > $(window).width()) {
				offset = $(window).width()*0.99;
			}
			aspectReset();
			resetDivSizes(offset);
			if (offset < $(window).width()*0.93) {
				g_separator_history = offset;
			} 
		}
	});
}
         
// Resizes plot and menu widths            
function resetDivSizes(width_left) {
	$('#plotToggle').width(width_left);
	$('#parallelPlotWrapper').width(width_left);
	$('#pcoaPlotWrapper').width(width_left);
	if(document.getElementById('parallel').checked) {
		togglePlots();
	}
	var width_right = $(window).width() - width_left - $('.separator').width()-1;                       
	$('#menu').width(width_right);
	g_separator_left = width_left/$(window).width();               
	if (g_separator_left > 1) {
		g_separator_left = 1;
	}
}

/*Builds the axes labels from ground up after changing the axes*/
function buildAxisLabels() {
	//build axis labels
	var axislabelhtml = "";
	var xcoords = toScreenXY(new THREE.Vector3(g_xMaximumValue, g_yMinimumValue, g_zMinimumValue),g_sceneCamera,$('#main_plot'));
	axislabelhtml += "<label id=\"pc1_label\" class=\"unselectable labels\" style=\"position:absolute; left:"+parseInt(xcoords['x'])+"px; top:"+parseInt(xcoords['y'])+"px;\">";
	axislabelhtml += g_pc1Label;
	axislabelhtml += "</label>";
	var ycoords = toScreenXY(new THREE.Vector3(g_xMinimumValue, g_yMaximumValue, g_zMinimumValue),g_sceneCamera,$('#main_plot'));
	axislabelhtml += "<label id=\"pc2_label\" class=\"unselectable labels\" style=\"position:absolute; left:"+parseInt(ycoords['x'])+"px; top:"+parseInt(ycoords['y'])+"px;\">";
	axislabelhtml += g_pc2Label;
	axislabelhtml += "</label>";
	var zcoords = toScreenXY(new THREE.Vector3(g_xMinimumValue, g_yMinimumValue, g_zMaximumValue),g_sceneCamera,$('#main_plot'));
	axislabelhtml += "<label id=\"pc3_label\" class=\"unselectable labels\" style=\"position:absolute; left:"+parseInt(zcoords['x'])+"px; top:"+parseInt(zcoords['y'])+"px;\">";
	axislabelhtml += g_pc3Label;
	axislabelhtml += "</label>";
	document.getElementById("axislabels").innerHTML = axislabelhtml;
}

//Unhides the info box if WebGL is disabled
function overlay() {
	overlay = document.getElementById("overlay");
	overlay.style.visibility = (overlay.style.visibility == "visible") ? "hidden" : "visible";
	parallel = document.getElementById("menu");
	parallel.style.visibility = (parallel.style.visibility == "invisible") ? "visible" : "hidden";
	separator = document.getElementById("separator");
	separator.style.visibility = (separator.style.visibility == "invisible") ? "visible" : "hidden";
	plotToggle = document.getElementById("plotToggle");
	plotToggle.style.visibility = (plotToggle.style.visibility == "invisible") ? "visible" : "hidden";
}

//Toggles fullscreen when double-clicking the separator
function separatorDoubleClick() {
	if (g_separator_left > 0.98) {
		if (g_separator_history/$(window).width() < .5) {
			g_separator_history = $(window).width()*.5;
			resetDivSizes(g_separator_history);
		}
		else {
			resetDivSizes(g_separator_history);
		}
	}
	else {
		resetDivSizes($(window).width()*0.99);
	}
	aspectReset();	
}

/*Setup and initialization function for the whole system

  This function will set all of the WebGL elements that are required to exist
  for the plot to work properly. This in turn will draw the ellipses, spheres
  and all the other elements that could be part of a plot.
*/
$(document).ready(function() {
	setJqueryUi()
	
	// Default sizes: g_separator_left is in percent and the others are in decimal
	g_separator_left = 0.73;
	g_separator_history = $(window).width()*0.73;
	containmentLeft = $(window).width()*0.5;
	containmentRight = $(window).width()*0.99;
	// Detecting that webgl is activated
	if ( ! Detector.webgl ) {
		overlay();
	}
	var main_plot = $('#main_plot');
	var particles, geometry, parameters, i, h, color;
	var mouseX = 0, mouseY = 0;

	var winWidth = Math.min(document.getElementById('pcoaPlotWrapper').offsetWidth,document.getElementById('pcoaPlotWrapper').offsetHeight), view_angle = 35, view_near = 0.0000001, view_far = 10000;
	var winAspect = document.getElementById('pcoaPlotWrapper').offsetWidth/document.getElementById('pcoaPlotWrapper').offsetHeight;

	$(window).resize(function() {	
		aspectReset();
		separator_draggable();
	});
	
	separator_draggable();
	
	// Validating the string for the saveas filename = taken from http://stackoverflow.com/questions/6741175/trim-input-field-value-to-only-alphanumeric-characters-separate-spaces-with-wi
    $('#saveas_name').keypress(function(event) {
        var code = (event.keyCode ? event.keyCode : event.which);
        if (g_validAsciiCodes.indexOf(code)==-1)
            event.preventDefault();
    });
    
    // Disables the enter key in the search bar
    $('#searchBox').keypress(function(event) {
    	if (event.keyCode == 13) {
        	event.preventDefault();
    	}
	});

	init();
	animate();

	function init() {
		// assign a position to the camera befor associating it with other
		// objects, else the original position will be lost and not make sense
		g_sceneCamera = new THREE.PerspectiveCamera(view_angle, winAspect, view_near, view_far);
		g_sceneCamera.position.set(0, 0, 0);

		$('#main_plot canvas').attr('width',document.getElementById('pcoaPlotWrapper').offsetWidth);
		$('#main_plot canvas').attr('height',document.getElementById('pcoaPlotWrapper').offsetHeight);

		g_mainScene = new THREE.Scene();
		g_mainScene.fog = new THREE.FogExp2( 0x000000, 0.0009);
		g_mainScene.add(g_sceneCamera);

		g_genericSphere = new THREE.SphereGeometry(g_radius, g_segments, g_rings);
		g_elementsGroup = new THREE.Object3D();
		g_mainScene.add(g_elementsGroup);

		drawSpheres();
		drawEllipses();
		drawTaxa();
		drawVectors();
		drawEdges(g_comparisonPositions);

		// set some of the scene properties
		g_plotIds = g_plotIds.sort();
		g_visiblePoints = g_plotIds.length;
		changePointCount(g_visiblePoints)

		// given that labels are turned off by default, leave a place holder
		var line = "";
		$("#labelcombo").append("<option>Select A Category...</option>");

		// this sorted list of headers is only used in the following loop
		// to create the 'color by', 'show by' and 'label by' drop-down menus
		sortedMappingFileHeaders = _splitAndSortNumericAndAlpha(g_mappingFileHeaders)
		for(var i in sortedMappingFileHeaders){
			var temp = [];
			for(var j in g_plotIds) {
				if(g_mappingFileData[g_plotIds[j]] == undefined){
					console.warning(g_plotIds[j] +" not in mapping")
					continue
				}
				temp.push(g_mappingFileData[g_plotIds[j]][i])
			}
			temp = dedupe(temp);
			
			// note that each category is added to all the dropdown menus in the
			// user interface, these are declared in _EMPEROR_FOOTER_HTML_STRING
			if (i==0) {
			    line = "<option selected value=\""+sortedMappingFileHeaders[i]+"\">"+sortedMappingFileHeaders[i]+"</option>"
			} else {
			    line = "<option value=\""+sortedMappingFileHeaders[i]+"\">"+sortedMappingFileHeaders[i]+"</option>"
			}
			$("#colorbycombo").append(line);
			$("#scalingbycombo").append(line);
			$("#showbycombo").append(line);
			$("#labelcombo").append(line);
		}

		setParallelPlots();

		colorByMenuChanged();
		showByMenuChanged();
		scalingByMenuChanged();

		togglePlots();

		// the light is attached to the camera to provide a 3d perspective
		g_sceneLight = new THREE.DirectionalLight(0x999999, 2);
		g_sceneLight.position.set(1,1,1).normalize();
		g_sceneCamera.add(g_sceneLight);

		// Adding camera
		g_sceneControl = new THREE.TrackballControls(g_sceneCamera, document.getElementById('main_plot'));
		g_sceneControl.rotateSpeed = 1.0;
		g_sceneControl.zoomSpeed = 1.2;
		g_sceneControl.panSpeed = 0.8;
		g_sceneControl.noZoom = false;
		g_sceneControl.noPan = false;
		g_sceneControl.staticMoving = true;
		g_sceneControl.dynamicDampingFactor = 0.3;
		g_sceneControl.keys = [ 65, 83, 68 ];

		// black is the default background color for the scene
		var rendererBackgroundColor = new THREE.Color();
		rendererBackgroundColor.setHex("0x000000");

		// renderer, the default background color is black
		g_mainRenderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
        
        // adding 'ctrl+p' to print screenshot
        g_screenshotBind = THREEx.Screenshot.bindKey(g_mainRenderer, {charCode: 16});
        
		g_mainRenderer.setClearColor(rendererBackgroundColor, 1);
		g_mainRenderer.setSize( document.getElementById('pcoaPlotWrapper').offsetWidth, document.getElementById('pcoaPlotWrapper').offsetHeight );
		g_mainRenderer.sortObjects = true;
		main_plot.append(g_mainRenderer.domElement);

		// build divs to hold point labels and position them
		var labelshtml = "";
		for(var i in g_plotIds) {
			var sid = g_plotIds[i];
			var divid = sid.replace(/\./g,'');
			mesh = g_plotSpheres[sid];
			var coords = toScreenXY(mesh.position,g_sceneCamera,$('#main_plot'));
			labelshtml += "<label id=\""+divid+"_label\" class=\"unselectable labels\" style=\"position:absolute; left:"+parseInt(coords['x'])+"px; top:"+parseInt(coords['y'])+"px;\">";
			labelshtml += sid;
			labelshtml += "</label>";
		}
		document.getElementById("labels").innerHTML = labelshtml;

		labelshtml = "";
		// add the labels with the taxonomic lineages to the taxalabels div
		for(var key in g_taxaPositions){

			// get the coordinate of this taxa sphere
			var coords = toScreenXY(g_plotTaxa[key].position,g_sceneCamera,$('#main_plot'));

			// labels are identified by the key they have in g_taxaPositions
			labelshtml += "<label id=\""+key+"_taxalabel\" class=\"unselectable labels\" style=\"position:absolute; left:"+parseInt(coords['x'])+"px; top:"+parseInt(coords['y'])+"px;\">";
			labelshtml += g_taxaPositions[key]['lineage'];
			labelshtml += "</label>";
		}
		document.getElementById("taxalabels").innerHTML = labelshtml
		
		// adding values for axes to display
		drawMenuAxesDisplayed();
		changeAxesDisplayed();
		drawAxisLines();

		buildAxisLabels();
	}

	function drawMenuAxesDisplayed() {
		if (!jQuery.isEmptyObject(g_vectorPositions) || !jQuery.isEmptyObject(g_taxaPositions) ||
			!jQuery.isEmptyObject(g_ellipsesDimensions) || g_number_of_custom_axes!=0) {
			text = '<table width="100%%">';
			text += '<tr><td><font color="red">This is disabled for custom axes, biplots, vectors, and jackknifed</font></td></tr>';
			text += '</table>';
			document.getElementById("axeslist").innerHTML = text;
			return;
		}

		text = '<table border="0" width="80%">';

		// Adding 1st axis
		text += '<tr>'
		text += '<td width="40px" class="unselectable lables">Axis 1:</td>'
		text += '<td><select id="pc1_axis" onchange="clean_label_refresh_axes();">';

		for (var i=1; i < g_fractionExplainedRounded.length + 1; i++) {
			if (i==1) {
				text += '<option selected value="P' + i + '">P' + i + " (" + g_fractionExplainedRounded[i-1] + "%)" + '</option>';
			} else {
				text += '<option value="P' + i + '">P' + i + " (" + g_fractionExplainedRounded[i-1] + "%)" + '</option>';
			}
		}
		text += '</select></td>'
		text += '<td>Negate values:<input id="flip_axes_1" class="checkbox" type="checkbox" style=""></td>';
		text += '</tr>';

		// Adding 2nd axis
		text += '<tr>'
		text += '<td width="40px" class="unselectable lables">Axis 2:</td>'
		text += '<td><select id="pc2_axis" onchange="clean_label_refresh_axes();">';
		for (var i=1; i < g_fractionExplained.length + 1; i++) {
			if (i==2) {
				text += '<option selected value="P' + i + '">P' + i + " (" + g_fractionExplainedRounded[i-1] + "%)" + '</option>';
			} else {
				text += '<option value="P' + i + '">P' + i + " (" + g_fractionExplainedRounded[i-1] + "%)" + '</option>';
			}
		}
		text += '</select></td>'
		text += '<td>Negate values:<input id="flip_axes_2" class="checkbox" type="checkbox" style=""></td>';
		text += '</tr>';

		// Adding 3rd axis
		text += '<tr>'
		text += '<td width="40px" class="unselectable lables">Axis 3:</td>'
		text += '<td><select id="pc3_axis" onchange="clean_label_refresh_axes();">';
		for (var i=1; i < g_fractionExplained.length + 1; i++) {
			if (i==3) {
				text += '<option selected value="P' + i + '">P' + i + " (" + g_fractionExplainedRounded[i-1] + "%)" + '</option>';
			} else {
				text += '<option value="P' + i + '">P' + i + " (" + g_fractionExplainedRounded[i-1] + "%)" + '</option>';
			}
		}
		text += '</select></td>'
		text += '<td>Negate values:<input id="flip_axes_3" class="checkbox" type="checkbox" style=""></td>';
		text += '</tr>';
		text += '</table>';

		// Adding button
		text += '<table width="100%%"><tr>';
		text += '<td width="20px"><input type="button" value="Refresh" onclick="changeAxesDisplayed();"></td>';
		text += '<td id="refresh_axes_label"></td></tr>';
		text += '</table>';
		document.getElementById("axeslist").innerHTML = text;
	}

	function animate() {
		requestAnimationFrame( animate );

		render();

		var labelCoordinates;

		// reposition the labels for the axes in the 3D plot
		labelCoordinates = toScreenXY(new THREE.Vector3(g_xMaximumValue, g_yMinimumValue, g_zMinimumValue), g_sceneCamera,$('#main_plot'));
		$("#pc1_label").css('left', labelCoordinates['x'])
		$("#pc1_label").css('top', labelCoordinates['y'])
		labelCoordinates = toScreenXY(new THREE.Vector3(g_xMinimumValue, g_yMaximumValue, g_zMinimumValue), g_sceneCamera,$('#main_plot'));
		$("#pc2_label").css('left', labelCoordinates['x'])
		$("#pc2_label").css('top', labelCoordinates['y'])
		labelCoordinates = toScreenXY(new THREE.Vector3(g_xMinimumValue, g_yMinimumValue, g_zMaximumValue), g_sceneCamera,$('#main_plot'));
		$("#pc3_label").css('left', labelCoordinates['x'])
		$("#pc3_label").css('top', labelCoordinates['y'])


		// move labels when the plot is moved
		if(document.plotoptions.elements[0].checked){
			for(var i in g_plotIds) {
				var sid = g_plotIds[i];
				mesh = g_plotSpheres[sid];
				var coords = toScreenXY(mesh.position, g_sceneCamera, $('#main_plot'));
				var divid = sid.replace(/\./g,'');
				$('#'+divid+"_label").css('left',coords['x']);
				$('#'+divid+"_label").css('top',coords['y']);
			}
		}
		// check if you have to reposition the taxa labels for each frame
		// this is something that will only happen when drawing biplots
		if(document.biplotoptions){
			if(document.biplotoptions.elements[0].checked){
				for(var key in g_taxaPositions) {
					// retrieve the position of the taxa on screen
					var coords = toScreenXY(g_plotTaxa[key].position,
						g_sceneCamera, $('#main_plot'));

					// add the label at the appropriate position
					$('#'+key+"_taxalabel").css('left',coords['x']);
					$('#'+key+"_taxalabel").css('top',coords['y']);
				}
			}
		}
		if(g_foundId) {
			var coords = toScreenXY(g_plotSpheres[g_foundId].position, g_sceneCamera, $('#main_plot'));
			$('#finder').css('left',coords['x']-15);
			$('#finder').css('top',coords['y']-5);
		}
	}
   
	function render() {
		g_sceneControl.update();
		g_mainRenderer.setSize( document.getElementById('pcoaPlotWrapper').offsetWidth, document.getElementById('pcoaPlotWrapper').offsetHeight );
		g_mainRenderer.render( g_mainScene, g_sceneCamera);
	}
	
});

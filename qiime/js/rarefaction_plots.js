//rarefaction_plots.js

window.onload=Main;
google.load("visualization", "1", {packages:["scatterchart"]});

setTimeout("make_collapsed_series()",0);

var colours = new Array();
var currentGraph = "";
var showHideOps = new Array();
var categories = new Array();
var categoryOps = new Array();
var sampleMat = new Array();
var sampleIDsarry = new Array();
var fileNames = new Array;
var collapsed_series_array = new Array();
var showonly_array = new Array();
var rareFilesData = new Array();
var contCalcs = new Array();

function Main() {
	loadExternalFiles();
	loadMapping(fileNames[0]);
    loadDiversityMeasures(fileNames.slice(1,fileNames.length));
    currentGraph = fileNames[1];
    makeDataTable();
    make_collapsed_series();
    makeOpsPanel();
    //make_showonly_array();
    //loadGraph(fileNames[1]);
}

function loadExternalFiles() {
	var datafiles = getFile("dataFilesforJS.txt");
	var filelines = datafiles.split("\n");
	for(var i = 0; i < filelines.length; i++)
		fileNames.push(filelines[i]);
	fileNames.pop()
}

//+ Carlos R. L. Rodrigues
//@ http://jsfromhell.com/array/average [rev. #1]

function average(a){
    var r = {mean: 0, variance: 0, deviation: 0}, t = a.length;
    for(var m, s = 0, l = t; l--; s += a[l]);
    for(m = r.mean = s / t, l = t, s = 0; l--; s += Math.pow(a[l] - m, 2));
    return r.deviation = Math.sqrt(r.variance = s / t), r;
}
// end

function makeDataTable() {
	//$('#raretable').dataTable();
}

function plotLines(series, seriestitles, xaxisvals, gtitle, legend){
    //colours = new Array();
    /*var c = ['#0000FF','#5500AA','#AA0055','#FF0000','#BF3F00','#7F7F00','#3FBF00','#00FF00']
	//var newcols = new Array();
	for(var i = 0; i < c.length; i++)
	{
	    for(var j = 0; j < 12; j++)
	        colours.push(c[i]);
	}*/
	    
	//console.log(series)
    var data = new google.visualization.DataTable();
	data.addColumn('number', 'Sequences Per Sample');
	
	for(var i = 0; i < series.length; i++)
	{
		data.addColumn('number',seriestitles[i]);
	}
	data.addRows(xaxisvals.length+1);
	for(var i = 0; i < xaxisvals.length; i++)
	{
		//console.log(xaxisvals[i])
		data.setValue(i, 0, xaxisvals[i]);
	}
	
	for(var i = 0; i < series.length; i++)
	{
		for(var j = 0; j < series[i].length; j++)
		{
			if(series[i][j] != 0)
				data.setValue(j, i+1, series[i][j]); // i+1 because 0 is for the xaxisvals
		}
	}

	document.getElementById('graph_container').innerHTML = "";
	var chart = new google.visualization.ScatterChart(document.getElementById('graph_container'));
    chart.draw(data, {width: 750, height: 400, lineSize: 1, pointSize: 1.5, legend: legend, colors: colours, title: gtitle, titleX: 'Sequences Per Sample', titleY: 'Diversity Measure'});
}

function loadDiversityMeasures(filenms) {
	//var filenms = ['PD_whole_tree_mouth.txt','chao1_mouth.txt'];
	var headers = ['SampleIDs'];
	var tempaave = new Array();
	for(f = 0; f < filenms.length; f++)
	{
		var tempData = new Array();
		var tempave = new Array();
		headers.push(filenms[f]);
		rareFilesData[filenms[f]] = loadFileData(filenms[f]);
		/*
		data.push(raremat); [0]
		data.push(rareIDs); [1]
		data.push(sampleIDs); [2]
		data.push(seqsPerSamp); [3]
		*/
		tempave = makeAverageSeries(rareFilesData[filenms[f]][0], rareFilesData[filenms[f]][1], rareFilesData[filenms[f]][2]);	
		sampleIDs = rareFilesData[filenms[f]][2]
		tempaave[f] = aveAverageSeries(tempave);
	}
	makeTable(tempaave, headers, sampleIDs);
}

function loadMapping(filenm) {
    sampleMat = new Array();
    sampleIDsarry = new Array();
    
	var mappingfl = getFile(filenm);
	//var mappingFileLines = mappingfl.split("\n");
	var mappingFileLines = mappingfl.split("\n");
    categories = mappingFileLines[0].split("\t");
    var categ = "";
    
	categoryOps = new Array();

    for(var i = 1; i < mappingFileLines.length; i++)
    {
        categ = mappingFileLines[i].split("\t");
        sampleIDsarry[i-1] = categ[0].toUpperCase();
        
        sampleMat[sampleIDsarry[i-1]] = new Array(); //categ[0] = sampleID
        for(var j = 0; j < categ.length; j++)
        {
			if(categoryOps[categories[j]] != null)
			{
				var seen = false;
				for(var k = 0; k < categoryOps[categories[j]].length; k++)
				{
					if(categoryOps[categories[j]][k] == categ[j])
						seen = true;
				}
				if(!seen)
					categoryOps[categories[j]].push(categ[j]);
			}
			else
				categoryOps[categories[j]] = [categ[j]]; //
            sampleMat[categ[0].toUpperCase()][j] = categ[j] ;
        }
    }
}

function contains(a, obj){
  for(var i = 0; i < a.length; i++) {
    if(a[i] == obj){
      return true;
    }
  }
  return false;
}

function makeErrorBarSeries(series){
	var errSeries = new Array();
	var aveSer = new Array();
	for(var col = 0; col < series[0].length; col++)
	{
		var stdDev = new Array();
		var values = new Array();
		
		for(var m = 0; m < series[0].length*2; m++) // init std devs all to zero
			stdDev[m] = 0;
		
		for(var row = 0; row < series.length; row++)
		{
			var value = series[row][col];
			values.push(value)
		}
		aveSer.push(average(values).mean)
		aveSer.push(average(values).mean)
		stdDev[col*2] = average(values).mean - average(values).deviation;
		stdDev[col*2+1] = average(values).mean + average(values).deviation;
		errSeries.push(stdDev)
	}
	errSeries.push(aveSer)
	return errSeries;
}

function get_showOnly(option, catOp) {
	var catToColor = new Array();
	colours = new Array();
	var data = rareFilesData[currentGraph];
	var ave = makeAverageSeries(data[0], data[1], data[2]);
	var series = new Array();
		
	for(var i = 0; i < categoryOps[option].length; i++)
		catToColor[categoryOps[option][i]] = makeColorGradient(.3,.3,.3,0,2,4,categoryOps[option].length)[i]
        
    
	var key = categories.indexOf(option)
	
	colours.push(catToColor[categoryOps[option][categoryOps[option].indexOf(catOp)]])
	
	for(var j = 0; j < sampleIDsarry.length; j++)
	{
		if(sampleMat[sampleIDsarry[j]][key] == catOp)
			series.push(ave[j])
	}
	
	var newSeries = makeErrorBarSeries(series);
	
	var newXaxis = new Array();
	for(var p = 0; p < data[3].length; p++)
	{
		newXaxis[p*2] = data[3][p]
		newXaxis[p*2+1] = data[3][p]
	}
	/*
	var c = ['#0000FF','#5500AA','#AA0055','#FF0000','#BF3F00','#7F7F00','#3FBF00','#00FF00']
	//var newcols = new Array();
	for(var i = 0; i < 12*8; i++)
	    colours.push(c[i/8]);
	//colours = newcols;
	*/
	
	var graphName = " Average Colored By " + option + " : " + catOp;
	//plotLines(newSeries, new Array(newSeries.length), newXaxis, currentGraph.split('.')[0]+ graphName, 'none')
    var result = new Array();
    result.push(newSeries);
    result.push(new Array(newSeries.length))
    result.push(newXaxis);
    result.push(graphName)
    result.push(colours)
    return result;
}

function showonly_runtime(option, catOp) {
    var ary = get_showOnly(option, catOp);
    colours = ary[4]
    plotLines(ary[0], ary[1], ary[2], ary[3], 'none')
}

function showOnly(option, catOp) {
    var ary = showonly_array[option][catOp];
    colours = ary[4]
    plotLines(ary[0], ary[1], ary[2], ary[3], 'none')
}

function collapseSeries(option) {
	colorBy(option); // make sure graph is set to this option
	var catToColor = new Array();
	colours = new Array();

	var categoryArray = new Array();
	var key = "";
	
	var data = rareFilesData[currentGraph];
	var ave = makeAverageSeries(data[0], data[1], data[2]);
		
	for(var j = 0; j < sampleIDsarry.length; j++)
	{
		var vals = sampleMat[sampleIDsarry[j]]; //all different options to color by (ie. sex)
		key = vals[categories.indexOf(option)]
		//alert(key)

		if(categoryArray[key] == null)
			categoryArray[key] = new Array();
		categoryArray[key].push(ave[j]); // (ie. categoryArray['M'].push(...))
	}
	
	var aveCatArray = new Array();
	
	for(var i = 0; i < categoryOps[option].length; i++)
		catToColor[categoryOps[option][i]] = makeColorGradient(.3,.3,.3,0,2,4,categoryOps[option].length)[i];
	
	var newSeries = new Array();
	var stdDevSeries = new Array();
	var seriesNames = new Array();//categoryOps[option].slice();
	//document.getElementById('debugging').innerHTML += "option: " + option;
	//document.getElementById('debugging').innerHTML += "catarray def " + categoryArray[categoryOps[option][0]][0];
	//document.getElementById('debugging').innerHTML += categoryOps[option]
	for(var i = 0; i < categoryOps[option].length; i++) // ie male or female so this would go two times
	{
		// need to do it for each first point of each series, then each second point, etc
		var currAve = new Array();
		colours.push(catToColor[categoryOps[option][i]]);
		seriesNames.push(categoryOps[option][i]);
		if(categoryArray[categoryOps[option][i]][0] == null)
		    document.getElementById('debuggingconsole').innerHTML += categoryOps[option][i]
		for(var l = 0; l < categoryArray[categoryOps[option][i]][0].length; l++) // for length of sequence line
		{
			var values = new Array();
			var stdDev = new Array();
			for(var m = 0; m < categoryArray[categoryOps[option][i]][0].length*2; m++)
				stdDev[m] = 0;
				
			for(var k = 0; k < categoryArray[categoryOps[option][i]].length; k++) // for each sequence in each option (ie M, F)
			{
				var value = categoryArray[categoryOps[option][i]][k][l];
				values.push(value);
			}
			// Have to add these values twice in order to get vertical
			// error bars
			currAve.push(average(values).mean);
			currAve.push(average(values).mean);
			stdDev[l*2] = average(values).mean + average(values).deviation;
			stdDev[l*2+1] = average(values).mean - average(values).deviation;
			colours.push(catToColor[categoryOps[option][i]]);
			newSeries.push(stdDev);
			seriesNames.push("");
		}
		newSeries.push(currAve);
	}
	
	var newXaxis = new Array();
	for(var p = 0; p < data[3].length; p++)
	{
		newXaxis[p*2] = data[3][p];
		newXaxis[p*2+1] = data[3][p];
	}
	/*
	var c = ['#0000FF','#5500AA','#AA0055','#FF0000','#BF3F00','#7F7F00','#3FBF00','#00FF00']
	//var newcols = new Array();
	for(var i = 0; i < 12*8; i++)
	    colours.push(c[i/8]);
	//colours = newcols;*/
	var graphName = currentGraph.split('.')[0]+ ' Average Colored By ' + option;
	//plotLines(newSeries, seriesNames, newXaxis, graphName, 'none')
	var result = new Array();
	result.push(newSeries);
	result.push(seriesNames);
	result.push(newXaxis);
	result.push(graphName);
	result.push(catToColor);
	return result
}

function make_collapsed_series() {
    for(var i = 0; i < categories.length; i++)
    {
        if(categoryOps[categories[i]].length == 1 || categoryOps[categories[i]].length == sampleIDsarry)
            continue;
        
        if(categoryOps[categories[i]].length > 500)
        {
            var r = confirm(categories[i] + " has > 500 categories to compare, would you like to continue calculations?");
            contCalcs[categories[i]] = r;
            if(r!=true)
                continue;
        }
        collapsed_series_array[categories[i]] = collapseSeries(categories[i]);
    }
}

function show_collapsed_series(option) {
    colorBy(option); // make sure graph is set to this option
    var ary = collapsed_series_array[option]
    colours = new Array();
    for(var i = 0; i < categoryOps[option].length; i++)
    {
        for(var j = 0; j < (ary[0][0].length+2)/2; j++)
            colours.push(ary[4][categoryOps[option][i]])
    }
    //colours = ary[4]
    plotLines(ary[0], ary[1], ary[2], ary[3], 'none')
}

function showHide(option) {
	for(var i = 0; i < showHideOps.length; i++)
	{
		document.getElementById(showHideOps[i]).style.display = 'none';
	}
	document.getElementById(option).style.display = 'block';
	show_collapsed_series(option);
}

function RGB2Color(r,g,b)
{
	return ' ' + byte2Hex(r) + byte2Hex(g) + byte2Hex(b);
}

function byte2Hex(n)

{
	var nybHexString = "0123456789ABCDEF";
	return String(nybHexString.substr((n >> 4) & 0x0F,1)) + nybHexString.substr(n & 0x0F,1);
}

function makeColorGradient(frequency1, frequency2, frequency3,phase1, phase2, phase3, center, width, len)
{
	var col = new Array();
	if (len == undefined) len = 50;
	if (center == undefined) center = 128;
	if (width == undefined) width = 127;

	for (var i = 0; i < len; ++i)
	{
		var red = Math.sin(frequency1*i + phase1) * width + center;
		var grn = Math.sin(frequency2*i + phase2) * width + center;
		var blu = Math.sin(frequency3*i + phase3) * width + center;
		col.push(RGB2Color(red,grn,blu));
	}
	var c = ['#0000FF','#5500AA','#AA0055','#FF0000','#BF3F00','#7F7F00','#3FBF00','#00FF00']
    var cols = new Array();
	for(var i = 0; i < len; i++)
	    cols.push(c[i%8])
	return cols;
}

function colorBy(option) {
	//console.log(option)
	var catToColor = new Array();
	colours = new Array();
	
	for(var i = 0; i < categoryOps[option].length; i++)//catToColor[categoryOps[option][i]] = '#'+Math.floor(Math.random()*16711680+255).toString(16);
		catToColor[categoryOps[option][i]] = makeColorGradient(.3,.3,.3,0,2,4,categoryOps[option].length)[i]
		
	for(var j = 0; j < sampleIDsarry.length; j++)
	{
		var vals = sampleMat[sampleIDsarry[j]]; //all different options to color by (ie. sex)
		var key = "";
		for(var k = 0; k < vals.length; k++)
		{
			if(categories[k] == option)
				key = vals[k] // colorby value for this specific sampleID (ie. male or female)
		}
		colours.push(catToColor[key]) // for each sample, push on the color corresponding to the colorby value
	}
	/*
	var c = ['#0000FF','#5500AA','#AA0055','#FF0000','#BF3F00','#7F7F00','#3FBF00','#00FF00']
	var newcols = new Array();
	for(var i = 0; i < 12*8; i++)
	    colours.push(c[i/8]);*/
}

function recolor(option) {
	colorBy(option)
	loadGraph(currentGraph)
}

function makeOpsPanel() {
	var htmllines = "<div class=\"ops clearfix\">";
	
	for(var i = 0; i < categories.length-1; i++)
	{
	    if(categoryOps[categories[i]].length == 1)
	    {
	        document.getElementById('debuggingconsole').innerHTML += "Ommited mapping file column "+categories[i]+" only one value to compare<br>";
            continue;
        }
        else if(categoryOps[categories[i]].length == sampleIDsarry.length)
	    {
	        document.getElementById('debuggingconsole').innerHTML += "Ommited mapping file column "+categories[i]+" values are equal to number of samples<br>";
            continue;
        }
        else if(categoryOps[categories[i]].length > 500 && !contCalcs[categories[i]])
        {
            document.getElementById('debuggingconsole').innerHTML += "Ommited mapping file column "+categories[i]+" > 500 categories<br>";
            continue;
	    }
		showHideOps.push(categories[i]);
		htmllines += "<a class=\"ops\" href=\"javascript:showHide(\'"+categories[i]+"\');\">"+categories[i]+"</a>";
		htmllines += "<a class=\"ops\" onMouseOver=\"javascript:show_collapsed_series(\'"+categories[i]+"\');\"> ave </a>";
		//htmllines += "<a class=\"ops\" onMouseOver=\"javascript:recolor(\'"+categories[i]+"\');\"> all </a>";
		htmllines += "<ul class=\"onoff\">";
		htmllines += "<span style=\"display:none\" id=\""+categories[i]+"\">";
		for(j = 0; j < categoryOps[categories[i]].length; j++)
			htmllines += "<li><a onMouseOver=\"javascript:showonly_runtime(\'"+categories[i]+"\',\'"+categoryOps[categories[i]][j]+"\');\">"+categoryOps[categories[i]][j]+"</a></li>"
		
		htmllines += "</span>";
		htmllines += "</ul>";
		htmllines += "<br />"
	}
	htmllines += "</div>"
	
	document.getElementById("colorops").innerHTML += htmllines;
	//colorBy(categories[1]); // 0 = sampleID
}

function loadGraph(filenm) {
	currentGraph = filenm;
	make_collapsed_series();
	//var data = loadFileData(filenm);
	//var ave = makeAverageSeries(data[0], data[1], data[2]);	
	//plotLines(ave, data[2], data[3], filenm.split('.')[0], 'right');
	//show_collapsed_series('#SampleID');
}

function loadFileData(filenm) {
    var rarefl = getFile(filenm);
    var rareFileLines = rarefl.split("\n");
    var headerln = rareFileLines[0];
    var headers = headerln.split("\t");
    var sampleIDs = headers.slice(); // copy headers array
    sampleIDs.shift(); // get rid of #
    sampleIDs.shift(); // get rid of seqs per sample
    sampleIDs.shift(); // get rid of iteration #
    for(var i = 0; i < sampleIDs.length; i++)
        sampleIDs[i] = sampleIDs[i].toUpperCase();
    rareFileLines.shift(); // get rid of header line so that for loop can start at 0
    var rareMat = new Array();
    var rareIDs = new Array();
	var seqsPerSamp = new Array();
    var vals = new Array();
    /*
    var cnt = 0;
    var iterNum = rareFileLines[cnt][2]; // iteration number
    var seen = new Array();
    while(!contains(seen, iterNum))
    {
        seen.push(iterNum);
        cnt += 1;
        iterNum = rareFileLines[cnt][2]
    }
    
    var maxIterations = cnt - 1;
    */
    
    var current = rareFileLines[0].split("\t")[1]; // second item on line = seqsPerSample
    var next;

    for(var i = 0; i < rareFileLines.length; i++)
    {
        vals = rareFileLines[i].split("\t");
        rareIDs[i] = vals[0]; // first item on the line should be rarefaction 
        next = vals[1];
        if(next != current)
        {
            seqsPerSamp.push(Number(current))
            current = next;
        }
        
		/*if(i%maxIterations == 0)
			seqsPerSamp[i/maxIteratons] = Number(vals[1]); // second item on line should be seqs per sample
            */
        
        rareMat[rareIDs[i]] = new Array();
        for(var j = 0; j < headers.length; j++)
        {
            if(isNaN(Number(vals[j])))
            {
				if(vals[j] == 'n/a')
				{
					rareMat[rareIDs[i]][j] = Number(0);
				}
				else
                	rareMat[rareIDs[i]][j] = vals[j];
            }
            else
            {
                rareMat[rareIDs[i]][j] = Number(vals[j]);
            }
        }    
    }
    seqsPerSamp.push(Number(next))
    
    /*
    var toRemove = new Array();
    // need to go through and get rid of rarefaction vals for seqIDs not found in mapping file
    for(var i = 0; i < sampleIDs.length; i++)
    {
        if(!contains(sampleIDsarry, sampleIDs[i]))
        {
            toRemove.push(i);
            document.getElementById('debuggingconsole').innerHTML += 'ID found in '+filenm +' not found in mapping file: ' + sampleIDs[i] +'<br>'
        }
    }
    
    for(j = 0; j < toRemove.length; j++)
    {
        for(var i = 0; i < rareIDs.length; i++)
        {
            rareMat[rareIDs[i]].splice(toRemove[j]-j,1);
            delete sampleMat[rareIDs[i]];
            sampleIDs.splice(toRemove[j]-j,1);
        }
    }*/

   	var data = new Array();
	data.push(rareMat);
	data.push(rareIDs);
	data.push(sampleIDs);
	data.push(seqsPerSamp);
	return data;
}

function makeTable(rarevals, headers, samIDs) {
    var tablelines = "<table id=\"raretable\"><thead><tr>";
    for(var i = 0; i < headers.length; i++)
    {
        tablelines += "<th>"
		if(i != 0)
		{
			tablelines += "<a href=\"javascript:loadGraph(\'"+headers[i]+"\');\">"+headers[i]+ "</a>"
		}
		else
		{
			tablelines += headers[i]
		}
		tablelines += "</th>";
    }
    tablelines += "</tr></thead><tbody>";
    for(var i = 0; i < samIDs.length; i++)
    {
        tablelines += "<tr><td>"+samIDs[i]+"</td>";
        for(var j = 0; j < rarevals.length; j++)
        {   
            try{
                tablelines += "<td>"+rarevals[j][i].toFixed(2)+"</td>";
            }
            catch(err){
                tablelines += "<td>"+rarevals[j][i]+"</td>";
            }
        }
        tablelines += "</tr>";
    }
    tablelines += "</tbody></table>";
    document.getElementById('raretablewrapper').innerHTML = tablelines;
}

function aveAverageSeries(ave) {
	var newmat = Array();
	var sum = 0;
	for(var i = 0; i < ave.length; i++)
	{
		for(var j = 0; j < ave[i].length; j++)
			sum += ave[i][j];
		newmat[i] = Number(sum);
		sum = 0;
	}
	return newmat;
}

function makeAverageSeries(matrix, IDs, samIDs) {
	var sseries = makeSeries(matrix, IDs, samIDs);
	var aveser = new Array();
	for(var k = 0; k < sseries.length; k++) 		// for each sample
	{
		aveser[k] = new Array(); 			// aveser[k] = []
		var curriter;
		var nextiter;
		var sum = 0;
		var num = 0;
		curriter = String(sseries[k][0][0]);
		for(var l = 0; l < sseries[k].length; l++) 	// for each iteration and number of seqs/sample
		{
			nextiter = String(sseries[k][l][0]);
			if(nextiter != curriter)
			{
				aveser[k][l/num-1] = sum/num;
				sum = 0;
				num = 0;
				curriter = nextiter;
			}
			else
			{
				sum += sseries[k][l][1];
			}
			num += 1;
		}
		aveser[k][sseries[k].length/num-1] = sum/num;
	}
	return aveser;
}

function makeAverageSeries1(matrix, IDs, samIDs) {
	var sseries = makeSeries(matrix, IDs, samIDs)
	//console.log(samIDs)
	var aveser = new Array();
	for(var k = 0; k < sseries.length; k++) 		// for each sample
	{
		aveser[k] = new Array(); 			// aveser[k] = []
		var curriter;
		var nextiter;
		var sum = 0;
		curriter = String(sseries[k][0][0])
		for(var l = 0; l < sseries[k].length; l++) 	// for each iteration and number of seqs/sample
			{
				nextiter = String(sseries[k][l][0]);
				if(nextiter != curriter)
				{
					aveser[k][l/10] = sum/10
					sum = 0
					curriter = nextiter
				}
				else
				{
					sum += sseries[k][l][1]
				}
			}
		aveser[k].shift()
		//console.log(aveser[k])
	}
	//console.log(aveser)
	return aveser;
}

//this function transposes the matrix
function makeSeries(matrix, IDs, samIDs) {
    var series = new Array(); 	
	var sseries = new Array();
	//document.getElementById('debugging').innerHTML += matrix[IDs[1]];
	for(var j = 3; j < matrix[IDs[0]].length; j++){
		sseries[j-3] = new Array();
		for(i = 0; i < IDs.length; i++)
			sseries[j-3][i] = new Array(matrix[IDs[i]][1], matrix[IDs[i]][j]);
	}
	//console.log(sseries[0])
    return sseries;
}

function getFile(fileName) {
    oxmlhttp = null;
    try{
        oxmlhttp = new XMLHttpRequest();
        oxmlhttp.overrideMimeType("text/xml");
    }
    catch(e){
        try{
            oxmlhttp = new ActiveXObject("Msxml2.XMLHTTP");
        }
        catch(e){
            return null;
        }
    }
    if(!oxmlhttp) return null;
    try{
       oxmlhttp.open("GET",fileName,false);
       oxmlhttp.send(null);
    }
    catch(e){
       return null;
    }
    return oxmlhttp.responseText;
}
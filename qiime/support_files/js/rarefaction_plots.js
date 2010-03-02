//rarefaction_plots.js

window.onload=Main;
var fileNames = new Array;
var graphNames = new Array;
var header = new Array;
var rarefaction_table = new Array;

function Main() {
	loadExternalFiles();
	makeOpsPanel();
	makeTable();
	$('#raretable').dataTable();
}

function loadExternalFiles() {
	var datafiles = getFile("graphNames.txt");
	var filelines = datafiles.split("\n");
	for(var i = 0; i < filelines.length; i++)
		graphNames.push(filelines[i]);
	graphNames.pop()
	
	var tablefile = getFile("rarefactionTable.txt");
	if(tablefile)
	{
	    var tablelines = tablefile.split("\n");
    	var current_header = "";
    	for(var i = 0; i < tablelines.length; i++)
    	{
    	    if(tablelines[i].match("^"+"#")=="#")
    	    {
    	        current_header = tablelines[i].split("#")[1];
    	        header.push(current_header)
    	        rarefaction_table[current_header] = new Array;
            }
            else
                rarefaction_table[current_header].push(tablelines[i])
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

function changeGraph(graph) {
    document.getElementById("graph_container").innerHTML = "<img src=\'"+graph+"\'>";
}

function makeOpsPanel() {
	var htmllines = "<div class=\"ops clearfix\">";
	
	htmllines += "<ul class=\"ops\">";
	for(var i = 0; i < graphNames.length; i++)
	{
	    htmllines +="<li>"
		htmllines += "<a class=\"ops\" onmouseover=\"javascript:changeGraph(\'"+graphNames[i]+"\')\">"+graphNames[i].split("/")[1].split(".")[0]+"</a>";
		htmllines += "<br /></li>";
	}
	htmllines+= "</ul>";
	htmllines += "</div>";
	
	document.getElementById("colorops").innerHTML += htmllines;
}

function makeTable() {
    var tablelines = "<table id=\"raretable\"><thead><tr>";
    for(var i = 0; i < header.length; i++)
    {
        tablelines += "<th>"
		if(i != 0)
		{
			tablelines += "<a href=\"javascript:loadGraph(\'"+header[i]+"\');\">"+header[i]+ "</a>"
		}
		else
		{
			tablelines += header[i]
		}
		tablelines += "</th>";
    }
    tablelines += "</tr></thead><tbody>";
    for(var i = 0; i < rarefaction_table[header[0]].length; i++)
    {
        tablelines += "<tr>";
        for(var j = 0; j < header.length; j++)
        {
            tablelines += "<td>"+rarefaction_table[header[j]][i]+"</td>";
        }
        tablelines += "</tr>";
    }
    tablelines += "</tbody></table>";
    document.getElementById('raretablewrapper').innerHTML = tablelines;
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
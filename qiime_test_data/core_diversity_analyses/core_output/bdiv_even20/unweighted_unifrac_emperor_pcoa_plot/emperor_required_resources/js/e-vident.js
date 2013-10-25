var slidervals = [];
var columnsIndex = 0;
var studydata = [];
var demodata = [];
var validColumns = [];

//jquery to activate the accordion info
$(function() {
	$("#accordion").accordion({
		collapsible: true,
		fillSpace: true
	});
});

//jquery to activate the inner accordion info
$(function() {
	$("#innerAccordion").accordion({
		collapsible: true,
	});
});

//jquery to activate the inner inner accordion info
$(function() {
	$("#innerAccordion2").accordion({
		collapsible: true,
	});
});

//jquery to activate the tabbed pane
$(function() {
	$("#tabs").tabs();
	$(".tabs-bottom .ui-tabs-nav, .tabs-bottom .ui-tabs-nav > *")
		.removeClass("ui-corner-all ui-corner-top")
		.addClass("ui-corner-bottom");
});

//jquery to make the tabs closeable
$("#tabs span.ui-icon-close").live("click", function() {
	var index = $("li", $("#tabs")).index($(this).parent());
	$("#tabs").tabs("remove", index);
	enableOptimize();
});

//function to add new tabs to the tabbed pane
function addTab(tid, title, content) {
    $('#tabs').tabs({
    tabTemplate: '<li><a id="#' + tid + '" href="#{href}">#{label}</a><span class=\'ui-icon ui-icon-close\'>Remove Tab</span></li>',
    add: function(event, ui) {
      var dataString = content;
      $(ui.panel).append(content);
    }
  });
  $('#tabs').tabs('add', '#'+tid, title);
  document.getElementById(tid).setAttribute('tabindex', 0);
  $(".tabs-bottom .ui-tabs-nav, .tabs-bottom .ui-tabs-nav > *")
		.removeClass("ui-corner-all ui-corner-top")
		.addClass("ui-corner-bottom");
}

//function called when html is loaded to populate the study select box
function load() {   
    $.ajax({ url: 'lib.psp',
                data: {fn: 'loadStudyData'},
                success: loadStudyData});
}

function loadStudyData(response) {
    rs = eval('('+response+')')
    demodata = rs['demo']
    studydata = rs['full']
    var html = '<option value=\"\">Select...</option>';
    for(var i in studydata)
        html += '<option value=\"'+i+'\">'+studydata[i][0]+'</option>';
    document.getElementById('studycombobox').innerHTML = html;
    disableInterface();
}

function enableOptimize() {
    document.getElementById('optimize').disabled = false;
}

function disableInterface() {
    $("#subjectslider").slider("disable");
        $("#sampleslider").slider("disable");
        $("#sequenceslider").slider("disable");
        $("#iterationslider").slider("disable");
        document.getElementById('optimize').disabled = true;
        $("#columns").css("color",'#ccc');
        for(i=0; i< document.visualizations.length; i++){
          document.visualizations[i].disabled = true;
       }
}

function enableInterface() {
        $("#subjectslider").slider("enable");
        $("#sampleslider").slider("enable");
        $("#sequenceslider").slider("enable");
        $("#iterationslider").slider("enable");
        document.getElementById('optimize').disabled = false;
        $("#columns").css("color",'#000');
        for(i=0; i< document.visualizations.length; i++){
          document.visualizations[i].disabled = false;
       }

       toggleSliders();
}

function toggleSliders() {
    var studyname = document.getElementById('studycombobox')[document.getElementById('studycombobox').selectedIndex].value;
    var lines = studydata[studyname][1].split('?');
    //first line is a header
    lines.shift();
    slidervals = [];

    for(var i=0; i < lines.length; i++)
    {
        slidervals.push(lines[i].split('\t'));
    }
    
    if(document.visualizations[0].checked)
    {        
        var seqs = parseInt(demodata[studyname]);     
        $("#sequenceslider").slider('option','max', parseInt(slidervals[slidervals.length-1][0]));
        $("#sequenceslider").slider('option','value', seqs);
        if(seqs == 0)
            document.getElementById('sequences').innerHTML = 'N/A'
            
        $("#subjectslider").slider("disable");
        $("#sampleslider").slider("disable");
        $("#sequenceslider").slider("disable");
        $("#iterationslider").slider("disable");
    }
    else
    {
        $("#subjectslider").slider("enable");
        $("#sampleslider").slider("enable");
        $("#sequenceslider").slider("enable");
        $("#iterationslider").slider("enable");
        
        $("#sequenceslider").slider('option','max', parseInt(slidervals[slidervals.length-1][0]));
        $("#sequenceslider").slider('option','value', parseInt(slidervals[0][0]));
    }
}

//function called to populate the interface sliders with corresponding
//info each time the study is changed
function studyChanged() {
    if(document.getElementById('studycombobox').selectedIndex == 0)
    {
        disableInterface();
        return;
    }
    
    enableInterface();
    
    //remove possible alert box around study selector
    document.getElementById('studycombop').className = document.getElementById("studycombobox").className.replace
      (/(?:^|\s)alert(?!\S)/ , '');
    
    if(document.visualizations[0].checked)
    {
        
    }
    else {    
        $("#sequenceslider").slider('option','max', parseInt(slidervals[slidervals.length-1][0]));
        $("#sequenceslider").slider('option','value', parseInt(slidervals[0][0]));
    }
    toggleSliders();
}

//remove alert box around visualizations
$(function(){
    $('input.checkbox').click(function(){
        document.getElementById('optimize').disabled = false;
        document.getElementById('opswrapper').className = document.getElementById('opswrapper').className.replace
      (/(?:^|\s)alert(?!\S)/ , '');
    });
});


//function to catch the optimize button and call python subsampling -> optimize
$(function(){
    $('input.optimize').click(function(){
        //remove alert box around study selector
        document.getElementById('studycombop').className = document.getElementById("studycombobox").className.replace
      (/(?:^|\s)alert(?!\S)/ , '');
        
        //if study selected index is 0, they have not selected a study
        if(document.getElementById('studycombobox').selectedIndex == 0)
        {
            //add css to put alert box around study selector
            document.getElementById('studycombop').className += " alert";
            document.getElementById('optimize').disabled = true;
            return;
        }
        
        //remove alert box around visualizations
      document.getElementById('opswrapper').className = document.getElementById('opswrapper').className.replace
      (/(?:^|\s)alert(?!\S)/ , '');
      //disable optimize button so it can only be clicked once
      document.getElementById('optimize').disabled = true;
      
      var checkedcnt = 0;
      
      //figure out how many frames need to be built
      for(i=0; i< document.visualizations.length; i++){
          if(document.visualizations[i].checked)
            checkedcnt += 1;
       }
      
      //if no visualizations have been checked, alert the user
        if(checkedcnt == 0)
        {
            document.getElementById('opswrapper').className += " alert";
            return;
        }
      
        //collapse the info accordion
        if($("#accordion").accordion('option','active') == 0)
            $("#accordion").accordion('option','active',false);
      
        //remove the blank home tab
        var len = $("#tabs").tabs("length");
        for(var i = 0; i < len; i++)
            $("#tabs").tabs("remove",0);
        
        var loadinghtml = '<div><img id=\"loading\" class=\"loading\" src=\"./img/loading.gif\"></div>'
        addTab('Loading','loading',loadinghtml);
        $.ajax({ url: 'lib.psp',
                data: { fn: 'subsample', 
                        study: document.getElementById('studycombobox')[document.getElementById('studycombobox').selectedIndex].value,
                        subjects: $("#subjectslider").slider("value"),
                        samples:  $("#sampleslider").slider("value"),
                        sequences: $("#sequenceslider").slider("value"),
                        iterations: $("#iterationslider").slider("value"),
                        demo: document.visualizations[0].checked,
                        pcoa: document.visualizations[1].checked,
                        alpha_stddev: document.visualizations[2].checked,
                        alpha_stderr: document.visualizations[3].checked
                      },
                success: process_optimize});
    });
});


//function to build the frame with relevant plot
function process_optimize() {
      var checkedcnt = 0;
      
      //figure out how many frames need to be built
      for(i=0; i< document.visualizations.length; i++){
          if(document.visualizations[i].checked)
            checkedcnt += 1;
       }
        
        //remove the loading tab
        var len = $("#tabs").tabs("length");
        for(var i = 0; i < len; i++)
            $("#tabs").tabs("remove",0);
        
        var studyname = document.getElementById('studycombobox')[document.getElementById('studycombobox').selectedIndex].value;
        //content holds the frame html
        var content = "";
        for(i=0; i< document.visualizations.length; i++){
            content = "";
            //if demo is checked don't build a frame for it
            if(document.visualizations[0] == document.visualizations[i])
                continue;
            
            //if current vis box is checked then build a frame for it
            if(document.visualizations[i].checked)
            {
                content += "<div><img id=\"loading\" class=\"loading\" src=\"./img/loading.gif\"></div>";
                content += "<iframe src=\"handler.psp?";
                content += "viz="+document.visualizations[i].name;
                content += "&iterations="+$("#iterationslider").slider("value");
                content += "&column=";
                for(var j = 0; j < validColumns.length; j++)
                   content += validColumns[j]+',';
                content = content.slice(0,-1)
                content += "\" onLoad=\"loaded()\" onabort=\"loaded()\"></iframe>";
                //add a tab with this frame in it 
                addTab(document.visualizations[i].name,document.visualizations[i].value,content);
            }
        }
        //jquery to activate new tabs
        $(function() {
            // $("#tabs").tabs();
         $(".tabs-bottom .ui-tabs-nav, .tabs-bottom .ui-tabs-nav > *")
             .removeClass("ui-corner-all ui-corner-top")
             .addClass("ui-corner-bottom viz");
        });
        document.getElementById('optimize').disabled = true;  
    }


//jquery to activate sample slider
$(function() {
		$("#sampleslider").slider({
			range: "max",
			min: 1,
			max: 50,
			value: 50,
			slide: function(event, ui) {
			    document.getElementById('optimize').disabled = false;
                document.getElementById('samples').innerHTML = ui.value+"/"+$("#sampleslider").slider('option','max');
			},
			change: function(event, ui) {
			    document.getElementById('optimize').disabled = false;
                document.getElementById('samples').innerHTML = ui.value+"/"+$("#sampleslider").slider('option','max');
			}
		});
		$("#sampleslider").slider("disable");
		document.getElementById('samples').innerHTML = $("#sampleslider").slider("value")+"/"+$("#sampleslider").slider('option','max');
});

//function called to update sliders once the number of sequences has been changed
function seqsChanged(ui) {
    document.getElementById('optimize').disabled = false;
	document.getElementById('sequences').innerHTML = ui.value+"/"+$("#sequenceslider").slider('option','max');
	
	var numseqs = ui.value;
	var seqsindex = 0;
	
	//figure out the highest index of seqs per sample
	//that the current number falls under
    for(var i = 0; i < slidervals.length; i++)
    {
        if(numseqs < slidervals[i][0])
        {
            seqsindex = Math.max(i-1,0); //make sure not to get -1
            break;
        }
        if(numseqs == slidervals[i][0])
        {
            seqsindex = i;
            break;
        }
    }
	
	//max number of subjects and samples corresponds to number of
	//seqs per sample
    var maxsubs = slidervals[seqsindex][1];
    var maxsamps = slidervals[seqsindex][2];
    
    //annoying hack to make the slider look right when the max
    //is 1
    if(maxsubs == 1)
    {
        $("#subjectslider").slider("disable");
        $("#subjectslider").slider('option','min', 0);
        $("#subjectslider").slider('option','max', 2);
        $("#subjectslider").slider('option','value', 2);
    }
    else
    {
        $("#subjectslider").slider("enable");
        $("#subjectslider").slider('option','min', 1);
    }
    
    $("#subjectslider").slider('option','max', maxsubs);
    $("#subjectslider").slider('option','value', maxsubs);
    
    //annoying hack to make the slider look right when the max
    //is 1
    if(maxsamps == 1)
    {
        $("#sampleslider").slider("disable");
        $("#sampleslider").slider('option','min', 0);
        $("#sampleslider").slider('option','max', 2);
        $("#sampleslider").slider('option','value', 2);
    }
    else
    {
        $("#sampleslider").slider("enable");
        $("#sampleslider").slider('option','min', 1);
    }
            
    $("#sampleslider").slider('option','max', maxsamps);
    $("#sampleslider").slider('option','value', maxsamps);
    
    //need to figure out the metadata columns for the current
    //number of seqs per sample
    columnsIndex = seqsindex;
    for(columnsIndex = seqsindex; columnsIndex >= 0; columnsIndex--)
    {
        if(slidervals[columnsIndex][3] != 'None')
            break;
    }
    
    validColumns = slidervals[columnsIndex][3].split(',').sort();
    //build a list of valid columns
    var columnsHTML = "<ul>";
    for(var i = 0; i < validColumns.length; i++)
        columnsHTML += "<li>"+validColumns[i]+"</li>";
    columnsHTML += "</ul>";
    document.getElementById("columns").innerHTML = columnsHTML;
}

//jquery to activate the sequences slider
$(function() {
		$("#sequenceslider").slider({
			range: "max",
			min: 1,
			max: 10000,
			value: 10000,
			slide: function(event, ui) {
			    seqsChanged(ui);
			},
			change: function(event, ui) {
			    seqsChanged(ui);
			}
		});
		$("#sequenceslider").slider("disable");
		document.getElementById('sequences').innerHTML = $("#sequenceslider").slider("value") +"/"+$("#sequenceslider").slider('option','max');
});

//jquery to activate the subjects slider
$(function() {
		$("#subjectslider").slider({
			range: "max",
			min: 1,
			max: 100,
			value: 100,
			slide: function(event, ui) {
			  document.getElementById('optimize').disabled = false;
			  document.getElementById('subjects').innerHTML = ui.value +"/"+$("#subjectslider").slider('option','max');  
			},
			change: function(event, ui) {
			    document.getElementById('optimize').disabled = false;
			  document.getElementById('subjects').innerHTML = ui.value +"/"+$("#subjectslider").slider('option','max');  
			}
		});
		$("#subjectslider").slider("disable");
		document.getElementById('subjects').innerHTML = $("#subjectslider").slider("value") +"/"+$("#subjectslider").slider('option','max');
});

//jquery to activate the iterations slider
$(function() {
		$("#iterationslider").slider({
			range: "max",
			min: 1,
			max: 10,
			value: 1,
			slide: function(event, ui) {
			    document.getElementById('optimize').disabled = false;
			    document.getElementById('iterations').innerHTML = ui.value +"/"+$("#iterationslider").slider('option','max');
			},
			change: function(event, ui) {
			    document.getElementById('optimize').disabled = false;
			    document.getElementById('iterations').innerHTML = ui.value +"/"+$("#iterationslider").slider('option','max');
			}
		});
		$("#iterationslider").slider("disable");
		document.getElementById('iterations').innerHTML = $("#iterationslider").slider("value") +"/"+$("#iterationslider").slider('option','max');
});
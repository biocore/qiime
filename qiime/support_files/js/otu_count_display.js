/*
#script type="javscript/text"
#file otu_count_display.js

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Rob Knight"
__email__ = "jesse.stombaugh@colorado.edu"

Author: Jesse Stombaugh (jesse.stombaugh@colorado.edu)
Status: Prototype

Tested on the following web browsers:
Mozilla Firefox 3.5.3
Safari 4.0.3

*/


var row_max=new Array();
var row_sums=new Array();
var col_max=new Array();
var col_sums=new Array();

var non_filtered_otu_table=OTU_table; //this stores the original table

function filter_by_otu_count(non_filtered_otu_table,otu_cutoff){
    /*Filter the otu table by number of otus per sample*/
	num_cols=non_filtered_otu_table[0].length;
	num_rows=non_filtered_otu_table.length-1;
	col_sums=get_col_sums(non_filtered_otu_table,num_cols,num_rows);
	
	//If no cutoff is given, then the default is 5.
	if (!otu_cutoff.match(/^\d+$/)){
		otu_cutoff=5;
	}
	
	var filtered_OTU_table=new Array();
	filtered_OTU_table[0]=new Array();
	filtered_OTU_table[0][0]=non_filtered_otu_table[0][0];

	for (var i=1;i<non_filtered_otu_table.length;i++){
		filtered_OTU_table[i]=new Array();
		filtered_OTU_table[i][0]=non_filtered_otu_table[i][0]
		var iterator=1;
		for (var j=1;j<non_filtered_otu_table[0].length;j++){
			if (col_sums[j]>otu_cutoff){
				filtered_OTU_table[0][iterator]=non_filtered_otu_table[0][j];
				filtered_OTU_table[i][iterator]=non_filtered_otu_table[i][j];
				iterator++;
			}
		}
	}
	return filtered_OTU_table;
}

var taxonomy=0; /*This is a global variable, so the script knows which table is 
                  being displayed*/

function create_OTU_intervals(){
    /*This function calls the script that generates the html heatmap based on
      SampleID's.*/
    
    otu_cutoff=document.getElementById("otu_count_cutoff").value;
    
    //First, we clear the table currently in the html.
	for (var i=document.all.otu_table_body.rows.length; i>=1; i--){
	         document.all.otu_table_body.deleteRow(i-1);
	}
	for (var i=document.all.otu_table_head.rows.length; i>=1; i--){
	         document.all.otu_table_head.deleteRow(i-1);
	}
    //Second, we filter by the otus per sample.
	OTU_table=filter_by_otu_count(non_filtered_otu_table,otu_cutoff);

	num_cols=OTU_table[0].length;
	num_rows=OTU_table.length-1;
	col_sums=get_col_sums(OTU_table,num_cols,num_rows);
	row_max=get_row_max(OTU_table,num_cols,num_rows);
	row_sums=get_row_sums(OTU_table,num_cols,num_rows);
	col_max=get_col_max(OTU_table,num_cols,num_rows);
	
	taxonomy=0;
	//This loop iterates through the table and generates each row.  Need to set
	//a timeout, or else the browser will stall when loading large tables.
	create_OTU_table(0,OTU_table);
	for (var i=1; i<(OTU_table.length-1); i++){
		setTimeout("create_OTU_table("+i+",OTU_table)",0);
	}
}

function write_taxon_heatmap(){
    /*This function calls the script that generates the html heatmap based on
      Taxons.*/
      
    otu_cutoff=document.getElementById("otu_count_cutoff").value;
    //First, we clear the table currently in the html.
	for (var i=document.all.otu_table_body.rows.length; i>=1; i--){
	         document.all.otu_table_body.deleteRow(i-1);
	}
	for (var i=document.all.otu_table_head.rows.length; i>=1; i--){
	         document.all.otu_table_head.deleteRow(i-1);
	}
    //Second, we filter by the otus per sample.
	OTU_table=filter_by_otu_count(non_filtered_otu_table,otu_cutoff);
	
	//Third, we need to transpose the original table and then sort by taxons.
	OTU_table=transpose_otu_table(OTU_table,otu_cutoff);
	//OTU_table=sort_by_taxonomy(OTU_table);
	
	num_cols=OTU_table[0].length-1;
	num_rows=OTU_table.length;

	col_sums=get_col_sums(OTU_table,num_cols,num_rows);
	row_max=get_row_max(OTU_table,num_cols,num_rows);
	row_sums=get_row_sums(OTU_table,num_cols,num_rows);
	col_max=get_col_max(OTU_table,num_cols,num_rows);
	taxonomy=1;
    
    //This loop iterates through the table and generates each row.  Need to set
	//a timeout, or else the browser will stall when loading large tables.
	create_OTU_table(0,OTU_table);
	for (var i=1; i<(OTU_table.length); i++){
		setTimeout("create_OTU_table("+i+",OTU_table)",0);
	}
}

function create_OTU_table(row_id,OTU_table){
    /*This function adds table rows to the heatmap using the DOM.*/
    
    //This is the rainbow set of colors used in the heatmap
	bkgrd_colors=new Array();
	bkgrd_colors[0]="#000080"
	bkgrd_colors[1]="#0000ff"
	bkgrd_colors[2]="#0063ff"
	bkgrd_colors[3]="#00d5ff"
	bkgrd_colors[4]="#4effa9"
	bkgrd_colors[5]="#a9ff4e"
	bkgrd_colors[6]="#ffe600"
	bkgrd_colors[7]="#ff7d00"
	bkgrd_colors[8]="#ff1400"
	bkgrd_colors[9]="#800000"
	
	var otu_table_html=document.getElementById("otu_table_html");
	var otu_table_body=document.getElementById("otu_table_body");
	var otu_table_head=document.getElementById("otu_table_head");
	var body = document.getElementsByTagName("body")[0];
	//When writing the rows, the header row is generated without color or 
	//mouseovers
	if (row_id==0){
		var row = document.createElement("tr");
		
		for (var i = 0; i < OTU_table[row_id].length; i++) {
		    var cell = document.createElement("th");
		    if (i==0 && taxonomy==1){
		        cell.setAttribute("class","lineage");
		    }else{
		        cell.setAttribute("class","rotate");
	        }

			var cellText = document.createTextNode(OTU_table[row_id][i]);
			cell.appendChild(cellText);
			row.appendChild(cell);
		}
		otu_table_head.appendChild(row);
		otu_table_html.appendChild(otu_table_head);
	}else{
		var row = document.createElement("tr");
		row_len=OTU_table[row_id].length;
       	for (var i = 0; i < row_len; i++) {
		   	var cell = document.createElement("td");
		   	//This handles the row headers separately from other cells
			if (i==0){
				var cellText= document.createTextNode(OTU_table[row_id][i]);
				cell.setAttribute("class","dragHandle");
				cell.appendChild(cellText);
			}else{
			    if (taxonomy==1 && i==(row_len-1)){
			        var cellText= document.createTextNode(OTU_table[row_id][i]);
 			        cell.appendChild(cellText);
 			        cell.setAttribute("bgcolor",'white');
			    }else if(OTU_table[row_id][i]==0){
			        var cellText= document.createTextNode('');
			        cell.appendChild(cellText);
			        cell.setAttribute("bgcolor",bkgrd_colors[0]);
			    }else{
    				var anchor=document.createElement("a");
    				var cellText= document.createTextNode(OTU_table[row_id][i]);
    				if(taxonomy==1){
					    var lineage=0;
					    normalized_value=(OTU_table[row_id][i]/col_max[i]);
					    anchor.setAttribute("onmouseover",create_mouseover(
					                        OTU_table,row_id,i,col_sums[i],
					                        row_sums[row_id],lineage,taxonomy));
    				}else{
    					normalized_value=(OTU_table[row_id][i]/row_max[row_id]);
    					var lineage=OTU_table.length-1;
    					anchor.setAttribute("onmouseover",create_mouseover(
    					                    OTU_table,row_id,i,row_sums[row_id],
    					                    col_sums[i],lineage,taxonomy));
    				}
				    anchor.setAttribute("onmouseout",'return nd();');
    				anchor.setAttribute("onclick","javascript:write_subset("+
    				                              "event,this,OTU_table);");
    				anchor.appendChild(cellText);
    				cell.appendChild(anchor);
				    var color_bin=0;
    				for (var c=0; c<bkgrd_colors.length;c++){
    					min=(c/10);
    					max=((c+1)/10);
    					if (min <= normalized_value && normalized_value <= max){
    						var color_bin=c;
    						break;
    					}
    				}
    				cell.setAttribute("bgcolor",bkgrd_colors[color_bin]);
                }
			}
			row.appendChild(cell);	
	
		}
		otu_table_body.appendChild(row);
        $('#otu_table_body').tableDnDUpdate();
		otu_table_html.appendChild(otu_table_body);
	}
	body.appendChild(otu_table_html);
	
}	

function create_mouseover(OTU_table,row_id,col_id,row_sum,col_sum,lineage,taxonomy){
    /*This creates the mouseover functionality of the cells*/
    var table_row_len=OTU_table[row_id].length-1;
    if (taxonomy==1){
        var split_lineage=OTU_table[row_id][lineage].split(";");
        var OTU_label=OTU_table[row_id][table_row_len];	
        var Sample_label=OTU_table[0][col_id];
    }else{
	    var split_lineage=OTU_table[lineage][col_id].split(";");
	    var OTU_label=OTU_table[0][col_id];	
	    var Sample_label=OTU_table[row_id][0];
    }
	mouseover_str="return overlib(\"";
	mouseover_str+="<p><b>OTU: ";
	mouseover_str+=OTU_label;
	mouseover_str+="</b><br>";
	mouseover_str+=OTU_table[row_id][col_id]
	mouseover_str+="/";
	mouseover_str+=col_sum;
	mouseover_str+=" (";
	mouseover_str+=(OTU_table[row_id][col_id]/col_sum*100).toFixed(2);
	mouseover_str+="%) Sequences";
	mouseover_str+="<br><br><b>SampleID: ";
	mouseover_str+=Sample_label;
	mouseover_str+="</b><br>";
	mouseover_str+=OTU_table[row_id][col_id];
	mouseover_str+="/";
	mouseover_str+=row_sum;
	mouseover_str+=" (";
	mouseover_str+=(OTU_table[row_id][col_id]/row_sum*100).toFixed(2);
	mouseover_str+="%) Displayed";
	mouseover_str+="<br><br><b>Lineage: ";
    for (var i=0;i<split_lineage.length-2;i++){
		mouseover_str+="<br>"+split_lineage[i];
	}
	var query_string=split_lineage[i].split(' ').join('');
	small_link='http://www.google.com/search?q='+query_string+'';

	open_window="javascript:window.open("+small_link+",'win2');";
	mouseover_str+="<br><a href="+"'"+small_link+"' target='_blank'"+">";
	mouseover_str+=split_lineage[i];
	mouseover_str+="</a>";
	
	mouseover_str+="</b></p>";
	mouseover_str+="\",STICKY,MOUSEOFF,RIGHT);"

	return mouseover_str;
}

var split_lineage;

function get_row_max(OTU_table,num_cols,num_rows){
    /*This gets the max value with the row, which is used for normalization*/
	row_max = new Array();
	row_max[0]=0;
	for (var b=1;b<num_rows;b++){
		row_max[b]=0;
		for (var i=1;i<num_cols;i++){
			if (OTU_table[b][i] > row_max[b]){
				row_max[b] = OTU_table[b][i];
			}
		}
	}
	return row_max;
};

function get_row_sums(OTU_table,num_cols,num_rows){
    /*This sums each row and is used for calculations*/
	row_sums = new Array();
	row_sums[0]=0;
	for (var b=1;b<num_rows;b++){
		row_sums[b]=0;
		for (var i=1;i<num_cols;i++){
			row_sums[b]+=OTU_table[b][i];
		}
	}
	return row_sums;
};

function get_col_max(OTU_table,num_cols,num_rows){
    /*This gets the max value for a each column*/
	col_max = new Array();
	col_max[0]=0;

	for (var i=1;i<num_cols;i++){
		col_max[i]=0;
		for (var b=1;b<num_rows;b++){
			if (OTU_table[b][i] > col_max[i]){
				col_max[i] = OTU_table[b][i];
			}
		}
	}
	return col_max;
};

function get_col_sums(OTU_table,num_cols,num_rows){
    /*This sums each column and is used for calculations*/
	col_sums = new Array();
	col_sums[0]=0;
	for (var i=1;i<num_cols;i++){
		col_sums[i]=0;
		for (var b=1;b<num_rows;b++){
			col_sums[i]+=OTU_table[b][i];
		}
	}
	return col_sums;
};

function write_subset(event,cell){
	/*This function creates the google heatmap popup window*/
	otu_count_dom=document.getElementById("otu_table_html");
	
	//These for loops extract the heatmap from the DOM.  The reason for this is
	//because the user can drag and drop rows, which is only updated in the DOM
	//and not in the javascript array.
	var otu_table_dom=new Array();	
	otu_table_dom[0]=new Array();
	for (var i=0;i<otu_count_dom.rows[0].cells.length;i++){
	    otu_table_dom[0][i]=otu_count_dom.rows[0].cells[i].innerHTML;
	}
	
	for (var j=1;j<otu_count_dom.rows.length;j++){
	    otu_table_dom[j]=new Array();
	    otu_table_dom[j][0]=otu_count_dom.rows[j].cells[0].innerHTML;
	}

	for (var i=1;i<otu_count_dom.rows.length;i++){
	    for (var j=1;j<otu_count_dom.rows[i].cells.length;j++){    
	        var otu_cell=otu_count_dom.rows[i].cells[j].innerHTML;
	        pattern=/>\d+</;
	        var parsestring=otu_cell.split(' ').join('').match(pattern);
            if( parsestring!=null ){
                str=parsestring.toString();
                str=str.replace(/>/,'');
                otu_table_dom[i][j]=str.replace(/</,'');
            }else{
                otu_table_dom[i][j]=0;
	        }
	    }
	}
	//In the heatmap based on taxonomy, this loop gets the last column which 
	//contains the OTU ID.
	if (taxonomy==1){
	    for (var j=0;j<otu_count_dom.rows.length;j++){
	        var rlen=(otu_count_dom.rows[j].cells.length-1);
	        otu_table_dom[j][rlen]=otu_count_dom.rows[j].cells[rlen].innerHTML;
	    }
	}
	OTU_table=otu_table_dom;

	cell_index=cell.parentNode.cellIndex;
	row_index=cell.parentNode.parentNode.rowIndex;
    //Depending on which cell is clicked, these statements take 20 cells in each
    //direction
	if (cell_index>21 && cell_index<(OTU_table[0].length-21)){
		x_start=cell_index-20;
		x_stop=cell_index+20;
	}else if(cell_index<21 && cell_index<(OTU_table[0].length-21)){
		x_start=1;
		x_stop=cell_index+20;
	}else if(cell_index>21 && cell_index>(OTU_table[0].length-21)){
		x_start=cell_index-20;
		if (taxonomy==1){
		    x_stop=OTU_table[0].length-2;
		}else{
		    x_stop=OTU_table[0].length-1;
		}
	}else{
		x_start=1;
		if (taxonomy==1){
		    x_stop=OTU_table[0].length-2;
		}else{
		    x_stop=OTU_table[0].length-1;
		}
	}
	if (row_index>21 && row_index<(OTU_table.length-22)){
		y_start=row_index-20;
		y_stop=row_index+20;
	}else if(row_index<21 && row_index<(OTU_table.length-22)){
		y_start=1;
		y_stop=row_index+20;
	}else if(row_index>21 && row_index>(OTU_table.length-22)){
		y_start=row_index-20;
		if (taxonomy==1){
		    y_stop=OTU_table.length-1;
		}else{
		    y_stop=OTU_table.length-1;
		}
	}else{
		y_start=1;
		if (taxonomy==1){
		    y_stop=OTU_table.length-1;
		}else{
		    y_stop=OTU_table.length-1;
		}
	}
    //This generates the data model, which the google heatmap uses to create
    //the heatmap.
	var setTable='';
	var iterator1=1;
	for (var j=y_start; j<=y_stop; j++) {
		var iterator2=1;
        for (var i=x_start; i<=x_stop; i++) {
            var row_length=OTU_table[j].length-1;
            if (taxonomy==1){
                setTable+="tableModel.setContentAt("+iterator1+",0,'OTU: "+
                                                    OTU_table[j][row_length]+
                                                    "');\n";
                setTable+="tableModel.setContentAt(0,"+iterator2+",'SampleID:"+
                                                    OTU_table[0][i]+"');\n";
            }else{
                setTable+="tableModel.setContentAt("+iterator1+",0,'SampleID: "+
                                                    OTU_table[j][0]+"');\n";
                setTable+="tableModel.setContentAt(0,"+iterator2+",'OTU:"+
                                                    OTU_table[0][i]+"');\n";			    	
			}
			setTable+="tableModel.setContentAt("+iterator1+", "+iterator2+", "+
			                                        OTU_table[j][i]+");\n";
			iterator2++;
		}
		iterator1++;
	}
	
	//Open a window and write the source code for that window.
	var win = window.open("", "win", "width=525,height=600"); // a window object

	win.document.writeln("<HTML><HEAD><TITLE>OTU Heatmap - Google API</TITLE><scr"+"ipt type='text/javascript' src='http://magic-table.googlecode.com/svn/trunk/magic-table/javascript/magic_table.js'></scr"+"ipt>"
	+"<\/HEAD><BODY onload='setTimeout(\"drawTable()\",5)';'><div id='tableTargetDiv'></div><scr"+"ipt type='text/javascript'>"
	+"func"+"tion drawTable(){var count=0;var targetElement = document.getElementById('tableTargetDiv'); var defaultRowHeight = 25; var defaultColumnWidth = 70; var tablePositionX = 0; var tablePositionY = 0; var tableHeight = 500; var tableWidth = 500; var rowHeaderCount = 1; var columnHeaderCount = 1; var rows = 41 + columnHeaderCount; var columns = 41 + rowHeaderCount;var tableModel = new greg.ross.visualisation.TableModel(rows, columns, defaultRowHeight, defaultColumnWidth, rowHeaderCount, columnHeaderCount);"
	+setTable
	+"tableModel.recalculateMinMaxValues();var fisheyeTable = new greg.ross.visualisation.FisheyeTable(tableModel, tablePositionX, tablePositionY, tableWidth, tableHeight, 'OTU Counts', targetElement);}</scr"+"ipt>"
	+"</BODY></HTML>");
	
	win.document.close();
}

function transpose_otu_table(tran_OTU_table,otu_cutoff){
    //This transposes the original otu count table.
	var table_len=tran_OTU_table.length-1;
    
    for (var j=0;j<tran_OTU_table[0].length;j++){
        var OTU_head=tran_OTU_table[0][j];
        tran_OTU_table[0][j]=tran_OTU_table[table_len][j];
        tran_OTU_table[table_len][j]=OTU_head;
    }
	
	var transpose_OTU_table=new Array();
	for (var i=0;i<tran_OTU_table[0].length;i++){
		transpose_OTU_table[i]=new Array();
		for (var j=0;j<tran_OTU_table.length;j++){
			transpose_OTU_table[i][j]=tran_OTU_table[j][i];
		}
	}
    return transpose_OTU_table;
}

function sort_by_taxonomy(tran_otu_table){
    //This sorts the taxonomy heatmap by the taxons.
    OTU_table_reduced=new Array();
	for (var i=1;i<tran_otu_table.length;i++){
		OTU_table_reduced[i-1]=tran_otu_table[i]
	}
	
	OTU_table_reduced=OTU_table_reduced.sort(comp_row_headers);
	for (var i=1;i<OTU_table_reduced.length;i++){
		tran_otu_table[i+1]=OTU_table_reduced[i];
	}
	return tran_otu_table;
}

function comp_row_headers(a, b) {
    //This is used to compare taxons and determines the order for sorting.
   // psudeo code.
	a=a[0];
    b=b[0];
    if (a < b) {
        return -1;
    }
    if (a > b) {
        return 1;
    }
    if (a == b) {
        return 0;
    }
}

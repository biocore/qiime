// DATA_TEMPLATE: js_data
oTest.fnStart( "fnServerData" );

/* Not testing anything beyond the default here. SSP will do that... */

$(document).ready( function () {
	/* Check the default */
	var oTable = $('#example').dataTable( {
		"aaData": gaaData
	} );
	var oSettings = oTable.fnSettings();
	var mPass;
	
	oTest.fnTest( 
		"Default should be $.getJSON",
		null,
		function () { return oSettings.fnServerData == $.getJSON; }
	);
	
	
	oTest.fnComplete();
} );
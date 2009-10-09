// DATA_TEMPLATE: empty_table
oTest.fnStart( "fnServerData" );

/* Not testing anything beyond the default here. SSP will do that... */

$(document).ready( function () {
	/* Check the default */
	var oTable = $('#example').dataTable( {
		"bServerSide": true,
		"sAjaxSource": "../../../examples/examples_support/server_processing.php"
	} );
	var oSettings = oTable.fnSettings();
	var mPass;
	
	oTest.fnWaitTest( 
		"Default should be $.getJSON",
		null,
		function () { return oSettings.fnServerData == $.getJSON; }
	);
	
	
	oTest.fnComplete();
} );
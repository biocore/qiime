/*
*  How to load a feed via the Feeds API.
*/
google.load("feeds", "1");
// var container = document.getElementById('news');
// Our callback function, for when a feed is loaded.
function newsFeedLoaded(result) {
    var container = document.getElementById("news");
  if (!result.error) {
	formatFeed(result, container);
  }
  else
  {
      container.innerHTML = "<p>Error loading feed.</p>";
  }
}

function citationFeedLoaded(result) {
    var container = document.getElementById("citations");
  if (!result.error) {
	formatFeed(result, container);
  }
  else
  {
      container.innerHTML = "<p>["+result.error.code+"] "+result.error.message+"</p>";
  }
}

function forumFeedLoaded(result) {
    var container = document.getElementById("forum");
  if (!result.error) {
	formatFeed(result, container);
  }
  else
  {
      container.innerHTML = "<p>Error loading feed.</p>";
  }
}

function formatFeed(result, container) {
	var html_string = '<ul class="feedlist">'
    for (var i = 0; i < result.feed.entries.length; i++) {
      var entry = result.feed.entries[i];
		html_string += '<li><a href=\"'+entry.link+'\" target="_blank">'+entry.title+'</a><p>'+entry.publishedDate+'<br>'+entry.contentSnippet+'</p></li>';
      // var div = document.createElement("div");
      // div.appendChild(document.createTextNode(entry.title));
      // container.appendChild(div);
    }
	html_string += '</ul>'
	container.innerHTML = html_string;
}

function load() {
	// var container = document.getElementById("news");
	var news_feed = new google.feeds.Feed("http://qiime.wordpress.com/feed/");

  // Calling load sends the request off.  It requires a callback function.
  news_feed.load(newsFeedLoaded);

  // container = document.getElementById("forum");
  var forums_feed = new google.feeds.Feed("https://groups.google.com/forum/feed/qiime-forum/topics/atom.xml?num=3");

  // Calling load sends the request off.  It requires a callback function.
  forums_feed.load(forumFeedLoaded);

  var citation_feed = new google.feeds.Feed("http://rss.webofknowledge.com/rss?e=9caba57f6a84aac7&c=7af8db36f4f5d8b8b377625e44613eda");

  // Calling load sends the request off.  It requires a callback function.
  citation_feed.load(citationFeedLoaded);
}

function loadCitations(selected_item)
{
    menuChange(selected_item);

    var citation_feed = new google.feeds.Feed("http://rss.webofknowledge.com/rss?e=9caba57f6a84aac7&c=7af8db36f4f5d8b8b377625e44613eda");
  citation_feed.setNumEntries(2000);
  citation_feed.includeHistoricalEntries();
  // Calling load sends the request off.  It requires a callback function.
  citation_feed.load(citationFeedLoaded);
}

// google.setOnLoadCallback(load);

function loadHeaderFooter(selected_page)
{
    if(selected_page == 'main')
        load();
    document.getElementById('header').innerHTML = getFile("/home_static/header.html");
    document.getElementById('leftcol').innerHTML = getFile("/home_static/leftcol.html");
    document.getElementById('footer').innerHTML = getFile("/home_static/footer.html");
}


// function menuChange(selected_item)
// {
//   // document.getElementById('main').className = "unselected";
//   // document.getElementById('dataFiles').className = "unselected";
//   // document.getElementById('calendar').className = "unselected";
//   // document.getElementById('citations').className = "unselected";
//
//   // for(var element in document.getElementsByName(selected_item))
//     // element.className = "selected";
//
//   document.getElementById('content').innerHTML = getFile(selected_item+".html");
//
//   if(selected_item == 'main')
//     load();
// }

function getFile(fileName){
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

function search() {
    var targets = document.forms["searchForm"]["target"];
    for(var i = 0; i<targets.length; i++) {
        if(targets[i].checked) {
           var target = targets[i].value;
        }
    }
    var query = document.forms["searchForm"]["query"].value;

    window.location = target+query;
    return false;
}

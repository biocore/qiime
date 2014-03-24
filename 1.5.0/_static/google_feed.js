google.load("feeds", "1");

function initialize() {
  var feed = new google.feeds.Feed("http://qiime.wordpress.com/feed/");
  feed.load(function(result) {
    if (!result.error) {
      var container = document.getElementById("feed");
      for (var i = 0; i < 3; i++) {
        var entry = result.feed.entries[i];
        var td =document.createElement('td')
        var link = document.createElement('a');
        link.setAttribute('href', entry.link);
        var dot = document.createElement('b');
        dot.setAttribute('style', 'color: red; font-size: 10pt');
        var dottext=document.createTextNode('• ')
        dot.appendChild(dottext)
        link.appendChild(dot)
        var title=document.createTextNode(entry.title)
        link.appendChild(title)
        td.appendChild(link)
        container.appendChild(td);
      }
    }
  });
}
google.setOnLoadCallback(initialize);
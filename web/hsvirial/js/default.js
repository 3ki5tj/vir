/* Copyright (c) Cheng Zhang 2010-2015 */
window.onload = init;



function init()
{
  //appendspaces(); // at spaces at the bottom of the page
  toc_init(); // initialize a toc
  dtbox_init(); // bar for data and time
  makeimgslinks(); // add links to the page
  makejava(); // substitute Java applet
}



// create a table-of-content sidebar
function toc_init()
{
  var nodes = document.getElementsByTagName("h2");
  var cnodes = document.getElementsByTagName("h3");

  var s = "";
  for (var i = 0; i < nodes.length; i++) {
    var cap = nodes[i].innerHTML;
    nodes[i].innerHTML = '<a name="h2_' + i + '"></a>' + cap; // add a marker
    s += '<li><a href="#h2_' + i + '">' + cap + '</a></li>'; // add to TOC
    var cs = "";
    // find h3 elements listed after this h2
    for (var j = 0; j < cnodes.length; j++) {
      var ccap = cnodes[j].innerHTML;
      for (var k = 0; k < nodes.length; k++) {
        if (cnodes[j].compareDocumentPosition(nodes[k]) == 4) {
          break; // cnode[j] is after nodes[k]
        }
      }
      // the parent <h2> is this one
      if ( k - 1 != i ) {
        continue;
      }
      console.log(j + " " + ccap + " " + nodes[i].innerHTML);
      cnodes[j].innerHTML = '<a name="h3_' + j + '"></a>' + ccap;
      cs += '<li><a href="#h3_' + j + '">' + ccap + '</a></li>';
    }
    if (cs !== "") {
      s += "<ul>" + cs + "</ul>";
    }
  }
  if (s === "") {
    return;
  } else {
    s = "<ul>" + s + "</ul>";
  }

  var env = document.createElement("div"); // envelope
  env.setAttribute("id", "tocbox");

  var title = document.createElement("div");
  title.setAttribute("id", "toctitle");
  title.innerHTML = document.getElementsByTagName('title')[0].innerHTML;
  env.appendChild(title);

  var list = document.createElement("div");
  list.setAttribute("id", "toclist");
  list.innerHTML = s;
  env.appendChild(list);

  var links = document.createElement("div");
  links.innerHTML = '<hr>'
    + '<div class="toclinkhead">Links:</div>'
    + '<ul class="toclinks">'
    + '<li><a href="/">Home</a>'
    + '</ul>';
  env.appendChild(links);

  env.onmouseover = function () { toc_show(1); };
  env.onmouseout = function () { toc_show(0); };

  nodes[0].parentNode.insertBefore(env, nodes[0]);
  toc_show(0);
}



// show or hide the table-of-contents sidebar
function toc_show(open)
{
  var box = document.getElementById("tocbox");
  var width = 200;
  box.style.width = "" + width  + "px";
  box.style.left =  (open ? "0px" : "" + (-width-40) + "px");
}



// the dynamic data/time bar at the bottom of the page
var dtbox;
function dtbox_init() // create a date/time box
{
  dtbox = document.createElement("div");
  dtbox.id = "dtbox";
  dtbox.setAttribute("style",
      "font: small-caps 100% sans-serif; position: fixed;" \
    + "bottom: 0; left: 0; right: 0; height: 20px; " \
    + "padding: 2px; text-align: center;");
  dtbox.onmouseover = function() { dtbox_update(true); };
  dtbox.onmouseout = function() { dtbox_update(false); };
  document.body.appendChild(dtbox);
  dtbox_update(false);
}



function dtbox_update(open)
{
  dtbox.style.backgroundColor = open ? "#F0F0F0" : "";
  dtbox.innerHTML = open ? new Date().toLocaleString() : "";
}


// remove the possible postfix of the file name `s'
function rmpostfix(s, pfx)
{
  var i = s.lastIndexOf(".");
  if (i == -1) {
    return s;
  }
  var m = pfx.length;
  if (i < m || s.substr(i - m, m) != pfx) {
    return s;
  }
  return s.substr(0, i - m) + s.substr(i, s.length - i);
}



// add links to all "demo" or "fig" class `img' elements
function makeimgslinks()
{
  var imgs = document.getElementsByTagName("img");

  for (var i = 0; i < imgs.length; i++) {
    if ( imgs[i].className == "demo" || imgs[i].className == "fig" ) {
      pr = imgs[i].parentNode;
      nb = document.createElement('a');
      // if the image name ends with 'sm', e.g. 'fig123sm.png'
      // then we remove the link
      imglink = rmpostfix(imgs[i].src, 'sm');
      nb.setAttribute('href', imglink);
      nb.setAttribute('title', imgs[i].getAttribute('alt'));
      nb.appendChild(imgs[i].cloneNode(true)); // destroy and create an img, so index is not messed up
      pr.replaceChild(nb, imgs[i]);
    }
  }
}



function appendspaces()
{
  document.body.innerHTML += "<p>&nbsp;</p><p>&nbsp;</p>";
}



function invokejava(img) {
  par = img.parentNode;
  p = document.createElement('div');
  p.innerHTML = '<applet code="' + img.alt + '" archive="' + img.alt + '.jar" ' \
      + 'width="' + img.width + '" height="' + img.height + '"></applet>';
  par.replaceChild(p, img);
}



function makejava() {
  imgs = document.getElementsByTagName("img");
  for (var i = 0; i < imgs.length; i++) {
    if (imgs[i].className == "javaimg") {
      imgs[i].setAttribute("onclick", "invokejava(this)");
    }
  }
}



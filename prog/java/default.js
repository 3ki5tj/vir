/* Copyright (c) Cheng Zhang 2010-2012 */
window.onload = init;

function init()
{
  appendspaces(); // at spaces at the bottom of the page
  toc_init(); // initialize a toc
  dtbox_init(); // bar for data and time
  makeimgslinks(); // add links to the page
  makejava(); // substitute Java applet
}

function toc_init() // create a table of content
{
  var nodes = document.getElementsByTagName("h2");
  for (s = "", i = 0; i < nodes.length; i++) {
    cap = nodes[i].innerHTML;
    nodes[i].innerHTML = '<a name="h2_' + i + '"></a>' + cap; // add a marker
    s += '<li><a href="#h2_' + i + '">' + cap + '</a></li>'; // add to TOC
  }
  if (s == "") return; else s = "<ul>" + s + "</ul>";

  var env = document.createElement("div");
  env.setAttribute("id", "tocbox");

  var title = document.createElement("div");
  title.setAttribute("id", "toctitle");
  title.innerHTML = document.getElementsByTagName('title')[0].innerHTML;
  env.appendChild(title);

  var list = document.createElement("div");
  list.setAttribute("id", "toclist");
  list.innerHTML = s;
  env.appendChild(list);

  env.onmouseover = function () { tocshow(1); };
  env.onmouseout = function () { tocshow(0); };

  nodes[0].parentNode.insertBefore(env, nodes[0]);
  tocshow(0);
}

function tocshow(open)
  { document.getElementById("tocbox").style.left = (open ? "0px" : "-200px"); }


var dtbox;
function dtbox_init() // create a date/time box
{
  dtbox = document.createElement("div");
  dtbox.id = "dtbox";
  dtbox.setAttribute("style", "font: small-caps 100% sans-serif; position: fixed;"
      + "bottom: 0; left: 0; right: 0; height: 20px; padding: 2px; text-align: center;")
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

function makeimgslinks() { // make all demo images links
  var imgs = document.getElementsByTagName("img");
  for (var i = 0; i < imgs.length; i++)
    if (imgs[i].className == "demo") {
      pr = imgs[i].parentNode;
      nb = document.createElement('a');
      nb.setAttribute('href', imgs[i].src);
      nb.setAttribute('title', imgs[i].getAttribute('alt'));
      nb.appendChild(imgs[i].cloneNode(true)); // destroy and create an img, so index is not messed up
      pr.replaceChild(nb, imgs[i]);
    }
}

function appendspaces() { document.body.innerHTML += "<p>&nbsp;<p>&nbsp;"; }

function invokejava(img) {
  par = img.parentNode;
  p = document.createElement('div');
  p.innerHTML = '<applet code="' + img.alt
      + '" width="' + img.width + '" height="' + img.height + '"></applet>';
  par.replaceChild(p, img);
}

function makejava() {
  imgs = document.getElementsByTagName("img");
  for (var i = 0; i < imgs.length; i++) {
    if (imgs[i].className == "javaimg")
      imgs[i].setAttribute("onclick", "invokejava(this)");
  }
}

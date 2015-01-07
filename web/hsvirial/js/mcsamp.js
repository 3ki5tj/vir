// Monte Carlo sampler

var radius = 16; // radius in pixels
var xy = [];
var amp = 16;
var canvas = document.getElementById("hsbox");
var ctx = canvas.getContext("2d");
var width = canvas.width;
var height = canvas.height;
var mcTot = 0, mcAcc = 0, step = 0;
var moveAtom = -1; // which atom is attempted
var moveAcc = false; // if the move is accepted


function dist(x, y) {
  x = Math.abs(x); if (x > width/2) x -= width;
  y = Math.abs(y); if (y > height/2) y -= height;
  return Math.sqrt(1.*x*x + 1.*y*y);
}

// draw a ball
function drawBall(ctx, x, y, r, color0, color1) {
  // the main blue ball
  var grd = ctx.createRadialGradient(x + r*.3, y - r*.4, r*.1, x, y, r);
  grd.addColorStop(0, color0);
  grd.addColorStop(1, color1);
  ctx.fillStyle = grd;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI);
  ctx.closePath();
  ctx.fill();
}

// initial positions
function initPos(n) {
  xy = [];
  for (var i = 0; i < 10000; i++) {
    var x = Math.random() * width;
    var y = Math.random() * height;
    for (var j = 0; j < xy.length; j++)
      if (dist(xy[j][0] - x, xy[j][1] - y) < 2*radius)
        break;
    if (j == xy.length) {
      xy.push([x, y]);
      if (xy.length >= n) break;
    }

  }
  return xy;
}

// draw all atoms in the box
function draw() {
  // draw the background
  ctx.fillStyle = "#f0f0f0";
  ctx.fillRect(0, 0, width, height);

  for (var i = 0; i < xy.length; i++) {
    var x = xy[i][0];
    var y = xy[i][1];
    var spotcolor = "#a0a0e0";
    var color = "#2040a0";
    if (i == moveAtom) {
      if (moveAcc) {
        spotcolor = "#a0e0a0";
        color = "#20c020";
      } else {
        spotcolor = "#e0a0a0";
        color = "#c02020";
      }
    }
    drawBall(ctx, x, y, radius, spotcolor, color);
  }
}

function mcStep() {
  var n = xy.length;
  moveAtom = Math.floor(Math.random() * n);
  var xi = xy[moveAtom][0] + (2*Math.random() - 1) * amp;
  if (xi < 0) xi += width; else if (xi > width) xi -= width;
  var yi = xy[moveAtom][1] + (2*Math.random() - 1) * amp;
  if (yi < 0) yi += height; else if (yi > height) yi -= height;

  step += 1;
  mcTot += 1;
  for (var j = 0; j < n; j++) {
    if (j == moveAtom) continue;
    var d = dist(xi - xy[j][0], yi - xy[j][1]);
    if (d < 2*radius) break;
  }
  if (j === n) {
    xy[moveAtom][0] = xi;
    xy[moveAtom][1] = yi;
    mcAcc += 1;
    moveAcc = true;
  } else {
    moveAcc = false;
  }
  accratio = 1.*mcAcc/mcTot*100;
  if (step % 100 === 0) {
    draw();
    document.getElementById("accratio").innerHTML
      = accratio.toPrecision(4) + "%";
  }
  return j === n;
}

initPos(55);
setInterval(mcStep, 10);



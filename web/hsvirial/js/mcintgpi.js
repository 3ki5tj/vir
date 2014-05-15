// Compute Pi from Monte Carlo integration

function mcgetPi() {
  var c = document.getElementById("mcpi");
  var ctx = c.getContext("2d");
  var w = c.width;
  var h = c.height;
  var r = Math.min(w/2, h/2) - 1;

  // draw the background
  ctx.fillStyle = "#f0f0f0";
  ctx.fillRect(0, 0, w, h);

  // fill the circle
  ctx.beginPath();
  ctx.arc(w/2, h/2, r, 0, Math.PI*2);
  ctx.closePath();
  ctx.fillStyle = "#c0c0c0";
  ctx.fill();

  // draw random dots
  ctx.fillStyle = "black";
  var s = document.getElementById("numSamps").value;
  if ( !isNaN(s) ) num = s; else num = 1;
  var maxNumDraw = w*h;
  var cnt = 0;
  document.getElementById("PiEstimate").innerHTML = "";
  document.getElementById("PiComment").innerHTML = num > maxNumDraw ? " (visualization disabled)" : "";

  for ( i = 0; i < num; i++ ) {
    var x = Math.random();
    var y = Math.random();
    if ((2*x - 1)*(2*x - 1) + (2*y - 1)*(2*y - 1) < 1)
      cnt += 1;
    // a dot === a one-by-one rectangle
    if ( num <= maxNumDraw )
      ctx.fillRect( x * w, y * h, 1, 1);
  }
  // write the output
  document.getElementById("PiEstimate").innerHTML = 4.0*cnt/num;
}

mcgetPi();


window.onload = function() {
    var S=100;
    var N=100;
    var phi=[new Array(N), new Array(N), new Array(N)];
    var w=new Array(N);
    var color=new Array(N);
    var wid=[0, 1.1, 2*Math.PI];
    var timer;
    var i, j;
    var interval=40;
    var canvas = [document.getElementById('cvs1'),document.getElementById('cvs2'),document.getElementById('cvs3')];
    if (!canvas[2]) {
	canvas.pop();
    }
    var ctx=[];
    for (i in canvas) {
	if (!canvas[i].getContext) {
	    return false;
	}
	canvas[i].width = 2*S;
	canvas[i].height = 2*S;
	ctx.push(canvas[i].getContext('2d'));
    }
    
    var hsv2rgb=function(H, S, V, alpha) {
	var Hi = Math.floor(H*6);
	var f = H*6-Hi;
	var p = V*(1-S);
	var q = V*(1-f*S);
	var t = V*(1-(1-f)*S);
	if (Hi==1) {
            R = q;
            G = V;
            B = p;
	} else if (Hi==2) {
            R = p;
            G = V;
            B = t;
	} else if (Hi==3) {
            R = p;
            G = q;
            B = V;
	} else if (Hi==4) {
            R = t;
            G = p;
            B = V;
	} else if (Hi==5) {
            R = V;
            G = p;
            B = q;
	} else {
            R = V;
            G = t;
            B = p;
	}
	return 'rgba('+Math.floor(256*R)+','+Math.floor(256*G)+','+Math.floor(256*B)+', '+alpha+')';
    };
    
    for (i in phi) {
	for (j=0 ; j<N ; j++) {
	    phi[i][j] = wid[i]*Math.PI*Math.random();
	}
    }
    for (j=0 ; j<N ; j++) {
	var U1=Math.random(), U2=Math.random();
	var mean=0.5;
	w[j] = mean+Math.sqrt(-2*Math.log(U1))*Math.cos(2*Math.PI*U2);
	color[j] = hsv2rgb(0.5+0.17*(w[j]-mean), 1, 1, 0.7);
    }

    var circle = function(x, y, r, color, i) {
	ctx[i].fillStyle = color;
	ctx[i].beginPath();
	ctx[i].arc(x, y, r, 0, Math.PI*2, false);
	ctx[i].fill();
    };

    var framecircle = function(x, y, r, color, i) {
	ctx[i].strokeStyle = color;
	ctx[i].beginPath();
	ctx[i].arc(x, y, r, 0, Math.PI*2, false);
	ctx[i].lineWidth = 1;
	ctx[i].stroke();
    };

    var strokeline = function(x1, y1, x2, y2, color, width, i) {
	ctx[i].strokeStyle = color;
	ctx[i].beginPath();
	ctx[i].moveTo(x1, y1);
	ctx[i].lineTo(x2, y2);
	ctx[i].lineWidth = width;
	ctx[i].stroke();
    };

    var move = function() {
	var h;
	if (!window.innerHeight) {
	    h = document.body.clientHeight;
	} else {
	    h = window.innerHeight;
	}
	if ((!document.hasFocus || document.hasFocus()) && canvas[0].getBoundingClientRect().top<h && canvas[1].getBoundingClientRect().bottom>0) {
	    var K=4;
	    var R=90, r=10;
	    var dt=0.1;
	    for (i in ctx) {
		ctx[i].clearRect(0, 0, 2*S, 2*S);
		framecircle(S, S, R, "rgba(0, 0, 0, 1)", i);
		var csum=0, ssum=0;
		for (j=0 ; j<N ; j++) {
		    csum += Math.cos(phi[i][j]);
		    ssum += Math.sin(phi[i][j]);
		}
		csum /= N;
		ssum /= N;
		for (j=0 ; j<N ; j++) {
		    var cosphi=Math.cos(phi[i][j]);
		    var sinphi=Math.sin(phi[i][j]);
		    circle(S+R*cosphi, S+R*sinphi, r, color[j], i);
		    phi[i][j] += dt*(w[j]+2*K*(ssum*cosphi-csum*sinphi)*(ssum*sinphi+csum*cosphi));
		    if (phi[i][j]>2*Math.PI) {
			phi[i][j] -= 2*Math.PI;
		    }
		}
		strokeline(S, S, S+R*csum, S+R*ssum, "rgba(255,0,0,1)", 3, i);
	    }
	}
	clearTimeout(timer);
	return timer=setTimeout(move, interval);
    };

    clearTimeout(timer);
    return timer=setTimeout(move, interval);
};

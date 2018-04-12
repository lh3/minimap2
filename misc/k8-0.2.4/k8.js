/*********************************************
 * getopt(): translated from the BSD version *
 *********************************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/* // getopt() example
var c;
while ((c = getopt(arguments, "i:xy")) != null) {
	switch (c) {
		case 'i': print(getopt.arg); break;
		case 'x': print("x"); break;
		case 'y': print("y"); break;
	}
}
*/

// print an object in a format similar to JSON
function obj2str(o)
{
	if (typeof(o) != 'object') {
		return o.toString();
	} else if (o == null) {
		return "null";
	} else if (Array.isArray(o)) {
		var s = "[";
		for (var i = 0; i < o.length; ++i) {
			if (i) s += ',';
			s += obj2str(o[i]);
		}
		return s + "]";
	} else {
		var i = 0, s = "{";
		for (var key in o) {
			if (i++) s += ',';
			s += key + ":";
			s += obj2str(o[key]);
		}
		return s + "}";
	}
}

/*********************
 * Special functions *
 *********************/

Math.lgamma = function(z) {
	var x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return Math.log(x) - 5.58106146679532777 - z + (z-0.5) * Math.log(z+6.5);
}

function _kf_gammap(s, z)
{
	var sum, x, k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < 1e-290) break;
	}
	return Math.exp(s * Math.log(z) - z - Math.lgamma(s + 1.) + Math.log(sum));
}

function _kf_gammaq(s, z)
{
	var C, D, f, KF_TINY = 1e-14;
	f = 1. + z - s; C = f; D = 0.;
	for (var j = 1; j < 100; ++j) {
		var a = j * (s - j), b = (j<<1) + 1 + z - s, d;
		D = b + a * D;
		if (D < KF_TINY) D = KF_TINY;
		C = b + a / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (Math.abs(d - 1.) < 1e-290) break;
	}
	return Math.exp(s * Math.log(z) - z - Math.lgamma(s) - Math.log(f));
}

Math.gammap = function(s, z) { return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z); }
Math.gammaq = function(s, z) { return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z); }
Math.chi2 = function(d, x2) { return Math.gammaq(.5 * d, .5 * x2); }

Math.erfc = function(x)
{
	var expntl, z, p;
	z = Math.abs(x) * 1.41421356237309504880168872420969808;
	if (z > 37.) return x > 0.? 0. : 2.;
	expntl = Math.exp(z * z * - .5);
	if (z < 7.07106781186547524400844362104849039) // for small z
	    p = expntl * ((((((.03526249659989109 * z + .7003830644436881) * z + 6.37396220353165) * z + 33.912866078383) * z + 112.0792914978709) * z + 221.2135961699311) * z + 220.2068679123761)
			/ (((((((.08838834764831844 * z + 1.755667163182642) * z + 16.06417757920695) * z + 86.78073220294608) * z + 296.5642487796737) * z + 637.3336333788311) * z + 793.8265125199484) * z + 440.4137358247522);
	else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
	return x > 0.? 2. * p : 2. * (1. - p);
}

Math.normald = function(x) { return .5 * Math.erfc(-x * 0.707106781186547524400844362104849039); }

Math.spearman = function(a)
{
	function aux_func(t) {
		return t == 1? 0 : (t * t - 1) * t / 12
	}
	var x = [], T = [], S = [];
	for (var i = 0; i < a.length; ++i)
		x[i] = [parseFloat(a[i][0]), parseFloat(a[i][1]), 0, 0]
	for (var k = 0; k < 2; ++k) {
		x.sort(function(a,b){return a[k]-b[k]});
		var same = 1;
		T[k] = 0;
		for (var i = 1; i <= x.length; ++i) {
			if (i < x.length && x[i-1][k] == x[i][k]) ++same;
			else {
				var rank = (i<<1) - same + 1;
				for (var j = i - same; j < i; ++j) x[j][k+2] = rank;
				if (same > 1) T[k] += aux_func(same), same = 1;
			}
		}
		S[k] = aux_func(x.length) - T[k];
	}
	var sum = 0.;
	for (var i = 0; i < x.length; ++i)
		sum += .25 * (x[i][2] - x[i][3]) * (x[i][2] - x[i][3]);
	return .5 * (S[0] + S[1] - sum) / Math.sqrt(S[0] * S[1]);
}

Math.fisher_exact = function(n11, n12, n21, n22)
{
	function lbinom(n, k) {
		if (k == 0 || n == k) return 0;
		return Math.lgamma(n+1) - Math.lgamma(k+1) - Math.lgamma(n-k+1);
	}

	function hypergeo(n11, n1_, n_1, n) {
		return Math.exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
	}

	function hypergeo_acc(n11, n1_, n_1, n, aux) {
		if (n1_ || n_1 || n) {
			aux.n11 = n11; aux.n1_ = n1_; aux.n_1 = n_1; aux.n = n;
		} else { // then only n11 changed; the rest fixed
			if (n11%11 && n11 + aux.n - aux.n1_ - aux.n_1) {
				if (n11 == aux.n11 + 1) { // incremental
					aux.p *= (aux.n1_ - aux.n11) / n11
						* (aux.n_1 - aux.n11) / (n11 + aux.n - aux.n1_ - aux.n_1);
					aux.n11 = n11;
					return aux.p;
				}
				if (n11 == aux.n11 - 1) { // incremental
					aux.p *= aux.n11 / (aux.n1_ - n11)
						* (aux.n11 + aux.n - aux.n1_ - aux.n_1) / (aux.n_1 - n11);
					aux.n11 = n11;
					return aux.p;
				}
			}
			aux.n11 = n11;
		}
		aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
		return aux.p;
	}

	var i, j, max, min;
	var p, q, left, right, two;
	var _aux = { n11:0, n1_:0, n_1:0, n:0, p:0. };
	var n1_, n_1, n;

	n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
	max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
	min = n1_ + n_1 - n;
	if (min < 0) min = 0; // min n11, for left tail
	if (min == max) return [1., 1., 1.]; // no need to do test
	q = hypergeo_acc(n11, n1_, n_1, n, _aux); // the probability of the current table
	// left tail
	p = hypergeo_acc(min, 0, 0, 0, _aux);
	for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) // loop until underflow
		left += p, p = hypergeo_acc(i, 0, 0, 0, _aux);
	--i;
	if (p < 1.00000001 * q) left += p;
	else --i;
	// right tail
	p = hypergeo_acc(max, 0, 0, 0, _aux);
	for (right = 0., j = max - 1; p < 0.99999999 * q; --j) // loop until underflow
		right += p, p = hypergeo_acc(j, 0, 0, 0, _aux);
	++j;
	if (p < 1.00000001 * q) right += p;
	else ++j;
	// two-tail
	two = left + right;
	if (two > 1.) two = 1.;
	// adjust left and right
	if (Math.abs(i - n11) < Math.abs(j - n11)) right = 1. - left + q;
	else left = 1.0 - right + q;
	return [two, left, right];
}

/*********************
 * Matrix operations *
 *********************/

Math.m = {};

Math.m.T = function(a) { // matrix transpose
	var b = [], m = a.length, n = a[0].length; // m rows and n cols
	for (var j = 0; j < n; ++j) b[j] = [];
	for (var i = 0; i < m; ++i)
		for (var j = 0; j < n; ++j)
			b[j].push(a[i][j]);
	return b;
}

Math.m.print = function(a) { // print a matrix to stdout
	var m = a.length, n = a[0].length;
	for (var i = 0; i < m; ++i) {
		var line = '';
		for (var j = 0; j < n; ++j)
			line += (j? "\t" : '') + a[i][j];
		print(line);
	}
}

Math.m.mul = function(a, b) { // matrix mul
	var m = a.length, n = a[0].length, s = b.length, t = b[0].length;
	if (n != s) return null;
	var x = [], c = Math.m.T(b);
	for (var i = 0; i < m; ++i) {
		x[i] = [];
		for (var j = 0; j < t; ++j) {
			var sum = 0;
			var ai = a[i], cj = c[j];
			for (var k = 0; k < n; ++k) sum += ai[k] * cj[k];
			x[i].push(sum);
		}
	}
	return x;
}

Math.m.add = function(a, b) { // matrix add
	var m = a.length, n = a[0].length, s = b.length, t = b[0].length;
	var x = [];
	if (m != s || n != t) return null; // different dimensions
	for (var i = 0; i < m; ++i) {
		x[i] = [];
		var ai = a[i], bi = b[i];
		for (var j = 0; j < n; ++j)
			x[i].push(ai[j] + bi[j]);
	}
	return x;
}

Math.m.solve = function(a, b) { // Gauss-Jordan elimination, translated from gaussj() in Numerical Recipes in C.
	// on return, a[n][n] is the inverse; b[n][m] is the solution
	var n = a.length, m = (b)? b[0].length : 0;
	if (a[0].length != n || (b && b.length != n)) return -1; // invalid input
	var xc = [], xr = [], ipiv = [];
	var i, ic, ir, j, l, tmp;

	for (j = 0; j < n; ++j) ipiv[j] = 0;
	for (i = 0; i < n; ++i) {
		var big = 0;
		for (j = 0; j < n; ++j) {
			if (ipiv[j] != 1) {
				for (k = 0; k < n; ++k) {
					if (ipiv[k] == 0) {
						if (Math.abs(a[j][k]) >= big) {
							big = Math.abs(a[j][k]);
							ir = j; ic = k;
						}
					} else if (ipiv[k] > 1) return -2; // singular matrix
				}
			}
		}
		++ipiv[ic];
		if (ir != ic) {
			for (l = 0; l < n; ++l) tmp = a[ir][l], a[ir][l] = a[ic][l], a[ic][l] = tmp;
			if (b) for (l = 0; l < m; ++l) tmp = b[ir][l], b[ir][l] = b[ic][l], b[ic][l] = tmp;
		}
		xr[i] = ir; xc[i] = ic;
		if (a[ic][ic] == 0) return -3; // singular matrix
		var pivinv = 1. / a[ic][ic];
		a[ic][ic] = 1.;
		for (l = 0; l < n; ++l) a[ic][l] *= pivinv;
		if (b) for (l = 0; l < m; ++l) b[ic][l] *= pivinv;
		for (var ll = 0; ll < n; ++ll) {
			if (ll != ic) {
				var dum = a[ll][ic];
				a[ll][ic] = 0;
				for (l = 0; l < n; ++l) a[ll][l] -= a[ic][l] * dum;
				if (b) for (l = 0; l < m; ++l) b[ll][l] -= b[ic][l] * dum;
			}
		}
	}
	for (l = n - 1; l >= 0; --l)
		if (xr[l] != xc[l])
			for (var k = 0; k < n; ++k)
				tmp = a[k][xr[l]], a[k][xr[l]] = a[k][xc[l]], a[k][xc[l]] = tmp;
	return 0;
}

Math.m.dup = function (a)
{
	var b = [];
	for (var i = 0; i < a.length; ++i)
		b[i] = a[i].slice(0);
	return b;
}

Math.m.eigen = function(a)
{
	function SQR(a) { return a * a; }
	function pythag(a, b)
	{
		var absa = Math.abs(a), absb = Math.abs(b);
		if (absa > absb) return absa * Math.sqrt(1.0 + SQR(absb / absa));
		else return (absb == 0.0 ? 0.0 : absb * Math.sqrt(1.0 + SQR(absa / absb)));
	}
	function tred2(a, d, e)
	{
		var n = a.length, l, i, j, k; // integers
		var scale, hh, h, g, f; // real numbers

		for (i = n - 1; i > 0; i--) {
			l = i - 1;
			h = scale = 0.0;
			if (l > 0) {
				for (k = 0; k < l + 1; k++)
					scale += Math.abs(a[i][k]);
				if (scale == 0.0)
					e[i] = a[i][l];
				else {
					for (k = 0; k < l + 1; k++) {
						a[i][k] /= scale;
						h += a[i][k] * a[i][k];
					}
					f = a[i][l];
					g = (f >= 0.0 ? -Math.sqrt(h) : Math.sqrt(h));
					e[i] = scale * g;
					h -= f * g;
					a[i][l] = f - g;
					f = 0.0;
					for (j = 0; j < l + 1; j++) {
						/* Next statement can be omitted if eigenvectors not wanted */
						a[j][i] = a[i][j] / h;
						g = 0.0;
						for (k = 0; k < j + 1; k++)
							g += a[j][k] * a[i][k];
						for (k = j + 1; k < l + 1; k++)
							g += a[k][j] * a[i][k];
						e[j] = g / h;
						f += e[j] * a[i][j];
					}
					hh = f / (h + h);
					for (j = 0; j < l + 1; j++) {
						f = a[i][j];
						e[j] = g = e[j] - hh * f;
						for (k = 0; k < j + 1; k++)
							a[j][k] -= (f * e[k] + g * a[i][k]);
					}
				}
			} else
				e[i] = a[i][l];
			d[i] = h;
		}
		/* Next statement can be omitted if eigenvectors not wanted */
		d[0] = 0.0;
		e[0] = 0.0;
		/* Contents of this loop can be omitted if eigenvectors not wanted except for statement d[i]=a[i][i]; */
		for (i = 0; i < n; i++) {
			l = i;
			if (d[i] != 0.0) {
				for (j = 0; j < l; j++) {
					g = 0.0;
					for (k = 0; k < l; k++)
						g += a[i][k] * a[k][j];
					for (k = 0; k < l; k++)
						a[k][j] -= g * a[k][i];
				}
			}
			d[i] = a[i][i];
			a[i][i] = 1.0;
			for (j = 0; j < l; j++)
				a[j][i] = a[i][j] = 0.0;
		}
	}
	function SIGN(a, b) { return b >= 0.? Math.abs(a) : -Math.abs(a); }
	function tqli(d, e, z)
	{
		var n = d.length, m, l, iter, i, k; // integers
		var s, r, p, g, f, dd, c, b; // real numbers

		for (i = 1; i < n; i++)
			e[i - 1] = e[i];
		e[n - 1] = 0.0;
		for (l = 0; l < n; l++) {
			iter = 0;
			do {
				for (m = l; m < n - 1; m++) {
					dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
					if (Math.abs(e[m]) + dd == dd)
						break;
				}
				if (m != l) {
					if (iter++ == 30) {
						warn("[tqli] Too many iterations in tqli.\n");
						break;
					}
					g = (d[l + 1] - d[l]) / (2.0 * e[l]);
					r = pythag(g, 1.0);
					g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
					s = c = 1.0;
					p = 0.0;
					for (i = m - 1; i >= l; i--) {
						f = s * e[i];
						b = c * e[i];
						e[i + 1] = (r = pythag(f, g));
						if (r == 0.0) {
							d[i + 1] -= p;
							e[m] = 0.0;
							break;
						}
						s = f / r;
						c = g / r;
						g = d[i + 1] - p;
						r = (d[i] - g) * s + 2.0 * c * b;
						d[i + 1] = g + (p = s * r);
						g = c * r - b;
						/* Next loop can be omitted if eigenvectors not wanted */
						for (k = 0; k < n; k++) {
							f = z[k][i + 1];
							z[k][i + 1] = s * z[k][i] + c * f;
							z[k][i] = c * z[k][i] - s * f;
						}
					}
					if (r == 0.0 && i >= l)
						continue;
					d[l] -= p;
					e[l] = g;
					e[m] = 0.0;
				}
			} while (m != l);
		}
	}

	var ev = [], e = [];
	tred2(a, ev, e);
	tqli(ev, e, a);
	return ev;
}

Math.m.Laplacian = function(a, is_norm)
{
	var d = [], n = a.length, L = [];
	for (var i = 0; i < n; ++i) {
		var s = 0;
		for (var j = 0; j < n; ++j)
			s += a[i][j];
		if (s == 0) return null;
		d[i] = s;
		L[i] = [];
	}
	if (is_norm) {
		for (var i = 0; i < n; ++i)
			d[i] = 1. / Math.sqrt(d[i]);
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < i; ++j)
				L[i][j] = -a[i][j] * d[i] * d[j];
			L[i][i] = 1 - a[i][i] * d[i] * d[i];
			for (var j = i + 1; j < n; ++j)
				L[i][j] = -a[i][j] * d[i] * d[j];
		}
	} else {
		for (var i = 0; i < n; ++i) {
			for (var j = 0; j < i; ++j)
				L[i][j] = -a[i][j];
			L[i][i] = d[i] - a[i][i];
			for (var j = i + 1; j < n; ++j)
				L[i][j] = -a[i][j];
		}
	}
	return L;
}

/* // Math.m example
x = [[1,2],[3,4]]; y = [[1,2],[3,4]];
Math.m.print(Math.m.add(Math.m.mul(x,y), x));
a = [[1,1],[1,-1]]; b = [[2],[0]];
a = [[1,2],[3,4]];
print(Math.m.solve(a));
Math.m.print(a);
*/
/* // eigen value example
var x = [[10, 1, 2, 3, 4],
		 [ 1, 9,-1, 2,-3],
		 [ 2,-1, 7, 3,-5],
		 [ 3, 2, 3,12,-1],
		 [ 4,-3,-5,-1,15]];
var ev = Math.m.eigen(x);
print(ev.join("\t")); // eigen values
Math.m.print(x); // eigen vectors
*/

/************************************
 * Fasta/Fastq reader in Javascript *
 ************************************/

Fastx = function(f) {
	this._file = f;
	this._last = 0;
	this._line = new Bytes();
	this._finished = false;
	this.s = new Bytes();
	this.q = new Bytes();
	this.n = new Bytes();
	this.c = new Bytes();
}

Fastx.prototype.read = function() {
	var c, f = this._file, line = this._line;
	if (this._last == 0) { // then jump to the next header line
		while ((c = f.read()) != -1 && c != 62 && c != 64);
		if (c == -1) return -1; // end of file
		this._last = c;
	} // else: the first header char has been read in the previous call
	this.c.length = this.s.length = this.q.length = 0;
	if ((c = f.readline(this.n, 0)) < 0) return -1; // normal exit: EOF
	if (c != 10) f.readline(this.c); // read FASTA/Q comment
	if (this.s.capacity == 0) this.s.capacity = 256;
	while ((c = f.read()) != -1 && c != 62 && c != 43 && c != 64) {
		if (c == 10) continue; // skip empty lines
		this.s.set(c);
		f.readline(this.s, 2, this.s.length); // read the rest of the line
	}
	if (c == 62 || c == 64) this._last = c; // the first header char has been read
	if (c != 43) return this.s.length; // FASTA
	this.q.capacity = this.s.capacity;
	c = f.readline(this._line); // skip the rest of '+' line
	if (c < 0) return -2; // error: no quality string
	var size = this.s.length;
	while (f.readline(this.q, 2, this.q.length) >= 0 && this.q.length < size);
	f._last = 0; // we have not come to the next header line
	if (this.q.length != size) return -2; // error: qual string is of a different length
	return size;
}

Fastx.prototype.destroy = function() {
	this.s.destroy(); this.q.destroy(); this.c.destroy(); this.n.destroy(); this._line.destroy();
	if (typeof(this._file.close) == 'object') this._file.close();
}

function intv_ovlp(intv, bits)
{
	if (typeof bits == "undefined") bits = 13;
	intv.sort(function(a,b) {return a[0]-b[0];});
	// merge overlapping regions
	var j = 0;
	for (var i = 1; i < intv.length; ++i) {
		if (intv[j][1] >= intv[i][0])
			intv[j][1] = intv[j][1] > intv[i][1]? intv[j][1] : intv[i][1];
		else intv[++j] = [intv[i][0], intv[i][1]];
	}
	intv.length = j + 1;
	// create the index
	var idx = [], max = 0;
	for (var i = 0; i < intv.length; ++i) {
		var b = intv[i][0]>>bits;
		var e = (intv[i][1]-1)>>bits;
		if (b != e) {
			for (var j = b; j <= e; ++j)
				if (idx[j] == null) idx[j] = i;
		} else if (idx[b] == null) idx[b] = i;
		max = max > e? max : e;
	}
	return function(_b, _e) { // closure
		var x = _b >> bits;
		if (x > max) return false;
		var off = idx[x];
		if (off == null) {
			var i;
			for (i = ((_e - 1) >> bits) - 1; i >= 0; --i)
				if (idx[i] != null) break;
			off = i < 0? 0 : idx[i];
		}
		for (var i = off; i < intv.length && intv[i][0] < _e; ++i)
			if (intv[i][1] > _b) return true;
		return false;
	}
}

/**************************
 * Bioinformatics related *
 **************************/
/*
Bio = {};
Bio.Seq.comp_pair = ['WSATUGCYRKMBDHVNwsatugcyrkmbdhvn', 'WSTAACGRYMKVHDBNwstaacgrymkvhdbn'];
Bio.Seq.ts_table = {'AG':1, 'GA':1, 'CT':2, 'TC':2};
*/
/* // Bio.Seq example
var s = new Bio.Seq();
var f = new File(arguments[0]);
while (s.next(function(){return f.next()}) != null) s.print();
*/

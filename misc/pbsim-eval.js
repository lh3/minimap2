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

var c, max_mapq = 60, mode = 0, err_out_q = 256, print_err = false;
while ((c = getopt(arguments, "Q:")) != null) {
	if (c == 'Q') err_out_q = parseInt(getopt.arg), print_err = true;
}

var file = arguments.length == getopt.ind? new File() : new File(arguments[getopt.ind]);
var buf = new Bytes();

var tot = [], err = [];
for (var q = 0; q <= max_mapq; ++q)
	tot[q] = err[q] = 0;

function is_correct(s, b)
{
	if (s[0] != b[0] || s[3] != b[3]) return false;
	var o, l;
	if (s[1] < b[1]) {
		if (s[2] <= b[1]) return false;
		o = (s[2] < b[2]? s[2] : b[2]) - b[1];
		l = (s[2] > b[2]? s[2] : b[2]) - s[1];
	} else {
		if (b[2] <= s[1]) return false;
		o = (s[2] < b[2]? s[2] : b[2]) - s[1];
		l = (s[2] > b[2]? s[2] : b[2]) - b[1];
	}
	return o/l > .333? true : false;
}

function count_err(qname, a, tot, err, mode)
{
	var s = qname.split("!");
	if (s.length < 5 || (s[4] != '+' && s[4] != '-'))
		throw Error("Failed to parse pbsim2fa read names '" + qname + "'");
	s[2] = parseInt(s[2]);
	s[3] = parseInt(s[3]);
	s.shift(); // skip pbsim orginal read name
	if (mode == 0) { // longest only
		var max = 0, max_i = -1;
		for (var i = 0; i < a.length; ++i)
			if (a[i][2] - a[i][1] > max)
				max = a[i][2] - a[i][1], max_i = i;
		var mapq = a[max_i][4];
		++tot[mapq];
		if (!is_correct(s, a[max_i])) {
			if (mapq >= err_out_q)
				print('E', qname, a[max_i].join("\t"));
			++err[mapq];
		}
	}
}

var lineno = 0, last = null, a = [];
while (file.readline(buf) >= 0) {
	var line = buf.toString();
	++lineno;
	if (line[0] != '@') {
		var t = line.split("\t");
		if (last != t[0]) {
			if (last != null) count_err(last, a, tot, err, mode);
			a = [], last = t[0];
		}
		if (t[4] == '+' || t[4] == '-') { // PAF
			if (/\ts1:i:\d+/.test(line) && !/\ts2:i:\d+/.test(line)) // secondary alignment in minimap2 PAF
				continue;
			var mapq = parseInt(t[11]);
			if (mapq > max_mapq) mapq = max_mapq;
			a.push([t[5], parseInt(t[7]), parseInt(t[8]), t[4], mapq]);
		}
	}
}
if (last != null) count_err(last, a, tot, err, mode);

buf.destroy();
file.close();

var sum_tot = 0, sum_err = 0, q_out = -1, sum_tot2 = 0, sum_err2 = 0;
for (var q = max_mapq; q >= 0; --q) {
	if (tot[q] == 0) continue;
	if (q_out < 0 || err[q] > 0) {
		if (q_out >= 0) print('Q', q_out, sum_tot, sum_err, (sum_err2/sum_tot2).toFixed(9));
		sum_tot = sum_err = 0, q_out = q;
	}
	sum_tot += tot[q], sum_err += err[q];
	sum_tot2 += tot[q], sum_err2 += err[q];
}
print('Q', q_out, sum_tot, sum_err, (sum_err2/sum_tot2).toFixed(9));

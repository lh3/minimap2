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

var re_cs = /([:=*+-])(\d+|[A-Za-z]+)/g;
var c, min_cov_len = 10000, min_var_len = 50000, gap_thres = 50;
while ((c = getopt(arguments, "l:L:g:")) != null) {
	if (c == 'l') min_cov_len = parseInt(getopt.arg);
	else if (c == 'L') min_var_len = parseInt(optarg.arg);
	else if (c == 'g') gap_thres = parseInt(optarg.arg);
}

if (arguments.length == getopt.ind) {
	print("Usage: k8 paf2diff.js [options] <with-cs.paf>");
	print("Options:");
	print("  -l INT    min alignment length to compute coverage ["+min_cov_len+"]");
	print("  -L INT    min alignment length to call variants ["+min_var_len+"]");
	print("  -g INT    short/long gap threshold ["+gap_thres+"]");
	exit(1);
}

var file = new File(arguments[getopt.ind]);
var buf = new Bytes();

var a = [], out = [];
var tot_len = 0, n_sub = 0, n_1_ins = 0, n_1_del = 0, n_short_ins = 0, n_short_del = 0, n_long_ins = 0, n_long_del = 0;
while (file.readline(buf) >= 0) {
	var line = buf.toString();
	if (!/\ts2:i:/.test(line)) continue; // skip secondary alignments
	var m, t = line.split("\t", 12);
	for (var i = 6; i <= 11; ++i)
		t[i] = parseInt(t[i]);
	if (t[10] < min_cov_len) continue;
	var ctg = t[5], x = t[7], end = t[8];
	// output variants ahead of this alignment
	while (out.length) {
		if (out[0][0] != ctg || out[0][2] <= x) {
			print(out[0].join("\t"));
			out.shift();
		} else break;
	}
	// update coverage
	for (var i = 0; i < out.length; ++i)
		if (out[i][1] >= x && out[i][2] <= end)
			++out[i][3];
	// drop alignments that don't overlap with the current one
	var k = 0;
	for (var i = 0; i < a.length; ++i)
		if (a[0][0] == ctg && a[0][2] > x)
			a[k++] = a[i];
	a.length = k;
	// core loop
	if (t[10] >= min_var_len) {
		if ((m = /\tcs:Z:(\S+)/.exec(line)) == null) continue; // no cs tag
		var cs = m[1];
		var blen = 0, n_diff = 0;
		tot_len += t[10];
		while ((m = re_cs.exec(cs)) != null) {
			var cov = 1;
			if (m[1] == '*' || m[1] == '+' || m[1] == '-')
				for (var i = 0; i < a.length; ++i)
					if (a[0][2] > x) ++cov;
			if (m[1] == '=' || m[1] == ':') {
				var l = m[1] == '='? m[2].length : parseInt(m[2]);
				x += l, blen += l;
			} else if (m[1] == '*') {
				++n_sub;
				out.push([t[5], x, x+1, cov, t[11], m[2].charAt(0), m[2].charAt(1)]);
				++x, ++blen, ++n_diff;
			} else if (m[1] == '+') {
				if (m[2].length == 1) ++n_1_ins;
				else if (m[2].length >= gap_thres) ++n_long_ins;
				else ++n_short_ins;
				out.push([t[5], x, x, cov, t[11], '-', m[2]]);
				++blen, ++n_diff;
			} else if (m[1] == '-') {
				if (m[2].length == 1) ++n_1_del;
				else if (m[2].length >= gap_thres) ++n_long_del;
				else ++n_short_del;
				out.push([t[5], x, x + m[2].length, cov, t[11], m[2], '-']);
				++blen, ++n_diff;
			}
		}
	}
	a.push([t[5], t[7], t[8]]);
}
while (out.length) {
	print(out[0].join("\t"));
	out.shift();
}

warn(tot_len + " contig bases considered in calling");
warn(n_sub + " substitutions");
warn(n_1_del + " 1bp deletions");
warn(n_1_ins + " 1bp insertions");
warn(n_short_del + " [2,"+gap_thres+") short deletions");
warn(n_short_ins + " [2,"+gap_thres+") short insertions");
warn(n_long_del + " >="+gap_thres+" long deletions");
warn(n_long_ins + " >="+gap_thres+" long insertions");

buf.destroy();
file.close();

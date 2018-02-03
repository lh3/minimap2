#!/usr/bin/env k8

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

var c, max_mapq = 60, mode = 0, err_out_q = 256, print_err = false, ovlp_ratio = 0.1, cap_short_mapq = false;
while ((c = getopt(arguments, "Q:r:m:c")) != null) {
	if (c == 'Q') err_out_q = parseInt(getopt.arg), print_err = true;
	else if (c == 'r') ovlp_ratio = parseFloat(getopt.arg);
	else if (c == 'm') mode = parseInt(getopt.arg);
	else if (c == 'c') cap_short_mapq = true;
}

if (arguments.length == getopt.ind) {
	warn("Usage: k8 sim-eval.js [options] <in.paf>|<in.sam>");
	warn("Options:");
	warn("  -r FLOAT   mapping correct if overlap_length/union_length>FLOAT [" + ovlp_ratio + "]");
	warn("  -Q INT     print wrong mappings with mapQ>INT [don't print]");
	warn("  -m INT     0: eval the longest aln only; 1: first aln only; 2: all primary aln [0]");
	exit(1);
}

var file = arguments[getopt.ind] == '-'? new File() : new File(arguments[getopt.ind]);
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
	return o/l > ovlp_ratio? true : false;
}

function count_err(qname, a, tot, err, mode)
{
	if (a.length == 0) return;

	var m, s;
	if ((m = /^(\S+)!(\S+)!(\d+)!(\d+)!([\+\-])$/.exec(qname)) != null) { // pbsim single-end reads
		s = [m[1], m[2], parseInt(m[3]), parseInt(m[4]), m[5]];
	} else if ((m = /^(\S+)!(\S+)!(\d+)_(\d+)!(\d+)_(\d+)!([\+\-])([\+\-])\/([12])$/.exec(qname)) != null) { // mason2 paired-end reads
		if (m[9] == '1') {
			s = [m[1], m[2], parseInt(m[3]), parseInt(m[5]), m[7]];
		} else {
			s = [m[1], m[2], parseInt(m[4]), parseInt(m[6]), m[8]];
		}
	} else throw Error("Failed to parse simulated read names '" + qname + "'");
	s.shift(); // skip the orginal read name

	if (mode == 0 || mode == 1) { // longest only or first only
		var max_i = 0;
		if (mode == 0) { // longest only
			var max = 0;
			for (var i = 0; i < a.length; ++i)
				if (a[i][5] > max)
					max = a[i][5], max_i = i;
		}
		var mapq = a[max_i][4];
		++tot[mapq];
		if (!is_correct(s, a[max_i])) {
			if (mapq >= err_out_q)
				print('E', qname, a[max_i].join("\t"));
			++err[mapq];
		}
	} else if (mode == 2) { // all primary mode
		var max_err_mapq = -1, max_mapq = 0, max_err_i = -1;
		if (cap_short_mapq) {
			var max = 0, max_q = 0;
			for (var i = 0; i < a.length; ++i)
				if (a[i][5] > max)
					max = a[i][5], max_q = a[i][4];
			for (var i = 0; i < a.length; ++i)
				a[i][4] = max_q < a[i][4]? max_q : a[i][4];
		}
		for (var i = 0; i < a.length; ++i) {
			max_mapq = max_mapq > a[i][4]? max_mapq : a[i][4];
			if (!is_correct(s, a[i]))
				if (a[i][4] > max_err_mapq)
					max_err_mapq = a[i][4], max_err_i = i;
		}
		if (max_err_mapq >= 0) {
			++tot[max_err_mapq], ++err[max_err_mapq];
			if (max_err_mapq >= err_out_q)
				print('E', qname, a[max_err_i].join("\t"));
		} else ++tot[max_mapq];
	}
}

var lineno = 0, last = null, a = [], n_unmapped = null;
var re_cigar = /(\d+)([MIDSHN])/g;
while (file.readline(buf) >= 0) {
	var m, line = buf.toString();
	++lineno;
	if (line[0] != '@') {
		var t = line.split("\t");
		if (t[4] == '+' || t[4] == '-') { // PAF
			if (last != t[0]) {
				if (last != null) count_err(last, a, tot, err, mode);
				a = [], last = t[0];
			}
			if (/\ts1:i:\d+/.test(line) && !/\ts2:i:\d+/.test(line)) // secondary alignment in minimap2 PAF
				continue;
			var mapq = parseInt(t[11]);
			if (mapq > max_mapq) mapq = max_mapq;
			a.push([t[5], parseInt(t[7]), parseInt(t[8]), t[4], mapq, parseInt(t[9])]);
		} else { // SAM
			var flag = parseInt(t[1]);
			var read_no = flag>>6&0x3;
			var qname = t[0];
			if (!/\/[12]$/.test(qname))
				qname = read_no == 1 || read_no == 2? t[0] + '/' + read_no : t[0];
			if (last != qname) {
				if (last != null) count_err(last, a, tot, err, mode);
				a = [], last = qname;
			}
			if (flag&0x100) continue; // secondary alignment
			if ((flag&0x4) || t[2] == '*') { // unmapped
				if (n_unmapped == null) n_unmapped = 0;
				++n_unmapped;
				continue;
			}
			var mapq = parseInt(t[4]);
			if (mapq > max_mapq) mapq = max_mapq;
			var pos = parseInt(t[3]) - 1, pos_end = pos;
			var n_gap = 0, mlen = 0;
			while ((m = re_cigar.exec(t[5])) != null) {
				var len = parseInt(m[1]);
				if (m[2] == 'M') pos_end += len, mlen += len;
				else if (m[2] == 'I') n_gap += len;
				else if (m[2] == 'D') n_gap += len, pos_end += len;
			}
			var score = pos_end - pos;
			if ((m = /\tNM:i:(\d+)/.exec(line)) != null) {
				var NM = parseInt(m[1]);
				if (NM >= n_gap) score = mlen - (NM - n_gap);
			}
			a.push([t[2], pos, pos_end, (flag&16)? '-' : '+', mapq, score]);
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
		if (q_out >= 0) print('Q', q_out, sum_tot, sum_err, (sum_err2/sum_tot2).toFixed(9), sum_tot2);
		sum_tot = sum_err = 0, q_out = q;
	}
	sum_tot += tot[q], sum_err += err[q];
	sum_tot2 += tot[q], sum_err2 += err[q];
}
print('Q', q_out, sum_tot, sum_err, (sum_err2/sum_tot2).toFixed(9), sum_tot2);
if (n_unmapped != null) print('U', n_unmapped);

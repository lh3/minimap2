#!/usr/bin/env k8

/*****************************
 ***** Library functions *****
 *****************************/

/*******************************
 * Command line option parsing *
 *******************************/

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

/***********************
 * Interval operations *
 ***********************/

Interval = {};

Interval.sort = function(a)
{
	if (typeof a[0] == 'number')
		a.sort(function(x, y) { return x - y });
	else a.sort(function(x, y) { return x[0] != y[0]? x[0] - y[0] : x[1] - y[1] });
}

Interval.merge = function(a, sorted)
{
	if (typeof sorted == 'undefined') sorted = true;
	if (!sorted) Interval.sort(a);
	var k = 0;
	for (var i = 1; i < a.length; ++i) {
		if (a[k][1] >= a[i][0])
			a[k][1] = a[k][1] > a[i][1]? a[k][1] : a[i][1];
		else a[++k] = a[i].slice(0);
	}
	a.length = k + 1;
}

Interval.index_end = function(a, sorted)
{
	if (a.length == 0) return;
	if (typeof sorted == 'undefined') sorted = true;
	if (!sorted) Interval.sort(a);
	a[0].push(0);
	var k = 0, k_en = a[0][1];
	for (var i = 1; i < a.length; ++i) {
		if (k_en <= a[i][0]) {
			for (++k; k < i; ++k)
				if (a[k][1] > a[i][0])
					break;
			k_en = a[k][1];
		}
		a[i].push(k);
	}
}

Interval.find_intv = function(a, x)
{
	var left = -1, right = a.length;
	if (typeof a[0] == 'number') {
		while (right - left > 1) {
			var mid = left + ((right - left) >> 1);
			if (a[mid] > x) right = mid;
			else if (a[mid] < x) left = mid;
			else return mid;
		}
	} else {
		while (right - left > 1) {
			var mid = left + ((right - left) >> 1);
			if (a[mid][0] > x) right = mid;
			else if (a[mid][0] < x) left = mid;
			else return mid;
		}
	}
	return left;
}

Interval.find_ovlp = function(a, st, en)
{
	if (a.length == 0 || st >= en) return [];
	var l = Interval.find_intv(a, st);
	var k = l < 0? 0 : a[l][a[l].length - 1];
	var b = [];
	for (var i = k; i < a.length; ++i) {
		if (a[i][0] >= en) break;
		else if (st < a[i][1])
			b.push(a[i]);
	}
	return b;
}

/**********************************
 * Reverse and reverse complement *
 **********************************/

Bytes.prototype.reverse = function()
{
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[i];
		this[i] = this[this.length - i - 1];
		this[this.length - i - 1] = tmp;
	}
}

// reverse complement a DNA string
Bytes.prototype.revcomp = function()
{
	if (Bytes.rctab == null) {
		var s1 = 'WSATUGCYRKMBDHVNwsatugcyrkmbdhvn';
		var s2 = 'WSTAACGRYMKVHDBNwstaacgrymkvhdbn';
		Bytes.rctab = [];
		for (var i = 0; i < 256; ++i) Bytes.rctab[i] = 0;
		for (var i = 0; i < s1.length; ++i)
			Bytes.rctab[s1.charCodeAt(i)] = s2.charCodeAt(i);
	}
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[this.length - i - 1];
		this[this.length - i - 1] = Bytes.rctab[this[i]];
		this[i] = Bytes.rctab[tmp];
	}
	if (this.length&1)
		this[this.length>>1] = Bytes.rctab[this[this.length>>1]];
}

/********************
 ***** paftools *****
 ********************/

// variant calling
function paf_call(args)
{
	var re_cs = /([:=*+-])(\d+|[A-Za-z]+)/g;
	var c, min_cov_len = 10000, min_var_len = 50000, gap_thres = 50, min_mapq = 5;
	while ((c = getopt(args, "l:L:g:q:B:")) != null) {
		if (c == 'l') min_cov_len = parseInt(getopt.arg);
		else if (c == 'L') min_var_len = parseInt(getopt.arg);
		else if (c == 'g') gap_thres = parseInt(getopt.arg);
		else if (c == 'q') min_mapq = parseInt(getopt.arg);
	}

	if (args.length == getopt.ind) {
		print("Usage: sort -k6,6 -k8,8n <with-cs.paf> | paftools.js call [options] -");
		print("Options:");
		print("  -l INT    min alignment length to compute coverage ["+min_cov_len+"]");
		print("  -L INT    min alignment length to call variants ["+min_var_len+"]");
		print("  -q INT    min mapping quality ["+min_mapq+"]");
		print("  -g INT    short/long gap threshold (for statistics only) ["+gap_thres+"]");
		exit(1);
	}

	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var buf = new Bytes();
	var tot_len = 0, n_sub = [0, 0, 0], n_ins = [0, 0, 0, 0], n_del = [0, 0, 0, 0];

	function count_var(o)
	{
		if (o[3] > 1) return;
		if (o[5] == '-' && o[6] == '-') return;
		if (o[5] == '-') { // insertion
			var l = o[6].length;
			if (l == 1) ++n_ins[0];
			else if (l == 2) ++n_ins[1];
			else if (l < gap_thres) ++n_ins[2];
			else ++n_ins[3];
		} else if (o[6] == '-') { // deletion
			var l = o[5].length;
			if (l == 1) ++n_del[0];
			else if (l == 2) ++n_del[1];
			else if (l < gap_thres) ++n_del[2];
			else ++n_del[3];
		} else {
			++n_sub[0];
			var s = o[5] + o[6];
			if (s == 'ag' || s == 'ga' || s == 'ct' || s == 'tc')
				++n_sub[1];
			else ++n_sub[2];
		}
	}

	var a = [], out = [];
	var c1_ctg = null, c1_start = 0, c1_end = 0, c1_counted = false, c1_len = 0;
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (!/\ts2:i:/.test(line)) continue; // skip secondary alignments
		var m, t = line.split("\t", 12);
		for (var i = 6; i <= 11; ++i)
			t[i] = parseInt(t[i]);
		if (t[10] < min_cov_len || t[11] < min_mapq) continue;
		for (var i = 1; i <= 3; ++i)
			t[i] = parseInt(t[i]);
		var ctg = t[5], x = t[7], end = t[8];
		var query = t[0], rev = (t[4] == '-'), y = rev? t[3] : t[2];
		// compute regions covered by 1 contig
		if (ctg != c1_ctg || x >= c1_end) {
			if (c1_counted && c1_end > c1_start) {
				c1_len += c1_end - c1_start;
				print('R', c1_ctg, c1_start, c1_end);
			}
			c1_ctg = ctg, c1_start = x, c1_end = end;
			c1_counted = (t[10] >= min_var_len);
		} else if (end > c1_end) { // overlap
			if (c1_counted && x > c1_start) {
				c1_len += x - c1_start;
				print('R', c1_ctg, c1_start, x);
			}
			c1_start = c1_end, c1_end = end;
			c1_counted = (t[10] >= min_var_len);
		} else { // contained
			if (c1_counted && x > c1_start) {
				c1_len += x - c1_start;
				print('R', c1_ctg, c1_start, x);
			}
			c1_start = end;
		}
		// output variants ahead of this alignment
		while (out.length) {
			if (out[0][0] != ctg || out[0][2] <= x) {
				count_var(out[0]);
				print('V', out[0].join("\t"));
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
				var qs, qe;
				if (m[1] == '=' || m[1] == ':') {
					var l = m[1] == '='? m[2].length : parseInt(m[2]);
					if (rev) y -= l;
					else y += l;
					x += l, blen += l;
				} else if (m[1] == '*') {
					if (rev) qs = y - 1, qe = y, --y;
					else qs = y, qe = y + 1, ++y;
					out.push([t[5], x, x+1, cov, t[11], m[2].charAt(0), m[2].charAt(1), query, qs, qe, rev? '-' : '+']);
					++x, ++blen, ++n_diff;
				} else if (m[1] == '+') {
					var l = m[2].length;
					if (rev) qs = y - l, qe = y, y -= l;
					else qs = y, qe = y + l, y += l;
					out.push([t[5], x, x, cov, t[11], '-', m[2], query, qs, qe, rev? '-' : '+']);
					++blen, ++n_diff;
				} else if (m[1] == '-') {
					var l = m[2].length;
					out.push([t[5], x, x + l, cov, t[11], m[2], '-', query, y, y, rev? '-' : '+']);
					x += l, ++blen, ++n_diff;
				}
			}
		}
		a.push([t[5], t[7], t[8]]);
	}
	if (c1_counted && c1_end > c1_start) {
		c1_len += c1_end - c1_start;
		print('R', c1_ctg, c1_start, c1_end);
	}
	while (out.length) {
		count_var(out[0]);
		print('V', out[0].join("\t"));
		out.shift();
	}

	//warn(tot_len + " alignment columns considered in calling");
	warn(c1_len + " reference bases covered by exactly one contig");
	warn(n_sub[0] + " substitutions; ts/tv = " + (n_sub[1]/n_sub[2]).toFixed(3));
	warn(n_del[0] + " 1bp deletions");
	warn(n_ins[0] + " 1bp insertions");
	warn(n_del[1] + " 2bp deletions");
	warn(n_ins[1] + " 2bp insertions");
	warn(n_del[2] + " [3,"+gap_thres+") deletions");
	warn(n_ins[2] + " [3,"+gap_thres+") insertions");
	warn(n_del[3] + " >="+gap_thres+" deletions");
	warn(n_ins[3] + " >="+gap_thres+" insertions");

	buf.destroy();
	file.close();
}

/**************************
 *** Evaluation related ***
 **************************/

// evaluate mapping accuracy
function paf_mapeval(args)
{
	var c, max_mapq = 60, mode = 0, err_out_q = 256, print_err = false, ovlp_ratio = 0.1, cap_short_mapq = false;
	while ((c = getopt(args, "Q:r:m:c")) != null) {
		if (c == 'Q') err_out_q = parseInt(getopt.arg), print_err = true;
		else if (c == 'r') ovlp_ratio = parseFloat(getopt.arg);
		else if (c == 'm') mode = parseInt(getopt.arg);
		else if (c == 'c') cap_short_mapq = true;
	}

	if (args.length == getopt.ind) {
		warn("Usage: paftools.js mapeval [options] <in.paf>|<in.sam>");
		warn("Options:");
		warn("  -r FLOAT   mapping correct if overlap_length/union_length>FLOAT [" + ovlp_ratio + "]");
		warn("  -Q INT     print wrong mappings with mapQ>INT [don't print]");
		warn("  -m INT     0: eval the longest aln only; 1: first aln only; 2: all primary aln [0]");
		exit(1);
	}

	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
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
}

// convert mason2 SAM to FASTQ
function paf_mason2fq(args)
{
	if (args.length == 0) {
		print("Usage: paftools.js mason2fq <mason.sam>");
		exit(1);
	}

	function print_se(a)
	{
		print('@' + a.slice(0, 5).join("!") + " " + a[8]);
		print(a[5]);
		print("+");
		print(a[6]);
	}

	var buf = new Bytes(), buf2 = new Bytes();
	var file = new File(args[0]);
	var re = /(\d+)([MIDSHN])/g;
	var last = null;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0].charAt(0) == '@') continue;
		var m, l_ref = 0;
		while ((m = re.exec(t[5])) != null)
			if (m[2] == 'D' || m[2] == 'M' || m[2] == 'N')
				l_ref += parseInt(m[1]);
		var flag = parseInt(t[1]);
		var rev = !!(flag&16);
		var seq, qual;
		if (rev) {
			buf2.length = 0;
			buf2.set(t[9], 0);
			buf2.revcomp();
			seq = buf2.toString();
			buf2.set(t[10], 0);
			buf2.reverse();
			qual = buf2.toString();
		} else seq = t[9], qual = t[10];
		var qname = t[0];
		qname = qname.replace(/^simulated./, "");
		var chr = t[2];
		var pos = parseInt(t[3]) - 1;
		var strand = (flag&16)? '-' : '+';
		var read_no = flag&0xc0;
		if (read_no == 0x40) read_no = 1;
		else if (read_no == 0x80) read_no = 2;
		else read_no = 0;
		var err = 0, snp = 0, indel = 0;
		for (var i = 11; i < t.length; ++i) {
			if ((m = /^XE:i:(\d+)/.exec(t[i])) != null) err = m[1];
			else if ((m = /^XS:i:(\d+)/.exec(t[i])) != null) snp = m[1];
			else if ((m = /^XI:i:(\d+)/.exec(t[i])) != null) indel = m[1];
		}
		var comment = [err, snp, indel].join(":");
		if (last == null) {
			last = [qname, chr, pos, pos + l_ref, strand, seq, qual, read_no, comment];
		} else if (last[0] != qname) {
			print_se(last);
			last = [qname, chr, pos, pos + l_ref, strand, seq, qual, read_no, comment];
		} else {
			if (read_no == 2) { // last[] is the first read
				if (last[7] != 1) throw Error("ERROR: can't find read1");
				var name = [qname, chr, last[2] + "_" + pos, last[3] + "_" + (pos + l_ref), last[4] + strand].join("!");
				print('@' + name + '/1' + ' ' + last[8]); print(last[5]); print("+"); print(last[6]);
				print('@' + name + '/2' + ' ' + comment); print(seq); print("+"); print(qual);
			} else {
				if (last[7] != 2) throw Error("ERROR: can't find read2");
				var name = [qname, chr, pos + "_" + last[2], (pos + l_ref) + "_" + last[3], strand + last[4]].join("!");
				print('@' + name + '/1' + ' ' + comment); print(seq); print("+"); print(qual);
				print('@' + name + '/2' + ' ' + last[8]); print(last[5]); print("+"); print(last[6]);
			}
			last = null;
		}
	}
	if (last != null) print_se(last);
	file.close();
	buf.destroy();
	buf2.destroy();
}

// convert pbsim MAF to FASTQ
function paf_pbsim2fq(args)
{
	if (args.length < 2) {
		print("Usage: paftools.js pbsim2fq <ref.fa.fai> <pbsim1.maf> [[pbsim2.maf] ...]");
		exit(1);
	}

	var file, buf = new Bytes(), buf2 = new Bytes();
	file = new File(args[0]);
	var chr_list = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split(/\s+/);
		chr_list.push(t[0]);
	}
	file.close();

	for (var k = 1; k < args.length; ++k) {
		var fn = args[k];
		file = new File(fn);
		var state = 0, reg;
		while (file.readline(buf) >= 0) {
			var line = buf.toString();
			if (state == 0 && line.charAt(0) == 'a') {
				state = 1;
			} else if (state == 1 && line.charAt(0) == 's') {
				var t = line.split(/\s+/);
				var st = parseInt(t[2]);
				reg = [st, st + parseInt(t[3])];
				state = 2;
			} else if (state == 2 && line.charAt(0) == 's') {
				var m, t = line.split(/\s+/);
				if ((m = /S(\d+)_\d+/.exec(t[1])) == null) throw Error("Failed to parse the read name");
				var chr_id = parseInt(m[1]) - 1;
				if (chr_id >= chr_list.length) throw Error("Index outside the chr list");
				var name = [t[1], chr_list[chr_id], reg[0], reg[1], t[4]].join("!");
				var seq = t[6].replace(/\-/g, "");
				if (seq.length != parseInt(t[5])) throw Error("Inconsistent read length");
				if (seq.indexOf("NN") < 0) {
					if (t[4] == '-') {
						buf2.set(seq, 0);
						buf2.length = seq.length;
						buf2.revcomp();
						seq = buf2.toString();
					}
					print(">" + name);
					print(seq);
				}
				state = 0;
			}
		}
		file.close();
	}
	buf.destroy();
	buf2.destroy();
}

function paf_ov_eval(args)
{
	var c, min_ovlp = 2000, min_frac = 0.95, min_mapq = 10;
	while ((c = getopt(args, "q:l:f:")) != null) {
		if (c == 'q') min_mapq = parseInt(getopt.arg);
		else if (c == 'l') min_ovlp = parseInt(getopt.arg);
		else if (c == 'f') min_frac = parseFloat(getopt.arg);
	}
	if (args.length - getopt.ind < 2) {
		print("Usage: sort -k6,6 -k8,8n to-ref.paf | paftools.js ov-eval [options] - <ovlp.paf>");
		print("Options:");
		print("  -l INT     min overlap length [2000]");
		print("  -q INT     min mapping quality [10]");
		print("  -f FLOAT   min fraction of mapped length [0.95]");
		exit(1);
	}

	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var a = [], h = {};
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var is_pri = false;
		if (parseInt(t[11]) < min_mapq) continue;
		for (var i = 12; i < t.length; ++i)
			if (t[i] == 'tp:A:P')
				is_pri = true;
		if (!is_pri) continue;
		for (var i = 1; i <= 3; ++i)
			t[i] = parseInt(t[i]);
		for (var i = 6; i <= 8; ++i)
			t[i] = parseInt(t[i]);
		if (t[3] - t[2] < min_ovlp || t[8] - t[7] < min_ovlp || (t[3] - t[2]) / t[1] < min_frac)
			continue;
		var ctg = t[5], st = t[7], en = t[8];
		while (a.length > 0) {
			if (a[0][0] == ctg && a[0][2] > st)
				break;
			else a.shift();
		}
		for (var j = 0; j < a.length; ++j) {
			if (a[j][3] == t[0]) continue;
			var len = (en > a[j][2]? a[j][2] : en) - st;
			if (len >= min_ovlp) {
				var key = a[j][3] < t[0]? a[j][3] + "\t" + t[0] : t[0] + "\t" + a[j][3];
				h[key] = len;
			}
		}
		a.push([ctg, st, en, t[0]]);
	}
	file.close();

	file = new File(args[getopt.ind + 1]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var key = t[0] < t[5]? t[0] + "\t" + t[5] : t[5] + "\t" + t[0];
		if (h[key] > 0) h[key] = -h[key];
	}
	file.close();
	buf.destroy();

	var n_ovlp = 0, n_missing = 0;
	for (var key in h) {
		++n_ovlp;
		if (h[key] > 0) ++n_missing;
	}
	print(n_ovlp + " overlaps inferred from the reference mapping");
	print(n_missing + " missed by the read overlapper");
	print((100 * (1 - n_missing / n_ovlp)).toFixed(2) + "% sensitivity");
}

/*************************
 ***** main function *****
 *************************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: paftools.js <command> [arguments]");
		print("Commands:");
		print("  call       call variants from asm-to-ref alignment with the cs tag");
		print("");
		print("  mapeval    evaluate mapping accuracy using mason2/PBSIM-simulated FASTQ");
		print("  mason2fq   convert mason2-simulated SAM to FASTQ");
		print("  pbsim2fq   convert PBSIM-simulated MAF to FASTQ");
		print("  ov-eval    evaluate read overlap sensitivity using read-to-ref mapping");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'call') paf_call(args);
	else if (cmd == 'mapeval') paf_mapeval(args);
	else if (cmd == 'mason2fq') paf_mason2fq(args);
	else if (cmd == 'pbsim2fq') paf_pbsim2fq(args);
	else if (cmd == 'ov-eval') paf_ov_eval(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);

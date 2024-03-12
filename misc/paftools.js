#!/usr/bin/env k8

var paftools_version = '2.27-r1193';

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

function fasta_read(fn)
{
	var h = {}, seqlen = [];
	var buf = new Bytes();
	var file = fn == '-'? new File() : new File(fn);
	if (typeof k8_version == "undefined") { // for k8-0.x
		var seq = null, name = null, gt = '>'.charCodeAt(0);
		while (file.readline(buf) >= 0) {
			if (buf[0] == gt) {
				if (seq != null && name != null) {
					seqlen.push([name, seq.length]);
					h[name] = seq;
					name = seq = null;
				}
				var m, line = buf.toString();
				if ((m = /^>(\S+)/.exec(line)) != null) {
					name = m[1];
					seq = new Bytes();
				}
			} else seq.set(buf);
		}
		if (seq != null && name != null) {
			seqlen.push([name, seq.length]);
			h[name] = seq;
		}
	} else { // for k8-1.x
		var seq = null, name = null;
		while (file.readline(buf) >= 0) {
			var line = buf.toString();
			if (line[0] == ">") {
				if (seq != null && name != null) {
					seqlen.push([name, seq.length]);
					h[name] = new Uint8Array(seq.buffer);
					name = seq = null;
				}
				var m;
				if ((m = /^>(\S+)/.exec(line)) != null) {
					name = m[1];
					seq = new Bytes();
				}
			} else seq.set(line);
		}
		if (seq != null && name != null) {
			seqlen.push([name, seq.length]);
			h[name] = new Uint8Array(seq.buffer);
		}
	}
	buf.destroy();
	file.close();
	return [h, seqlen];
}

function fasta_free(fa)
{
	if (typeof k8_version == "undefined")
		for (var name in fa)
			fa[name].destroy();
	// FIXME: for k8-1.0, sequences are not freed. This is ok for now but not general.
}

Bytes.prototype.reverse = function()
{
	if (typeof k8_version === "undefined") { // k8-0.x
		for (var i = 0; i < this.length>>1; ++i) {
			var tmp = this[i];
			this[i] = this[this.length - i - 1];
			this[this.length - i - 1] = tmp;
		}
	} else { // k8-1.x
		var buf = new Uint8Array(this.buffer);
		for (var i = 0; i < buf.length>>1; ++i) {
			var tmp = buf[i];
			buf[i] = buf[buf.length - i - 1];
			buf[buf.length - i - 1] = tmp;
		}
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
	if (typeof k8_version === "undefined") { // k8-0.x
		for (var i = 0; i < this.length>>1; ++i) {
			var tmp = this[this.length - i - 1];
			this[this.length - i - 1] = Bytes.rctab[this[i]];
			this[i] = Bytes.rctab[tmp];
		}
		if (this.length&1)
			this[this.length>>1] = Bytes.rctab[this[this.length>>1]];
	} else { // k8-1.x
		var buf = new Uint8Array(this.buffer);
		for (var i = 0; i < buf.length>>1; ++i) {
			var tmp = buf[buf.length - i - 1];
			buf[buf.length - i - 1] = Bytes.rctab[buf[i]];
			buf[i] = Bytes.rctab[tmp];
		}
		if (buf.length&1)
			buf[buf.length>>1] = Bytes.rctab[buf[buf.length>>1]];
	}
}

/********************
 ***** paftools *****
 ********************/

/*****************
 * Miscellaneous *
 *****************/

// liftover
function paf_liftover(args)
{
	function read_bed(fn, to_merge)
	{
		if (fn == null) return null;
		if (typeof to_merge == 'undefined') to_merge = true;
		var file = fn == '-'? new File() : new File(fn);
		var buf = new Bytes();
		var bed = {};
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (bed[t[0]] == null) bed[t[0]] = [];
			bed[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
		}
		buf.destroy();
		file.close();

		for (var chr in bed) {
			Interval.sort(bed[chr]);
			if (to_merge)
				Interval.merge(bed[chr], true);
			Interval.index_end(bed[chr], true);
		}
		return bed;
	}

	var re_cigar = /(\d+)([MID])/g, re_tag = /^(\S\S):([AZif]):(\S+)$/;
	var c, to_merge = false, min_mapq = 5, min_len = 50000, max_div = 2.0;
	var re = /(\d+)([MID])/g;
	while ((c = getopt(args, "mq:l:d:")) != null) {
		if (c == 'm') to_merge = true;
		else if (c == 'q') min_mapq = parseInt(getopt.arg);
		else if (c == 'l') min_len = parseInt(getopt.arg);
		else if (c == 'd') max_div = parseFloat(getopt.arg);
	}
	if (args.length - getopt.ind < 2) {
		print("Usage: paftools.js liftover [options] <aln.paf> <query.bed>");
		print("Options:");
		print("  -q INT    min mapping quality [" + min_mapq + "]");
		print("  -l INT    min alignment length [" + min_len + "]");
		print("  -d FLOAT  max sequence divergence (>=1 to disable) [1]");
		exit(1);
	}
	var bed = read_bed(args[getopt.ind+1], to_merge);

	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");

		if (bed[t[0]] == null) continue; // sequence not present in BED; skip

		// parse tp and cg tags
		var m, tp = null, cg = null;
		for (var i = 12; i < t.length; ++i) {
			if ((m = re_tag.exec(t[i])) != null) {
				if (m[1] == 'tp') tp = m[3];
				else if (m[1] == 'cg') cg = m[3];
			}
		}
		if (tp != 'P' && tp != 'I') continue; // only process primary alignments
		if (cg == null) throw Error("unable to find the 'cg' tag");

		// filter out bad alignments and check overlaps
		for (var i = 1; i <= 3; ++i)
			t[i] = parseInt(t[i]);
		for (var i = 6; i <= 11; ++i)
			t[i] = parseInt(t[i]);
		if (t[11] < min_mapq || t[10] < min_len) continue;
		var regs = Interval.find_ovlp(bed[t[0]], t[2], t[3]);
		if (regs.length == 0) continue; // not overlapping any regions in input BED
		if (max_div >= 0.0 && max_div < 1.0) {
			var n_gaps = 0, n_opens = 0;
			while ((m = re_cigar.exec(cg)) != null)
				if (m[2] == 'I' || m[2] == 'D')
					n_gaps += parseInt(m[1]), ++n_opens;
			var n_mm = t[10] - t[9] - n_gaps;
			var n_diff2 = n_mm + n_opens;
			if (n_diff2 / (n_diff2 + t[9]) > max_div)
				continue;
		}

		// extract start and end positions
		var a = [], r = [], strand = t[4];
		for (var i = 0; i < regs.length; ++i) {
			var s = regs[i][0], e = regs[i][1];
			if (strand == '+') {
				a.push([s, 0, i, -2]);
				a.push([e - 1, 1, i, -2]);
			} else {
				a.push([t[1] - e, 0, i, -2]);
				a.push([t[1] - s - 1, 1, i, -2]);
			}
			r.push([-2, -2]);
		}
		a.sort(function(x, y) { return x[0] - y[0] });

		// lift start/end positions
		var k = 0, x = t[7], y = strand == '+'? t[2] : t[1] - t[3];
		while ((m = re_cigar.exec(cg)) != null) { // TODO: be more careful about edge cases
			var len = parseInt(m[1]);
			if (m[2] == 'D') { // do nothing for D
				x += len;
				continue;
			}
			while (k < a.length && a[k][0] < y) ++k; // skip out-of-range positions
			for (var i = k; i < a.length; ++i) {
				if (y <= a[i][0] && a[i][0] < y + len)
					a[i][3] = m[2] == 'M'? x + (a[i][0] - y) : x;
				else break;
			}
			y += len;
			if (m[2] == 'M') x += len;
		}
		if (x != t[8] || (strand == '+' && y != t[3]) || (strand == '-' && y != t[1] - t[2]))
			throw Error("CIGAR is inconsistent with mapping coordinates");

		// generate result
		for (var i = 0; i < a.length; ++i) {
			if (a[i][1] == 0) r[a[i][2]][0] = a[i][3];
			else r[a[i][2]][1] = a[i][3] + 1; // change to half-close-half-open
		}
		for (var i = 0; i < r.length; ++i) {
			var name = [t[0], regs[i][0], regs[i][1]].join("_");
			if (r[i][0] < 0) name += "_t5", r[i][0] = t[7];
			if (r[i][1] < 0) name += "_t3", r[i][1] = t[8];
			print(t[5], r[i][0], r[i][1], name, 0, strand);
		}
	}
	buf.destroy();
	file.close();
}

// variant calling
function paf_call(args)
{
	var re_cs = /([:=*+-])(\d+|[A-Za-z]+)/g, re_tag = /\t(\S\S:[AZif]):(\S+)/g;
	var c, min_cov_len = 10000, min_var_len = 50000, gap_thres = 50, gap_thres_long = 1000, min_mapq = 5;
	var fa_tmp = null, fa, fa_lens, is_vcf = false, sample_name = "sample";
	while ((c = getopt(args, "l:L:g:q:B:f:s:")) != null) {
		if (c == 'l') min_cov_len = parseInt(getopt.arg);
		else if (c == 'L') min_var_len = parseInt(getopt.arg);
		else if (c == 'g') gap_thres = parseInt(getopt.arg);
		else if (c == 'G') gap_thres_long = parseInt(getopt.arg);
		else if (c == 'q') min_mapq = parseInt(getopt.arg);
		else if (c == 'f') fa_tmp = fasta_read(getopt.arg, fa_lens);
		else if (c == 's') sample_name = getopt.arg;
	}
	if (fa_tmp != null) fa = fa_tmp[0], fa_lens = fa_tmp[1], is_vcf = true;

	if (args.length == getopt.ind) {
		print("Usage: sort -k6,6 -k8,8n <with-cs.paf> | paftools.js call [options] -");
		print("Options:");
		print("  -l INT    min alignment length to compute coverage ["+min_cov_len+"]");
		print("  -L INT    min alignment length to call variants ["+min_var_len+"]");
		print("  -q INT    min mapping quality ["+min_mapq+"]");
		print("  -g INT    short/long gap threshold (for statistics only) ["+gap_thres+"]");
		print("  -f FILE   reference sequences (enabling VCF output) [null]");
		print("  -s NAME   sample name in VCF header ["+sample_name+"]");
		exit(1);
	}

	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var buf = new Bytes();
	var tot_len = 0, n_sub = [0, 0, 0], n_ins = [0, 0, 0, 0, 0], n_del = [0, 0, 0, 0, 0];

	function print_vcf(o, fa)
	{
		var v = null;
		if (o[3] != 1) return; // coverage is one; skip
		if (o[5] == '-' && o[6] == '-') return;
		if (o[5] != '-' && o[6] != '-') { // snp
			v = [o[0], o[1] + 1, '.', o[5].toUpperCase(), o[6].toUpperCase()];
		} else if (o[1] > 0) { // shouldn't happen in theory
			if (fa[o[0]] == null) throw Error('sequence "' + o[0] + '" is absent from the reference FASTA');
			if (o[1] >= fa[o[0]].length) throw Error('position ' + o[1] + ' exceeds the length of sequence "' + o[0] + '"');
			var ref = String.fromCharCode(fa[o[0]][o[1]-1]).toUpperCase();
			if (o[5] == '-') // insertion
				v = [o[0], o[1], '.', ref, ref + o[6].toUpperCase()];
			else // deletion
				v = [o[0], o[1], '.', ref + o[5].toUpperCase(), ref];
		}
		v.push(o[4], '.', 'QNAME=' + o[7] + ';QSTART=' + (o[8]+1) + ';QSTRAND=' + (rev? '-' : '+'), 'GT', '1/1');
		if (v == null) throw Error("unexpected variant: [" + o.join(",") + "]");
		print(v.join("\t"));
	}

	function count_var(o)
	{
		if (o[3] > 1) return;
		if (o[5] == '-' && o[6] == '-') return;
		if (o[5] == '-') { // insertion
			var l = o[6].length;
			if (l == 1) ++n_ins[0];
			else if (l == 2) ++n_ins[1];
			else if (l < gap_thres) ++n_ins[2];
			else if (l < gap_thres_long) ++n_ins[3];
			else ++n_ins[4];
		} else if (o[6] == '-') { // deletion
			var l = o[5].length;
			if (l == 1) ++n_del[0];
			else if (l == 2) ++n_del[1];
			else if (l < gap_thres) ++n_del[2];
			else if (l < gap_thres_long) ++n_del[3];
			else ++n_del[4];
		} else {
			++n_sub[0];
			var s = (o[5] + o[6]).toLowerCase();
			if (s == 'ag' || s == 'ga' || s == 'ct' || s == 'tc')
				++n_sub[1];
			else ++n_sub[2];
		}
	}

	if (is_vcf) {
		print('##fileformat=VCFv4.1');
		for (var i = 0; i < fa_lens.length; ++i)
			print('##contig=<ID=' + fa_lens[i][0] + ',length=' + fa_lens[i][1] + '>');
		print('##INFO=<ID=QNAME,Number=1,Type=String,Description="Query name">');
		print('##INFO=<ID=QSTART,Number=1,Type=Integer,Description="Query start">');
		print('##INFO=<ID=QSTRAND,Number=1,Type=String,Description="Query strand">');
		print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">');
		print('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'+sample_name);
	}

	var a = [], out = [];
	var c1_ctg = null, c1_start = 0, c1_end = 0, c1_counted = false, c1_len = 0;
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var m, t = line.split("\t", 12);
		if (t.length < 12 || t[5] == '*') continue; // unmapped
		for (var i = 6; i <= 11; ++i)
			t[i] = parseInt(t[i]);
		if (t[10] < min_cov_len || t[11] < min_mapq) continue;
		//print(t[0], t[7], t[8], c1_start, c1_end);
		for (var i = 1; i <= 3; ++i)
			t[i] = parseInt(t[i]);
		var ctg = t[5], x = t[7], end = t[8];
		var query = t[0], rev = (t[4] == '-'), y = rev? t[3] : t[2];
		// collect tags
		var cs = null, tp = null, have_s1 = false, have_s2 = false;
		while ((m = re_tag.exec(line)) != null) {
			if (m[1] == 'cs:Z') cs = m[2];
			else if (m[1] == 'tp:A') tp = m[2];
			else if (m[1] == 's1:i') have_s1 = true;
			else if (m[1] == 's2:i') have_s2 = true;
		}
		if (have_s1 && !have_s2) continue;
		if (tp != null && (tp == 'S' || tp == 'i')) continue;
		// compute regions covered by 1 contig
		if (ctg != c1_ctg || x >= c1_end) {
			if (c1_counted && c1_end > c1_start) {
				c1_len += c1_end - c1_start;
				if (!is_vcf) print('R', c1_ctg, c1_start, c1_end);
			}
			c1_ctg = ctg, c1_start = x, c1_end = end;
			c1_counted = (t[10] >= min_var_len);
		} else if (end > c1_end) { // overlap
			if (c1_counted && x > c1_start) {
				c1_len += x - c1_start;
				if (!is_vcf) print('R', c1_ctg, c1_start, x);
			}
			c1_start = c1_end, c1_end = end;
			c1_counted = (t[10] >= min_var_len);
		} else if (end > c1_start) { // contained
			if (c1_counted && x > c1_start) {
				c1_len += x - c1_start;
				if (!is_vcf) print('R', c1_ctg, c1_start, x);
			}
			c1_start = end;
		} // else, the alignment precedes the cov1 region; do nothing
		// output variants ahead of this alignment
		while (out.length) {
			if (out[0][0] != ctg || out[0][2] <= x) {
				count_var(out[0]);
				if (is_vcf) print_vcf(out[0], fa);
				else print('V', out[0].join("\t"));
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
			if (a[i][0] == ctg && a[i][2] > x)
				a[k++] = a[i];
		a.length = k;
		// core loop
		if (t[10] >= min_var_len) {
			if (cs == null) continue; // no cs tag
			var blen = 0, n_diff = 0;
			tot_len += t[10];
			while ((m = re_cs.exec(cs)) != null) {
				var cov = 1;
				if (m[1] == '*' || m[1] == '+' || m[1] == '-')
					for (var i = 0; i < a.length; ++i)
						if (a[i][2] > x) ++cov;
				var qs, qe;
				if (m[1] == '=' || m[1] == ':') {
					var l = m[1] == '='? m[2].length : parseInt(m[2]);
					if (rev) y -= l;
					else y += l;
					x += l, blen += l;
				} else if (m[1] == '*') {
					if (rev) qs = y - 1, qe = y, --y;
					else qs = y, qe = y + 1, ++y;
					var br = m[2].charAt(0), bq = m[2].charAt(1);
					if (br != 'n' && bq != 'n') { // don't call a SNP if there is an ambiguous base
						out.push([t[5], x, x+1, cov, t[11], br, bq, query, qs, qe, rev? '-' : '+']);
						++n_diff;
					}
					++x, ++blen;
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
		if (!is_vcf) print('R', c1_ctg, c1_start, c1_end);
	}
	while (out.length) {
		count_var(out[0]);
		if (is_vcf) print_vcf(out[0], fa);
		else print('V', out[0].join("\t"));
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
	warn(n_del[3] + " ["+gap_thres+","+gap_thres_long+") deletions");
	warn(n_ins[3] + " ["+gap_thres+","+gap_thres_long+") insertions");
	warn(n_del[4] + " >=" + gap_thres_long + " deletions");
	warn(n_ins[4] + " >=" + gap_thres_long + " insertions");

	buf.destroy();
	file.close();
	if (fa != null) fasta_free(fa);
}

function paf_asmstat(args)
{
	var c, min_query_len = 0, min_seg_len = 10000, max_diff = 0.01, bp_flank_len = 0, bp_gap_len = 0;
	while ((c = getopt(args, "l:d:b:g:q:")) != null) {
		if (c == 'l') min_seg_len = parseInt(getopt.arg);
		else if (c == 'd') max_diff = parseFloat(getopt.arg);
		else if (c == 'b') bp_flank_len = parseInt(getopt.arg);
		else if (c == 'g') bp_gap_len = parseInt(getopt.arg);
		else if (c == 'q') min_query_len = parseInt(getopt.arg);
	}
	if (getopt.ind == args.length) {
		print("Usage: paftools.js asmstat [options] <ref.fa.fai> <asm1.paf> [...]");
		print("Options:");
		print("  -q INT     ignore query shorter than INT [0]");
		print("  -l INT     min alignment block length [" + min_seg_len + "]");
		print("  -d FLOAT   max gap-compressed sequence divergence [" + max_diff + "]");
		exit(1);
	}

	var file, buf = new Bytes();

	var ref_len = 0;
	file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		ref_len += parseInt(t[1]);
	}
	file.close();

	function process_query(qblocks, qblock_len, bp, qi) {
		qblocks.sort(function(a,b) { return a[0]-b[0]; });
		var last_k = null, last_blen = null, st = -1, en = -1, qcov = 0;
		for (var k = 0; k < qblocks.length; ++k) {
			var blen = qblocks[k][1] - qblocks[k][0];
			if (k > 0 && qblocks[k][0] < qblocks[k-1][1]) {
				if (qblocks[k][1] < qblocks[k-1][1]) continue;
				blen = qblocks[k][1] - qblocks[k-1][1];
			}
			qblock_len.push(blen);
			if (qblocks[k][0] > en) {
				qcov += en - st;
				st = qblocks[k][0];
				en = qblocks[k][1];
			} else en = en > qblocks[k][1]? en : qblocks[k][1];
			if (last_k != null) {
				var gap = 1000000000;
				if (qblocks[k][2] == qblocks[last_k][2] && qblocks[k][3] == qblocks[last_k][3]) { // same chr and strand
					var g1 = qblocks[k][0] - qblocks[last_k][1];
					var g2 = qblocks[k][2] == '+'? qblocks[k][4] - qblocks[last_k][5] : qblocks[last_k][4] - qblocks[k][5];
					gap = g1 > g2? g1 - g2 : g2 - g1;
				}
				var min = blen < last_blen? blen : last_blen;
				var flank = k == 0? min : blen;
				bp.push([flank, gap]);
				qi.bp.push([flank, gap]);
			}
			last_k = k, last_blen = blen;
		}
		qcov += en - st;
		return qcov;
	}

	function N50(lens, tot, quantile) {
		lens.sort(function(a,b) { return b - a; });
		if (tot == null) {
			tot = 0;
			for (var k = 0; k < lens.length; ++k)
				tot += lens[k];
		}
		var sum = 0;
		for (var k = 0; k < lens.length; ++k) {
			if (sum <= quantile * tot && sum + lens[k] > quantile * tot)
				return lens[k];
			sum += lens[k];
		}
	}

	function AUN(lens, tot) {
		lens.sort(function(a,b) { return b - a; });
		if (tot == null) {
			tot = 0;
			for (var k = 0; k < lens.length; ++k)
				tot += lens[k];
		}
		var x = 0, y = 0;
		for (var k = 0; k < lens.length; ++k) {
			var l = x + lens[k] <= tot? lens[k] : tot - x;
			x += lens[k];
			y += l * (l / tot);
			if (x >= tot) break;
		}
		return y.toFixed(0);
	}

	function count_bp(bp, min_blen, min_gap) {
		var n_bp = 0;
		for (var k = 0; k < bp.length; ++k)
			if (bp[k][0] >= min_blen && bp[k][1] >= min_gap)
				++n_bp;
		return n_bp;
	}

	function compute_diff(cigar, NM) {
		var m, re = /(\d+)([MID])/g;
		var n_M = 0, n_gapo = 0, n_gaps = 0;
		while ((m = re.exec(cigar)) != null) {
			var len = parseInt(m[1]);
			if (m[2] == 'M') n_M += len;
			else ++n_gapo, n_gaps += len;
		}
		if (NM < n_gaps) throw Error('NM is smaller the number of gaps');
		return (NM - n_gaps + n_gapo) / (n_M + n_gapo);
	}

	var labels = ['Length', 'l_cov', 'Rcov', 'Rdup', 'Qcov', 'NG75', 'NG50', 'NGA50', 'AUNGA', '#breaks', 'bp(' + min_seg_len + ',0)', 'bp(' + min_seg_len + ',10k)'];
	var rst = [];
	for (var i = 0; i < labels.length; ++i)
		rst[i] = [];

	var n_asm = args.length - (getopt.ind + 1);
	var header = ["Metric"];
	for (var i = 0; i < n_asm; ++i) {
		var n_breaks = 0, qcov = 0;
		var fn = args[getopt.ind + 1 + i];
		var label = fn.replace(/.paf(.gz)?$/, "");
		header.push(label);
		var ref_blocks = [], qblock_len = [], qblocks = [], bp = [];
		var query = {}, qinfo = {};
		var last_qname = null;
		file = new File(fn);
		while (file.readline(buf) >= 0) {
			var m, line = buf.toString();
			var t = line.split("\t");
			t[1] = parseInt(t[1]);
			if (t[1] < min_query_len) continue;
			if (t.length < 2) continue;
			query[t[0]] = t[1];
			if (qinfo[t[0]] == null) qinfo[t[0]] = {};
			qinfo[t[0]].len = t[1];
			qinfo[t[0]].bp = [];
			if (t.length < 9 || t[5] == "*") continue;
			if (!/\ttp:A:[PI]/.test(line)) continue;
			var cigar = (m = /\tcg:Z:(\S+)/.exec(line)) != null? m[1] : null;
			var NM = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;
			var diff = cigar != null && NM != null? compute_diff(cigar, NM) : 0;
			t[2] = parseInt(t[2]);
			t[3] = parseInt(t[3]);
			t[7] = parseInt(t[7]);
			t[8] = parseInt(t[8]);
			if (t[0] == last_qname) ++n_breaks;
			if (diff > max_diff) continue;
			if (t[3] - t[2] < min_seg_len) continue;
			if (t[0] != last_qname) {
				if (last_qname != null)
					qcov += process_query(qblocks, qblock_len, bp, qinfo[last_qname]);
				qblocks = [];
				last_qname = t[0];
			}
			ref_blocks.push([t[5], t[7], t[8]]);
			qblocks.push([t[2], t[3], t[4], t[5], t[7], t[8]]);
		}
		if (last_qname != null)
			qcov += process_query(qblocks, qblock_len, bp, qinfo[last_qname]);
		file.close();

		// compute NG50
		var asm_len = 0, asm_lens = []
		for (var ctg in query) {
			asm_len += query[ctg];
			asm_lens.push(query[ctg]);
		}
		rst[0][i] = asm_len;
		rst[5][i] = N50(asm_lens, ref_len, 0.75);
		rst[6][i] = N50(asm_lens, ref_len, 0.5);

		// compute coverage
		var l_cov = 0;
		ref_blocks.sort(function(a, b) { return a[0] > b[0]? 1 : a[0] < b[0]? -1 : a[1] - b[1]; });
		var last_ref = null, st = -1, en = -1;
		for (var j = 0; j < ref_blocks.length; ++j) {
			if (ref_blocks[j][0] != last_ref || ref_blocks[j][1] > en) {
				l_cov += en - st;
				last_ref = ref_blocks[j][0];
				st = ref_blocks[j][1];
				en = ref_blocks[j][2];
			} else en = en > ref_blocks[j][2]? en : ref_blocks[j][2];
		}
		l_cov += en - st;
		rst[1][i] = l_cov;
		rst[2][i] = (100.0 * (l_cov / ref_len)).toFixed(2) + '%';
		rst[4][i] = (100.0 * (qcov  / asm_len)).toFixed(2) + '%';

		// compute cov1 and cov2+ lengths; see paf_call() for details
		var c1_ctg = null, c1_start = 0, c1_end = 0, c1_len = 0;
		for (var j = 0; j < ref_blocks.length; ++j) {
			if (ref_blocks[j][0] != c1_ctg || ref_blocks[j][1] >= c1_end) {
				if (c1_end > c1_start)
					c1_len += c1_end - c1_start;
				c1_ctg = ref_blocks[j][0], c1_start = ref_blocks[j][1], c1_end = ref_blocks[j][2];
			} else if (ref_blocks[j][2] > c1_end) { // overlap
				if (ref_blocks[j][1] > c1_start)
					c1_len += ref_blocks[j][1] - c1_start;
				c1_start = c1_end, c1_end = ref_blocks[j][2];
			} else if (ref_blocks[j][2] > c1_start) { // contained
				if (ref_blocks[j][1] > c1_start)
					c1_len += ref_blocks[j][1] - c1_start;
				c1_start = ref_blocks[j][2];
			}
			//print(ref_blocks[j][0], ref_blocks[j][1], ref_blocks[j][2], c1_start, c1_end, c1_len);
		}
		if (c1_end > c1_start)
			c1_len += c1_end - c1_start;
		rst[3][i] = (100 * (l_cov - c1_len) / l_cov).toFixed(2) + '%';

		// compute NGA50
		rst[7][i] = N50(qblock_len, ref_len, 0.5);

		// compute AUNGA
		rst[8][i] = AUN(qblock_len, ref_len);

		// compute break points
		rst[9][i] = n_breaks;
		rst[10][i] = count_bp(bp, 500, 0);
		rst[11][i] = count_bp(bp, 500, 10000);

		// nb-plot; NOT USED
		/*
		var qa = [];
		for (var qn in qinfo)
			qa.push([qinfo[qn].len, qinfo[qn].bp]);
		qa = qa.sort(function(a, b) { return b[0] - a[0] });
		var sum = 0, n_bp = 0, next_quantile = 0.1;
		for (var j = 0; j < qa.length; ++j) {
			sum += qa[j][0];
			for (var k = 0; k < qa[j][1].length; ++k)
				if (qa[j][1][k][0] >= bp_flank_len && qa[j][1][k][1] >= bp_gap_len)
					++n_bp;
			if (sum >= ref_len * next_quantile) {
				print(label, Math.floor(next_quantile * 100 + .5), qa[j][0], (sum / n_bp).toFixed(0), n_bp);
				next_quantile += 0.1;
				if (next_quantile >= 1.0) break;
			}
		}
		*/
	}
	buf.destroy();

	if (bp_flank_len <= 0) {
		print(header.join("\t"));
		for (var i = 0; i < labels.length; ++i)
			print(labels[i], rst[i].join("\t"));
	}
}

function paf_asmgene(args)
{
	var c, opt = { min_cov:0.99, min_iden:0.99 }, print_err = false, auto_only = false;
	while ((c = getopt(args, "i:c:ea")) != null)
		if (c == 'i') opt.min_iden = parseFloat(getopt.arg);
		else if (c == 'c') opt.min_cov = parseFloat(getopt.arg);
		else if (c == 'e') print_err = true;
		else if (c == 'a') auto_only = true;

	var n_fn = args.length - getopt.ind;
	if (n_fn < 2) {
		print("Usage: paftools.js asmgene [options] <ref-splice.paf> <asm-splice.paf> [...]");
		print("Options:");
		print("  -i FLOAT     min identity [" + opt.min_iden + "]");
		print("  -c FLOAT     min coverage [" + opt.min_cov + "]");
		print("  -a           only evaluate genes mapped to the autosomes");
		print("  -e           print fragmented/missing genes");
		exit(1);
	}

	function process_query(opt, a) {
		var b = [], cnt = [0, 0, 0];
		for (var j = 0; j < a.length; ++j) {
			if (a[j][4] < a[j][5] * opt.min_iden)
				continue;
			b.push(a[j].slice(0));
		}
		if (b.length == 0) return cnt;
		// count full
		var n_full = 0;
		for (var j = 0; j < b.length; ++j)
			if (b[j][3] - b[j][2] >= b[j][1] * opt.min_cov)
				++n_full;
		cnt[0] = n_full;
		// compute coverage
		b = b.sort(function(x, y) { return x[2] - y[2] });
		var l_cov = 0, st = b[0][2], en = b[0][3];
		for (var j = 1; j < b.length; ++j) {
			if (b[j][2] <= en)
				en = b[j][3] > en? b[j][3] : en;
			else l_cov += en - st;
		}
		l_cov += en - st;
		cnt[1] = l_cov / b[0][1];
		cnt[2] = b.length;
		return cnt;
	}

	var buf = new Bytes();
	var gene = {}, header = [], refpos = {};
	for (var i = getopt.ind; i < args.length; ++i) {
		var fn = args[i];
		var label = fn.replace(/.paf(.gz)?$/, "");
		header.push(label);
		var file = new File(fn), a = [];
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var ql = parseInt(t[1]), qs = parseInt(t[2]), qe = parseInt(t[3]), mlen = parseInt(t[9]), blen = parseInt(t[10]), mapq = parseInt(t[11]);
			if (i == getopt.ind) refpos[t[0]] = [t[0], t[1], t[5], t[7], t[8]];
			if (gene[t[0]] == null) gene[t[0]] = [];
			if (a.length && t[0] != a[0][0]) {
				gene[a[0][0]][i - getopt.ind] = process_query(opt, a);
				a = [];
			}
			a.push([t[0], ql, qs, qe, mlen, blen]);
		}
		if (a.length)
			gene[t[0]][i - getopt.ind] = process_query(opt, a);
		file.close();
	}

	// select the longest genes (not optimal, but should be good enough)
	var gene_list = [], gene_nr = {};
	for (var g in refpos)
		gene_list.push(refpos[g]);
	gene_list = gene_list.sort(function(a, b) { return a[2] < b[2]? -1 : a[2] > b[2]? 1 : a[3] - b[3] });
	var last = 0;
	for (var j = 1; j < gene_list.length; ++j) {
		if (gene_list[j][2] != gene_list[last][2] || gene_list[j][3] >= gene_list[last][4]) {
			gene_nr[gene_list[last][0]] = 1;
			last = j;
		} else if (gene_list[j][1] > gene_list[last][1]) {
			last = j;
		}
	}
	gene_nr[gene_list[last][0]] = 1;

	// count and print
	var col1 = ["full_sgl", "full_dup", "frag", "part50+", "part10+", "part10-", "dup_cnt", "dup_sum"];
	var rst = [];
	for (var k = 0; k < col1.length; ++k) {
		rst[k] = [];
		for (var i = 0; i < n_fn; ++i)
			rst[k][i] = 0;
	}
	for (var g in gene) { // count single-copy genes
		if (gene[g][0] == null || gene[g][0][0] != 1) continue;
		if (gene_nr[g] == null) continue;
		if (auto_only && /^(chr)?[XY]$/.test(refpos[g][2])) continue;
		for (var i = 0; i < n_fn; ++i) {
			if (gene[g][i] == null) {
				rst[5][i]++;
				if (print_err) print('M', header[i], refpos[g].join("\t"));
			} else if (gene[g][i][0] == 1) {
				rst[0][i]++;
			} else if (gene[g][i][0] > 1) {
				rst[1][i]++;
				if (print_err) print('D', header[i], refpos[g].join("\t"));
			} else if (gene[g][i][1] >= opt.min_cov) {
				rst[2][i]++;
				if (print_err) print('F', header[i], refpos[g].join("\t"));
			} else if (gene[g][i][1] >= 0.5) {
				rst[3][i]++;
				if (print_err) print('5', header[i], refpos[g].join("\t"));
			} else if (gene[g][i][1] >= 0.1) {
				rst[4][i]++;
				if (print_err) print('1', header[i], refpos[g].join("\t"));
			} else {
				rst[5][i]++;
				if (print_err) print('0', header[i], refpos[g].join("\t")); // TODO: reduce code duplicates...
			}
		}
	}
	for (var g in gene) { // count multi-copy genes
		if (gene[g][0] == null || gene[g][0][0] <= 1) continue;
		if (gene_nr[g] == null) continue;
		if (auto_only && /^(chr)?[XY]$/.test(refpos[g][2])) continue;
		for (var i = 0; i < n_fn; ++i) {
			if (gene[g][i] != null) rst[7][i] += gene[g][i][0];
			if (gene[g][i] != null && gene[g][i][0] > 1) {
				rst[6][i]++;
			} else if (print_err) {
				print('d', header[i], gene[g][0][0], refpos[g].join("\t"));
			}
		}
	}
	print('H', 'Metric', header.join("\t"));
	for (var k = 0; k < rst.length; ++k) {
		print('X', col1[k], rst[k].join("\t"));
	}
	buf.destroy();
}

function paf_stat(args)
{
	var c, gap_out_len = null, count_err = false;
	while ((c = getopt(args, "cl:")) != null)
		if (c == 'l') gap_out_len = parseInt(getopt.arg);
		else if (c == 'c') count_err = true;

	if (getopt.ind == args.length) {
		print("Usage: paftools.js stat [-c] [-l gapOutLen] <in.sam>|<in.paf>");
		exit(1);
	}

	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var re = /(\d+)([MIDSHNX=])/g;

	var lineno = 0, n_pri = 0, n_2nd = 0, n_seq = 0, n_cigar_64k = 0, l_tot = 0, l_cov = 0;
	var n_gap = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]], n_sub = 0;

	function cov_len(regs)
	{
		regs.sort(function(a,b) {return a[0]-b[0]});
		var st = regs[0][0], en = regs[0][1], l = 0;
		for (var i = 1; i < regs.length; ++i) {
			if (regs[i][0] < en)
				en = en > regs[i][1]? en : regs[i][1];
			else l += en - st, st = regs[i][0], en = regs[i][1];
		}
		l += en - st;
		return l;
	}

	var last = null, last_qlen = null, regs = [];
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		++lineno;
		if (line.charAt(0) != '@') {
			var t = line.split("\t", 12);
			var m, rs, cigar = null, is_pri = false, is_sam = false, is_rev = false, tname = null;
			var atlen = null, aqlen, qs, qe, mapq, ori_qlen, NM = null, nn = 0;
			if (t.length < 2) continue;
			if (t[4] == '+' || t[4] == '-' || t[4] == '*') { // PAF
				if (t[4] == '*') continue; // unmapped
				if (!/\ts2:i:\d+/.test(line)) {
					++n_2nd;
					continue;
				}
				if ((m = /\tNM:i:(\d+)/.exec(line)) != null)
					NM = parseInt(m[1]);
				if ((m = /\tnn:i:(\d+)/.exec(line)) != null)
					nn = parseInt(m[1]);
				if ((m = /\tcg:Z:(\S+)/.exec(line)) != null)
					cigar = m[1];
				if (cigar == null) {
					warn("WARNING: no CIGAR at line " + lineno);
					continue;
				}
				tname = t[5];
				qs = parseInt(t[2]), qe = parseInt(t[3]);
				aqlen = qe - qs;
				is_rev = t[4] == '+'? false : true;
				rs = parseInt(t[7]);
				atlen = parseInt(t[8]) - rs;
				mapq = parseInt(t[11]);
				ori_qlen = parseInt(t[1]);
			} else { // SAM
				var flag = parseInt(t[1]);
				if ((flag & 4) || t[2] == '*' || t[5] == '*') continue;
				if (flag & 0x100) {
					++n_2nd;
					continue;
				}
				if ((m = /\tNM:i:(\d+)/.exec(line)) != null)
					NM = parseInt(m[1]);
				if ((m = /\tnn:i:(\d+)/.exec(line)) != null)
					nn = parseInt(m[1]);
				cigar = t[5];
				tname = t[2];
				rs = parseInt(t[3]) - 1;
				mapq = parseInt(t[4]);
				aqlen = t[9].length;
				is_sam = true;
				is_rev = !!(flag&0x10);
			}
			++n_pri;
			if (last != t[0]) {
				if (last != null) {
					l_tot += last_qlen;
					l_cov += cov_len(regs);
				}
				regs = [];
				++n_seq, last = t[0];
			}
			var M = 0, tl = 0, ql = 0, clip = [0, 0], n_cigar = 0, sclip = 0;
			var n_gapo = 0, n_gap_all = 0, l_match = 0;
			while ((m = re.exec(cigar)) != null) {
				var l = parseInt(m[1]);
				++n_cigar;
				if (m[2] == 'M' || m[2] == '=' || m[2] == 'X') {
					tl += l, ql += l, M += l;
					l_match += l;
				} else if (m[2] == 'I' || m[2] == 'D') {
					var type;
					if (l < 50) type = 0;
					else if (l < 100) type = 1;
					else if (l < 300) type = 2;
					else if (l < 400) type = 3;
					else if (l < 1000) type = 4;
					else type = 5;
					if (m[2] == 'I') ql += l, ++n_gap[0][type];
					else tl += l, ++n_gap[1][type];
					if (gap_out_len != null && l >= gap_out_len)
						print(t[0], ql, is_rev? '-' : '+', tname, rs + tl, m[2], l);
					++n_gapo, n_gap_all += l;
				} else if (m[2] == 'N') {
					tl += l;
				} else if (m[2] == 'S') {
					clip[M == 0? 0 : 1] = l, sclip += l;
				} else if (m[2] == 'H') {
					clip[M == 0? 0 : 1] = l;
				}
			}
			if (NM != null) {
				var tmp = NM - n_gap_all - nn;
				if (tmp < 0 && nn == 0) warn("WARNING: NM is smaller than the number of gaps at line " + lineno + ": NM=" + NM + ", nn=" + nn + ", G=" + n_gap_all);
				if (tmp < 0) tmp = 0;
				n_sub += tmp;
			}
			if (n_cigar > 65535) ++n_cigar_64k;
			if (ql + sclip != aqlen)
				warn("WARNING: aligned query length is inconsistent with CIGAR at line " + lineno + " (" + (ql+sclip) + " != " + aqlen + ")");
			if (atlen != null && atlen != tl)
				warn("WARNING: aligned reference length is inconsistent with CIGAR at line " + lineno);
			if (is_sam) {
				qs = clip[is_rev? 1 : 0], qe = qs + ql;
				ori_qlen = clip[0] + ql + clip[1];
			}
			if (count_err && NM != null) {
				var n_mm = NM - n_gap_all;
				if (n_mm < 0) warn("WARNING: NM is smaller than the number of gaps at line " + lineno);
				if (n_mm < 0) n_mm = 0;
				print(t[0], ori_qlen, t[11], ori_qlen - (qe - qs), NM, l_match + n_gap_all, n_mm + n_gapo, l_match + n_gapo);
			}
			regs.push([qs, qe]);
			last_qlen = ori_qlen;
		}
	}
	if (regs.length) {
		l_tot += last_qlen;
		l_cov += cov_len(regs);
	}

	file.close();
	buf.destroy();

	if (gap_out_len == null && !count_err) {
		print("Number of mapped sequences: " + n_seq);
		print("Number of primary alignments: " + n_pri);
		print("Number of secondary alignments: " + n_2nd);
		print("Number of primary alignments with >65535 CIGAR operations: " + n_cigar_64k);
		print("Number of bases in mapped sequences: " + l_tot);
		print("Number of mapped bases: " + l_cov);
		print("Number of substitutions: " + n_sub);
		print("Number of insertions in [0,50): " + n_gap[0][0]);
		print("Number of insertions in [50,100): " + n_gap[0][1]);
		print("Number of insertions in [100,300): " + n_gap[0][2]);
		print("Number of insertions in [300,400): " + n_gap[0][3]);
		print("Number of insertions in [400,1000): " + n_gap[0][4]);
		print("Number of insertions in [1000,inf): " + n_gap[0][5]);
		print("Number of deletions in [0,50): " + n_gap[1][0]);
		print("Number of deletions in [50,100): " + n_gap[1][1]);
		print("Number of deletions in [100,300): " + n_gap[1][2]);
		print("Number of deletions in [300,400): " + n_gap[1][3]);
		print("Number of deletions in [400,1000): " + n_gap[1][4]);
		print("Number of deletions in [1000,inf): " + n_gap[1][5]);
	}
}

function paf_bedcov(args)
{
	function read_bed(fn, to_merge, to_dedup)
	{
		var file = new File(fn);
		var buf = new Bytes();
		var h = {};
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (h[t[0]] == null)
				h[t[0]] = [];
			var bst = parseInt(t[1]);
			var ben = parseInt(t[2]);
			if (t.length >= 12 && /^\d+$/.test(t[9])) {
				t[9] = parseInt(t[9]);
				var sz = t[10].split(",");
				var st = t[11].split(",");
				for (var i = 0; i < t[9]; ++i) {
					st[i] = parseInt(st[i]);
					sz[i] = parseInt(sz[i]);
					h[t[0]].push([bst + st[i], bst + st[i] + sz[i], 0, 0, 0]);
				}
			} else {
				h[t[0]].push([bst, ben, 0, 0, 0]);
			}
		}
		buf.destroy();
		file.close();
		for (var chr in h) {
			if (to_merge) Interval.merge(h[chr], false);
			else if (to_dedup) Interval.dedup(h[chr], false);
			else Interval.sort(h[chr]);
			Interval.index_end(h[chr]);
		}
		return h;
	}

	var c, print_len = false, to_merge = true, to_dedup = false, fn_excl = null;
	while ((c = getopt(args, "pde:")) != null) {
		if (c == 'p') print_len = true;
		else if (c == 'd') to_dedup = true, to_merge = false;
		else if (c == 'e') fn_excl = getopt.arg;
	}

	if (args.length - getopt.ind < 2) {
		print("Usage: paftools.js bedcov [options] <regions.bed> <target.bed>");
		print("Options:");
		print("  -e FILE    exclude target regions (2nd file) overlapping BED FILE []");
		print("  -p         print number of covered bases for each target");
		exit(1);
	}

	var excl = fn_excl != null? read_bed(fn_excl, true, false) : null;
	var target = read_bed(args[getopt.ind], to_merge, to_dedup);

	var file, buf = new Bytes();
	var tot_len = 0, hit_len = 0;
	file = args[getopt.ind+1] != '-'? new File(args[getopt.ind+1]) : new File();
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var a = [];
		var bst = parseInt(t[1]);
		var ben = parseInt(t[2]);
		if (t.length >= 12 && /^\d+$/.test(t[9])) { // BED12
			t[9] = parseInt(t[9]);
			var sz = t[10].split(",");
			var st = t[11].split(",");
			for (var i = 0; i < t[9]; ++i) {
				st[i] = parseInt(st[i]);
				sz[i] = parseInt(sz[i]);
				a.push([bst + st[i], bst + st[i] + sz[i], false]);
			}
		} else a.push([bst, ben, false]); // 3-column BED
		var feat_len = 0;
		for (var i = 0; i < a.length; ++i) {
			if (excl != null && excl[t[0]] != null) {
				var oe = Interval.find_ovlp(excl[t[0]], a[i][0], a[i][1]);
				if (oe.length > 0)
					continue;
			}
			a[i][2] = true;
			feat_len += a[i][1] - a[i][0];
		}
		tot_len += feat_len;
		if (target[t[0]] == null) continue;
		var b = [];
		for (var i = 0; i < a.length; ++i) {
			if (!a[i][2]) continue;
			var o = Interval.find_ovlp(target[t[0]], a[i][0], a[i][1]);
			for (var j = 0; j < o.length; ++j) {
				var max_st = o[j][0] > a[i][0]? o[j][0] : a[i][0];
				var min_en = o[j][1] < a[i][1]? o[j][1] : a[i][1];
				b.push([max_st, min_en]);
				o[j][2] += min_en - max_st;
				++o[j][3];
				if (max_st == o[j][0] && min_en == o[j][1])
					++o[j][4];
			}
		}
		// find the length covered
		var feat_hit_len = 0;
		if (b.length > 0) {
			b.sort(function(a,b) {return a[0]-b[0]});
			var st = b[0][0], en = b[0][1];
			for (var i = 1; i < b.length; ++i) {
				if (b[i][0] <= en) en = en > b[i][1]? en : b[i][1];
				else feat_hit_len += en - st, st = b[i][0], en = b[i][1];
			}
			feat_hit_len += en - st;
		}
		hit_len += feat_hit_len;
		if (print_len) print('F', t.slice(0, 4).join("\t"), feat_len, feat_hit_len);
	}
	file.close();

	buf.destroy();

	warn("# target bases: " + tot_len);
	warn("# target bases overlapping regions: " + hit_len + ' (' + (100.0 * hit_len / tot_len).toFixed(2) + '%)');
}

function paf_vcfpair(args)
{
	var c, is_male = false, sample = 'syndip', hgver = null;
	var PAR = { '37':[[0, 2699520], [154931043, 155260560]] };
	while ((c = getopt(args, "ms:g:")) != null) {
		if (c == 'm') is_male = true;
		else if (c == 's') sample = getopt.arg;
		else if (c == 'g') hgver = getopt.arg;
	}
	if (is_male && (hgver == null || PAR[hgver] == null))
		throw("for a male, -g must be specified to properly handle PARs on chrX");

	if (getopt.ind == args.length) {
		print("Usage: paftools.js vcfpair [options] <in.pair.vcf>");
		print("Options:");
		print("  -m       the sample is male");
		print("  -g STR   human genome version '37' []");
		print("  -s STR   sample name [" + sample + "]");
		exit(1);
	}

	var re_ctg = is_male? /^(chr)?([0-9]+|X|Y)$/ : /^(chr)?([0-9]+|X)$/;
	var label = ['1', '2'];
	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if (line.charAt(0) == '#') {
			if (/^##(source|reference)=/.test(line)) continue;
			if ((m = /^##contig=.*ID=([^\s,]+)/.exec(line)) != null) {
				if (!re_ctg.test(m[1])) continue;
			} else if (/^#CHROM/.test(line)) {
				var t = line.split("\t");
				--t.length;
				t[t.length-1] = sample;
				line = t.join("\t");
				print('##FILTER=<ID=HET1,Description="Heterozygous in the first haplotype">');
				print('##FILTER=<ID=HET2,Description="Heterozygous in the second haplotype">');
				print('##FILTER=<ID=GAP1,Description="Uncalled in the first haplotype">');
				print('##FILTER=<ID=GAP2,Description="Uncalled in the second haplotype">');
			}
			print(line);
			continue;
		}
		var t = line.split("\t");
		if (!re_ctg.test(t[0])) continue;
		var GT = null, AD = null, FILTER = [], HT = [null, null];
		for (var i = 0; i < 2; ++i) {
			if ((m = /^(\.|[0-9]+)\/(\.|[0-9]+):(\S+)/.exec(t[9+i])) == null) {
				warn(line);
				throw Error("malformatted VCF");
			}
			var s = m[3].split(",");
			if (AD == null) {
				AD = [];
				for (var j = 0; j < s.length; ++j)
					AD[j] = 0;
			}
			for (var j = 0; j < s.length; ++j)
				AD[j] += parseInt(s[j]);
			if (m[1] == '.') {
				FILTER.push('GAP' + label[i]);
				HT[i] = '.';
			} else if (m[1] != m[2]) {
				FILTER.push('HET' + label[i]);
				HT[i] = '.';
			} else HT[i] = m[1];
		}
		--t.length;
		// test if this is in a haploid region
		var hap = 0, st = parseInt(t[1]), en = st + t[3].length;
		if (is_male) {
			if (/^(chr)?X/.test(t[0])) {
				if (hgver != null && PAR[hgver] != null) {
					var r = PAR[hgver], in_par = false;
					for (var i = 0; i < r.length; ++i)
						if (r[i][0] <= st && en <= r[i][1])
							in_par = true;
					hap = in_par? 0 : 2;
				}
			} else if (/^(chr)?Y/.test(t[0])) {
				hap = 1;
			}
		}
		// special treatment for haploid regions
		if (hap > 0 && FILTER.length == 1) {
			if ((hap == 2 && FILTER[0] == "GAP1") || (hap == 1 && FILTER[0] == "GAP2"))
				FILTER.length = 0;
		}
		// update VCF
		t[5] = 30; // fake QUAL
		t[6] = FILTER.length? FILTER.join(";") : ".";
		t[9] = HT.join("|") + ":" + AD.join(",");
		print(t.join("\t"));
	}
	file.close();
	buf.destroy();
}

/**********************
 * Conversion related *
 **********************/

function paf_view(args)
{
	var c, line_len = 80, fmt = "aln";
	while ((c = getopt(args, "f:l:")) != null) {
		if (c == 'f') {
			fmt = getopt.arg;
			if (fmt != "aln" && fmt != "lastz-cigar" && fmt != "maf")
				throw Error("format must be one of aln, lastz-cigar and maf");
		} else if (c == 'l') line_len = parseInt(getopt.arg);
	}
	if (line_len == 0) line_len = 0x7fffffff;

	if (getopt.ind == args.length) {
		print("Usage: paftools.js view [options] <in.paf>");
		print("Options:");
		print("  -f STR    output format: aln (BLAST-like), maf or lastz-cigar [aln]");
		print("  -l INT    line length in BLAST-like output [80]");
		exit(1);
	}

	function padding_str(x, len, right)
	{
		var s = x.toString();
		if (s.length < len) {
			if (right) s += Array(len - s.length + 1).join(" ");
			else s = Array(len - s.length + 1).join(" ") + s;
		}
		return s;
	}

	function update_aln(s_ref, s_qry, s_mid, type, seq, slen)
	{
		var l = type == '*'? 1 : seq.length;
		if (type == '=' || type == ':') {
			s_ref.set(seq);
			s_qry.set(seq);
			s_mid.set(Array(l+1).join("|"));
			slen[0] += l, slen[1] += l;
		} else if (type == '*') {
			s_ref.set(seq.charAt(0));
			s_qry.set(seq.charAt(1));
			s_mid.set(' ');
			slen[0] += 1, slen[1] += 1;
		} else if (type == '+') {
			s_ref.set(Array(l+1).join("-"));
			s_qry.set(seq);
			s_mid.set(Array(l+1).join(" "));
			slen[1] += l;
		} else if (type == '-') {
			s_ref.set(seq);
			s_qry.set(Array(l+1).join("-"));
			s_mid.set(Array(l+1).join(" "));
			slen[0] += l;
		}
	}

	function print_aln(rs, qs, strand, slen, elen, s_ref, s_qry, s_mid)
	{
		print(["Ref+:", padding_str(rs + slen[0] + 1, 10, false), s_ref.toString(), padding_str(rs + elen[0], 10, true)].join(" "));
		print("                 " + s_mid.toString());
		var st, en;
		if (strand == '+') st = qs + slen[1] + 1, en = qs + elen[1];
		else st = qs - slen[1], en = qs - elen[1] + 1;
		print(["Qry" + strand + ":", padding_str(st, 10, false), s_qry.toString(), padding_str(en, 10, true)].join(" "));
	}

	var s_ref = new Bytes(), s_qry = new Bytes(), s_mid = new Bytes(); // these are used to show padded alignment
	var re_cs = /([:=\-\+\*])(\d+|[A-Za-z]+)/g;
	var re_cg = /(\d+)([MIDNSHP=X])/g;

	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var lineno = 0;
	if (fmt == "maf") print("##maf version=1\n");
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		var t = line.split("\t", 12);
		++lineno;
		s_ref.length = s_qry.length = s_mid.length = 0;
		var slen = [0, 0], elen = [0, 0];
		if (fmt == "lastz-cigar") { // LASTZ-cigar output
			var cg = (m = /\tcg:Z:(\S+)/.exec(line)) != null? m[1] : null;
			if (cg == null) {
				warn("WARNING: converting to LASTZ-cigar format requires the 'cg' tag, which is absent on line " + lineno);
				continue;
			}
			var score = (m = /\tAS:i:(\d+)/.exec(line)) != null? m[1] : 0;
			var out = ['cigar:', t[0], t[2], t[3], t[4], t[5], t[7], t[8], '+', score];
			while ((m = re_cg.exec(cg)) != null)
				out.push(m[2], m[1]);
			print(out.join(" "));
		} else if (fmt == "maf") { // MAF output
			var cs = (m = /\tcs:Z:(\S+)/.exec(line)) != null? m[1] : null;
			if (cs == null) {
				warn("WARNING: converting to MAF requires the 'cs' tag, which is absent on line " + lineno);
				continue;
			}
			while ((m = re_cs.exec(cs)) != null) {
				if (m[1] == ':')
					throw Error("converting to MAF only works with 'minimap2 --cs=long'");
				update_aln(s_ref, s_qry, s_mid, m[1], m[2], elen);
			}
			var score = (m = /\tAS:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : 0;
			var len = t[0].length > t[5].length? t[0].length : t[5].length;
			print("a " + score);
			print(["s", padding_str(t[5], len, true), padding_str(t[7], 10, false), padding_str(parseInt(t[8]) - parseInt(t[7]), 10, false),
					"+", padding_str(t[6], 10, false), s_ref.toString()].join(" "));
			var qs, qe, ql = parseInt(t[1]);
			if (t[4] == '+') {
				qs = parseInt(t[2]);
				qe = parseInt(t[3]);
			} else {
				qs = ql - parseInt(t[3]);
				qe = ql - parseInt(t[2]);
			}
			print(["s", padding_str(t[0], len, true), padding_str(qs, 10, false), padding_str(qe - qs, 10, false),
					t[4], padding_str(ql, 10, false), s_qry.toString()].join(" "));
			print("");
		} else { // BLAST-like output
			var cs = (m = /\tcs:Z:(\S+)/.exec(line)) != null? m[1] : null;
			if (cs == null) {
				warn("WARNING: converting to BLAST-like alignment requires the 'cs' tag, which is absent on line " + lineno);
				continue;
			}
			var n_mm = 0, n_oi = 0, n_od = 0, n_ei = 0, n_ed = 0;
			while ((m = re_cs.exec(cs)) != null) {
				if (m[1] == '*') ++n_mm;
				else if (m[1] == '+') ++n_oi, n_ei += m[2].length;
				else if (m[1] == '-') ++n_od, n_ed += m[2].length;
			}
			line = line.replace(/\tc[sg]:Z:\S+/g, ""); // get rid of cs or cg tags
			print('>' + line + "\tmm:i:"+n_mm + "\toi:i:"+n_oi + "\tei:i:"+n_ei + "\tod:i:"+n_od + "\ted:i:"+n_ed);
			var rs = parseInt(t[7]), qs = t[4] == '+'? parseInt(t[2]) : parseInt(t[3]);
			var n_blocks = 0;
			while ((m = re_cs.exec(cs)) != null) {
				if (m[1] == ':') m[2] = Array(parseInt(m[2]) + 1).join("=");
				var start = 0, rest = m[1] == '*'? 1 : m[2].length;
				while (rest > 0) {
					var l_proc;
					if (s_ref.length + rest >= line_len) {
						l_proc = line_len - s_ref.length;
						update_aln(s_ref, s_qry, s_mid, m[1], m[1] == '*'? m[2] : m[2].substr(start, l_proc), elen);
						if (n_blocks > 0) print("");
						print_aln(rs, qs, t[4], slen, elen, s_ref, s_qry, s_mid);
						++n_blocks;
						s_ref.length = s_qry.length = s_mid.length = 0;
						slen[0] = elen[0], slen[1] = elen[1];
					} else {
						l_proc = rest;
						update_aln(s_ref, s_qry, s_mid, m[1], m[1] == '*'? m[2] : m[2].substr(start, l_proc), elen);
					}
					rest -= l_proc, start += l_proc;
				}
			}
			if (s_ref.length > 0) {
				if (n_blocks > 0) print("");
				print_aln(rs, qs, t[4], slen, elen, s_ref, s_qry, s_mid);
				++n_blocks;
			}
			print("//");
		}
	}
	file.close();
	buf.destroy();

	s_ref.destroy(); s_qry.destroy(); s_mid.destroy();
}

function paf_gff2bed(args)
{
	var c, fn_ucsc_fai = null, is_short = false, keep_gff = false, print_junc = false, output_gene = false, ens_canon_only = false;
	while ((c = getopt(args, "u:sgjGe")) != null) {
		if (c == 'u') fn_ucsc_fai = getopt.arg;
		else if (c == 's') is_short = true;
		else if (c == 'g') keep_gff = true;
		else if (c == 'j') print_junc = true;
		else if (c == 'G') output_gene = true;
		else if (c == 'e') ens_canon_only = true;
	}

	if (getopt.ind == args.length) {
		print("Usage: paftools.js gff2bed [options] <in.gff>");
		print("Options:");
		print("  -j       output junction BED");
		print("  -s       print names in the short form");
		print("  -u FILE  hg38.fa.fai for chr name conversion");
		print("  -e       only show transcript tagged with 'Ensembl_canonical'");
		print("  -g       output GFF (used with -u)");
		exit(1);
	}

	var ens2ucsc = {};
	if (fn_ucsc_fai != null) {
		var buf = new Bytes();
		var file = new File(fn_ucsc_fai);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var s = t[0];
			if (/_(random|alt|decoy)$/.test(s)) {
				s = s.replace(/_(random|alt|decoy)$/, '');
				s = s.replace(/^chr\S+_/, '');
			} else {
				s = s.replace(/^chrUn_/, '');
			}
			s = s.replace(/v(\d+)/, ".$1");
			if (s != t[0]) ens2ucsc[s] = t[0];
		}
		file.close();
		buf.destroy();
	}

	var colors = {
		'protein_coding':'0,128,255',
		'mRNA':'0,128,255',
		'lincRNA':'0,192,0',
		'snRNA':'0,192,0',
		'miRNA':'0,192,0',
		'misc_RNA':'0,192,0'
	};

	function print_bed12(exons, cds_st, cds_en, is_short, print_junc)
	{
		if (exons.length == 0) return;
		var name = is_short? exons[0][7] + "|" + exons[0][5] : exons[0].slice(4, 7).join("|");
		var a = exons.sort(function(a,b) {return a[1]-b[1]});
		if (print_junc) {
		for (var i = 1; i < a.length; ++i)
			print(a[i][0], a[i-1][2], a[i][1], name, 1000, a[i][3]);
			return;
		}
		var sizes = [], starts = [], st, en;
		st = a[0][1];
		en = a[a.length - 1][2];
		if (cds_st == 1<<30) cds_st = st;
		if (cds_en == 0) cds_en = en;
		if (cds_st < st || cds_en > en)
			throw Error("inconsistent thick start or end for transcript " + a[0][4]);
		for (var i = 0; i < a.length; ++i) {
			sizes.push(a[i][2] - a[i][1]);
			starts.push(a[i][1] - st);
		}
		var color = colors[a[0][5]];
		if (color == null) color = '196,196,196';
		print(a[0][0], st, en, name, 1000, a[0][3], cds_st, cds_en, color, a.length, sizes.join(",") + ",", starts.join(",") + ",");
	}

	var re_gtf  = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name|tag) "([^"]+)";/g;
	var re_gff3 = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name)=([^;]+)/g;
	var re_gtf_gene  = /\b(gene_id|gene_type|gene_name) "([^;]+)";/g;
	var re_gff3_gene = /\b(gene_id|gene_type|source_gene|gene_biotype|gene_name)=([^;]+);/g;
	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);

	var exons = [], cds_st = 1<<30, cds_en = 0, last_id = null;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (keep_gff) {
			if (t[0].charAt(0) != '#' && ens2ucsc[t[0]] != null)
				t[0] = ens2ucsc[t[0]];
			print(t.join("\t"));
			continue;
		}
		if (t[0].charAt(0) == '#') continue;
		if (output_gene) {
			var id = null, src = null, biotype = null, type = "", name = "N/A";
			if (t[2] != "gene") continue;
			while ((m = re_gtf_gene.exec(t[8])) != null) {
				if (m[1] == "gene_id") id = m[2];
				else if (m[1] == "gene_type") type = m[2];
				else if (m[1] == "gene_name") name = m[2];
			}
			while ((m = re_gff3_gene.exec(t[8])) != null) {
				if (m[1] == "gene_id") id = m[2];
				else if (m[1] == "source_gene") src = m[2];
				else if (m[1] == "gene_type") type = m[2];
				else if (m[1] == "gene_biotype") biotype = m[2];
				else if (m[1] == "gene_name") name = m[2];
			}
			if (src != null) id = src;
			if (type == "" && biotype != null) type = biotype;
			print(t[0], parseInt(t[3]) - 1, t[4], [id, type, name].join("|"), 1000, t[6]);
			continue;
		}
		if (t[2] != "CDS" && t[2] != "exon") continue;
		t[3] = parseInt(t[3]) - 1;
		t[4] = parseInt(t[4]);
		var id = null, type = "", name = "N/A", biotype = "", m, tname = "N/A", ens_canonical = false;
		while ((m = re_gtf.exec(t[8])) != null) {
			if (m[1] == "transcript_id") id = m[2];
			else if (m[1] == "transcript_type") type = m[2];
			else if (m[1] == "transcript_biotype" || m[1] == "gbkey") biotype = m[2];
			else if (m[1] == "gene_name" || m[1] == "gene_id") name = m[2];
			else if (m[1] == "transcript_name") tname = m[2];
			else if (m[1] == "tag" && m[2] == "Ensembl_canonical") ens_canonical = true;
		}
		while ((m = re_gff3.exec(t[8])) != null) {
			if (m[1] == "transcript_id") id = m[2];
			else if (m[1] == "transcript_type") type = m[2];
			else if (m[1] == "transcript_biotype" || m[1] == "gbkey") biotype = m[2];
			else if (m[1] == "gene_name" || m[1] == "gene_id") name = m[2];
			else if (m[1] == "transcript_name") tname = m[2];
		}
		if (ens_canon_only && !ens_canonical) continue;
		if (type == "" && biotype != "") type = biotype;
		if (id == null) throw Error("No transcript_id");
		if (id != last_id) {
			print_bed12(exons, cds_st, cds_en, is_short, print_junc);
			exons = [], cds_st = 1<<30, cds_en = 0;
			last_id = id;
		}
		if (t[2] == "CDS") {
			cds_st = cds_st < t[3]? cds_st : t[3];
			cds_en = cds_en > t[4]? cds_en : t[4];
		} else if (t[2] == "exon") {
			if (fn_ucsc_fai != null) {
				if (ens2ucsc[t[0]] != null)
					t[0] = ens2ucsc[t[0]];
				else if (/^[A-Z]+\d+\.\d+$/.test(t[0]))
					t[0] = t[0].replace(/([A-Z]+\d+)\.(\d+)/, "chrUn_$1v$2");
			}
			exons.push([t[0], t[3], t[4], t[6], id, type, name, tname]);
		}
	}
	if (last_id != null)
		print_bed12(exons, cds_st, cds_en, is_short, print_junc);

	file.close();
	buf.destroy();
}

function paf_sam2paf(args)
{
	var c, pri_only = false, long_cs = false;
	while ((c = getopt(args, "pL")) != null) {
		if (c == 'p') pri_only = true;
		else if (c == 'L') long_cs = true;
	}
	if (args.length == getopt.ind) {
		print("Usage: paftools.js sam2paf [options] <in.sam>");
		print("Options:");
		print("  -p      convert primary or supplementary alignments only");
		print("  -L      output the cs tag in the long form");
		exit(1);
	}

	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var buf = new Bytes();
	var re = /(\d+)([MIDSHNX=])/g, re_MD = /(\d+)|(\^[A-Za-z]+)|([A-Za-z])/g, re_tag = /\t(\S\S:[AZif]):(\S+)/g;

	var ctg_len = {}, lineno = 0;
	while (file.readline(buf) >= 0) {
		var m, n_cigar = 0, line = buf.toString();
		++lineno;
		if (line.charAt(0) == '@') {
			if (/^@SQ/.test(line)) {
				var name = (m = /\tSN:(\S+)/.exec(line)) != null? m[1] : null;
				var l = (m = /\tLN:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;
				if (name != null && l != null) ctg_len[name] = l;
			}
			continue;
		}
		var t = line.split("\t", 11);
		var flag = parseInt(t[1]);
		if (t[9] != '*' && t[10] != '*' && t[9].length != t[10].length)
			throw Error("at line " + lineno + ": inconsistent SEQ and QUAL lengths - " + t[9].length + " != " + t[10].length);
		if (t[2] == '*' || (flag&4) || t[5] == '*') continue;
		if (pri_only && (flag&0x100)) continue;
		var tlen = ctg_len[t[2]];
		if (tlen == null) throw Error("at line " + lineno + ": can't find the length of contig " + t[2]);
		// find tags
		var nn = 0, NM = null, MD = null, cs_str = null, md_list = [];
		while ((m = re_tag.exec(line)) != null) {
			if (m[1] == "NM:i") NM = parseInt(m[2]);
			else if (m[1] == "nn:i") nn = parseInt(m[2]);
			else if (m[1] == "MD:Z") MD = m[2];
			else if (m[1] == "cs:Z") cs_str = m[2];
		}
		if (t[9] == '*') MD = cs_str = null;
		// infer various lengths from CIGAR
		var clip = [0, 0], soft_clip = 0, I = [0, 0], D = [0, 0], M = 0, N = 0, mm = 0, have_M = false, have_ext = false, cigar = [];
		while ((m = re.exec(t[5])) != null) {
			var l = parseInt(m[1]), op = m[2];
			if (op == 'M') M += l, have_M = true;
			else if (op == 'I') ++I[0], I[1] += l;
			else if (op == 'D') ++D[0], D[1] += l;
			else if (op == 'N') N += l;
			else if (op == 'S') clip[n_cigar == 0? 0 : 1] = l, soft_clip += l;
			else if (op == 'H') clip[n_cigar == 0? 0 : 1] = l;
			else if (op == '=') M += l, have_ext = true, op = 'M';
			else if (op == 'X') M += l, mm += l, have_ext = true, op = 'M';
			++n_cigar;
			if (MD != null && op != 'H') {
				if (cigar.length > 0 && cigar[cigar.length-1][1] == op)
					cigar[cigar.length-1][0] += l;
				else cigar.push([l, op]);
			}
		}
		var ql = M + I[1] + soft_clip;
		var tl = M + D[1] + N;
		var ts = parseInt(t[3]) - 1, te = ts + tl;
		// checking coordinate and length consistencies
		if (n_cigar > 65535)
			warn("WARNING at line " + lineno + ": " + n_cigar + " CIGAR operations");
		if (te > tlen) {
			warn("WARNING at line " + lineno + ": alignment end position larger than ref length; skipped");
			continue;
		}
		if (t[9] != '*' && t[9].length != ql) {
			warn("WARNING at line " + lineno + ": SEQ length inconsistent with CIGAR (" + t[9].length + " != " + ql + "); skipped");
			continue;
		}
		// parse MD
		var cs = [];
		if (MD != null && cs_str == null && t[9] != "*") {
			var k = 0, cx = 0, cy = 0, mx = 0, my = 0; // cx: cigar ref position; cy: cigar query; mx: MD ref; my: MD query
			while ((m = re_MD.exec(MD)) != null) {
				if (m[2] != null) { // deletion from the reference
					var len = m[2].length - 1;
					cs.push('-', m[2].substr(1));
					mx += len, cx += len, ++k;
				} else { // copy or mismatch
					var ml = m[1] != null? parseInt(m[1]) : 1;
					while (k < cigar.length && cigar[k][1] != 'D') {
						var cl = cigar[k][0], op = cigar[k][1];
						if (op == 'M') {
							if (my + ml < cy + cl) {
								if (ml > 0) {
									if (m[3] != null) cs.push('*', m[3], t[9][my]);
									else if (long_cs) cs.push('=', t[9].substr(my, ml));
									else cs.push(':', ml);
								}
								mx += ml, my += ml, ml = 0;
								break;
							} else {
								var dl = cy + cl - my;
								if (long_cs) cs.push('=', t[9].substr(my, dl));
								else cs.push(':', dl);
								cx += cl, cy += cl, ++k;
								mx += dl, my += dl, ml -= dl;
							}
						} else if (op == 'I') {
							cs.push('+', t[9].substr(cy, cl));
							cy += cl, my += cl, ++k;
						} else if (op == 'S') {
							cy += cl, my += cl, ++k;
						} else throw Error("at line " + lineno + ": inconsistent MD tag");
					}
					if (ml != 0) throw Error("at line " + lineno + ": inconsistent MD tag");
				}
			}
			if (cx != mx || cy != my) throw Error("at line " + lineno + ": inconsistent MD tag");
		}
		// compute matching length, block length and calibrate NM
		if (have_ext && !have_M) { // extended CIGAR
			if (NM != null && NM != I[1] + D[1] + mm)
				warn("WARNING at line " + lineno + ": NM is different from sum of gaps and mismatches");
			NM = I[1] + D[1] + mm;
		} else if (NM != null) { // standard CIGAR; NM present
			if (NM < I[1] + D[1]) {
				warn("WARNING at line " + lineno + ": NM is less than the total number of gaps (" + NM + " < " + (I[1]+D[1]) + ")");
				NM = I[1] + D[1];
			}
			mm = NM - (I[1] + D[1]);
		} else { // no way to compute mm
			warn("WARNING at line " + lineno + ": unable to find the number of mismatches; assuming zero");
			mm = 0;
		}
		var mlen = M - mm;
		var blen = M + I[1] + D[1];
		// find query name, start and end
		var qlen = M + I[1] + clip[0] + clip[1];
		var qname = t[0], qs, qe;
		if ((flag&1) && (flag&0x40)) qname += '/1';
		if ((flag&1) && (flag&0x80)) qname += '/2';
		if (flag&16) qs = clip[1], qe = qlen - clip[0];
		else qs = clip[0], qe = qlen - clip[1];
		// optional tags
		var type = flag&0x100? 'S' : 'P';
		var tags = ["tp:A:" + type];
		if (NM != null) tags.push("mm:i:"+mm);
		tags.push("gn:i:"+(I[1]+D[1]), "go:i:"+(I[0]+D[0]), "cg:Z:" + t[5].replace(/\d+[SH]/g, ''));
		if (cs_str != null) tags.push("cs:Z:" + cs_str);
		else if (cs.length > 0) tags.push("cs:Z:" + cs.join(""));
		// print out
		var a = [qname, qlen, qs, qe, flag&16? '-' : '+', t[2], tlen, ts, te, mlen, blen, t[4]];
		print(a.join("\t"), tags.join("\t"));
	}

	buf.destroy();
	file.close();
}

function paf_delta2paf(args)
{
	if (args.length == 0) {
		print("Usage: paftools.js delta2paf <in.delta>");
		exit(1);
	}

	var buf = new Bytes();
	var file = args[0] == '-'? new File() : new File(args[0]);

	var rname, qname, rlen, qlen, qs, qe, rs, re, strand, NM, cigar, x, y, seen_gt = false;
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if ((m = /^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/.exec(line)) != null) {
			rname = m[1], qname = m[2], rlen = parseInt(m[3]), qlen = parseInt(m[4]);
			seen_gt = true;
			continue;
		}
		if (!seen_gt) continue;
		var t = line.split(" ");
		if (t.length == 7) {
			for (var i = 0; i < 5; ++i)
				t[i] = parseInt(t[i]);
			strand = ((t[0] < t[1] && t[2] < t[3]) || (t[0] > t[1] && t[2] > t[3]))? 1 : -1;
			rs = (t[0] < t[1]? t[0] : t[1]) - 1;
			re = t[1] > t[0]? t[1] : t[0];
			qs = (t[2] < t[3]? t[2] : t[3]) - 1;
			qe = t[3] > t[2]? t[3] : t[2];
			x = y = 0;
			NM = parseInt(t[4]);
			cigar = [];
		} else if (t.length == 1) {
			var d = parseInt(t[0]);
			if (d == 0) {
				var blen = 0, cigar_str = [];
				if (re - rs - x != qe - qs - y) throw Error("inconsisnt alignment");
				cigar.push((re - rs - x) << 4);
				for (var i = 0; i < cigar.length; ++i) {
					blen += cigar[i] >> 4;
					cigar_str.push((cigar[i]>>4) + "MID".charAt(cigar[i]&0xf));
				}
				print([qname, qlen, qs, qe, strand > 0? '+' : '-', rname, rlen, rs, re, blen - NM, blen, 0, "NM:i:" + NM, "cg:Z:" + cigar_str.join("")].join("\t"));
			} else if (d > 0) {
				var l = d - 1;
				x += l + 1, y += l;
				if (l > 0) cigar.push(l<<4);
				if (cigar.length > 0 && (cigar[cigar.length-1]&0xf) == 2)
					cigar[cigar.length-1] += 1<<4;
				else cigar.push(1<<4|2); // deletion
			} else {
				var l = -d - 1;
				x += l, y += l + 1;
				if (l > 0) cigar.push(l<<4);
				if (cigar.length > 0 && (cigar[cigar.length-1]&0xf) == 1)
					cigar[cigar.length-1] += 1<<4;
				else cigar.push(1<<4|1); // insertion
			}
		}
	}
	file.close();
	buf.destroy();
}

function paf_splice2bed(args)
{
	var colors = ["0,128,255", "255,0,0", "0,192,0"];

	function print_lines(a, fmt, keep_multi)
	{
		if (a.length == 0) return;
		if (fmt == "bed") {
			var n_pri = 0;
			for (var i = 0; i < a.length; ++i)
				if (a[i][8] == 0) ++n_pri;
			if (n_pri > 1) {
				for (var i = 0; i < a.length; ++i)
					if (a[i][8] == 0) a[i][8] = 1;
			} else if (n_pri == 0) {
				warn("Warning: " + a[0][3] + " doesn't have a primary alignment");
			}
			for (var i = 0; i < a.length; ++i) {
				if (!keep_multi && a[i][8] == 2) continue;
				a[i][8] = colors[a[i][8]];
				print(a[i].join("\t"));
			}
		}
		a.length = 0;
	}

	var re = /(\d+)([MIDNSHP=X])/g;
	var c, fmt = "bed", fn_name_conv = null, keep_multi = false;
	while ((c = getopt(args, "f:n:m")) != null) {
		if (c == 'f') fmt = getopt.arg;
		else if (c == 'n') fn_name_conv = getopt.arg;
		else if (c == 'm') keep_multi = true;
	}
	if (getopt.ind == args.length) {
		print("Usage: paftools.js splice2bed [options] <in.paf>|<in.sam>");
		print("Options:");
		print("  -m      keep multiple mappings (SAM flag 0x100)");
		exit(1);
	}

	var conv = null;
	if (fn_name_conv != null) {
		conv = new Map();
		var file = new File(fn_name_conv);
		var buf = new Bytes();
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			conv.put(t[0], t[1]);
		}
		buf.destroy();
		file.close();
	}

	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var buf = new Bytes();
	var a = [];
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '@') continue; // skip SAM header lines
		var t = line.split("\t");
		var is_pri = false, cigar = null, a1;
		var qname = conv != null? conv.get(t[0]) : null;
		if (qname != null) t[0] = qname;
		if (t.length >= 10 && t[4] != '+' && t[4] != '-' && /^\d+/.test(t[1])) { // SAM
			var flag = parseInt(t[1]);
			if (flag&1) t[0] += '/' + (flag>>6&3);
		}
		if (a.length && a[0][3] != t[0]) {
			print_lines(a, fmt, keep_multi);
			a = [];
		}
		if (t.length >= 12 && (t[4] == '+' || t[4] == '-')) { // PAF
			for (var i = 12; i < t.length; ++i) {
				if (t[i].substr(0, 5) == 'cg:Z:') {
					cigar = t[i].substr(5);
				} else if (t[i].substr(0, 5) == 's2:i:') {
					is_pri = true;
				}
			}
			a1 = [t[5], t[7], t[8], t[0], Math.floor(t[9]/t[10]*1000), t[4]];
		} else if (t.length >= 10) { // SAM
			var flag = parseInt(t[1]);
			if ((flag&4) || a[2] == '*') continue;
			cigar = t[5];
			is_pri = (flag&0x100)? false : true;
			a1 = [t[2], parseInt(t[3])-1, null, t[0], 1000, (flag&16)? '-' : '+'];
		} else {
			throw Error("unrecognized input format");
		}
		if (cigar == null) throw Error("missing CIGAR");
		var m, x0 = 0, x = 0, bs = [], bl = [];
		while ((m = re.exec(cigar)) != null) {
			if (m[2] == 'M' || m[2] == 'D') {
				x += parseInt(m[1]);
			} else if (m[2] == 'N') {
				bs.push(x0);
				bl.push(x - x0);
				x += parseInt(m[1]);
				x0 = x;
			}
		}
		bs.push(x0);
		bl.push(x - x0);
		// write the BED12 line
		if (a1[2] == null) a1[2] = a1[1] + x;
		a1.push(a1[1], a1[2]); // thick start/end is the same as start/end
		a1.push(is_pri? 0 : 2, bs.length, bl.join(",")+",", bs.join(",")+",");
		a.push(a1);
	}
	print_lines(a, fmt, keep_multi);
	buf.destroy();
	file.close();
	if (conv != null) conv.destroy();
}

/**********************
 * Evaluation related *
 **********************/

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
		warn("  -Q INT     print wrong mappings with mapQ>=INT [don't print]");
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

function paf_junceval(args)
{
	var c, l_fuzzy = 0, print_ovlp = false, print_err_only = false, first_only = false, chr_only = false, aa = false, is_bed = false;
	while ((c = getopt(args, "l:epcab1")) != null) {
		if (c == 'l') l_fuzzy = parseInt(getopt.arg);
		else if (c == 'e') print_err_only = print_ovlp = true;
		else if (c == 'p') print_ovlp = true;
		else if (c == 'c') chr_only = true;
		else if (c == 'a') aa = true;
		else if (c == 'b') is_bed = true;
		else if (c == '1') first_only = true;
	}

	if (args.length - getopt.ind < 1) {
		print("Usage: paftools.js junceval [options] <gene.gtf> <aln.sam>");
		print("Options:");
		print("  -l INT    tolerance of junction positions (0 for exact) [0]");
		print("  -p        print overlapping introns");
		print("  -e        print erroreous overlapping introns");
		print("  -c        only consider alignments to /^(chr)?([0-9]+|X|Y)$/");
		print("  -a        miniprot PAF as input");
		print("  -b        BED as input");
		print("  -1        only process the first alignment of each query");
		exit(1);
	}

	var file, buf = new Bytes();

	var tr = {};
	file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, t = buf.toString().split("\t");
		if (t[0].charAt(0) == '#') continue;
		if (t[2] != 'exon') continue;
		var st = parseInt(t[3]) - 1;
		var en = parseInt(t[4]);
		if ((m = /transcript_id "(\S+)"/.exec(t[8])) == null) continue;
		var tid = m[1];
		if (tr[tid] == null) tr[tid] = [t[0], t[6], 0, 0, []];
		tr[tid][4].push([st, en]);
	}
	file.close();

	var anno = {};
	for (var tid in tr) {
		var t = tr[tid];
		Interval.sort(t[4]);
		t[2] = t[4][0][0];
		t[3] = t[4][t[4].length - 1][1];
		if (anno[t[0]] == null) anno[t[0]] = [];
		var s = t[4];
		for (var i = 0; i < s.length - 1; ++i) {
			if (s[i][1] >= s[i+1][0])
				warn("WARNING: incorrect annotation for transcript "+tid+" ("+s[i][1]+" >= "+s[i+1][0]+")")
					anno[t[0]].push([s[i][1], s[i+1][0]]);
		}
	}
	tr = null;

	for (var chr in anno) {
		var e = anno[chr];
		if (e.length == 0) continue;
		Interval.sort(e);
		var k = 0;
		for (var i = 1; i < e.length; ++i) // dedup
			if (e[i][0] != e[k][0] || e[i][1] != e[k][1])
				e[++k] = e[i].slice(0);
		e.length = k + 1;
		Interval.index_end(e);
	}

	var n_pri = 0, n_unmapped = 0, n_mapped = 0;
	var n_sgl = 0, n_splice = 0, n_splice_hit = 0, n_splice_novel = 0;

	file = getopt.ind+1 >= args.length || args[getopt.ind+1] == '-'? new File() : new File(args[getopt.ind+1]);
	var last_qname = null;
	var re_cigar = /(\d+)([MIDNSHP=XFGUV])/g;
	while (file.readline(buf) >= 0) {
		var m, t = buf.toString().split("\t");
		var ctg_name = null, cigar = null, pos = null, qname;

		if (t[0].charAt(0) == '@') continue;
		if (t[0] == "##PAF") t.shift();
		qname = t[0];
		if (is_bed) {
			ctg_name = t[0], pos = parseInt(t[1]), cigar == null;
		} else if (t[4] == '+' || t[4] == '-' || t[4] == '*') { // PAF
			ctg_name = t[5], pos = parseInt(t[7]);
			var type = 'P';
			for (i = 12; i < t.length; ++i) {
				if ((m = /^(tp:A|cg:Z):(\S+)/.exec(t[i])) != null) {
					if (m[1] == 'tp:A') type = m[2];
					else cigar = m[2];
				}
			}
			if (type == 'S') continue; // secondary
		} else { // SAM
			ctg_name = t[2], pos = parseInt(t[3]) - 1, cigar = t[5];
			var flag = parseInt(t[1]);
			if (flag&0x100) continue; // secondary
		}

		if (chr_only && !/^(chr)?([0-9]+|X|Y)$/.test(ctg_name)) continue;
		if (first_only && last_qname == qname) continue;
		if (ctg_name == '*') { // unmapped
			++n_unmapped;
			continue;
		} else {
			++n_pri;
			if (last_qname != qname) {
				++n_mapped;
				last_qname = qname;
			}
		}

		var intron = [];
		if (is_bed) {
			intron.push([pos, parseInt(t[2])]);
		} else if (aa) {
			var tmp_junc = [], tmp = 0;
			while ((m = re_cigar.exec(cigar)) != null) {
				var len = parseInt(m[1]), op = m[2];
				if (op == 'N') {
					tmp_junc.push([tmp, tmp + len]);
					tmp += len;
				} else if (op == 'U') {
					tmp_junc.push([tmp + 1, tmp + len - 2]);
					tmp += len;
				} else if (op == 'V') {
					tmp_junc.push([tmp + 2, tmp + len - 1]);
					tmp += len;
				} else if (op == 'M' || op == 'X' || op == '=' || op == 'D') {
					tmp += len * 3;
				} else if (op == 'F' || op == 'G') {
					tmp += len;
				}
			}
			if (t[4] == '+') {
				for (var i = 0; i < tmp_junc.length; ++i)
					intron.push([pos + tmp_junc[i][0], pos + tmp_junc[i][1]]);
			} else if (t[4] == '-') {
				var glen = parseInt(t[8]) - parseInt(t[7]);
				for (var i = tmp_junc.length - 1; i >= 0; --i)
					intron.push([pos + (glen - tmp_junc[i][1]), pos + (glen - tmp_junc[i][0])]);
			}
		} else {
			while ((m = re_cigar.exec(cigar)) != null) {
				var len = parseInt(m[1]), op = m[2];
				if (op == 'N') {
					intron.push([pos, pos + len]);
					pos += len;
				} else if (op == 'M' || op == 'X' || op == '=' || op == 'D') pos += len;
			}
		}
		if (intron.length == 0) {
			++n_sgl;
			continue;
		}
		n_splice += intron.length;

		var chr = anno[ctg_name];
		if (chr != null) {
			for (var i = 0; i < intron.length; ++i) {
				var o = Interval.find_ovlp(chr, intron[i][0], intron[i][1]);
				if (o.length > 0) {
					var hit = false;
					for (var j = 0; j < o.length; ++j) {
						var st_diff = intron[i][0] - o[j][0];
						var en_diff = intron[i][1] - o[j][1];
						if (st_diff < 0) st_diff = -st_diff;
						if (en_diff < 0) en_diff = -en_diff;
						if (st_diff <= l_fuzzy && en_diff <= l_fuzzy)
							++n_splice_hit, hit = true;
						if (hit) break;
					}
					if (print_ovlp) {
						var type = hit? 'C' : 'P';
						if (hit && print_err_only) continue;
						var x = '[';
						for (var j = 0; j < o.length; ++j) {
							if (j) x += ', ';
							x += '(' + o[j][0] + "," + o[j][1] + ')';
						}
						x += ']';
						print(type, qname, i+1, ctg_name, intron[i][0], intron[i][1], x);
					}
				} else {
					++n_splice_novel;
					if (print_ovlp)
						print('N',  qname, i+1, ctg_name, intron[i][0], intron[i][1]);
				}
			}
		} else {
			n_splice_novel += intron.length;
		}
	}
	file.close();

	buf.destroy();

	if (!print_ovlp) {
		print("# unmapped reads: " + n_unmapped);
		print("# mapped reads: " + n_mapped);
		print("# primary alignments: " + n_pri);
		print("# singletons: " + n_sgl);
		print("# predicted introns: " + n_splice);
		print("# non-overlapping introns: " + n_splice_novel);
		print("# correct introns: " + n_splice_hit + " (" + (n_splice_hit / n_splice * 100).toFixed(2) + "%)");
	}
}

function paf_exoneval(args) // adapted from paf_junceval()
{
	var c, l_fuzzy = 0, print_ovlp = false, print_err_only = false, first_only = false, chr_only = false, aa = false, is_bed = false, use_cds = false, eval_base = false;
	while ((c = getopt(args, "l:epcab1ds")) != null) {
		if (c == 'l') l_fuzzy = parseInt(getopt.arg);
		else if (c == 'e') print_err_only = print_ovlp = true;
		else if (c == 'p') print_ovlp = true;
		else if (c == 'c') chr_only = true;
		else if (c == 'a') aa = true, use_cds = true;
		else if (c == 'b') is_bed = true;
		else if (c == '1') first_only = true;
		else if (c == 'd') use_cds = true;
		else if (c == 's') eval_base = true;
	}

	if (args.length - getopt.ind < 1) {
		print("Usage: paftools.js exoneval [options] <gene.gtf> <aln.sam>");
		print("Options:");
		print("  -l INT    tolerance of junction positions (0 for exact) [0]");
		print("  -d        evaluate coding regions only (exon regions by default)");
		print("  -a        miniprot PAF as input (force -d)");
		print("  -p        print overlapping exons");
		print("  -e        print erroreous overlapping exons");
		print("  -c        only consider alignments to /^(chr)?([0-9]+|X|Y)$/");
		print("  -1        only process the first alignment of each query");
		print("  -b        BED as input");
		print("  -s        compute base Sn and Sp (more memory)");
		exit(1);
	}

	var file, buf = new Bytes();

	warn("Reading reference GTF...");
	var tr = {};
	file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, t = buf.toString().split("\t");
		if (t[0].charAt(0) == '#') continue;
		if (use_cds) {
			if (t[2] != "cds" && t[2] != "CDS") continue;
		} else {
			if (t[2] != 'exon') continue;
		}
		var st = parseInt(t[3]) - 1;
		var en = parseInt(t[4]);
		if ((m = /transcript_id "(\S+)"/.exec(t[8])) == null) continue;
		var tid = m[1];
		if (tr[tid] == null) tr[tid] = [t[0], t[6], 0, 0, []];
		tr[tid][4].push([st, en]); // this keeps transcript
	}
	file.close();

	var anno = {};
	for (var tid in tr) { // traverse each transcript
		var t = tr[tid];
		Interval.sort(t[4]);
		t[2] = t[4][0][0];
		t[3] = t[4][t[4].length - 1][1];
		if (anno[t[0]] == null) anno[t[0]] = [];
		var s = t[4];
		for (var i = 0; i < s.length; ++i) // traverse each exon
			anno[t[0]].push([s[i][0], s[i][1]]);
	}
	tr = null;

	for (var chr in anno) { // index exons
		var e = anno[chr];
		if (e.length == 0) continue;
		Interval.sort(e);
		var k = 0;
		for (var i = 1; i < e.length; ++i) // dedup
			if (e[i][0] != e[k][0] || e[i][1] != e[k][1])
				e[++k] = e[i].slice(0);
		e.length = k + 1;
		Interval.index_end(e);
	}

	var n_pri = 0, n_unmapped = 0, n_mapped = 0;
	var n_exon = 0, n_exon_hit = 0, n_exon_novel = 0;

	file = getopt.ind+1 >= args.length || args[getopt.ind+1] == '-'? new File() : new File(args[getopt.ind+1]);
	var last_qname = null, qexon = {};
	var re_cigar = /(\d+)([MIDNSHP=XFGUV])/g;

	warn("Evaluating alignments...");
	while (file.readline(buf) >= 0) {
		var m, t = buf.toString().split("\t");
		var ctg_name = null, cigar = null, pos = null, qname;

		if (t[0].charAt(0) == '@') continue;
		if (t[0] == "##PAF") t.shift();
		qname = t[0];
		if (is_bed) {
			ctg_name = t[0], pos = parseInt(t[1]), cigar == null;
		} else if (t[4] == '+' || t[4] == '-' || t[4] == '*') { // PAF
			ctg_name = t[5], pos = parseInt(t[7]);
			var type = 'P';
			for (i = 12; i < t.length; ++i) {
				if ((m = /^(tp:A|cg:Z):(\S+)/.exec(t[i])) != null) {
					if (m[1] == 'tp:A') type = m[2];
					else cigar = m[2];
				}
			}
			if (type == 'S') continue; // secondary
		} else { // SAM
			ctg_name = t[2], pos = parseInt(t[3]) - 1, cigar = t[5];
			var flag = parseInt(t[1]);
			if (flag&0x100) continue; // secondary
		}

		if (chr_only && !/^(chr)?([0-9]+|X|Y)$/.test(ctg_name)) continue;
		if (first_only && last_qname == qname) continue;
		if (ctg_name == '*') { // unmapped
			++n_unmapped;
			continue;
		} else {
			++n_pri;
			if (last_qname != qname) {
				++n_mapped;
				last_qname = qname;
			}
		}

		var exon = [];
		if (is_bed) { // BED
			exon.push([pos, parseInt(t[2])]);
		} else if (aa) {
			var tmp_exon = [], tmp = 0, tmp_st = 0;
			while ((m = re_cigar.exec(cigar)) != null) {
				var len = parseInt(m[1]), op = m[2];
				if (op == 'N') {
					tmp_exon.push([tmp_st, tmp]);
					tmp_st = tmp + len, tmp += len;
				} else if (op == 'U') {
					tmp_exon.push([tmp_st, tmp + 1]);
					tmp_st = tmp + len - 2, tmp += len;
				} else if (op == 'V') {
					tmp_exon.push([tmp_st, tmp + 2]);
					tmp_st = tmp + len - 1, tmp += len;
				} else if (op == 'M' || op == 'X' || op == '=' || op == 'D') {
					tmp += len * 3;
				} else if (op == 'F' || op == 'G') {
					tmp += len;
				}
			}
			tmp_exon.push([tmp_st, tmp]);
			if (t[4] == '+') {
				for (var i = 0; i < tmp_exon.length; ++i)
					exon.push([pos + tmp_exon[i][0], pos + tmp_exon[i][1]]);
			} else if (t[4] == '-') { // For protein-to-genome alignment, the coordinates are on the query strand. Need to flip them.
				var glen = parseInt(t[8]) - parseInt(t[7]);
				for (var i = tmp_exon.length - 1; i >= 0; --i)
					exon.push([pos + (glen - tmp_exon[i][1]), pos + (glen - tmp_exon[i][0])]);
			}
		} else {
			var tmp_st = pos;
			while ((m = re_cigar.exec(cigar)) != null) {
				var len = parseInt(m[1]), op = m[2];
				if (op == 'N') {
					exon.push([tmp_st, pos]);
					tmp_st = pos + len, pos += len;
				} else if (op == 'M' || op == 'X' || op == '=' || op == 'D') pos += len;
			}
			exon.push([tmp_st, pos]);
		}
		n_exon += exon.length;

		var chr = anno[ctg_name];
		if (chr != null) {
			for (var i = 0; i < exon.length; ++i) {
				if (eval_base) {
					if (qexon[ctg_name] == null) qexon[ctg_name] = [];
					qexon[ctg_name].push([exon[i][0], exon[i][1]]);
				}
				var o = Interval.find_ovlp(chr, exon[i][0], exon[i][1]);
				if (o.length > 0) {
					var hit = false;
					for (var j = 0; j < o.length; ++j) {
						var st_diff = exon[i][0] - o[j][0];
						var en_diff = exon[i][1] - o[j][1];
						if (st_diff < 0) st_diff = -st_diff;
						if (en_diff < 0) en_diff = -en_diff;
						if (st_diff <= l_fuzzy && en_diff <= l_fuzzy)
							++n_exon_hit, hit = true;
						if (hit) break;
					}
					if (print_ovlp) {
						var type = hit? 'C' : 'P';
						if (hit && print_err_only) continue;
						var x = '[';
						for (var j = 0; j < o.length; ++j) {
							if (j) x += ', ';
							x += '(' + o[j][0] + "," + o[j][1] + ')';
						}
						x += ']';
						print(type, qname, i+1, ctg_name, exon[i][0], exon[i][1], x);
					}
				} else {
					++n_exon_novel;
					if (print_ovlp)
						print('N',  qname, i+1, ctg_name, exon[i][0], exon[i][1]);
				}
			}
		} else {
			n_exon_novel += exon.length;
		}
	}
	file.close();

	buf.destroy();

	if (!print_ovlp) {
		print("# unmapped reads: " + n_unmapped);
		print("# mapped reads: " + n_mapped);
		print("# primary alignments: " + n_pri);
		print("# predicted exons: " + n_exon);
		print("# non-overlapping exons: " + n_exon_novel);
		print("# correct exons: " + n_exon_hit + " (" + (n_exon_hit / n_exon * 100).toFixed(2) + "%)");
	}

	function merge_and_index(ex) {
		for (var chr in ex) {
			var a = [];
			e = ex[chr];
			Interval.sort(e);
			var st = e[0][0], en = e[0][1];
			for (var i = 1; i < e.length; ++i) { // merge
				if (e[i][0] > en) {
					a.push([st, en]);
					st = e[i][0], en = e[i][1];
				} else {
					en = en > e[i][1]? en : e[i][1];
				}
			}
			a.push([st, en]);
			Interval.index_end(a);
			ex[chr] = a;
		}
	}

	function cal_sn(a0, a1) {
		var tot = 0, cov = 0;
		for (var chr in a1) {
			var e0 = a0[chr], e1 = a1[chr];
			for (var i = 0; i < e1.length; ++i)
				tot += e1[i][1] - e1[i][0];
			if (e0 == null) continue;
			for (var i = 0; i < e1.length; ++i) {
				var o = Interval.find_ovlp(e0, e1[i][0], e1[i][1]);
				for (var j = 0; j < o.length; ++j) { // this only works when there are no overlaps between intervals
					var st = e1[i][0] > o[j][0]? e1[i][0] : o[j][0];
					var en = e1[i][1] < o[j][1]? e1[i][1] : o[j][1];
					cov += en - st;
				}
			}
		}
		return [tot, cov];
	}

	if (eval_base) {
		warn("Computing base Sn and Sp...");
		merge_and_index(qexon);
		merge_and_index(anno);
		var sn = cal_sn(qexon, anno);
		var sp = cal_sn(anno, qexon);
		print("Base Sn: " + sn[1] + " / " + sn[0] + " = " + (sn[1] / sn[0] * 100).toFixed(2) + "%");
		print("Base Sp: " + sp[1] + " / " + sp[0] + " = " + (sp[1] / sp[0] * 100).toFixed(2) + "%");
	}
}

// evaluate overlap sensitivity
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

function paf_vcfstat(args)
{
	var c, ts = { "AG":1, "GA":1, "CT":1, "TC":1 };
	while ((c = getopt(args, "")) != null) {
	}
	var buf = new Bytes();
	var file = args.length == getopt.ind? new File() : new File(args[getopt.ind]);
	var x = { sub:0, ts:0, tv:0, ins:0, del:0, ins1:0, del1:0, ins2:0, del2:0, ins50:0, del50:0, ins1k:0, del1k:0, ins7k:0, del7k:0, insinf:0, delinf:0 };
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0][0] == '#') continue;
		var alt = t[4].split(",");
		var ref = t[3];
		for (var i = 0; i < alt.length; ++i) {
			var a = alt[i];
			if (a[0] == '<' || a[1] == '>') continue;
			var l = ref.length < a.length? ref.length : a.length;
			for (var j = 0; j < l; ++j) {
				if (ref[j] != a[j]) {
					++x.sub;
					if (ts[ref[j] + a[j]]) ++x.ts;
					else ++x.tv;
				}
			}
			var d = a.length - ref.length;
			if (d > 0) {
				++x.ins;
				if (d == 1) ++x.ins1;
				else if (d == 2) ++x.ins2;
				else if (d < 50) ++x.ins50;
				else if (d < 1000) ++x.ins1k;
				else if (d < 7000) ++x.ins7k;
				else ++x.insinf;
			} else if (d < 0) {
				d = -d;
				++x.del;
				if (d == 1) ++x.del1;
				else if (d == 2) ++x.del2;
				else if (d < 50) ++x.del50;
				else if (d < 1000) ++x.del1k;
				else if (d < 7000) ++x.del7k;
				else ++x.delinf;
			}
		}
	}
	file.close();
	buf.destroy();
	print("# substitutions: " + x.sub);
	print("ts/tv: " + (x.ts / x.tv).toFixed(3));
	print("# insertions: " + x.ins);
	print("# 1bp insertions: " + x.ins1);
	print("# 2bp insertions: " + x.ins2);
	print("# [3,50) insertions: " + x.ins50);
	print("# [50,1000) insertions: " + x.ins1k);
	print("# [1000,7000) insertions: " + x.ins7k);
	print("# >=7000 insertions: " + x.insinf);
	print("# deletions: " + x.del);
	print("# 1bp deletions: " + x.del1);
	print("# 2bp deletions: " + x.del2);
	print("# [3,50) deletions: " + x.del50);
	print("# [50,1000) deletions: " + x.del1k);
	print("# [1000,7000) deletions: " + x.del7k);
	print("# >=7000 deletions: " + x.delinf);
}

function paf_parseNum(s) {
	var m, x = null;
	if ((m = /^(\d*\.?\d*)([mMgGkK]?)/.exec(s)) != null) {
		x = parseFloat(m[1]);
		if (m[2] == 'k' || m[2] == 'K') x *= 1000;
		else if (m[2] == 'm' || m[2] == 'M') x *= 1000000;
		else if (m[2] == 'g' || m[2] == 'G') x *= 1000000000;
	}
	return Math.floor(x + .499);
}

function paf_misjoin(args)
{
	var c, min_seg_len = 1000000, max_gap = 1000000, fn_cen = null, show_long = false, show_err = false, cen_ratio = 0.5;
	var n_diff = [0, 0], n_gap = [0, 0], n_inv = [0, 0], n_inv_end = [0, 0];
	while ((c = getopt(args, "l:g:c:per:")) != null) {
		if (c == 'l') min_seg_len = paf_parseNum(getopt.arg);
		else if (c == 'g') max_gap = paf_parseNum(getopt.arg);
		else if (c == 'c') fn_cen = getopt.arg;
		else if (c == 'r') cen_ratio = parseFloat(getopt.arg);
		else if (c == 'p') show_long = true;
		else if (c == 'e') show_err = true;
	}
	if (args.length == getopt.ind) {
		print("Usage: paftools.js misjoin [options] <in.paf>");
		print("Options:");
		print("  -c FILE   BED for centromeres []");
		print("  -r FLOAT  count a centromeric event if overlap ratio > FLOAT [" + cen_ratio + "]");
		print("  -l NUM    min alignment block length [1m]");
		print("  -g NUM    max gap size [1m]");
		print("  -e        output misjoins not involving centromeres");
		print("  -p        output long alignment blocks for debugging");
		return;
	}
	var cen = {};
	var file, buf = new Bytes();
	if (fn_cen != null) {
		file = new File(fn_cen);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (cen[t[0]] == null) cen[t[0]] = [];
			cen[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
		}
		file.close();
	}

	function test_cen(cen, chr, st, en) {
		var b = cen[chr], len = 0;
		if (b == null) return false;
		for (var j = 0; j < b.length; ++j)
			if (b[j][0] < en && b[j][1] > st) {
				var s = b[j][0] > st? b[j][0] : st;
				var e = b[j][1] < en? b[j][1] : en;
				len += e - s;
			}
		return len < (en - st) * cen_ratio? false : true;
	}

	function test_cen_point(cen, chr, x) {
		var b = cen[chr];
		if (b == null) return false;
		for (var j = 0; j < b.length; ++j)
			if (x >= b[j][0] && x < b[j][1])
				return true;
		return false;
	}

	if (show_err || show_long) {
		print("C\tJ  inter-chromosomal misjoin");
		print("C\tj  inter-chromosomal misjoin with both breakpoints ending in centromeres");
		print("C\tG  long gap on the reference genome");
		print("C\tg  long gap on the reference genome with both breakpoints ending in centromeres");
		print("C\tM  closed inversion");
		print("C");
	}
	function process(a) {
		var k = 0;
		for (var i = 0; i < a.length; ++i) {
			for (var j = 1; j <= 3;  ++j) a[i][j] = parseInt(a[i][j]);
			for (var j = 6; j <= 11; ++j) a[i][j] = parseInt(a[i][j]);
			if (a[i][10] >= min_seg_len) a[k++] = a[i];
		}
		a.length = k;
		if (a.length == 1) return;
		a = a.sort(function(x,y){return x[2]-y[2]});
		if (show_long) for (var i = 0; i < a.length; ++i) print(a[i].join("\t"));
		for (var i = 1; i < a.length; ++i) {
			var ov = [false, false], end_cen = [false, false];
			ov[0] = test_cen(cen, a[i-1][5], a[i-1][7], a[i-1][8]);
			ov[1] = test_cen(cen, a[i][5], a[i][7], a[i][8]);
			end_cen[0] = test_cen_point(cen, a[i-1][5], a[i-1][4] == '+'? a[i-1][8] : a[i-1][7]);
			end_cen[1] = test_cen_point(cen, a[i][5],   a[i][4] == '+'?   a[i][7]   : a[i][8]);
			if (a[i-1][5] != a[i][5]) { // different chr
				if (ov[0] || ov[1]) ++n_diff[1];
				else if (show_err) {
					var label = end_cen[0] && end_cen[1]? 'j' : 'J';
					print(label, a[i-1].slice(0, 12).join("\t"));
					print(label, a[i].slice(0, 12).join("\t"));
				}
				++n_diff[0];
			} else if (a[i-1][4] == a[i][4]) { // a gap
				var dq = a[i][2] - a[i-1][3];
				var dr = a[i][4] == '+'? a[i][7] - a[i-1][8] : a[i-1][7] - a[i][8];
				var gap = dr > dq? dr - dq : dq - dr;
				if (gap > max_gap) {
					if (ov[0] || ov[1]) ++n_gap[1];
					else if (show_err) {
						var label = end_cen[0] && end_cen[1]? 'g' : 'G';
						print(label, a[i-1].slice(0, 12).join("\t"));
						print(label, a[i].slice(0, 12).join("\t"));
					}
					++n_gap[0];
				}
			} else if (i + 1 < a.length && a[i+1][4] == a[i-1][4]) { // bracketed inversion
				if (ov[0] || ov[1]) ++n_inv[1];
				else if (show_err) {
					print("M", a[i-1].slice(0, 12).join("\t"));
					print("M", a[i].slice(0, 12).join("\t"));
					print("M", a[i+1].slice(0, 12).join("\t"));
				}
				++n_inv[0];
				++i;
			} else { // hanging inversion
				if (ov[0] || ov[1]) ++n_inv_end[1];
				++n_inv_end[0];
			}
		}
	}

	file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var a = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (a.length > 0 && a[0][0] != t[0]) {
			process(a);
			a.length = 0;
		}
		a.push(t);
	}
	if (a.length > 0) process(a);
	file.close();
	buf.destroy();
	print("# inter-chromosomal misjoins: " + n_diff.join(","));
	print("# intra-chromosomal gaps: " + n_gap.join(","));
	print("# candidate inversions in the middle: " + n_inv.join(","));
	print("# candidate inversions at contig ends: " + n_inv_end.join(","));
}

function _paf_get_alen(t)
{
	var svlen = null, alen = null;
	if ((m = /(^|;)SVLEN=(-?\d+)/.exec(t[7])) != null)
		svlen = parseInt(m[2]);
	var s = t[4].split(",");
	var min_abs_diff = 1<<30, max_abs_diff = 0;
	if (svlen != null && svlen != 0)
		alen = svlen, min_abs_diff = max_abs_diff = svlen > 0? svlen : -svlen;
	var rlen = t[3].length;
	for (var i = 0; i < s.length; ++i) {
		if (/^<\S+>$/.test(s[i])) continue;
		var diff = s[i].length - rlen;
		var abs_diff = diff > 0? diff : -diff;
		min_abs_diff = min_abs_diff < abs_diff? min_abs_diff : abs_diff;
		if (max_abs_diff < abs_diff)
			max_abs_diff = abs_diff, alen = diff;
	}
	return [alen, min_abs_diff, max_abs_diff];
}

function paf_sveval(args)
{
	var c, min_flt = 30, min_size = 50, max_size = 100000, win_size = 500, print_err = false, print_match = false, bed_fn = null;
	var len_diff_ratio = 0.5;
	while ((c = getopt(args, "f:i:x:w:er:pd:")) != null) {
		if (c == 'f') min_flt = paf_parseNum(getopt.arg);
		else if (c == 'i') min_size = paf_parseNum(getopt.arg);
		else if (c == 'x') max_size = paf_parseNum(getopt.arg);
		else if (c == 'w') win_size = paf_parseNum(getopt.arg);
		else if (c == 'd') len_diff_ratio = parseFloat(getopt.arg);
		else if (c == 'r') bed_fn = getopt.arg;
		else if (c == 'e') print_err = true;
		else if (c == 'p') print_match = true;
	}
	if (args.length - getopt.ind < 2) {
		print("Usage: paftools.js sveval [options] <base.vcf> <call.vcf>");
		print("Options:");
		print("  -r FILE    confident region in BED []");
		print("  -f INT     min length to discard [" + min_flt + "]");
		print("  -i INT     min SV length [" + min_size + "]");
		print("  -x INT     max SV length [" + max_size + "]");
		print("  -w INT     fuzzy windown size [" + win_size + "]");
		print("  -d FLOAT   max allele diff if there is a single allele in the window [" + len_diff_ratio + "]");
		print("  -e         print errors");
		return;
	}

	function read_bed(fn) {
		var buf = new Bytes();
		var file = new File(fn);
		var bed = {};
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (bed[t[0]] == null) bed[t[0]] = [];
			bed[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
		}
		file.close();
		buf.destroy();
		for (var x in bed) {
			Interval.sort(bed[x]);
			Interval.merge(bed[x]);
			Interval.index_end(bed[x]);
		}
		return bed;
	}

	var bed = bed_fn != null? read_bed(bed_fn) : null;

	function read_vcf(fn, bed) {
		var buf = new Bytes();
		var file = new File(fn);
		var v = {};
		while (file.readline(buf) >= 0) {
			var m, t = buf.toString().split("\t");
			if (t[0][0] == '#') continue;
			if (bed != null && bed[t[0]] == null) continue;
			if (t[4] == '<INV>' || t[4] == '<INVDUP>') continue; // no inversion
			if (/[\[\]]/.test(t[4])) continue; // no break points
			var st = parseInt(t[1]) - 1, en = st + t[3].length;
			// parse svlen
			var b = _paf_get_alen(t), svlen = b[0];
			var abslen = svlen == null? 0 : svlen > 0? svlen : -svlen;
			if (abslen < min_flt || abslen > max_size) continue;
			// update end
			if ((m = /(^|;)END=(\d+)/.exec(t[7])) != null)
				en = parseInt(m[2]);
			else if (svlen != null && svlen < 0)
				en = st + (-svlen);
			if (en < st) en = st;
			if (st == en) --st, ++en;
			if (bed != null && Interval.find_ovlp(bed[t[0]], st, en).length == 0) continue;
			// insert
			if (v[t[0]] == null) v[t[0]] = [];
			v[t[0]].push([st, en, svlen, abslen]);
		}
		file.close();
		buf.destroy();
		for (var x in v) {
			Interval.sort(v[x]);
			Interval.index_end(v[x]);
		}
		return v;
	}

	function compare_vcf(v0, v1, label) {
		var m = 0, n = 0;
		for (var x in v1) {
			var a1 = v1[x], a0 = v0[x];
			for (var i = 0; i < a1.length; ++i) {
				if (a1[i][3] < min_size) continue;
				++n;
				if (a0 == null) continue;
				var ws = win_size + (a1[i][3]>>1);
				var st = a1[i][0] > ws? a1[i][0] - ws : 0;
				b = Interval.find_ovlp(a0, st, a1[i][1] + ws);
				var n_ins = 0, n_del = 0, sv_del = null, sv_ins = null;
				for (var j = 0; j < b.length; ++j) {
					if (b[j][2] < 0) ++n_del, sv_del = -b[j][2];
					else if (b[j][2] > 0) ++n_ins, sv_ins = b[j][2];
					if (print_match)
						print("MA", x, a1[i].slice(0, 3).join("\t"), b[j].slice(0, 3).join("\t"));
				}
				var match = false;
				if (a1[i][2] > 0) { // insertion
					if (n_ins == 1) {
						var diff = sv_ins - a1[i][3];
						if (diff < 0) diff = -diff;
						if (diff < min_size || diff / a1[i][3] < len_diff_ratio)
							match = true;
					} else if (n_ins > 1) match = true; // multiple insertions; ambiguous
				} else if (a1[i][2] < 0) {
					if (n_del == 1) { // deletion
						var diff = sv_del - a1[i][3];
						if (diff < 0) diff = -diff;
						if (diff < min_size || diff / a1[i][3] < len_diff_ratio)
							match = true;
					} else if (n_del > 1) match = true; // multiple deletions; ambiguous
				}
				if (match) ++m;
				else if (print_err) {
					if ((a1[i][2] > 0 && n_ins > 0) || (a1[i][2] < 0 && n_del > 0))
						print("MM", x, a1[i].slice(0, 3).join("\t"));
					print(label, x, a1[i].slice(0, 3).join("\t"));
				}
			}
		}
		return [n, m];
	}

	var v_base = read_vcf(args[getopt.ind+0], bed);
	var v_call = read_vcf(args[getopt.ind+1], bed);
	var fn = compare_vcf(v_call, v_base, 'FN');
	var fp = compare_vcf(v_base, v_call, 'FP');
	print('SN', fn[0], fn[1], (fn[1] / fn[0]).toFixed(6));
	print('PC', fp[0], fp[1], (fp[1] / fp[0]).toFixed(6));
	print('F1', ((fn[1] / fn[0] + fp[1] / fp[0]) / 2).toFixed(6));
}

function paf_vcfsel(args)
{
	var c, min_l = 0, max_l = 1<<30;
	while ((c = getopt(args, "l:L:")) != null) {
		if (c == 'l') min_l = parseInt(getopt.arg);
		else if (c == 'L') max_l = parseInt(getopt.arg);
	}

	var buf = new Bytes();
	if (getopt.ind == args.length) {
		print("Usage: paftools.js vcfsel [options] <in.vcf>");
		return 1;
	}
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if (line[0] == '#') {
			print(line);
			continue;
		}
		var t = line.split("\t");
		var st = parseInt(t[1]), en = st + t[3].length - 1;
		if ((m = /(^|;)END=(\d+)/.exec(t[7])) != null)
			en = parseInt(m[2]);
		if (en < st) {
			warn("END is smaller than POS: " + en + " < " + st);
			en = st;
		}
		var b = _paf_get_alen(t);
		var alen = b[0], min_abs_diff = b[1], max_abs_diff = b[2];
		if (max_abs_diff < min_l || min_abs_diff > max_l)
			continue;
		print(line);
	}
	file.close();
	buf.destroy();
}

function paf_pafcmp(args)
{
	var c, opt = { min_len:5000, min_mapq:10, min_ovlp:0.5 };
	while ((c = getopt(args, "q:")) != null) {
		if (c == 'q') opt.min_mapq = parseInt(getopt.arg);
	}

	var buf = new Bytes();
	if (args.length - getopt.ind < 2) {
		print("Usage: paftools.js pafcmp [options] <base.paf> <test.paf>");
		print("Options:");
		print("  -q INT    min mapping quality [" + opt.min_mapq + "]");
		return 1;
	}

	var eval = { n_base:0, n_test:0, n_out_high:0, n_out_low:0, n_hit:0, n_wrong:0, n_miss:0 };

	function process_base(base, a) {
		if (a.length != 1) return;
		for (var i = 1; i < 4; ++i)
			a[0][i] = parseInt(a[0][i]);
		for (var i = 6; i < 12; ++i)
			a[0][i] = parseInt(a[0][i]);
		if (a[0][1] < opt.min_len) return;
		if (a[0][11] >= opt.min_mapq) ++eval.n_base;
		base[a[0][0]] = [a[0][5], a[0][7], a[0][8], a[0][11], 0, 0];
	}

	var file = new File(args[getopt.ind]);
	warn("Reading " + args[getopt.ind] + "...");
	var a = [], base = {};
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t");
		if (/\ttp:A:S/.test(line)) continue;
		if (a.length > 0 && a[0][0] != t[0]) {
			process_base(base, a);
			a = [];
		}
		a.push(t);
	}
	process_base(base, a);
	file.close();

	function process_test(base, a) {
		for (var i = 1; i < 4; ++i)
			a[0][i] = parseInt(a[0][i]);
		for (var i = 6; i < 12; ++i)
			a[0][i] = parseInt(a[0][i]);
		if (a[0][1] < opt.min_len) return;
		if (a[0][11] >= opt.min_mapq) ++eval.n_test;
		var c = [a[0][5], a[0][7], a[0][8], a[0][11]];
		if (base[a[0][0]] == null) {
			if (c[3] >= opt.min_mapq) ++opt.n_out_high;
			else ++opt.n_out_low;
		} else {
			var b = base[a[0][0]];
			var inter = 0, union = (b[2] - b[1]) + (c[2] - c[1]);
			if (b[0] == c[0]) { // same chr
				if (b[1] < c[1]) {
					if (b[2] > c[1])
						inter = b[2] - c[1], union = c[2] - b[1];
				} else { // c[1] < b[1]
					if (c[2] > b[1])
						inter = c[2] - b[1], union = b[2] - c[1];
				}
			}
			if (inter >= union * opt.min_ovlp) {
				if (b[3] >= opt.min_mapq) ++eval.n_hit;
				++b[4];
			} else {
				if (b[3] >= opt.min_mapq) {
					print("W", a[0][0], b.slice(0, 4).join("\t"), c.join("\t"));
					++eval.n_wrong;
				}
				++b[5];
			}
		}
	}

	file = new File(args[getopt.ind+1]);
	warn("Reading " + args[getopt.ind+1] + "...");
	a = [];
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t");
		if (/\ttp:A:S/.test(line)) continue;
		if (a.length > 0 && a[0][0] != t[0]) {
			process_test(base, a);
			a = [];
		}
		a.push(t);
	}
	process_test(base, a);
	file.close();

	for (var r in base) {
		var b = base[r];
		if (b[3] >= opt.min_mapq && b[4] == 0 && b[5] == 0) {
			++eval.n_miss;
			print("M", r, b.slice(0, 4).join("\t"));
		}
	}

	print("X", eval.n_base     + " base alignments with mapQ>=" + opt.min_mapq);
//	print("X", eval.n_test     + " test alignments with mapQ>=" + opt.min_mapq);
	print("X", eval.n_hit      + " base alignments correctly mapped by test");
	print("X", eval.n_wrong    + " wrong test alignment");
	print("X", eval.n_miss     + " base alignments missing");
	print("X", eval.n_out_high + " additional test alignments with mapQ>=" + opt.min_mapq);

	buf.destroy();
}

function paf_longcs2seq(args) {
	var c, opt = { query:false };
	while ((c = getopt(args, "q")) != null)
		if (c == 'q') opt.query = true;
	if (args.length == getopt.ind) {
		print("Usage: paftools.js longcs2seq [-q] <long-cs.paf>");
		return;
	}
	var re_cs = /([:=*+-])(\d+|[A-Za-z]+)/g
	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, cs = null, t = buf.toString().split("\t");
		for (var i = 12; i < t.length; ++i)
			if ((m = /^cs:Z:(\S+)/.exec(t[i])) != null) {
				cs = m[1];
				break;
			}
		if (cs == null) continue;
		var ts = "", qs = "";
		while ((m = re_cs.exec(cs)) != null) {
			if (m[1] == "=") ts += m[2], qs += m[2];
			else if (m[1] == "+") qs += m[2].toUpperCase();
			else if (m[1] == "-") ts += m[2].toUpperCase();
			else if (m[1] == "*") ts += m[2][0].toUpperCase(), qs += m[2][1].toUpperCase();
			else if (m[1] == ":") throw Error("Long cs is required");
		}
		if (opt.query) {
			print(">" + t[0] + "_" + t[2] + "_" + t[3]);
			print(qs);
		} else {
			print(">" + t[5] + "_" + t[7] + "_" + t[8]);
			print(ts);
		}
	}
	file.close();
	buf.destroy();
}

function paf_paf2gff(args) {
	var c, opt = { aa:false };
	var re_cigar = /(\d+)([A-Z=])/g;
	while ((c = getopt(args, "a")) != null) {
		if (c == 'a') opt.aa = true;
	}
	if (args.length == getopt.ind) {
		print("Usage: paftools.js paf2gff [-a] <in.paf>");
		return;
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var hid = 1, last_name = null;
	while (file.readline(buf) >= 0) {
		var m, t = buf.toString().split("\t");
		if (t[5] == '*') continue; // skip unmapped lines

		if (t[0] != last_name) last_name = t[0], hid = 1;
		else ++hid;
		for (var i = 1; i <= 3; ++i) t[i] = parseInt(t[i]);
		for (var i = 6; i <= 11; ++i) t[i] = parseInt(t[i]);
		var cigar = null, score = null, np = null, dist_stop = null, dist_start = null;
		for (var i = 12; i < t.length; ++i) {
			if ((m = /^(cg:Z|AS:i|np:i|da:i|do:i):(\S+)/.exec(t[i])) != null) {
				if (m[1] == 'cg:Z') cigar = m[2];
				else if (m[1] == 'AS:i') score = parseInt(m[2]);
				else if (m[1] == 'np:i') np = parseInt(m[2]);
				else if (m[1] == 'do:i') dist_stop = parseInt(m[2]);
				else if (m[1] == 'da:i') dist_start = parseInt(m[2]);
			}
		}
		if (cigar == null) throw Error("failed to find the cg:Z tag");
		if (score == null) throw Error("failed to find the AS:i tag");

		var st = 0, en = 0, phase = 0, pseudo = false, fs = 0, a = [];
		if (dist_start != null && dist_start == 0)
			a.push([t[5], 'paf2gff', 'start_codon', 0, 3, 0, t[4], '.', 0]);
		while ((m = re_cigar.exec(cigar)) != null) {
			var len = parseInt(m[1]);
			if (m[2] == 'M' || m[2] == 'D') {
				en += opt.aa? len * 3 : len;
			} else if (m[2] == 'F' || m[2] == 'G' || m[2] == 'R') {
				en += len, pseudo = true, fs = 1;
			} else if (m[2] == 'N') {
				a.push([t[5], 'paf2gff', 'exon', st, en, 0, t[4], phase, fs]);
				st = en + len, en += len, phase = 0, fs = 0;
			} else if (m[2] == 'U') { // ...xGT...AGxx...
				a.push([t[5], 'paf2gff', 'exon', st, en + 1, 0, t[4], phase, fs]);
				st = en + len - 2, en += len, phase = 2, fs = 0;
			} else if (m[2] == 'V') { // ...xxGT...AGx...
				a.push([t[5], 'paf2gff', 'exon', st, en + 2, 0, t[4], phase, fs]);
				st = en + len - 1, en += len, phase = 1, fs = 0;
			}
		}
		a.push([t[5], 'paf2gff', 'exon', st, en, 0, t[4], phase, fs]);
		if (en != t[8] - t[7]) throw Error("inconsistent cigar");
		if (dist_stop != null && dist_stop == 0)
			a.push([t[5], 'paf2gff', 'stop_codon', en, en + 3, 0, t[4], '.', 0]);
		var type = pseudo? 'pseudogene' : 'protein_coding';
		var attr = ['transcript_id=' + t[0] + '#' + hid, 'transcript_type=' + type].join(";");
		var trans_attr = 'identity=' + (t[9] / t[10]).toFixed(4);
		if (np != null) trans_attr += ';positive=' + (np * 3 / t[10]).toFixed(4);
		trans_attr += ';aa_start=' + t[2];
		trans_attr += ';aa_end=' + (t[1] - t[3]);
		if (dist_start != null && dist_start >= 0) trans_attr += ';dist_start_codon=' + dist_start;
		if (dist_stop != null && dist_stop >= 0) trans_attr += ';dist_stop_codon=' + dist_stop;
		var trans_st = t[7], trans_en = t[8];
		if (dist_stop != null && dist_stop == 0) {
			if (t[4] == '-') trans_st -= 3;
			else trans_en += 3;
		}
		print([t[5], 'paf2gff', 'transcript', trans_st + 1, trans_en, score, t[4], '.', attr + ';' + trans_attr].join("\t"));
		if (opt.aa && t[4] == '-') {
			var b = [], len = t[8] - t[7];
			for (var i = a.length - 1; i >= 0; --i) {
				var x = len - a[i][3];
				a[i][3] = len - a[i][4];
				a[i][4] = x;
				//a[i][7] = a[i][7] == 0? 0 : 3 - a[i][7]; // not sure if this line is needed
				b.push(a[i]);
			}
			a = b;
		}
		for (var i = 0; i < a.length; ++i) {
			if (!pseudo && a[i][2] == "exon") a[i][2] = "CDS";
			a[i][3] += t[7] + 1;
			a[i][4] += t[7];
			a[i][8] = attr + ";frameshift=" + a[i][8];
			print(a[i].join("\t"));
		}
	}
	file.close();
	buf.destroy();
}

function paf_gff2junc(args) {
	var c, feat = "CDS";
	while ((c = getopt(args, "f:")) != null) {
		if (c == 'f') feat = getopt.arg;
	}
	if (getopt.ind == args.length) {
		print("Usage: paftools.js gff2junc [-f feature] <in.gff3>");
		return;
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);

	function process_a(a) {
		if (a.length < 2) return;
		a = a.sort(function(x, y) { return x[4] - y[4] });
		for (var i = 1; i < a.length; ++i)
			print([a[i][1], a[i-1][5], a[i][4], a[i][0], 0, a[i][7]].join("\t"));
	}

	var a = [];
	while (file.readline(buf) >= 0) {
		var m, t = buf.toString().split("\t");
		if (t[0][0] == '#') continue;
		if (t[2].toLowerCase() != feat.toLowerCase()) continue;
		//print(t.join("\t"));
		if ((m = /\bParent=([^;]+)/.exec(t[8])) == null) {
			warn("Can't find Parent");
			continue;
		}
		t[3] = parseInt(t[3]) - 1;
		t[4] = parseInt(t[4]);
		t.unshift(m[1]);
		if (a.length > 0 && a[0][0] != m[1]) {
			process_a(a);
			a.length = 0;
			a.push(t);
		} else a.push(t);
	}
	process_a(a);
	file.close();
	buf.destroy();
}

/*************************
 ***** main function *****
 *************************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: paftools.js <command> [arguments]");
		print("Commands:");
		print("  view       convert PAF to BLAST-like (for eyeballing) or MAF");
		print("  splice2bed convert spliced alignment in PAF/SAM to BED12");
		print("  sam2paf    convert SAM to PAF");
		print("  delta2paf  convert MUMmer's delta to PAF");
		print("  gff2bed    convert GTF/GFF3 to BED12");
		print("  gff2junc   convert GFF3 to junction BED");
		print("  longcs2seq convert long-cs PAF to sequences");
//		print("  paf2gff    convert PAF to GFF3 (tested for miniprot only)");
		print("");
		print("  stat       collect basic mapping information in PAF/SAM");
		print("  asmstat    collect basic assembly information");
		print("  asmgene    evaluate gene completeness");
		print("  misjoin    evaluate large-scale misjoins");
		print("  liftover   simplistic liftOver");
		print("  call       call variants from asm-to-ref alignment with the cs tag");
		print("  bedcov     compute the number of bases covered");
		print("  vcfstat    VCF statistics");
		print("  sveval     compare two SV callsets in VCF");
		print("  version    print paftools.js version");
		print("");
		print("  mapeval    evaluate mapping accuracy using mason2/PBSIM-simulated FASTQ");
		print("  pafcmp     compare two PAF files");
		print("  mason2fq   convert mason2-simulated SAM to FASTQ");
		print("  pbsim2fq   convert PBSIM-simulated MAF to FASTQ");
		print("  junceval   evaluate splice junction consistency with known annotations");
		print("  exoneval   evaluate exon-level consistency with known annotations");
		print("  ov-eval    evaluate read overlap sensitivity using read-to-ref mapping");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'view') paf_view(args);
	else if (cmd == 'sam2paf') paf_sam2paf(args);
	else if (cmd == 'delta2paf') paf_delta2paf(args);
	else if (cmd == 'splice2bed') paf_splice2bed(args);
	else if (cmd == 'gff2bed') paf_gff2bed(args);
	else if (cmd == 'gff2junc') paf_gff2junc(args);
	else if (cmd == 'stat') paf_stat(args);
	else if (cmd == 'asmstat') paf_asmstat(args);
	else if (cmd == 'asmgene') paf_asmgene(args);
	else if (cmd == 'misjoin') paf_misjoin(args);
	else if (cmd == 'liftover' || cmd == 'liftOver') paf_liftover(args);
	else if (cmd == 'vcfpair') paf_vcfpair(args);
	else if (cmd == 'call') paf_call(args);
	else if (cmd == 'mapeval') paf_mapeval(args);
	else if (cmd == 'pafcmp') paf_pafcmp(args);
	else if (cmd == 'bedcov') paf_bedcov(args);
	else if (cmd == 'mason2fq') paf_mason2fq(args);
	else if (cmd == 'pbsim2fq') paf_pbsim2fq(args);
	else if (cmd == 'junceval') paf_junceval(args);
	else if (cmd == 'exoneval') paf_exoneval(args);
	else if (cmd == 'ov-eval') paf_ov_eval(args);
	else if (cmd == 'vcfstat') paf_vcfstat(args);
	else if (cmd == 'sveval') paf_sveval(args);
	else if (cmd == 'vcfsel') paf_vcfsel(args);
	else if (cmd == 'longcs2seq') paf_longcs2seq(args);
	else if (cmd == 'paf2gff') paf_paf2gff(args);
	else if (cmd == 'version') print(paftools_version);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);

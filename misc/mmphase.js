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

function read_fastx(file, buf)
{
	if (file.readline(buf) < 0) return null;
	var m, line = buf.toString();
	if ((m = /^([>@])(\S+)/.exec(line)) == null)
		throw Error("wrong fastx format");
	var is_fq = (m[1] == '@');
	var name = m[2];
	if (file.readline(buf) < 0)
		throw Error("missing sequence line");
	var seq = buf.toString();
	if (is_fq) { // skip quality
		file.readline(buf);
		file.readline(buf);
	}
	return [name, seq];
}

function filter_paf(a, opt)
{
	if (a.length == 0) return;
	var k = 0;
	for (var i = 0; i < a.length; ++i) {
		var ai = a[i];
		if (ai[10] < opt.min_blen) continue;
		if (ai[9] < ai[10] * opt.min_iden) continue;
		var clip = [0, 0];
		if (ai[4] == '+') {
			clip[0] = ai[2] < ai[7]? ai[2] : ai[7];
			clip[1] = ai[1] - ai[3] < ai[6] - ai[8]? ai[1] - ai[3] : ai[6] - ai[8];
		} else {
			clip[0] = ai[2] < ai[6] - ai[8]? ai[2] : ai[6] - ai[8];
			clip[1] = ai[1] - ai[3] < ai[7]? ai[1] - ai[3] : ai[7];
		}
		if (clip[0] > opt.max_clip_len || clip[1] > opt.max_clip_len) continue;
		a[k++] = ai;
	}
	a.length = k;
}

function parse_events(t, ev, id, buf)
{
	var re = /(:(\d+))|(([\+\-\*])([a-z]+))/g;
	var m, cs = null;
	for (var j = 12; j < t.length; ++j) {
		if ((m = /^cs:Z:(\S+)/.exec(t[j])) != null) {
			cs = m[1].toLowerCase();
			break;
		}
	}
	if (cs == null) {
		warn("Warning: no cs tag for read '" + t[0] + "'");
		return;
	}
	var st = t[2], en = t[3];
	var x = st;
	while ((m = re.exec(cs)) != null) {
		var l;
		if (m[2] != null) { // an identitcal match ":\d+"
			l = parseInt(m[2]);
			// [start, end, type, index, changed_base]
			ev.push([x, x + l, 0, id]);
		} else {
			if (m[4] == '*') {
				l = 1;
				ev.push([x, x + 1, 1, id, m[5][0]]);
			} else if (m[4] == '+') {
				l = m[5].length;
				ev.push([x, x + l, 2, id]);
			} else if (m[4] == '-') {
				l = 0;
				ev.push([x, x, -1, id, m[5]]);
			}
		}
		x += l;
	}
	if (x != en)
		throw Error("inconsistent cs for read '" + t[0] + "'");
}

function find_het_sub(ev, a, opt)
{
	var n = a.length, last0_i = -1, h = [], d = [];
	for (var i = 0; i < n; ++i) h[i] = [], d[i] = [];
	for (var i = 0; i < ev.length; ++i) {
		if (ev[i][2] == 0) {
			if (last0_i < 0 || ev[i][0] != ev[last0_i][0]) last0_i = i;
			else if (ev[i][1] > ev[last0_i][1])
				last0_i = i;
		} else if (ev[i][2] == 1 && last0_i >= 0 && ev[i][0] < ev[last0_i][1]) {
			if (ev[last0_i][1] - ev[last0_i][0] >= opt.min_mlen) {
				if (opt.dbg_ev) print("EV", ev[last0_i].join("\t"), "|", ev[i].join("\t"));
				var e0 = ev[last0_i], hl = h[e0[3]];
				if (hl.length == 0 || hl[hl.length-1][0] != e0[0])
					hl.push([e0[0], e0[1]]);
				d[ev[i][3]].push([ev[i][0], e0[1] - e0[0]]);
			}
		}
	}
	var b = [];
	for (var i = 0; i < n; ++i) {
		var sh = 0, dh = 0;
		for (var j = 0; j < h[i].length; ++j)
			sh += h[i][j][1] - h[i][j][0];
		for (var j = 0; j < d[i].length; ++j)
			dh += d[i][j][1];
		// [start, end, index, #consistent, lenConsistent, #conflictive, lenConflictive, identity, mlen]
		b[i] = [a[i][2], a[i][3], i, h[i].length, sh, d[i].length, dh, a[i][9] / a[i][10], a[i][9]];
	}
	return b;
}

function flt_utg_for_ec(b, opt)
{
	var k = 0;
	for (var i = 0; i < b.length; ++i) {
		var bi = b[i];
		if (bi[4] == 0 && bi[6] == 0) b[k++] = bi; // entirely ambiguous
		else if (bi[6] < (bi[4] + bi[6]) * opt.max_ratio0) b[k++] = bi;
	}
	b.length = k;
	if (b.length == 0) return;
	// find the longest contiguous segment
	b.sort(function(x,y) { return x[0]-y[0] });
	var st = b[0][0], en = b[0][1], max_st = 0, max_en = 0, max_max_en = en;
	for (var i = 1; i < b.length; ++i) {
		if (b[i][0] > en) {
			if (en - st > max_en - max_st)
				max_st = st, max_en = en;
			st = b[i][0], en = b[i][1];
		} else {
			en = en > b[i][1]? en : b[i][1];
		}
		max_max_en = max_max_en > b[i][1]? max_max_en : b[i][1];
	}
	if (en - st > max_en - max_st)
		max_st = st, max_en = en;
	if (max_max_en != en || st != b[0][0]) {
		var k = 0;
		for (var i = 0; i < b.length; ++i)
			if (b[i][0] < max_en && b[i][1] > max_st)
				b[k++] = b[i];
		b.length = k;
	}
}

function flt_utg_for_bin(b, opt) // filter out alignments clearly on the wrong phase
{
	var k = 0;
	for (var i = 0; i < b.length; ++i) {
		var bi = b[i];
		if (bi[4] + bi[6] == 0 || bi[4] >= (bi[4] + bi[6]) * opt.max_ratio0) b[k++] = bi;
	}
	b.length = k;
}

function ec_core(b, n_a, ev, buf, ecb) // error correction
{
	var intv = [];
	for (var i = 0; i < n_a; ++i)
		intv[i] = null;
	intv[b[0][2]] = [b[0][0], b[0][1]];
	var en = b[0][1];
	for (var i = 1; i < b.length; ++i) {
		if (b[i][1] <= en) continue;
		intv[b[i][2]] = [en, b[i][1]];
		en = b[i][1];
	}
	var k = 0;
	ecb.capacity = buf.capacity;
	ecb.length = 0;
	for (var i = 0; i < ev.length; ++i) {
		var e = ev[i], I = intv[e[3]];
		if (I == null) continue;
		if (e[0] >= I[0] && e[0] < I[1]) { // this is to reduce duplicated events around junctions
			//print("X", e.join("\t"));
			if (e[2] == 0) {
				ecb.length += e[1] - e[0];
				for (var j = e[0]; j < e[1]; ++j)
					ecb[k++] = buf[j];
			} else if (e[2] == 1) {
				++ecb.length;
				ecb[k++] = e[4].charCodeAt(0);
			} else if (e[2] < 0) {
				ecb.length += e[4].length;
				for (var j = 0; j < e[4].length; ++j)
					ecb[k++] = e[4].charCodeAt(j);
			} // else, skip e[2] == 2
		}
	}
	if (ecb.length != k) throw Error("BUG!");
}

function process_paf(a, opt, fp_seq, buf, ecb)
{
	if (a.length == 0) return;
	var len = a[0][1], name = a[0][0], seq = null;
	if (len < opt.min_rlen) return;
	if (fp_seq) {
		var ret;
		while ((ret = read_fastx(fp_seq, buf)) != null)
			if (ret[0] == a[0][0])
				break;
		if (ret == null)
			throw Error("failed to find sequence for read '" + a[0][0] + "'");
		name = ret[0], seq = ret[1];
		if (seq.length != len)
			throw Error("inconsistent length for read '" + name + "'");
	}
	filter_paf(a, opt);
	if (a.length == 0) return;
	var ev = [];
	for (var i = 0; i < a.length; ++i)
		parse_events(a[i], ev, i, buf);
	ev.sort(function(x,y) { return x[0]!=y[0]? x[0]-y[0] : x[2]-y[2] });
	if (seq == null) print("SQ", name, a[0][1], a.length);
	var b = find_het_sub(ev, a, opt);
	if (opt.ec) flt_utg_for_ec(b, opt);
	else flt_utg_for_bin(b, opt);
	if (seq == null) {
		for (var i = 0; i < b.length; ++i) {
			var m, ai = a[b[i][2]], score = 0;
			for (var j = 10; j < ai.length; ++j)
				if ((m = /^AS:i:(\d+)/.exec(ai[j])) != null)
					score = m[1];
			print("TS", b[i][2], b[i][0], b[i][1], ai.slice(5, 9).join("\t"), b[i].slice(3, 7).join("\t"), score);
		}
		print("//");
	} else { // error correction
		if (b.length == 0) return;
		buf.set(seq, 0);
		ec_core(b, a.length, ev, buf, ecb);
		print(">" + name);
		print(ecb);
	}
}

function main(args)
{
	var c, opt = { min_rlen:5000, min_blen:5000, min_iden:0.8, min_mlen:5, max_clip_len:500, max_ratio0:0.25, dbg_ev:false };
	while ((c = getopt(args, "l:b:d:m:c:r:E")) != null) {
		if (c == 'l') opt.min_rlen = parseInt(getopt.arg);
		else if (c == 'b') opt.min_blen = parseInt(getopt.arg);
		else if (c == 'd') opt.min_iden = parseFloat(getopt.arg);
		else if (c == 'm') opt.min_slen = parseInt(getopt.arg);
		else if (c == 'c') opt.max_clip_len = parseInt(getopt.arg);
		else if (c == 'r') opt.max_ratio0 = parseFloat(getopt.arg);
		else if (c == 'E') opt.dbg_ev = true;
	}
	if (args.length - getopt.ind < 1) {
		print("Usage: mmphase.js [options] <map-with-cs.paf> [reads.fa]");
		print("Options:");
		print("  -l INT    min read length [" + opt.min_rlen + "]");
		print("  -b INT    min alignment length [" + opt.min_blen + "]");
		print("  -d FLOAT  min identity [" + opt.min_iden + "]");
		print("  -s INT    min match length [" + opt.min_mlen + "]");
		print("  -c INT    max clip length [" + opt.max_clip_len + "]");
		print("  -r FLOAT  initial ratio for haplotype filtering [" + opt.max_ratio0 + "]");
		return 0;
	}

	opt.ec = args.length - getopt.ind < 2? false : true;
	if (!opt.ec) {
		print("CC");
		print("CC", "SQ  qName  qLen    nHits");
		print("CC", "TS  index  qStart  qEnd  tName  tLen  tStart  tEnd  nConsistent  lCons  nConflictive  lConf  score");
		print("CC");
	}

	var buf = new Bytes(), ecb = new Bytes();
	var fp_paf = new File(args[getopt.ind]);
	var fp_seq = args.length - getopt.ind >= 2? new File(args[getopt.ind+1]) : null;
	var a = [];
	while (fp_paf.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (a.length > 0 && a[0][0] != t[0]) {
			process_paf(a, opt, fp_seq, buf, ecb);
			a.length = 0;
		}
		for (var i = 1; i <= 3; ++i) t[i] = parseInt(t[i]);
		if (t[1] < opt.min_rlen) continue;
		for (var i = 6; i <= 10; ++i) t[i] = parseInt(t[i]);
		if (t[10] < opt.min_blen) continue;
		a.push(t);
	}
	if (a.length >= 0)
		process_paf(a, opt, fp_seq, buf, ecb);
	if (fp_seq) fp_seq.close();
	fp_paf.close();
	ecb.destroy();
	buf.destroy();
}

var ret = main(arguments)
exit(ret)

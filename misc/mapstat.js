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

var c, gap_out_len = null;
while ((c = getopt(arguments, "l:")) != null)
	if (c == 'l') gap_out_len = parseInt(getopt.arg);

if (getopt.ind == arguments.length) {
	print("Usage: k8 mapstat.js [-l gapOutLen] <in.sam>|<in.paf>");
	exit(1);
}

var buf = new Bytes();
var file = new File(arguments[getopt.ind]);
var re = /(\d+)([MIDSHNX=])/g;

var lineno = 0, n_pri = 0, n_2nd = 0, n_seq = 0, n_cigar_64k = 0, l_tot = 0, l_cov = 0;
var n_gap = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]];

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
		var atlen = null, aqlen, qs, qe, mapq, ori_qlen;
		if (t[4] == '+' || t[4] == '-') { // PAF
			if (!/\ts2:i:\d+/.test(line)) {
				++n_2nd;
				continue;
			}
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
		while ((m = re.exec(cigar)) != null) {
			var l = parseInt(m[1]);
			++n_cigar;
			if (m[2] == 'M' || m[2] == '=' || m[2] == 'X') {
				tl += l, ql += l, M += l;
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
			} else if (m[2] == 'N') {
				tl += l;
			} else if (m[2] == 'S') {
				clip[M == 0? 0 : 1] = l, sclip += l;
			} else if (m[2] == 'H') {
				clip[M == 0? 0 : 1] = l;
			}
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
		regs.push([qs, qe]);
		last_qlen = ori_qlen;
	}
}
l_tot += last_qlen;
l_cov += cov_len(regs);

file.close();
buf.destroy();

if (gap_out_len == null) {
	print("Number of mapped sequences: " + n_seq);
	print("Number of primary alignments: " + n_pri);
	print("Number of secondary alignments: " + n_2nd);
	print("Number of primary alignments with >65535 CIGAR operations: " + n_cigar_64k);
	print("Number of bases in mapped sequences: " + l_tot);
	print("Number of mapped bases: " + l_cov);
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

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

var c, pri_only = false;
while ((c = getopt(arguments, "p")) != null)
	if (c == 'p') pri_only = true;

var file = arguments.length == getopt.ind? new File() : new File(arguments[getopt.ind]);
var buf = new Bytes();
var re = /(\d+)([MIDSHNX=])/g;

var len = {}, lineno = 0;
while (file.readline(buf) >= 0) {
	var m, n_cigar = 0, line = buf.toString();
	++lineno;
	if (line.charAt(0) == '@') {
		if (/^@SQ/.test(line)) {
			var name = (m = /\tSN:(\S+)/.exec(line)) != null? m[1] : null;
			var l = (m = /\tLN:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;
			if (name != null && l != null) len[name] = l;
		}
		continue;
	}
	var t = line.split("\t");
	var flag = parseInt(t[1]);
	if (t[9] != '*' && t[10] != '*' && t[9].length != t[10].length) throw Error("ERROR at line " + lineno + ": inconsistent SEQ and QUAL lengths - " + t[9].length + " != " + t[10].length);
	if (t[2] == '*' || (flag&4)) continue;
	if (pri_only && (flag&0x100)) continue;
	var tlen = len[t[2]];
	if (tlen == null) throw Error("ERROR at line " + lineno + ": can't find the length of contig " + t[2]);
	var nn = (m = /\tnn:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : 0;
	var NM = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;
	var have_NM = NM == null? false : true;
	NM += nn;
	var clip = [0, 0], I = [0, 0], D = [0, 0], M = 0, N = 0, ql = 0, tl = 0, mm = 0, ext_cigar = false;
	while ((m = re.exec(t[5])) != null) {
		var l = parseInt(m[1]);
		if (m[2] == 'M') M += l, ql += l, tl += l, ext_cigar = false;
		else if (m[2] == 'I') ++I[0], I[1] += l, ql += l;
		else if (m[2] == 'D') ++D[0], D[1] += l, tl += l;
		else if (m[2] == 'N') N += l, tl += l;
		else if (m[2] == 'S') clip[M == 0? 0 : 1] = l, ql += l;
		else if (m[2] == 'H') clip[M == 0? 0 : 1] = l;
		else if (m[2] == '=') M += l, ql += l, tl += l, ext_cigar = true;
		else if (m[2] == 'X') M += l, ql += l, tl += l, mm += l, ext_cigar = true;
		++n_cigar;
	}
	if (n_cigar > 65535)
		warn("WARNING at line " + lineno + ": " + n_cigar + " CIGAR operations");
	if (tl + parseInt(t[3]) - 1 > tlen) {
		warn("WARNING at line " + lineno + ": alignment end position larger than ref length; skipped");
		continue;
	}
	if (t[9] != '*' && t[9].length != ql) {
		warn("WARNING at line " + lineno + ": SEQ length inconsistent with CIGAR (" + t[9].length + " != " + ql + "); skipped");
		continue;
	}
	if (!have_NM || ext_cigar) NM = I[1] + D[1] + mm;
	if (NM < I[1] + D[1] + mm) {
		warn("WARNING at line " + lineno + ": NM is less than the total number of gaps (" + NM + " < " + (I[1]+D[1]+mm) + ")");
		NM = I[1] + D[1] + mm;
	}
	var extra = ["mm:i:"+(NM-I[1]-D[1]), "io:i:"+I[0], "in:i:"+I[1], "do:i:"+D[0], "dn:i:"+D[1]];
	var match = M - (NM - I[1] - D[1]);
	var blen = M + I[1] + D[1];
	var qlen = M + I[1] + clip[0] + clip[1];
	var qs, qe;
	if (flag&16) qs = clip[1], qe = qlen - clip[0];
	else qs = clip[0], qe = qlen - clip[1];
	var ts = parseInt(t[3]) - 1, te = ts + M + D[1] + N;
	var a = [t[0], qlen, qs, qe, flag&16? '-' : '+', t[2], tlen, ts, te, match, blen, t[4]];
	print(a.join("\t"), extra.join("\t"));
}

buf.destroy();
file.close();

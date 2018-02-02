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

var colors = ["0,128,255", "255,0,0", "0,192,0"];

function print_lines(a, fmt) {
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
			a[i][8] = colors[a[i][8]];
			print(a[i].join("\t"));
		}
	}
	a.length = 0;
}

function main(args) {
	var re = /(\d+)([MIDNSH])/g;
	var c, fmt = "bed", fn_name_conv = null;
	while ((c = getopt(args, "f:n:")) != null) {
		if (c == 'f') fmt = getopt.arg;
		else if (c == 'n') fn_name_conv = getopt.arg;
	}
	if (getopt.ind == args.length) {
		warn("Usage: k8 splice2bed.js <in.paf>|<in.sam>");
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
			print_lines(a, fmt);
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
	print_lines(a, fmt);
	buf.destroy();
	file.close();
	if (conv != null) conv.destroy();
}

main(arguments);

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
			if (a[i][9] == 0) ++n_pri;
		if (n_pri > 1) {
			for (var i = 0; i < a.length; ++i)
				if (a[i][9] == 0) a[i][9] = 1;
		} else if (n_pri == 0) {
			warn("Warning: " + a[0][0] + " doesn't have a primary alignment");
		}
		for (var i = 0; i < a.length; ++i) {
			a[i][9] = colors[a[i][9]];
			a[i].shift();
			print(a[i].join("\t"));
		}
	}
	a.length = 0;
}

function main(args) {
	var re = /(\d+)([MIDNSH])/g;
	var c, fmt = "bed", with_hdr = false, hdr_only = false, name = 'splice2bed', pos = null;
	while ((c = getopt(args, "Hhn:p:")) != null) {
		if (c == 'h') with_hdr = true;
		else if (c == 'p') pos = getopt.arg;
		else if (c == 'H') with_hdr = hdr_only = true;
		else if (c == 'n') name = getopt.arg;
	}
	if (getopt.ind == args.length && !hdr_only) {
		warn("Usage: k8 splice2bed.js [options] <in.paf>");
		exit(1);
	}
	if (with_hdr) {
		if (pos != null)
			print('browser position ' + pos);
		print('track name=' + name + ' useScore=1 visibility=2 itemRgb="On"');
	}
	if (hdr_only) return;
	var file = new File(args[getopt.ind]);
	var buf = new Bytes();
	var a = [];
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '@') continue; // skip SAM header lines
		var t = line.split("\t");
		if (t.length >= 12 && (t[4] == '+' || t[4] == '-')) {
			if (a.length && a[0][0] != t[0]) {
				print_lines(a, fmt);
				a = [];
			}
			var is_pri = false, cigar = null;
			for (var i = 12; i < t.length; ++i) {
				if (t[i].substr(0, 5) == 'cg:Z:') {
					cigar = t[i].substr(5);
				} else if (t[i].substr(0, 5) == 's2:i:') {
					is_pri = true;
				}
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
			a.push([t[0], t[5], t[7], t[8], t[0], Math.floor(t[9]/t[10]*1000), t[4], t[7], t[8], is_pri? 0:2, bs.length, bl.join(",")+",", bs.join(",")+","]);
		} else {
			throw Error("unrecognized input format");
		}
	}
	print_lines(a, fmt);
	buf.destroy();
	file.close();
}

main(arguments);

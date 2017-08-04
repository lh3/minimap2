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

var c, maf_out = false, line_len = 0;
while ((c = getopt(arguments, "ml:")) != null) {
	if (c == 'm') maf_out = true;
	else if (c == 'l') line_len = parseInt(getopt.arg); // TODO: not implemented yet
}

if (getopt.ind == arguments.length) {
	print("Usage: k8 paf2aln.js [options] <with-cs.paf>");
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

var s_ref = new Bytes(), s_qry = new Bytes(), s_mid = new Bytes();
var re = /([=\-\+\*])([A-Za-z]+)/g;

var buf = new Bytes();
var file = new File(arguments[getopt.ind]);
if (maf_out) print("##maf version=1\n");
while (file.readline(buf) >= 0) {
	var m, line = buf.toString();
	var t = line.split("\t", 12);
	if ((m = /\tcs:Z:(\S+)/.exec(line)) == null) continue;
	var cs = m[1];
	s_ref.length = s_qry.length = s_mid.length = 0;
	if (line_len == 0 || maf_out) {
		var ql = 0, tl = 0, al = 0;
		var n_mm = 0, n_gap = 0;
		while ((m = re.exec(cs)) != null) {
			var l = m[2].length;
			if (m[1] == '=') {
				s_ref.set(m[2]);
				s_qry.set(m[2]);
				s_mid.set(Array(l+1).join("|"));
			} else if (m[1] == '*') {
				s_ref.set(m[2].charAt(0));
				s_qry.set(m[2].charAt(1));
				s_mid.set(' ');
				++n_mm;
			} else if (m[1] == '+') {
				s_ref.set(Array(l+1).join("-"));
				s_qry.set(m[2]);
				s_mid.set(Array(l+1).join(" "));
				n_gap += l;
			} else if (m[1] == '-') {
				s_ref.set(m[2]);
				s_qry.set(Array(l+1).join("-"));
				s_mid.set(Array(l+1).join(" "));
				n_gap += l;
			}
		}
		if (maf_out) {
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
		} else {
			print(line);
			print(n_mm, n_gap);
			print(s_ref); print(s_mid); print(s_qry);
		}
	}
}
file.close();
buf.destroy();

s_ref.destroy(); s_qry.destroy(); s_mid.destroy();

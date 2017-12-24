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

var c, min_ovlp = 2000, min_frac = 0.95, min_mapq = 10;
while ((c = getopt(arguments, "q:l:f:")) != null) {
	if (c == 'q') min_mapq = parseInt(getopt.arg);
	else if (c == 'l') min_ovlp = parseInt(getopt.arg);
	else if (c == 'f') min_frac = parseFloat(getopt.arg);
}
if (getopt.ind == arguments.length) {
	print("Usage: sort -k6,6 -k8,8n to-ref.paf | k8 paf2ovlp.js [options] - > out.ovlp");
	print("Options:");
	print("  -l INT     min overlap length [2000]");
	print("  -q INT     min mapping quality [10]");
	print("  -f FLOAT   min fraction of mapped length [0.95]");
	exit(1);
}

var buf = new Bytes();
var file = arguments[getopt.ind] == '-'? new File() : new File(arguments[getopt.ind]);
var a = [];
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
		var len = en > a[j][2]? a[j][2] - st : en - st;
		if (len >= min_ovlp) {
			if (a[j][3] < t[0]) print(a[j][3], t[0], len);
			else print(t[0], a[j][3], len);
		}
	}
	a.push([ctg, st, en, t[0]]);
}
file.close();
buf.destroy();

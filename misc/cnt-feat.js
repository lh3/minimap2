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

Interval.dedup = function(a, sorted)
{
	if (typeof sorted == 'undefined') sorted = true;
	if (!sorted) Interval.sort(a);
	var k = 0;
	for (var i = 1; i < a.length; ++i)
		if (a[k][0] != a[i][0] || a[k][1] != a[i][1])
			a[++k] = a[i].slice(0);
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

/*****************
 * Main function *
 *****************/

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

function main(args)
{
	var c, print_len = false, to_merge = true, to_dedup = false, fn_excl = null;
	while ((c = getopt(args, "pde:")) != null) {
		if (c == 'p') print_len = true;
		else if (c == 'd') to_dedup = true, to_merge = false;
		else if (c == 'e') fn_excl = getopt.arg;
	}

	if (args.length - getopt.ind < 2) {
		print("Usage: k8 cnt-feat.js [options] <target.bed> <feature.bed>");
		print("Options:");
		print("  -e FILE    exclude features overlapping regions in BED FILE []");
		print("  -p         print number of covered bases for each feature");
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

	warn("# feature bases: " + tot_len);
	warn("# feature bases overlapping targets: " + hit_len + ' (' + (100.0 * hit_len / tot_len).toFixed(2) + '%)');
}

main(arguments);

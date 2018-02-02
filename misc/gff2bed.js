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

var c, fn_ucsc_fai = null, is_short = false;
while ((c = getopt(arguments, "u:s")) != null) {
	if (c == 'u') fn_ucsc_fai = getopt.arg;
	else if (c == 's') is_short = true;
}

if (getopt.ind == arguments.length) {
	print("Usage: k8 gff2bed.js [-u ucsc-genome.fa.fai] <in.gff>");
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
		'lincRNA':'0,192,0',
		'snRNA':'0,192,0',
		'miRNA':'0,192,0',
		'misc_RNA':'0,192,0'
	};

function print_bed12(exons, cds_st, cds_en, is_short)
{
	if (exons.length == 0) return;
	var name = is_short? exons[0][7] + "|" + exons[0][5] : exons[0].slice(4, 7).join("|");
	var a = exons.sort(function(a,b) {return a[1]-b[1]});
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

var re_gtf = /(transcript_id|transcript_type|transcript_biotype|gene_name|transcript_name) "([^"]+)";/g;
var re_gff3 = /(transcript_id|transcript_type|transcript_biotype|gene_name|transcript_name)=([^;]+)/g;
var buf = new Bytes();
var file = arguments[getopt.ind] == '-'? new File() : new File(arguments[getopt.ind]);

var exons = [], cds_st = 1<<30, cds_en = 0, last_id = null;
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (t[0].charAt(0) == '#') continue;
	if (t[2] != "CDS" && t[2] != "exon") continue;
	t[3] = parseInt(t[3]) - 1;
	t[4] = parseInt(t[4]);
	var id = null, type = "", gname = "N/A", biotype = "", m, tname = "N/A";
	while ((m = re_gtf.exec(t[8])) != null) {
		if (m[1] == "transcript_id") id = m[2];
		else if (m[1] == "transcript_type") type = m[2];
		else if (m[1] == "transcript_biotype") biotype = m[2];
		else if (m[1] == "gene_name") name = m[2];
		else if (m[1] == "transcript_name") tname = m[2];
	}
	while ((m = re_gff3.exec(t[8])) != null) {
		if (m[1] == "transcript_id") id = m[2];
		else if (m[1] == "transcript_type") type = m[2];
		else if (m[1] == "transcript_biotype") biotype = m[2];
		else if (m[1] == "gene_name") name = m[2];
		else if (m[1] == "transcript_name") tname = m[2];
	}
	if (type == "" && biotype != "") type = biotype;
	if (id == null) throw Error("No transcript_id");
	if (id != last_id) {
		print_bed12(exons, cds_st, cds_en, is_short);
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
	print_bed12(exons, cds_st, cds_en, is_short);

file.close();
buf.destroy();

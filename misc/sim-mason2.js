#!/usr/bin/env k8

Bytes.prototype.reverse = function()
{
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[i];
		this[i] = this[this.length - i - 1];
		this[this.length - i - 1] = tmp;
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
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[this.length - i - 1];
		this[this.length - i - 1] = Bytes.rctab[this[i]];
		this[i] = Bytes.rctab[tmp];
	}
	if (this.length&1)
		this[this.length>>1] = Bytes.rctab[this[this.length>>1]];
}

if (arguments.length == 0) {
	print("Usage: k8 sim-mason2.js <mason.sam>");
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
var file = new File(arguments[0]);
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

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

if (arguments.length < 2) {
	print("Usage: k8 sim-pbsim.js <ref.fa.fai> <pbsim1.maf> [[pbsim2.maf] ...]");
	exit(1);
}

var file, buf = new Bytes(), buf2 = new Bytes();
file = new File(arguments[0]);
var chr_list = [];
while (file.readline(buf) >= 0) {
	var t = buf.toString().split(/\s+/);
	chr_list.push(t[0]);
}
file.close();

for (var k = 1; k < arguments.length; ++k) {
	var fn = arguments[k];
	file = new File(fn);
	var state = 0, reg;
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (state == 0 && line.charAt(0) == 'a') {
			state = 1;
		} else if (state == 1 && line.charAt(0) == 's') {
			var t = line.split(/\s+/);
			var st = parseInt(t[2]);
			reg = [st, st + parseInt(t[3])];
			state = 2;
		} else if (state == 2 && line.charAt(0) == 's') {
			var m, t = line.split(/\s+/);
			if ((m = /S(\d+)_\d+/.exec(t[1])) == null) throw Error("Failed to parse the read name");
			var chr_id = parseInt(m[1]) - 1;
			if (chr_id >= chr_list.length) throw Error("Index outside the chr list");
			var name = [t[1], chr_list[chr_id], reg[0], reg[1], t[4]].join("!");
			var seq = t[6].replace(/\-/g, "");
			if (seq.length != parseInt(t[5])) throw Error("Inconsistent read length");
			if (seq.indexOf("NN") < 0) {
				if (t[4] == '-') {
					buf2.set(seq, 0);
					buf2.length = seq.length;
					buf2.revcomp();
					seq = buf2.toString();
				}
				print(">" + name);
				print(seq);
			}
			state = 0;
		}
	}
	file.close();
}
buf.destroy();
buf2.destroy();

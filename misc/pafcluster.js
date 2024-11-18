#!/usr/bin/env k8

"use strict";

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

function merge_hits(b) {
	if (b.length == 1)
		return { name1:b[0].name1, name2:b[0].name2, len1:b[0].len1, len2:b[0].len2, min_cov:b[0].min_cov, max_cov:b[0].max_cov, cov1:b[0].cov1, cov2:b[0].cov2, s1:b[0].s1, dv:b[0].dv };
	b.sort(function(x, y) { return x.st1 - y.st1 });
	let f = [], bt = [];
	for (let i = 0; i < b.length; ++i)
		f[i] = b[i].s1, bt[i] = -1;
	for (let i = 0; i < b.length; ++i) {
		for (let j = 0; j < i; ++j) {
			if (b[j].st2 < b[i].st2) {
				if (b[j].en1 >= b[i].en1) continue;
				if (b[j].en2 >= b[i].en2) continue;
				const ov1 = b[j].en1 <= b[i].st1? 0 : b[i].st1 - b[j].en1;
				const li1 = b[i].en1 - b[i].st1;
				const s11 = b[i].s1 / li1 * (li1 - ov1);
				const ov2 = b[j].en2 <= b[i].st2? 0 : b[i].st2 - b[j].en2;
				const li2 = b[i].en2 - b[i].st2;
				const s12 = b[i].s1 / li2 * (li2 - ov2);
				const s1 = s11 < s12? s11 : s12;
				if (f[i] < f[j] + s1)
					f[i] = f[j] + s1, bt[i] = j;
			}
		}
	}
	let max_i = -1, max_f = 0, d = [];
	for (let i = 0; i < b.length; ++i)
		if (max_f < f[i])
			max_f = f[i], max_i = i;
	for (let k = max_i; k >= 0; k = bt[k])
		d.push(k);
	d = d.reverse();
	let dv = 0, tot = 0, cov1 = 0, cov2 = 0, st1 = 0, en1 = 0, st2 = 0, en2 = 0;
	for (let k = 0; k < d.length; ++k) {
		const i = d[k];
		tot += b[i].blen;
		dv += b[i].dv * b[i].blen;
		if (b[i].st1 > en1) {
			cov1 += en1 - st1;
			st1 = b[i].st1, en1 = b[i].en1;
		} else en1 = en1 > b[i].en1? en1 : b[i].en1;
		if (b[i].st2 > en2) {
			cov2 += en2 - st2;
			st2 = b[i].st2, en2 = b[i].en2;
		} else en2 = en2 > b[i].en2? en2 : b[i].en2;
	}
	dv /= tot;
	cov1 = (cov1 + (en1 - st1)) / b[0].len1;
	cov2 = (cov2 + (en2 - st2)) / b[0].len2;
	const min_cov = cov1 < cov2? cov1 : cov2;
	const max_cov = cov1 > cov2? cov1 : cov2;
	//warn(d.length, b[0].name1, b[0].name2, min_cov, max_cov);
	return { name1:b[0].name1, name2:b[0].name2, len1:b[0].len1, len2:b[0].len2, min_cov:min_cov, max_cov:max_cov, cov1:cov1, cov2:cov2, s1:max_f, dv:dv };
}

function main(args) {
	let opt = { min_cov:.9, max_dv:.015, max_diff:20000 };
	for (const o of getopt(args, "c:d:e:", [])) {
		if (o.opt == '-c') opt.min_cov = parseFloat(o.arg);
		else if (o.opt == '-d') opt.max_dv = parseFloat(o.arg);
		else if (o.opt == '-e') opt.max_diff = parseFloat(o.arg);
	}
	if (args.length == 0) {
		print("Usage: pafcluster.js [options] <ava.paf>");
		print("Options:");
		print(`  -c FLOAT     min coverage [${opt.min_cov}]`);
		print(`  -d FLOAT     max divergence [${opt.max_dv}]`);
		print(`  -e FLOAT     max difference [${opt.max_diff}]`);
		return;
	}

	// read
	let a = [], len = {}, name2len = {};
	for (const line of k8_readline(args[0])) {
		let m, t = line.split("\t");
		if (t[4] != "+") continue;
		for (let i = 1; i < 4;  ++i) t[i] = parseInt(t[i]);
		for (let i = 6; i < 11; ++i) t[i] = parseInt(t[i]);
		const len1 = t[1], len2 = t[6];
		let s1 = -1, dv = -1.0;
		for (let i = 12; i < t.length; ++i) {
			if ((m = /^(s1|dv):\S:(\S+)/.exec(t[i])) != null) {
				if (m[1] == "s1") s1 = parseInt(m[2]);
				else if (m[1] == "dv") dv = parseFloat(m[2]);
			}
		}
		if (s1 < 0 || dv < 0) continue;
		const cov1 = (parseInt(t[3]) - parseInt(t[2])) / len1;
		const cov2 = (parseInt(t[8]) - parseInt(t[7])) / len2;
		const min_cov = cov1 < cov2? cov1 : cov2;
		const max_cov = cov1 > cov2? cov1 : cov2;
		name2len[t[0]] = len1;
		name2len[t[5]] = len2;
		a.push({ name1:t[0], name2:t[5], len1:len1, len2:len2, min_cov:min_cov, max_cov:max_cov, s1:s1, dv:dv, cov1:cov1, cov2:cov2, st1:t[2], en1:t[3], st2:t[7], en2:t[8], blen:t[10] });
		len[t[0]] = len1, len[t[5]] = len2;
	}
	warn(`Read ${a.length} hits`);

	// merge duplicated hits
	let h = {};
	for (let i = 0; i < a.length; ++i) {
		const key = `${a[i].name1}\t${a[i].name2}`;
		if (h[key] == null) h[key] = [];
		h[key].push(a[i]);
	}
	a = [];
	for (const key in h)
		a.push(merge_hits(h[key]));

	// core loop
	while (a.length > 1) {
		// select the sequence with the highest sum of s1
		let h = {};
		for (let i = 0; i < a.length; ++i) {
			if (h[a[i].name1] == null) h[a[i].name1] = 0;
			h[a[i].name1] += a[i].s1;
		}
		let max_s1 = 0, max_name = "";
		for (const name in h)
			if (max_s1 < h[name])
				max_s1 = h[name], max_name = name;
		// find contigs in the same group
		h = {};
		h[max_name] = 1;
		for (let i = 0; i < a.length; ++i) {
			if (a[i].name1 != max_name && a[i].name2 != max_name)
				continue;
			const diff1 = a[i].len1 * (1.0 - a[i].cov1);
			const diff2 = a[i].len2 * (1.0 - a[i].cov2);
			if (a[i].min_cov >= opt.min_cov && a[i].dv <= opt.max_dv && diff1 <= opt.max_diff && diff2 <= opt.max_diff)
				h[a[i].name1] = h[a[i].name2] = 1;
		}
		let n = 0;
		for (const key in h) {
			++n;
			delete name2len[key];
		}
		print(`SD\t${max_name}\t${n}`);
		for (const key in h) print(`CL\t${key}\t${len[key]}`);
		print("//");
		// filter out redundant hits
		let b = [];
		for (let i = 0; i < a.length; ++i)
			if (h[a[i].name1] == null && h[a[i].name2] == null)
				b.push(a[i]);
		warn(`Reduced the number of hits from ${a.length} to ${b.length}`);
		a = b;
	}

	// output remaining singletons
	for (const key in name2len) {
		print(`SD\t${key}\t1`);
		print(`CL\t${key}\t${name2len[key]}`);
		print(`//`);
	}
}

main(arguments);

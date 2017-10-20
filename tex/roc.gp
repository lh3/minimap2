set t po eps enh co so "Helvetica,26"

set style line 1 lt 1 pt 1 lc rgb "#e41a1c" lw 2;
set style line 2 lt 1 pt 2 lc rgb "#377eb8" lw 2;
set style line 3 lt 1 pt 3 lc rgb "#4daf4a" lw 2;
set style line 4 lt 1 pt 4 lc rgb "#984ea3" lw 2;
set style line 5 lt 1 pt 6 lc rgb "#ff7f00" lw 2;
set style line 6 lt 1 pt 8 lc rgb "#f781bf" lw 2;

set out "roc-color.eps"

set pointsize 2.0
set size 1.59,1.04
set multiplot layout 1,2

set label "(a)" at graph -0.245,1.06 font "Helvetica-bold,40"
set xlab "Error rate of mapped PacBio reads"
set ylab "Fraction of mapped reads" off +1.8
set ytics 0.02
set yran [0.9:1]

set size 0.8,1
set log x
set format x "10^{%L}"
set key bot right
plot "<./eval2roc.pl blasr-mc.eval" u 2:3 t "blasr-mc" w lp ls 4, \
     "<./eval2roc.pl bwa.eval" u 2:3 t "bwa-mem" w lp ls 2, \
     "<./eval2roc.pl graphmap.eval" u 2:3 t "graphmap" w lp ls 3, \
     "<./eval2roc.pl minialign.eval" u 2:3 t "minialign" w lp ls 1, \
     "<./eval2roc.pl mm2.eval" u 2:3 t "minimap2" w lp ls 6, \
     "<./eval2roc.pl ngmlr.eval" u 2:3 t "ngm-lr" w lp ls 5
unset label

set origin 0.8,0
set size 0.79,1
set label "(b)" at graph -0.245,1.06 font "Helvetica-bold,40"
set xlab "Error rate of mapped short reads"

set key top left
plot "<./eval2roc.pl -n2e7 bowtie2-s3.sam.eval" u 2:3 t "bowtie2" w lp ls 5, \
	 "<./eval2roc.pl -n2e7 bwa-s3.sam.eval" u 2:3 t "bwa-mem" w lp ls 2, \
	 "<./eval2roc.pl -n2e7 mm2-s3.sam.eval" u 2:3 t "minimap2" w lp ls 6, \
	 "<./eval2roc.pl -n2e7 snap-s3.sam.eval" u 2:3 t "snap" w lp ls 3

#unset log
#unset format
#unset key
#set log y
#set ylab "Accumulative mapping error rate" off +0
#set xlab "Mapping quality"
#set yran [1e-5:0.1]
#set ytics 1e-5,0.1
#set format y "10^{%L}"
#set xran [60:0] reverse
#plot "<./eval2roc.pl blasr-mc.eval" u 1:2 w lp ls 4, \
#     "<./eval2roc.pl bwa.eval" u 1:2 t "bwa-mem" w lp ls 2, \
#     "<./eval2roc.pl graphmap.eval" u 1:2 t "graphmap" w lp ls 3, \
#     "<./eval2roc.pl minialign.eval" u 1:2 t "minialign" w lp ls 1, \
#     "<./eval2roc.pl mm2.eval" u 1:2 t "minimap2" w lp ls 6, \
#     "<./eval2roc.pl ngmlr.eval" u 1:2 t "ngm-lr" w lp ls 5

set t po eps enh co so "Helvetica,18"

set style line 1 lt 1 pt 1 lc rgb "#FF0000" lw 2;
set style line 2 lt 1 pt 2 lc rgb "#00C000" lw 2;
set style line 3 lt 1 pt 3 lc rgb "#0080FF" lw 2;
set style line 4 lt 1 pt 4 lc rgb "#C000FF" lw 2;
set style line 5 lt 1 pt 5 lc rgb "#00EEEE" lw 2;
set style line 6 lt 1 pt 6 lc rgb "#C04000" lw 2;
set style line 7 lt 1 lc rgb "#C8C800" lw 2;
set style line 8 lt 1 lc rgb "#FF80FF" lw 2;
set style line 9 lt 1 lc rgb "#4E642E" lw 2;
set style line 10 lt 1 lc rgb "#800000" lw 2;
set style line 11 lt 1 lc rgb "#67B7F7" lw 2;
set style line 12 lt 1 lc rgb "#FFC127" lw 2;

set xlab "False positive rate"
set ylab "Sensitivity"
set yran [0.9:1]

set out "roc-color.eps"
set log x
set format x "10^{%L}"
set key bot right
plot "<./eval2roc.pl blasr-mc.eval" u 2:3 t "blasr-mc" w lp ls 1, \
     "<./eval2roc.pl bwa.eval" u 2:3 t "bwa-mem" w lp ls 2, \
     "<./eval2roc.pl minialign.eval" u 2:3 t "minialign" w lp ls 3, \
     "<./eval2roc.pl mm2.eval" u 2:3 t "minimap2" w lp ls 4, \
     "<./eval2roc.pl ngmlr.eval" u 2:3 t "ngm-lr" w lp ls 5

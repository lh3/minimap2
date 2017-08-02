.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:minimap2.pdf

roc-color.eps:roc.gp
		gnuplot roc.gp

minimap2.pdf:minimap2.tex minimap2.bib roc-color.pdf
		pdflatex minimap2; bibtex minimap2; pdflatex minimap2; pdflatex minimap2;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ minimap2.pdf

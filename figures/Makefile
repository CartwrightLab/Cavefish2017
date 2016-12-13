FIGURES= \
	fig-deltaq-many.pdf \
	fig-hlines.pdf \
	fig-equ.pdf \
	fig-newgrid.pdf \
	fig-Qlines.pdf \
	fig-ulines.pdf \
	fig-newgrid-k.pdf

RSCRIPT=Rscript --vanilla

.PHONY: clean

.SUFFIXES:
.SUFFIXES: .R .pdf

all: $(FIGURES) fig

.R.pdf:
	$(RSCRIPT) $<

clean:
	@$(RM) $(FIGURES)

fig: figure-1.pdf figure-2.pdf figure-3.pdf figure-4.pdf figure-5.pdf figure-6.pdf figure-7.pdf

figure-1.pdf : fig-deltaq-many.pdf
	cp $< $@

figure-2.pdf : fig-hlines.pdf
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=1 -dLastPage=1 -sOutputFile=$@ $<

figure-3.pdf : fig-Qlines.pdf
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=1 -dLastPage=1 -sOutputFile=$@ $<

figure-4.pdf : fig-ulines.pdf
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=1 -dLastPage=1 -sOutputFile=$@ $<

figure-5.pdf : fig-equ.pdf
	cp $< $@

figure-6.pdf : fig-newgrid.pdf
	cp $< $@

figure-7.pdf : fig-newgrid-k.pdf
	cp $< $@

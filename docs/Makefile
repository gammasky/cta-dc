TARGET=cta_1dc_gps_sky_model
LATEX=pdflatex
BIBTEX=bibtex

all: $(TARGET).pdf

$(TARGET).pdf: *.tex $(TARGET).bib

%.pdf: %.tex
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<
	$(LATEX) $<

make clean:
	- rm -f $(TARGET).pdf $(TARGET).aux $(TARGET).log $(TARGET).bbl $(TARGET).blg *.aux

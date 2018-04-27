FILENAME=user_manual

.SUFFIXES:
.SUFFIXES: .bib .pdf .tex
.PHONY: clean

run: $(FILENAME).pdf

$(FILENAME).pdf: $(FILENAME).bbl $(FILENAME).tex
	pdflatex $(FILENAME).tex -draftmode
	pdflatex $(FILENAME).tex 

$(FILENAME).bbl: $(FILENAME).aux
	bibtex $(FILENAME).aux

$(FILENAME).aux: ../../NLP.bib
	pdflatex $(FILENAME).tex -draftmode
	pdflatex $(FILENAME).tex -draftmode

clean:
	rm -rf *.aux *.out *.lof *.log *.lot *.toc *.bbl *.blg *pdf *bak *~


#Makefile für TUM-PH-Thesis
#
# Dieses makefile geht davon aus, dass
# - die Hauptdatei den Namen $(BASEFILENAME).tex trägt und
# - alle TeX-Dateien, die inkludiert werden im aktuellen Verzeichnis liegen,
# - Abbildungen (PDF, PNG, JPG, EPS, PS) im Unterverzeichnis ./Abbildungen/ liegen und
#
# Wenn Abbildungen vermittels PST-PDF  eingebunden werden, die Variable PSTRICKS=1 setzen.

BASEFILENAME = thesis
PSTRICKS = 0

FOLDERSEP = /
PATHSEP = :

ifneq (,$(findstring Windows,$(OS)))
	FOLDERSEP = \\
	PATHSEP = ;
endif

ifeq ($(shell kpsewhich PH.pdf | grep -c ''),0)
	override TEXINPUTS := ./Abbildungen/$(PATHSEP)$(TEXINPUTS)
	export TEXINPUTS
endif

TEXDATEIEN = $(wildcard *.tex) $(BASEFILENAME).bib

PSABBILDUNGEN = $(wildcard ./Abbildungen/*.eps) $(wildcard ./Abbildungen/*.ps) $(wildcard ./Abbildungen/*.pfx) $(wildcard ./Abbildungen/*.ovp)

ifeq ($(PSTRICKS),0)
	PDFABBILDUNGEN = $(wildcard ./Abbildungen/*.png) $(wildcard ./Abbildungen/*.pdf) $(wildcard ./Abbildungen/*.jpg)
else
	PDFABBILDUNGEN = $(wildcard ./Abbildungen/*.png) $(wildcard ./Abbildungen/*.pdf) $(wildcard ./Abbildungen/*.jpg) $(BASEFILENAME)-pics.pdf
endif

CLEANDATEIEN = $(wildcard *.aux) $(BASEFILENAME).out $(BASEFILENAME).log $(wildcard tmp/*) \
	$(BASEFILENAME).toc $(BASEFILENAME).bbl $(BASEFILENAME).blg $(wildcard bu?.bbl) $(wildcard bu?.blg) $(BASEFILENAME)-pics.pdf

pdf: $(BASEFILENAME).pdf

%.bbl: %.aux
	biber $*

$(BASEFILENAME)-pics.pdf: $(PSABBILDUNGEN)
	mkdir -p tmp/
	latex -output-directory=tmp/ $(BASEFILENAME).tex
	dvips -Ppdf -o tmp/$(BASEFILENAME)-pics.ps tmp/$(BASEFILENAME).dvi
	ps2pdf tmp/$(BASEFILENAME)-pics.ps tmp/$(BASEFILENAME)-pics.pdf
	pdfcrop tmp/$(BASEFILENAME)-pics.pdf $(BASEFILENAME)-pics.pdf

$(BASEFILENAME).pdf: $(TEXDATEIEN) $(PDFABBILDUNGEN)
	pdflatex -draftmode "\PassOptionsToPackage{draft}{pst-pdf}\input{$(BASEFILENAME)}"
	biber $(BASEFILENAME)
	pdflatex $(BASEFILENAME).tex
	@if [ "`grep -e "Rerun to get .* right" $(BASEFILENAME).log`" ]; then pdflatex $(BASEFILENAME).tex; fi
	@if [ "`grep -e "Rerun to get .* right" $(BASEFILENAME).log`" ]; then pdflatex $(BASEFILENAME).tex; fi

clean:
	rm -f $(CLEANDATEIEN)

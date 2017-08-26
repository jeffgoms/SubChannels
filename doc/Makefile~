#! /bin/bash

default:	CoupledThermodynamicModel.pdf README.pdf

CoupledThermodynamicModel.pdf:  CoupledThermodynamicModel.tex
CoupledThermodynamicModel.idx:  CoupledThermodynamicModel.idx
README.pdf: README.tex
README.idx: README.idx


%.pdf:	%.tex
	rm -f $*.pdf
	pdflatex $<
	makeindex $*.idx
	pdflatex $<
	if fgrep "Rerun to" $*.log; then $(MAKE) --assume-new $< $@;fi
ifndef NOBIBTEX
	if fgrep "There were undefined"  $*.log;\
	then bibtex $*; \
	$(MAKE) NOBIBTEX=true --assume-new $< $@;fi
endif

%.idx:	%.idx
	makeindex $<

%.aux:	%.aux
	bibtex $<

.PHONY: clean clean-pdf

# We don't rm *.ps because that might clobber image files. Instead, put any
# ps files which should go in $(PSFILES)
clean:	*.pdf
	rm $(TARGETS)*.aux *.log *.bbl *.blg *.toc *.lof *ind *idx *ilg *nav *snm \
*.out *pdf *.pdftex *_tex.pdf *.pdftex_t *.pstex *.pstex_t *lot\
$(addsuffix .png, $(XFIG_IMAGES)) 2>/dev/null||true

htclean_%:
	rm $*.4ct $*.4tc $*.css $*.idv $*.lg $*.tmp $*.xref *.html 2>/dev/null||true
clean-aux:
		rm $(TARGETS)*.aux *.log *.bbl *.blg *.toc *.lof *ind *idx *ilg *nav *snm \
*.out 
clean-pdf: 
		rm *.pdf 

ACLOCAL_AMFLAGS = -I m4
PAPER_SRCS = tutorialBiogeme.tex
PAPER_PDFS = $(PAPER_SRCS:%.tex=%.pdf)

pdf: $(PAPER_PDFS)

clean:
	$(LATEXMK) -c $(PAPER_SRCS)

# This would be clean-local, but we distribute the PDFs
veryclean:
	$(LATEXMK) -C $(PAPER_SRCS)

EXTRA_DIST = $(PAPER_SRCS) $(PAPER_PDFS)

.tex.pdf:
	$(LATEXMK) -pdf $<

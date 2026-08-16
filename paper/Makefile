# Paper build. Targets mirror the verification gates in the plan.
#
#   make numbers   — regenerate frozen numbers (pass FREEZE=v1|v2)
#   make figs      — build all figures into figs/ (never solves; lists
#                    missing runs as paper/runs/ commands on a cache miss)
#   make pdf       — compile main.tex (RevTeX 4.2)
#   make additivity— gate G9: the paper must compile with the equal-time
#                    cross-check subsection deleted, zero dangling refs
#   make check     — gates that are scriptable so far

FREEZE ?= v1
PY ?= python3

.PHONY: numbers figs pdf check additivity clean

numbers:
	$(PY) make_numbers.py --freeze $(FREEZE)

figs:
	$(PY) make_figures.py

pdf:
	pdflatex -interaction=nonstopmode main.tex
	bibtex main || true
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

# G9: compile with the equal-time slot excluded. The subsection lives in its
# own file precisely so exclusion is a file-level operation.
additivity:
	@test ! -f sections/10e-equal-time-crosscheck.tex || \
	  ( mv sections/10e-equal-time-crosscheck.tex /tmp/10e.tex.bak && \
	    $(MAKE) pdf; status=$$?; \
	    mv /tmp/10e.tex.bak sections/10e-equal-time-crosscheck.tex; \
	    exit $$status )

check:
	$(PY) check_gates.py

clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.pdf

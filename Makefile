# closure_quotient_landscape paper build
# Standard latexmk pipeline (pdflatex), with `make clean` and
# `make watch` (continuous rebuild on save).

TEX     = closure_quotient_landscape.tex
PDF     = $(TEX:.tex=.pdf)
LATEXMK = latexmk -pdf -interaction=nonstopmode -halt-on-error

.PHONY: all clean watch verify-scripts

all: $(PDF)

$(PDF): $(TEX)
	$(LATEXMK) $(TEX)

watch:
	$(LATEXMK) -pvc $(TEX)

clean:
	$(LATEXMK) -C $(TEX)
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.fdb_latexmk *.fls *.synctex.gz

# Re-run every probe script and check exit codes. Uses the closure-v5
# venv (must have been bootstrapped via `python3 -m venv ../.venv`).
VENV_PY = ../.venv/bin/python
verify-scripts:
	@for s in scripts/*.py; do \
		echo "==> $$s"; \
		$(VENV_PY) "$$s" > /dev/null || { echo "FAIL: $$s"; exit 1; }; \
	done
	@echo "All probe scripts ran successfully."

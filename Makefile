# One-command reproduction of every number, table and figure of the manuscript.
#   make all        results -> tables -> figures -> article
# See REPRODUCING.md for data sources, seeds and run times.

PYTHON ?= python3
RUN = PYTHONPATH=. $(PYTHON)

ANALYSES = analyze_email_eu_core analyze_openflights analyze_celegans \
           analyze_turbulence network_topology compute_article_numbers

.PHONY: all results tables figures article submission site publish-site test-fast

all: results tables figures article

# All analysis scripts -> results/*.json (fixed seeds; ~13 min, dominated by
# compute_article_numbers and analyze_turbulence).
results:
	@mkdir -p results
	@for s in $(ANALYSES); do \
		echo "== scripts/$$s.py =="; \
		$(RUN) scripts/$$s.py || exit 1; \
	done

# LaTeX tables -> article/generated/*.tex (from results/*.json only).
tables:
	$(RUN) scripts/make_tables.py

# Figures 1-6 -> figures/*.pdf (+ 300 dpi .png previews).
figures:
	$(RUN) scripts/make_theory_figures.py
	$(RUN) scripts/make_data_figures.py

article:
	cd article && latexmk -pdf main.tex

# Flat bundle for Editorial Manager -> article/submission/ (paths rewritten,
# tables inlined, highlights.txt, cover-letter.pdf). main.tex is not touched.
submission:
	$(RUN) scripts/make_submission.py

# pkgdown site of htheoryR -> htheoryR/docs/ (published on the gh-pages branch).
site:
	cd htheoryR && Rscript -e 'pkgdown::build_site(preview = FALSE, install = TRUE, new_process = TRUE)'
	$(PYTHON) scripts/fix_site.py htheoryR/docs

# Replace the gh-pages branch with the current htheoryR/docs (run make site first).
publish-site:
	scripts/publish_site.sh

test-fast:
	$(RUN) -m pytest tests/ -q -m "not slow"

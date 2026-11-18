# Reproducibility

## Environment

Python 3.12 suffices for all exact checks and reference tests. These steps use the standard library and read only repository-relative inputs. Keep assertions enabled for both mathematical checkers. The rendering environment is separate. Install the pinned Python packages in an isolated environment for the figures.

| Command from the repository root | What it verifies |
|---|---|
| python -I -B scripts/check_artifacts.py | SHA-256 of the tape, checker source, manuscript source and reviewed PDF; decoded tape length and digest |
| python -I -B certificate/verify_s07.py certificate/three-star-tape.txt | Sign chamber geometry by cofactors, all nondominant chamber/matching trees, rational/algebraic leaf tests and tight directions |
| python -I -B certificate/verify_s08.py certificate/three-star-tape.txt | The same plan with an alternative clipping and arithmetic implementation; exact four-star fixtures |
| python -I -B certificate/verify_s08.py certificate/three-star-tape.txt --geometry | The preceding checks plus eight-coordinate partition counts: 2,470 maximal chambers, 206 rays, maximum 37 rays in one chamber |
| python -B -m unittest discover -s tests -v | Exact rational variance and admissibility fixtures, numerical Gram values, invalid-input boundaries and the plotted witness |

The eight-coordinate count concerns geometry. A uniform variance comparison over that arrangement is open. The two exact implementations consume the same certificate plan. The run does not search for new witnesses.

CI runs the artifact checks, the reference tests and both standard certificate replays. The optional eight-coordinate count is separate because it is more expensive.

## Reference implementation

[filtered_variance.py](filtered_variance.py) uses rational masses and rational weights for the canonical law, with integer clock endpoints. Deterministic sign flips reduce signed weights to magnitudes. The law recurrence convolves independent child totals, then adds the signed root magnitude. A fresh fair sign is used at zero. The Gram norm is computed numerically by orthogonal child sums and root alignment. The support of the exact law can grow exponentially.

The separate sign-LP verifier checks a supplied sign distribution. It checks every observed matured-history balance and permits crossing schedules, as does the LP itself. It does not optimize the LP or claim the laminar comparison for a crossing schedule.

~~~python
from filtered_variance import Score, laminar_optima

result = laminar_optima([Score(0, 2, 1), Score(0, 1, 1), Score(1, 2, 1)])
print(result.variance)        # 5, exactly
print(result.covariance_sdp) # 3 + 2 sqrt(2), numerically
~~~

## Figures

~~~sh
python -m venv .venv-figures
# Activate .venv-figures using the command appropriate for your shell.
python -m pip install -r requirements-figures.txt
python -B -m scripts.render_figures
~~~

The three files in [figures/](figures/) are SVG with letterforms drawn as paths. Each carries a title. The generator omits dates and fixes the Matplotlib SVG hash salt. Use --out-dir to render to a different directory, and --preview for additional PNG inspection copies. Numeric sources are the interval examples, the proved constants, and the Rademacher moment and Gram formulas in the manuscript. The slice minimizes the quotient over nonnegative outer weight at each heavy coefficient. It illustrates the analytic result and does not prove it numerically. Different font engines or renderers may change raw SVG bytes.

## Manuscript

Enter [manuscript/](manuscript/) before compiling. Build dependencies are a regular TeX installation with AMS packages, geometry, hyperref and BibTeX.

~~~sh
cd manuscript
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
bibtex manuscript
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
~~~

Alternatively use Tectonic 0.17.0 with the compatible resource bundle:

~~~sh
tectonic --untrusted manuscript.tex
~~~

Use a separate output directory for a verification build to preserve the included reviewed PDF. Bibliography page and link numbering are verified against the review build. PDF metadata and timestamps may vary.

## Exact certificate format

[certificate/three-star-tape.txt](certificate/three-star-tape.txt) preserves the frozen base64 tape between explicit markers. Its decoded 1,854 bytes encode the preorder traversal of 300 binary subdivision trees on nondominant chamber–matching cases. The manuscript appendix specifies the ordering and the split or active-set interpretation of every byte. The two unchanged exact verifiers consume that minimal tape directly. The artifact manifest verifies the wrapper and the decoded sequence separately.

## Citation and permissions

The technical paper's complete bibliography is supplied alongside the LaTeX. [CITATION.bib](CITATION.bib) identifies the repository materials. Cite a commit for source-specific results. See the [license scope and file inventory](LICENSE_SCOPE.md) and [third-party notices](THIRD_PARTY_NOTICES.md) for component-specific terms.

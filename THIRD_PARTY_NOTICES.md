# Third-party notices

## Manuscript PDF

The PDF embeds subsets of 19 Type 1 fonts. Latin Modern uses the GUST Font License, which references the LaTeX Project Public License 1.3c or later. Latin Modern is by Bogusław Jackowski and Janusz M. Nowacki on behalf of TeX Users Groups, based on Donald E. Knuth's Computer Modern. The Computer Modern and AMS mathematical Type 1 fonts used here are distributed by the American Mathematical Society under SIL OFL 1.1. The document uses these fonts for rendering. The repository supplies no separate font files or font package, and the manuscript's CC BY terms apply to its original expression. [Latin Modern at CTAN](https://ctan.org/pkg/lm) · [GUST Font License](LICENSES/third-party/Latin-Modern-GUST.txt) · [LPPL 1.3c](LICENSES/third-party/LPPL-1.3c.txt) · [AMSFonts at CTAN](https://ctan.org/pkg/amsfonts) · [AMSFonts copyright and OFL terms](LICENSES/third-party/AMSFonts-OFL-1.1.txt)

## Figures

The three SVGs use DejaVu Sans letterforms as paths. The seven-weight slice also uses a glyph from STIX Size One Symbols. The figures include glyph outlines used for rendering, with no standalone font file. DejaVu retains the applicable Bitstream Vera and Arev notices and the public-domain statement for DejaVu changes. STIX retains the copyright notices and SIL OFL 1.1 terms in the Matplotlib distribution. These notices do not expand the project's grant over third-party font material. [DejaVu notice](LICENSES/third-party/DejaVu.txt) · [STIX notice](LICENSES/third-party/STIX.txt) · [OFL guidance on embedding and artwork](https://openfontlicense.org/ofl-faq/)

## Build and CI dependencies

Python, the figure-rendering dependencies and the TeX toolchain are obtained separately. GitHub Actions fetches the two pinned workflow actions from their respective repositories. Those components retain their upstream terms and are not bundled with the project. [Figure requirements](requirements-figures.txt) · [Matplotlib terms](https://matplotlib.org/stable/project/license.html) · [CI workflow](.github/workflows/verify.yml) · [actions/checkout](https://github.com/actions/checkout) · [actions/setup-python](https://github.com/actions/setup-python)

## Cited work

The manuscript credits the relevant mathematical literature. The repository includes bibliographic metadata and links to cited works. The licenses here grant no rights in the cited works themselves. [Manuscript bibliography](manuscript/references.bib)

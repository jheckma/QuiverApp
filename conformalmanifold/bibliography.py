"""Bibliography of the papers behind the physics in this app.

Citation data pulled verbatim from INSPIRE-HEP (https://inspirehep.net,
`/api/arxiv/<eprint>` and its `?format=bibtex` serializer) on 2026-07-06.
Entries are ordered roughly by the app cards they support: the orbifold tab
(quiver construction, conformal manifold, superconformal data), the toric tab
(webs/diagrams, named geometries, R-charges, brane tilings and the inverse
algorithm, Seiberg duality), and the AdS6/CFT5 card (5d fixed points,
geometric engineering, defect groups, SymTFT, the AdS6 duals).

Each entry carries paper-style display fields, the INSPIRE BibTeX record, and
helpers producing a LaTeX `thebibliography` block.  Served via
/api/bibliography.

To add a paper: fetch https://inspirehep.net/api/arxiv/<id>?format=bibtex
(and the JSON for the display fields) and append an entry below.
"""

from __future__ import annotations

ENTRIES = [
    {
        "key": 'Douglas:1996sw',
        "authors": 'M. R. Douglas and G. W. Moore',
        "title": 'D-branes, quivers, and ALE instantons',
        "journal": '',
        "volume": '',
        "year": '',
        "pages": '',
        "eprint": 'hep-th/9603167',
        "category": 'hep-th',
        "doi": '',
        "role": 'D-branes at orbifold singularities and their quiver gauge theories',
        "bibtex": '@article{Douglas:1996sw,\n    author = "Douglas, Michael R. and Moore, Gregory W.",\n    title = "{D-branes, quivers, and ALE instantons}",\n    eprint = "hep-th/9603167",\n    archivePrefix = "arXiv",\n    reportNumber = "RU-96-15, YCTP-P5-96",\n    month = "3",\n    year = "1996"\n}',
    },
    {
        "key": 'Kachru:1998ys',
        "authors": 'S. Kachru and E. Silverstein',
        "title": '4-D conformal theories and strings on orbifolds',
        "journal": 'Phys.Rev.Lett.',
        "volume": '80',
        "year": '1998',
        "pages": '4855',
        "eprint": 'hep-th/9802183',
        "category": 'hep-th',
        "doi": '10.1103/PhysRevLett.80.4855',
        "role": '4d conformal field theories from string orbifolds',
        "bibtex": '@article{Kachru:1998ys,\n    author = "Kachru, Shamit and Silverstein, Eva",\n    title = "{4-D conformal theories and strings on orbifolds}",\n    eprint = "hep-th/9802183",\n    archivePrefix = "arXiv",\n    reportNumber = "SLAC-PUB-7756, LBL-41440, LBNL-41440, UCB-PTH-98-12",\n    doi = "10.1103/PhysRevLett.80.4855",\n    journal = "Phys. Rev. Lett.",\n    volume = "80",\n    pages = "4855--4858",\n    year = "1998"\n}',
    },
    {
        "key": 'Lawrence:1998ja',
        "authors": 'A. E. Lawrence, N. Nekrasov and C. Vafa',
        "title": 'On conformal field theories in four-dimensions',
        "journal": 'Nucl.Phys.B',
        "volume": '533',
        "year": '1998',
        "pages": '199',
        "eprint": 'hep-th/9803015',
        "category": 'hep-th',
        "doi": '10.1016/S0550-3213(98)00495-7',
        "role": 'large-N orbifold projections and quiver conformal field theories',
        "bibtex": '@article{Lawrence:1998ja,\n    author = "Lawrence, Albion E. and Nekrasov, Nikita and Vafa, Cumrun",\n    title = "{On conformal field theories in four-dimensions}",\n    eprint = "hep-th/9803015",\n    archivePrefix = "arXiv",\n    reportNumber = "HUTP-98-A015, ITEP-TH-15-98",\n    doi = "10.1016/S0550-3213(98)00495-7",\n    journal = "Nucl. Phys. B",\n    volume = "533",\n    pages = "199--209",\n    year = "1998"\n}',
    },
    {
        "key": 'Leigh:1995ep',
        "authors": 'R. G. Leigh and M. J. Strassler',
        "title": 'Exactly marginal operators and duality in four-dimensional N=1 supersymmetric gauge theory',
        "journal": 'Nucl.Phys.B',
        "volume": '447',
        "year": '1995',
        "pages": '95',
        "eprint": 'hep-th/9503121',
        "category": 'hep-th',
        "doi": '10.1016/0550-3213(95)00261-P',
        "role": 'exactly marginal operators in 4d N=1 gauge theories',
        "bibtex": '@article{Leigh:1995ep,\n    author = "Leigh, Robert G. and Strassler, Matthew J.",\n    title = "{Exactly marginal operators and duality in four-dimensional N=1 supersymmetric gauge theory}",\n    eprint = "hep-th/9503121",\n    archivePrefix = "arXiv",\n    reportNumber = "RU-95-2",\n    doi = "10.1016/0550-3213(95)00261-P",\n    journal = "Nucl. Phys. B",\n    volume = "447",\n    pages = "95--136",\n    year = "1995"\n}',
    },
    {
        "key": 'Green:2010da',
        "authors": 'D. Green, Z. Komargodski, N. Seiberg, Y. Tachikawa and B. Wecht',
        "title": 'Exactly Marginal Deformations and Global Symmetries',
        "journal": 'JHEP',
        "volume": '06',
        "year": '2010',
        "pages": '106',
        "eprint": '1005.3546',
        "category": 'hep-th',
        "doi": '10.1007/JHEP06(2010)106',
        "role": 'the general counting of exactly marginal couplings (conformal manifolds)',
        "bibtex": '@article{Green:2010da,\n    author = "Green, Daniel and Komargodski, Zohar and Seiberg, Nathan and Tachikawa, Yuji and Wecht, Brian",\n    title = "{Exactly Marginal Deformations and Global Symmetries}",\n    eprint = "1005.3546",\n    archivePrefix = "arXiv",\n    primaryClass = "hep-th",\n    doi = "10.1007/JHEP06(2010)106",\n    journal = "JHEP",\n    volume = "06",\n    pages = "106",\n    year = "2010"\n}',
    },
    {
        "key": 'Intriligator:2003jj',
        "authors": 'K. A. Intriligator and B. Wecht',
        "title": 'The Exact superconformal R symmetry maximizes a',
        "journal": 'Nucl.Phys.B',
        "volume": '667',
        "year": '2003',
        "pages": '183',
        "eprint": 'hep-th/0304128',
        "category": 'hep-th',
        "doi": '10.1016/S0550-3213(03)00459-0',
        "role": 'a-maximization: the exact superconformal R-symmetry',
        "bibtex": '@article{Intriligator:2003jj,\n    author = "Intriligator, Kenneth A. and Wecht, Brian",\n    title = "{The Exact superconformal R symmetry maximizes a}",\n    eprint = "hep-th/0304128",\n    archivePrefix = "arXiv",\n    reportNumber = "UCSD-PTH-03-02",\n    doi = "10.1016/S0550-3213(03)00459-0",\n    journal = "Nucl. Phys. B",\n    volume = "667",\n    pages = "183--200",\n    year = "2003"\n}',
    },
    {
        "key": 'Aharony:1997bh',
        "authors": 'O. Aharony, A. Hanany and B. Kol',
        "title": 'Webs of (p,q) five-branes, five-dimensional field theories and grid diagrams',
        "journal": 'JHEP',
        "volume": '01',
        "year": '1998',
        "pages": '002',
        "eprint": 'hep-th/9710116',
        "category": 'hep-th',
        "doi": '10.1088/1126-6708/1998/01/002',
        "role": '(p,q) 5-brane webs, 5d field theories, and their grid (toric) diagrams',
        "bibtex": '@article{Aharony:1997bh,\n    author = "Aharony, Ofer and Hanany, Amihay and Kol, Barak",\n    title = "{Webs of (p,q) five-branes, five-dimensional field theories and grid diagrams}",\n    eprint = "hep-th/9710116",\n    archivePrefix = "arXiv",\n    reportNumber = "IASSNS-HEP-97-113, RU-97-81, SU-ITP-97-40",\n    doi = "10.1088/1126-6708/1998/01/002",\n    journal = "JHEP",\n    volume = "01",\n    pages = "002",\n    year = "1998"\n}',
    },
    {
        "key": 'Klebanov:1998hh',
        "authors": 'I. R. Klebanov and E. Witten',
        "title": 'Superconformal field theory on three-branes at a Calabi-Yau singularity',
        "journal": 'Nucl.Phys.B',
        "volume": '536',
        "year": '1998',
        "pages": '199',
        "eprint": 'hep-th/9807080',
        "category": 'hep-th',
        "doi": '10.1016/S0550-3213(98)00654-3',
        "role": 'the conifold gauge theory: D3-branes at a Calabi-Yau singularity',
        "bibtex": '@article{Klebanov:1998hh,\n    author = "Klebanov, Igor R. and Witten, Edward",\n    title = "{Superconformal field theory on three-branes at a Calabi-Yau singularity}",\n    eprint = "hep-th/9807080",\n    archivePrefix = "arXiv",\n    reportNumber = "IASSNS-HEP-98-64, PUPT-1804",\n    doi = "10.1016/S0550-3213(98)00654-3",\n    journal = "Nucl. Phys. B",\n    volume = "536",\n    pages = "199--218",\n    year = "1998"\n}',
    },
    {
        "key": 'Butti:2005ps',
        "authors": 'A. Butti and A. Zaffaroni',
        "title": 'From toric geometry to quiver gauge theory: The Equivalence of a-maximization and Z-minimization',
        "journal": 'Fortsch.Phys.',
        "volume": '54',
        "year": '2006',
        "pages": '309',
        "eprint": 'hep-th/0512240',
        "category": 'hep-th',
        "doi": '10.1002/prop.200510276',
        "role": 'R-charges from toric data: a-maximization equals Z-minimization',
        "bibtex": '@article{Butti:2005ps,\n    author = "Butti, Agostino and Zaffaroni, Alberto",\n    editor = "Bakas, I. and Lust, D.",\n    title = "{From toric geometry to quiver gauge theory: The Equivalence of a-maximization and Z-minimization}",\n    eprint = "hep-th/0512240",\n    archivePrefix = "arXiv",\n    reportNumber = "BICOCCA-FT-05-27",\n    doi = "10.1002/prop.200510276",\n    journal = "Fortsch. Phys.",\n    volume = "54",\n    pages = "309--316",\n    year = "2006"\n}',
    },
    {
        "key": 'Hanany:2005ve',
        "authors": 'A. Hanany and K. D. Kennaway',
        "title": 'Dimer models and toric diagrams',
        "journal": '',
        "volume": '',
        "year": '',
        "pages": '',
        "eprint": 'hep-th/0503149',
        "category": 'hep-th',
        "doi": '',
        "role": 'dimer models for toric quiver gauge theories',
        "bibtex": '@article{Hanany:2005ve,\n    author = "Hanany, Amihay and Kennaway, Kristian D.",\n    title = "{Dimer models and toric diagrams}",\n    eprint = "hep-th/0503149",\n    archivePrefix = "arXiv",\n    reportNumber = "MIT-CTP-3613",\n    month = "3",\n    year = "2005"\n}',
    },
    {
        "key": 'Franco:2005rj',
        "authors": 'S. Franco, A. Hanany, K. D. Kennaway, D. Vegh and B. Wecht',
        "title": 'Brane dimers and quiver gauge theories',
        "journal": 'JHEP',
        "volume": '01',
        "year": '2006',
        "pages": '096',
        "eprint": 'hep-th/0504110',
        "category": 'hep-th',
        "doi": '10.1088/1126-6708/2006/01/096',
        "role": 'brane tilings: the dimer / quiver / toric dictionary',
        "bibtex": '@article{Franco:2005rj,\n    author = "Franco, Sebastian and Hanany, Amihay and Kennaway, Kristian D. and Vegh, David and Wecht, Brian",\n    title = "{Brane dimers and quiver gauge theories}",\n    eprint = "hep-th/0504110",\n    archivePrefix = "arXiv",\n    reportNumber = "MIT-CTP-3619",\n    doi = "10.1088/1126-6708/2006/01/096",\n    journal = "JHEP",\n    volume = "01",\n    pages = "096",\n    year = "2006"\n}',
    },
    {
        "key": 'Gulotta:2008ef',
        "authors": 'D. R. Gulotta',
        "title": 'Properly ordered dimers, R-charges, and an efficient inverse algorithm',
        "journal": 'JHEP',
        "volume": '10',
        "year": '2008',
        "pages": '014',
        "eprint": '0807.3012',
        "category": 'hep-th',
        "doi": '10.1088/1126-6708/2008/10/014',
        "role": 'an efficient inverse algorithm: from toric diagram to brane tiling',
        "bibtex": '@article{Gulotta:2008ef,\n    author = "Gulotta, Daniel R.",\n    title = "{Properly ordered dimers, R-charges, and an efficient inverse algorithm}",\n    eprint = "0807.3012",\n    archivePrefix = "arXiv",\n    primaryClass = "hep-th",\n    reportNumber = "PUPT-2273",\n    doi = "10.1088/1126-6708/2008/10/014",\n    journal = "JHEP",\n    volume = "10",\n    pages = "014",\n    year = "2008"\n}',
    },
    {
        "key": 'Seiberg:1994pq',
        "authors": 'N. Seiberg',
        "title": 'Electric - magnetic duality in supersymmetric nonAbelian gauge theories',
        "journal": 'Nucl.Phys.B',
        "volume": '435',
        "year": '1995',
        "pages": '129',
        "eprint": 'hep-th/9411149',
        "category": 'hep-th',
        "doi": '10.1016/0550-3213(94)00023-8',
        "role": 'Seiberg duality of N=1 gauge theories',
        "bibtex": '@article{Seiberg:1994pq,\n    author = "Seiberg, N.",\n    title = "{Electric - magnetic duality in supersymmetric nonAbelian gauge theories}",\n    eprint = "hep-th/9411149",\n    archivePrefix = "arXiv",\n    reportNumber = "RU-94-82, IASSNS-HEP-94-98",\n    doi = "10.1016/0550-3213(94)00023-8",\n    journal = "Nucl. Phys. B",\n    volume = "435",\n    pages = "129--146",\n    year = "1995"\n}',
    },
    {
        "key": 'Seiberg:1996bd',
        "authors": 'N. Seiberg',
        "title": 'Five-dimensional SUSY field theories, nontrivial fixed points and string dynamics',
        "journal": 'Phys.Lett.B',
        "volume": '388',
        "year": '1996',
        "pages": '753',
        "eprint": 'hep-th/9608111',
        "category": 'hep-th',
        "doi": '10.1016/S0370-2693(96)01215-4',
        "role": 'interacting 5d fixed points with E_n global symmetry',
        "bibtex": '@article{Seiberg:1996bd,\n    author = "Seiberg, Nathan",\n    title = "{Five-dimensional SUSY field theories, nontrivial fixed points and string dynamics}",\n    eprint = "hep-th/9608111",\n    archivePrefix = "arXiv",\n    reportNumber = "RU-96-69",\n    doi = "10.1016/S0370-2693(96)01215-4",\n    journal = "Phys. Lett. B",\n    volume = "388",\n    pages = "753--760",\n    year = "1996"\n}',
    },
    {
        "key": 'Intriligator:1997pq',
        "authors": 'K. A. Intriligator, D. R. Morrison and N. Seiberg',
        "title": 'Five-dimensional supersymmetric gauge theories and degenerations of Calabi-Yau spaces',
        "journal": 'Nucl.Phys.B',
        "volume": '497',
        "year": '1997',
        "pages": '56',
        "eprint": 'hep-th/9702198',
        "category": 'hep-th',
        "doi": '10.1016/S0550-3213(97)00279-4',
        "role": '5d gauge theories from degenerations of Calabi-Yau spaces',
        "bibtex": '@article{Intriligator:1997pq,\n    author = "Intriligator, Kenneth A. and Morrison, David R. and Seiberg, Nathan",\n    title = "{Five-dimensional supersymmetric gauge theories and degenerations of Calabi-Yau spaces}",\n    eprint = "hep-th/9702198",\n    archivePrefix = "arXiv",\n    reportNumber = "RU-96-99, IASSNS-HEP-96-112",\n    doi = "10.1016/S0550-3213(97)00279-4",\n    journal = "Nucl. Phys. B",\n    volume = "497",\n    pages = "56--100",\n    year = "1997"\n}',
    },
    {
        "key": 'Morrison:2020ool',
        "authors": 'D. R. Morrison, S. Schafer-Nameki and B. Willett',
        "title": 'Higher-Form Symmetries in 5d',
        "journal": 'JHEP',
        "volume": '09',
        "year": '2020',
        "pages": '024',
        "eprint": '2005.12296',
        "category": 'hep-th',
        "doi": '10.1007/JHEP09(2020)024',
        "role": 'higher-form symmetries of 5d theories from Calabi-Yau intersection data',
        "bibtex": '@article{Morrison:2020ool,\n    author = "Morrison, David R. and Schafer-Nameki, Sakura and Willett, Brian",\n    title = "{Higher-Form Symmetries in 5d}",\n    eprint = "2005.12296",\n    archivePrefix = "arXiv",\n    primaryClass = "hep-th",\n    doi = "10.1007/JHEP09(2020)024",\n    journal = "JHEP",\n    volume = "09",\n    pages = "024",\n    year = "2020"\n}',
    },
    {
        "key": 'Albertini:2020mdx',
        "authors": 'F. Albertini, M. Del Zotto, I. García Etxebarria and S. S. Hosseini',
        "title": 'Higher Form Symmetries and M-theory',
        "journal": 'JHEP',
        "volume": '12',
        "year": '2020',
        "pages": '203',
        "eprint": '2005.12831',
        "category": 'hep-th',
        "doi": '10.1007/JHEP12(2020)203',
        "role": 'defect groups and linking pairings of 5d theories from M-theory',
        "bibtex": '@article{Albertini:2020mdx,\n    author = "Albertini, Federica and Del Zotto, Michele and Garc{\\\'\\i}a Etxebarria, I{\\~n}aki and Hosseini, Saghar S.",\n    title = "{Higher Form Symmetries and M-theory}",\n    eprint = "2005.12831",\n    archivePrefix = "arXiv",\n    primaryClass = "hep-th",\n    doi = "10.1007/JHEP12(2020)203",\n    journal = "JHEP",\n    volume = "12",\n    pages = "203",\n    year = "2020"\n}',
    },
    {
        "key": 'Apruzzi:2021nmk',
        "authors": 'F. Apruzzi, F. Bonetti, I. García Etxebarria, S. S. Hosseini and S. Schafer-Nameki',
        "title": 'Symmetry TFTs from String Theory',
        "journal": 'Commun.Math.Phys.',
        "volume": '402',
        "year": '2023',
        "pages": '895',
        "eprint": '2112.02092',
        "category": 'hep-th',
        "doi": '10.1007/s00220-023-04737-2',
        "role": 'symmetry TFTs from reducing string theory on the boundary geometry',
        "bibtex": '@article{Apruzzi:2021nmk,\n    author = "Apruzzi, Fabio and Bonetti, Federico and Garc{\\\'\\i}a Etxebarria, I{\\~n}aki and Hosseini, Saghar S. and Schafer-Nameki, Sakura",\n    title = "{Symmetry TFTs from String Theory}",\n    eprint = "2112.02092",\n    archivePrefix = "arXiv",\n    primaryClass = "hep-th",\n    doi = "10.1007/s00220-023-04737-2",\n    journal = "Commun. Math. Phys.",\n    volume = "402",\n    number = "1",\n    pages = "895--949",\n    year = "2023"\n}',
    },
    {
        "key": 'DHoker:2016ujz',
        "authors": "E. D'Hoker, M. Gutperle, A. Karch and C. F. Uhlemann",
        "title": 'Warped $AdS_6\\times S^2$ in Type IIB supergravity I: Local solutions',
        "journal": 'JHEP',
        "volume": '08',
        "year": '2016',
        "pages": '046',
        "eprint": '1606.01254',
        "category": 'hep-th',
        "doi": '10.1007/JHEP08(2016)046',
        "role": 'warped AdS_6 x S^2 solutions of Type IIB supergravity: gravity duals of 5d SCFTs',
        "bibtex": '@article{DHoker:2016ujz,\n    author = "D\'Hoker, Eric and Gutperle, Michael and Karch, Andreas and Uhlemann, Christoph F.",\n    title = "{Warped $AdS_6\\times S^2$ in Type IIB supergravity I: Local solutions}",\n    eprint = "1606.01254",\n    archivePrefix = "arXiv",\n    primaryClass = "hep-th",\n    doi = "10.1007/JHEP08(2016)046",\n    journal = "JHEP",\n    volume = "08",\n    pages = "046",\n    year = "2016"\n}',
    },
]


def _tex_authors(display: str) -> str:
    """'E. D'Hoker, M. Gutperle and A. Karch' -> LaTeX with non-breaking ~."""
    return display.replace(". ", ".~")


def latex_bibitem(e: dict) -> str:
    """One paper-style \\bibitem, JHEP-review style (arXiv-only if unpublished)."""
    pub = (f"{e['journal']} \\textbf{{{e['volume']}}} ({e['year']}) {e['pages']}, "
           if e["journal"] else "")
    return (f"\\bibitem{{{e['key']}}}\n"
            f"{_tex_authors(e['authors'])},\n"
            f"``{e['title']},''\n"
            f"{pub}[arXiv:{e['eprint']} [{e['category']}]].")


def latex_bibliography() -> str:
    """A complete paste-ready thebibliography block."""
    items = "\n\n".join(latex_bibitem(e) for e in ENTRIES)
    return (f"\\begin{{thebibliography}}{{{len(ENTRIES)}}}\n\n"
            f"{items}\n\n\\end{{thebibliography}}")


def bibtex_all() -> str:
    """All INSPIRE BibTeX records, concatenated."""
    return "\n\n".join(e["bibtex"] for e in ENTRIES)


def entries_json() -> dict:
    """JSON payload for the web UI."""
    return {
        "fetched_from": "INSPIRE-HEP, 2026-07-06",
        "entries": [{**e, "latex": latex_bibitem(e)} for e in ENTRIES],
        "bibtex_all": bibtex_all(),
        "latex_all": latex_bibliography(),
    }

"""Bibliography: INSPIRE-sourced entries, paper-style LaTeX, JSON payload."""

from conformalmanifold import bibliography as B


def test_entries_well_formed():
    assert len(B.ENTRIES) == 21
    keys = [e["key"] for e in B.ENTRIES]
    assert len(set(keys)) == len(keys)                     # unique texkeys
    for e in B.ENTRIES:
        for field in ("key", "authors", "title", "eprint", "category",
                      "bibtex", "role"):
            assert e[field], f"{e['key']}: empty {field}"
        assert e["key"] in e["bibtex"]                     # bibtex matches entry
        assert e["eprint"] in e["bibtex"]


def test_card_ordering_and_captions():
    keys = [e["key"] for e in B.ENTRIES]
    # orbifold-tab construction papers first, AdS6 duals last
    assert keys[0] == "Douglas:1996sw"
    assert keys[-1] == "DHoker:2016ujz"
    # the 5d/AdS6 block sits together at the end, in concept order
    tail = keys[-8:]
    assert tail == ["Seiberg:1996bd", "Intriligator:1997pq", "Morrison:2020ool",
                    "Albertini:2020mdx", "Apruzzi:2021nmk", "Gukov:2020btk",
                    "BenettiGenolini:2020doj", "DHoker:2016ujz"]
    # captions are physicist-facing: no internal project language
    for e in B.ENTRIES:
        low = e["role"].lower()
        for banned in ("fived", "beat", "anchor paper", "app", "module", "card"):
            assert banned not in low, f"{e['key']}: internal language in caption"


def test_unpublished_entries_render_arxiv_only():
    dm = next(e for e in B.ENTRIES if e["key"] == "Douglas:1996sw")
    assert dm["journal"] == ""
    item = B.latex_bibitem(dm)
    assert "textbf" not in item                            # no fake journal ref
    assert "[arXiv:hep-th/9603167 [hep-th]]" in item


def test_anchor_entries_present():
    by_key = {e["key"]: e for e in B.ENTRIES}
    dgku = by_key["DHoker:2016ujz"]                        # the beat's anchor
    assert (dgku["journal"], dgku["volume"], dgku["year"], dgku["pages"]) == \
        ("JHEP", "08", "2016", "046")
    assert "Albertini:2020mdx" in by_key                   # defect-group basis
    assert by_key["Apruzzi:2021nmk"]["journal"] == "Commun.Math.Phys."


def test_latex_bibitem_paper_style():
    e = next(x for x in B.ENTRIES if x["key"] == "DHoker:2016ujz")
    item = B.latex_bibitem(e)
    assert item.startswith("\\bibitem{DHoker:2016ujz}")
    assert "``" in item and ",''" in item                  # quoted title
    assert "\\textbf{08} (2016) 046" in item               # journal ref
    assert "[arXiv:1606.01254 [hep-th]]" in item


def test_latex_bibliography_block():
    block = B.latex_bibliography()
    assert block.startswith("\\begin{thebibliography}")
    assert block.rstrip().endswith("\\end{thebibliography}")
    for e in B.ENTRIES:
        assert f"\\bibitem{{{e['key']}}}" in block


def test_entries_json_payload():
    j = B.entries_json()
    assert len(j["entries"]) == 21
    assert all("latex" in e and "bibtex" in e for e in j["entries"])
    assert j["bibtex_all"].count("@article") == 21
    assert "INSPIRE" in j["fetched_from"]

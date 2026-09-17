"""Render a narrative's markdown-subset text (headings, bold/italic, bullet lists,
GitHub-flavored pipe tables) to PDF bytes, for the Insights tab's export buttons.

Deliberately NOT a general markdown-to-PDF converter: it only understands the
constructs our own narrative text actually produces (see agents/narrative.py's
system prompt and deterministic_narrative()/comparative_table_markdown()), so it
stays small and dependency-light (fpdf2 only -- no system libraries like Cairo/
Pango, unlike weasyprint) rather than pulling in a second markdown-parsing
dependency for a bounded, self-controlled vocabulary.
"""
from __future__ import annotations

import re
from datetime import datetime, timezone
from io import BytesIO

from fpdf import FPDF

_HEADING_RE = re.compile(r"^(#{1,3})\s+(.*)$")
_BULLET_RE = re.compile(r"^[-*]\s+(.*)$")
_TABLE_ROW_RE = re.compile(r"^\|(.+)\|\s*$")
_TABLE_SEP_RE = re.compile(r"^\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)+\|?$")

_NAVY = (0x17, 0x32, 0x4F)
_CRIMSON = (0xD1, 0x19, 0x47)
_MUTED = (0x5B, 0x6B, 0x75)
_LINE = (0xD7, 0xE0, 0xE6)

# fpdf2's core Helvetica font is Latin-1 only; our narrative text (from the LLM and our
# own deterministic formatting) routinely uses characters outside that range (em/en
# dashes, Greek delta, arrows, curly quotes, comparison operators, bullets). Rather than
# bundling a Unicode TTF as a new asset/dependency, map the common cases to a clean ASCII
# equivalent, then let a `latin-1` replace-pass catch anything unexpected instead of
# crashing the export.
_UNICODE_MAP = {
    "\u2014": "--", "\u2013": "-", "\u2212": "-",
    "\u2018": "'", "\u2019": "'", "\u201c": '"', "\u201d": '"',
    "\u2026": "...", "\u2022": "-", "\u00b7": "-",
    "\u2192": "->", "\u2190": "<-", "\u2194": "<->",
    "\u2265": ">=", "\u2264": "<=", "\u2260": "!=",
    "\u0394": "Delta", "\u03b1": "alpha", "\u03b2": "beta", "\u03b3": "gamma",
    "\u00b1": "+/-", "\u00d7": "x", "\u00b0": "deg",
    "\U0001f916": "[LLM]",
}


def _sanitize(text: str) -> str:
    for uni, ascii_ in _UNICODE_MAP.items():
        text = text.replace(uni, ascii_)
    return text.encode("latin-1", errors="replace").decode("latin-1")


def _split_inline(text: str) -> list[tuple[str, bool, bool, bool]]:
    """Split a line into (text, bold, italic, highlight) runs for **bold**,
    *italic*/_italic_, and ==highlight== (a bounded, app-defined syntax -- not
    standard markdown -- for the LLM to mark its single most important
    takeaway per gene; see gene_narrative_llm()'s system prompt). The LLM never
    emits raw HTML/color codes; this is the only "emphasis" channel beyond
    bold/italic, and it is parsed and rendered entirely by this app's own code.
    """
    text = _sanitize(text)
    runs: list[tuple[str, bool, bool, bool]] = []
    pattern = re.compile(r"(==.+?==|\*\*.+?\*\*|\*.+?\*|_.+?_)")
    pos = 0
    for m in pattern.finditer(text):
        if m.start() > pos:
            runs.append((text[pos:m.start()], False, False, False))
        tok = m.group(0)
        if tok.startswith("=="):
            runs.append((tok[2:-2], False, False, True))
        elif tok.startswith("**"):
            runs.append((tok[2:-2], True, False, False))
        elif tok.startswith("*"):
            runs.append((tok[1:-1], False, True, False))
        else:
            runs.append((tok[1:-1], False, True, False))
        pos = m.end()
    if pos < len(text):
        runs.append((text[pos:], False, False, False))
    return runs or [(text, False, False, False)]


class _NarrativePDF(FPDF):
    def footer(self) -> None:
        self.set_y(-12)
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(*_MUTED)
        self.cell(0, 8, _sanitize(f"SCRAP-AI \u00b7 page {self.page_no()}"), align="C")


def _write_inline(pdf: FPDF, line: str, size: float = 10.5) -> None:
    pdf.set_font("Helvetica", "", size)
    pdf.set_text_color(0x1A, 0x2A, 0x33)
    for chunk, bold, italic, highlight in _split_inline(line):
        style = ("B" if bold else "") + ("I" if italic else "")
        pdf.set_font("Helvetica", style, size)
        pdf.set_text_color(*_CRIMSON) if highlight else pdf.set_text_color(0x1A, 0x2A, 0x33)
        pdf.write(6, chunk)
    pdf.ln(7)


def _parse_table_block(lines: list[str], i: int) -> tuple[list[str], list[list[str]], int]:
    """lines[i] is a header row, lines[i+1] the GFM separator. Returns (headers, rows, next_i)."""
    def cells(row: str) -> list[str]:
        inner = row.strip()
        if inner.startswith("|"):
            inner = inner[1:]
        if inner.endswith("|"):
            inner = inner[:-1]
        return [c.strip().replace("\\|", "|") for c in inner.split("|")]

    headers = cells(lines[i])
    j = i + 2
    rows: list[list[str]] = []
    while j < len(lines) and _TABLE_ROW_RE.match(lines[j].strip()):
        rows.append(cells(lines[j]))
        j += 1
    return headers, rows, j


def _render_table(pdf: FPDF, headers: list[str], rows: list[list[str]]) -> None:
    pdf.ln(2)
    n_cols = len(headers)
    # size relative column widths by the longest content in each column (header or any
    # cell), with a floor -- avoids the default equal-width split producing awkward
    # mid-word wraps on wide reports like the 11-column comparative summary table.
    widths = []
    for i, h in enumerate(headers):
        longest = len(str(h))
        for r in rows:
            if i < len(r):
                longest = max(longest, len(str(r[i])))
        widths.append(max(10, longest))
    font_size = 9.5 if n_cols <= 6 else 7.0
    prev_size = pdf.font_size_pt
    pdf.set_font("Helvetica", "", font_size)
    with pdf.table(
        text_align="LEFT",
        line_height=5.0,
        v_align="TOP",
        borders_layout="ALL",
        col_widths=tuple(widths),
    ) as table:
        header_row = table.row()
        for h in headers:
            header_row.cell(_sanitize(str(h)))
        for r in rows:
            row = table.row()
            for c in r:
                row.cell(_sanitize(str(c)))
    pdf.set_font("Helvetica", "", prev_size)
    pdf.ln(3)


def markdown_to_pdf_bytes(title: str, subtitle: str, body_markdown: str,
                          images: list[tuple[str, bytes]] | None = None) -> bytes:
    """Render `body_markdown` (our narrative subset) into a titled PDF, return bytes.
    `images`, if given, is a list of (caption, png_bytes) pairs -- rasterized Plotly
    figures reused from other tabs (see streamlit_app.py's _fig_to_png_bytes()) --
    appended as a captioned gallery after the main body.
    """
    pdf = _NarrativePDF(format="Letter")
    pdf.set_auto_page_break(auto=True, margin=16)
    pdf.set_margins(16, 14, 16)
    pdf.add_page()

    pdf.set_font("Helvetica", "B", 17)
    pdf.set_text_color(*_CRIMSON)
    pdf.write(9, "SCRAP-AI")
    pdf.ln(9)
    pdf.set_font("Helvetica", "B", 14)
    pdf.set_text_color(*_NAVY)
    pdf.multi_cell(0, 7, _sanitize(title), new_x="LMARGIN", new_y="NEXT")
    if subtitle:
        pdf.set_font("Helvetica", "", 10)
        pdf.set_text_color(*_MUTED)
        pdf.multi_cell(0, 5, _sanitize(subtitle), new_x="LMARGIN", new_y="NEXT")
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    pdf.set_font("Helvetica", "I", 8.5)
    pdf.set_text_color(*_MUTED)
    pdf.multi_cell(0, 5, _sanitize(f"Generated {ts}"), new_x="LMARGIN", new_y="NEXT")
    pdf.set_draw_color(*_LINE)
    pdf.ln(2)
    pdf.line(pdf.l_margin, pdf.get_y(), pdf.w - pdf.r_margin, pdf.get_y())
    pdf.ln(4)

    lines = body_markdown.replace("\r\n", "\n").split("\n")
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        stripped = line.strip()

        if not stripped:
            pdf.ln(2)
            i += 1
            continue

        if (_TABLE_ROW_RE.match(stripped) and i + 1 < n
                and _TABLE_SEP_RE.match(lines[i + 1].strip())):
            headers, rows, nxt = _parse_table_block(lines, i)
            _render_table(pdf, headers, rows)
            i = nxt
            continue

        m = _HEADING_RE.match(stripped)
        if m:
            level = len(m.group(1))
            text = m.group(2).strip()
            size = {1: 15, 2: 13, 3: 11.5}.get(level, 11)
            pdf.ln(2)
            pdf.set_font("Helvetica", "B", size)
            pdf.set_text_color(*_NAVY)
            pdf.multi_cell(0, 7, _sanitize(text), new_x="LMARGIN", new_y="NEXT")
            pdf.ln(1)
            i += 1
            continue

        m = _BULLET_RE.match(stripped)
        if m:
            pdf.set_x(pdf.l_margin + 4)
            pdf.set_font("Helvetica", "", 10.5)
            pdf.set_text_color(0x1A, 0x2A, 0x33)
            pdf.write(6, "-  ")
            for chunk, bold, italic, highlight in _split_inline(m.group(1)):
                style = ("B" if bold else "") + ("I" if italic else "")
                pdf.set_font("Helvetica", style, 10.5)
                pdf.set_text_color(*_CRIMSON) if highlight else pdf.set_text_color(0x1A, 0x2A, 0x33)
                pdf.write(6, chunk)
            pdf.ln(6.5)
            i += 1
            continue

        # italic-only caption line, e.g. "*(deterministic ...)*"
        if stripped.startswith("*(") and stripped.endswith(")*"):
            pdf.set_font("Helvetica", "I", 9)
            pdf.set_text_color(*_MUTED)
            pdf.multi_cell(0, 5, _sanitize(stripped[1:-1]), new_x="LMARGIN", new_y="NEXT")
            i += 1
            continue

        _write_inline(pdf, stripped)
        i += 1

    if images:
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 13)
        pdf.set_text_color(*_NAVY)
        pdf.multi_cell(0, 7, "Plots from other tabs", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(1)
        for caption, png_bytes in images:
            if pdf.get_y() > pdf.h - pdf.b_margin - 90:
                pdf.add_page()
            pdf.set_font("Helvetica", "B", 10.5)
            pdf.set_text_color(0x1A, 0x2A, 0x33)
            pdf.multi_cell(0, 6, _sanitize(caption), new_x="LMARGIN", new_y="NEXT")
            pdf.ln(1)
            try:
                img_w = pdf.w - pdf.l_margin - pdf.r_margin
                pdf.image(BytesIO(png_bytes), w=img_w)
            except Exception as exc:  # never let one bad image sink the whole PDF
                pdf.set_font("Helvetica", "I", 9)
                pdf.set_text_color(*_MUTED)
                pdf.multi_cell(0, 5, _sanitize(f"(could not embed image: {exc})"),
                              new_x="LMARGIN", new_y="NEXT")
            pdf.ln(4)

    return bytes(pdf.output())

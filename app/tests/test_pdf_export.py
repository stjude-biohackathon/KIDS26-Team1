"""Tests for the PDF export helper. Run: uv run pytest app/tests -q"""
from __future__ import annotations

import re
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from pdf_export import markdown_to_pdf_bytes, _split_inline, _sanitize  # noqa: E402


def _font_names(pdf_bytes: bytes) -> set[str]:
    return {m.decode() for m in re.findall(rb"/BaseFont\s*/([A-Za-z,-]+)", pdf_bytes)}


def test_markdown_to_pdf_bytes_is_a_valid_pdf():
    data = markdown_to_pdf_bytes("Title", "Subtitle", "Some plain text.")
    assert data.startswith(b"%PDF-")
    assert b"%%EOF" in data


def test_pdf_renders_headings_bullets_and_bold_italic_as_distinct_fonts():
    md = (
        "## A Heading\n\n"
        "Plain paragraph with **bold** and *italic* text.\n\n"
        "- bullet one\n"
        "- bullet two\n"
    )
    data = markdown_to_pdf_bytes("Report", "Sub", md)
    fonts = _font_names(data)
    # bold/italic must be genuinely distinct embedded font resources, not just
    # visually-identical text -- this is what makes emphasis real in the PDF.
    assert "Helvetica" in fonts
    assert "Helvetica-Bold" in fonts
    assert "Helvetica-Italic" in fonts or "Helvetica-Oblique" in fonts


def test_pdf_renders_a_pipe_table():
    md = (
        "## Comparative Summary Table\n\n"
        "| Gene | Rank |\n"
        "|---|---|\n"
        "| BRD7 | 1/19 |\n"
        "| DNAJC11 | 4/19 |\n"
    )
    data = markdown_to_pdf_bytes("Report", "Sub", md)
    assert data.startswith(b"%PDF-")
    # a real table draws border/cell content streams -- length is a crude but
    # useful signal that *something* beyond bare paragraphs was rendered
    assert len(data) > 1500


def test_unicode_characters_do_not_crash_export():
    """Regression test: fpdf2's core Helvetica font is Latin-1 only. Real narrative
    text (ours and the LLM's) routinely contains em-dashes, Greek delta, arrows,
    etc. -- these must be sanitized, never raise FPDFUnicodeEncodingException."""
    md = "Title with \u2014 em-dash, \u0394 delta, \u2192 arrow, \u2265 gte, \u2022 bullet char inline."
    data = markdown_to_pdf_bytes("Title \u2014 em dash", "Sub \u00b7 dot", md)
    assert data.startswith(b"%PDF-")


def test_sanitize_maps_common_unicode_to_ascii():
    assert _sanitize("a \u2014 b") == "a -- b"
    assert _sanitize("\u0394x") == "Deltax"
    assert "?" not in _sanitize("normal text")


def test_split_inline_bold_and_italic():
    runs = _split_inline("**bold** then *italic* then plain")
    assert runs[0] == ("bold", True, False)
    assert any(r == ("italic", False, True) for r in runs)
    assert runs[-1] == (" then plain", False, False)


def test_empty_body_still_produces_valid_pdf():
    data = markdown_to_pdf_bytes("Title", "", "")
    assert data.startswith(b"%PDF-")

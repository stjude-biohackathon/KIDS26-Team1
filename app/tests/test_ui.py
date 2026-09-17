"""Tests for app/ui.py's presentational helpers. Run: uv run pytest app/tests -q"""
from __future__ import annotations

import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

import ui  # noqa: E402


def test_render_highlighted_html_wraps_double_equals_in_mark():
    out = ui.render_highlighted_html("plain ==important fact== plain")
    assert "<mark" in out
    assert "important fact</mark>" in out
    assert out.startswith("plain ")
    assert out.endswith(" plain")


def test_render_highlighted_html_multiple_spans():
    out = ui.render_highlighted_html("==first== and ==second==")
    assert out.count("<mark") == 2
    assert "first</mark>" in out and "second</mark>" in out


def test_render_highlighted_html_escapes_stray_html_first():
    """The model/deterministic text must never be able to inject arbitrary
    HTML/JS -- only this function's own <mark> substitution should ever
    appear as real HTML."""
    out = ui.render_highlighted_html("<script>alert(1)</script> ==ok==")
    assert "<script>" not in out
    assert "&lt;script&gt;" in out
    assert "<mark" in out  # our own substitution still works alongside the escape


def test_render_highlighted_html_no_highlight_syntax_passthrough():
    out = ui.render_highlighted_html("plain text with **bold** and *italic*")
    assert "<mark" not in out
    assert "**bold**" in out  # markdown bold left untouched for st.markdown to render


def test_render_highlighted_html_escapes_ampersand():
    out = ui.render_highlighted_html("R&D ==finding==")
    assert "R&amp;D" in out

"""Presentational helpers: brand theme (light/dark), header, metric cards,
status chips, sidebar sections, and unified Plotly templates.

Pure UI — no analysis logic. Import and call from streamlit_app.py.
"""
from __future__ import annotations

import base64
import re
from pathlib import Path

import plotly.graph_objects as go
import plotly.io as pio
import streamlit as st

_ASSETS = Path(__file__).resolve().parent.parent / "assets"


def _b64(path: Path) -> str:
    try:
        return base64.b64encode(path.read_bytes()).decode()
    except Exception:
        return ""

# ---- palette ----
NAVY = "#0B2E4F"
BLUE = "#1C5D99"
TEAL = "#1C7293"
SEA = "#00A896"
INK = "#1A2A33"
MUTE = "#5B6B75"

PALETTE = [BLUE, SEA, TEAL, "#E29578", "#8172B3", "#55A868", NAVY, "#C44E52"]

THEMES = {
    "Light": {
        "bg": "#FFFFFF", "surface": "#F4F8FA", "surface2": "#ECF3F7",
        "text": INK, "muted": MUTE, "border": "#D7E0E6",
        "header_from": "#D11947", "header_to": "#D11947", "header_text": "#FFFFFF",
        "header_sub": "#FBD9E1", "header_accent": "#FFE3EA",
        "plot": "kady_light", "plot_bg": "#FFFFFF",
        "grid": "#E6EDF1",
    },
    "Dark": {
        "bg": "#0E1B24", "surface": "#152632", "surface2": "#1B303E",
        "text": "#E6EEF2", "muted": "#9FB3C0", "border": "#243B49",
        "header_from": "#D11947", "header_to": "#D11947", "header_text": "#FFFFFF",
        "header_sub": "#FBD9E1", "header_accent": "#FFE3EA",
        "plot": "kady_dark", "plot_bg": "#152632",
        "grid": "#243B49",
    },
}


def register_plotly_templates() -> None:
    for name, tvals in [("kady_light", THEMES["Light"]), ("kady_dark", THEMES["Dark"])]:
        pio.templates[name] = go.layout.Template(
            layout=dict(
                colorway=PALETTE,
                font=dict(family="Inter, Segoe UI, Helvetica, Arial, sans-serif",
                          color=tvals["text"], size=13),
                paper_bgcolor=tvals["plot_bg"],
                plot_bgcolor=tvals["plot_bg"],
                title=dict(font=dict(size=15, color=tvals["text"])),
                xaxis=dict(gridcolor=tvals["grid"], zerolinecolor=tvals["grid"],
                           linecolor=tvals["border"]),
                yaxis=dict(gridcolor=tvals["grid"], zerolinecolor=tvals["grid"],
                           linecolor=tvals["border"]),
                hoverlabel=dict(font_size=12, font_family="Inter, sans-serif"),
                margin=dict(l=60, r=20, t=48, b=48),
                legend=dict(bgcolor="rgba(0,0,0,0)"),
            )
        )


def current_theme() -> dict:
    return THEMES[st.session_state.get("ui_theme", "Light")]


def force_sidebar_expanded_on_load() -> None:
    """Streamlit persists the sidebar collapse state in the browser's
    localStorage (key `stSidebarCollapsed-`), which silently overrides
    `initial_sidebar_state="expanded"` on every future load in that same
    browser once a user has ever manually collapsed it once. This clears that
    stale flag and, if the sidebar is currently rendered collapsed, clicks the
    same toggle a user would to re-expand it -- going through Streamlit's own
    expand logic rather than fighting its layout with CSS (forcing the sidebar
    open via CSS alone leaves the main content area still sized for a
    collapsed sidebar, causing an overlap).

    Guarded to run only once per real page load (a `window` flag), so it never
    fights a deliberate in-session collapse triggered by later Streamlit
    reruns (e.g. clicking any other widget).
    """
    st.components.v1.html("""
    <script>
    (function() {
      try {
        var w = window.parent;
        if (w.__scrapaiSidebarForced) return;
        w.__scrapaiSidebarForced = true;
        var doc = w.document;
        var sidebar = doc.querySelector('section[data-testid="stSidebar"]');
        if (sidebar && sidebar.getAttribute('aria-expanded') === 'false') {
          var btn = doc.querySelector('[data-testid="stSidebarCollapseButton"] button');
          if (btn) { btn.click(); }
        }
        w.localStorage.removeItem('stSidebarCollapsed-');
      } catch (e) { /* no-op if blocked (e.g. cross-origin embed) */ }
    })();
    </script>
    """, height=0)


def force_dark_mode_text_colors() -> None:
    """Several native Streamlit components hardcode text in the STATIC
    textColor from .streamlit/config.toml (#1A2A33, this app's light-mode ink)
    regardless of this app's own Light/Dark toggle -- confirmed via computed
    styles on stWidgetLabel, stCaptionContainer, stExpander summaries,
    stRadioOption, stTooltipIcon (SVG stroke), stMetricLabel/stMetricValue,
    and sidebar headings. A plain CSS override wins for most of these, but at
    least stMetricLabel's text kept the static color even with a maximally
    specific `!important` rule (an unresolved specificity/ordering fight with
    Streamlit's own injected styles) -- so this uses direct inline styling via
    JS instead, which always wins regardless of stylesheet specificity.

    Runs once immediately, then keeps re-applying via a MutationObserver so it
    also catches elements that render later (switching tabs, expanding a
    section, a fresh Streamlit rerun replacing DOM nodes) -- a plain one-shot
    fix would only catch what already existed at the moment it ran. Only
    active in Dark mode; a full page rerun when switching themes disconnects
    the previous observer before installing a fresh one.
    """
    t = current_theme()
    is_dark = st.session_state.get("ui_theme", "Light") == "Dark"
    st.components.v1.html(f"""
    <script>
    (function() {{
      var w = window.parent;
      if (w.__scrapaiDarkObserver) {{ w.__scrapaiDarkObserver.disconnect(); w.__scrapaiDarkObserver = null; }}
      if (!{str(is_dark).lower()}) return;
      var COLOR = {t['text']!r};
      var COLOR_SELECTORS = [
        '[data-testid="stWidgetLabel"] p', '[data-testid="stWidgetLabel"] span',
        '[data-testid="stCaptionContainer"] p', '[data-testid="stCaptionContainer"] span',
        '[data-testid="stExpander"] summary', '[data-testid="stExpander"] summary p',
        '[data-testid="stExpander"] summary span',
        '[data-testid="stRadioOption"] p', '[data-testid="stRadioOption"] span',
        '[data-testid="stMetricLabel"] p', '[data-testid="stMetricLabel"] span',
        '[data-testid="stMetricValue"]',
        'section[data-testid="stSidebar"] h1', 'section[data-testid="stSidebar"] h2',
        'section[data-testid="stSidebar"] h3', 'section[data-testid="stSidebar"] h4'
      ].join(',');
      var STROKE_SELECTORS = '[data-testid="stTooltipIcon"] svg, [data-testid="stTooltipIcon"] path';
      try {{
        var doc = w.document;
        function fix() {{
          doc.querySelectorAll(COLOR_SELECTORS).forEach(function(el) {{
            if (el.style.getPropertyValue('color') !== COLOR) {{
              el.style.setProperty('color', COLOR, 'important');
            }}
          }});
          doc.querySelectorAll(STROKE_SELECTORS).forEach(function(el) {{
            if (el.style.getPropertyValue('stroke') !== COLOR) {{
              el.style.setProperty('stroke', COLOR, 'important');
            }}
          }});
        }}
        fix();
        var obs = new MutationObserver(fix);
        obs.observe(doc.body, {{childList: true, subtree: true, attributes: true,
                              attributeFilter: ['style', 'class']}});
        w.__scrapaiDarkObserver = obs;
      }} catch (e) {{ /* no-op if blocked */ }}
    }})();
    </script>
    """, height=0)


def apply_theme() -> None:
    """Inject CSS for the active theme and set the default Plotly template."""
    t = current_theme()
    pio.templates.default = t["plot"]
    st.markdown(f"""
    <style>
      @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
      @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@500;600;700&display=swap');
      html, body, [class*="css"] {{ font-family: 'Inter', system-ui, sans-serif; }}
      .stApp {{ background: {t['bg']}; color: {t['text']}; }}
      section[data-testid="stSidebar"] {{ background: {t['surface']}; border-right: 1px solid {t['border']}; }}
      /* Streamlit's native widget labels/captions/expander titles use the
         STATIC textColor from .streamlit/config.toml -- that file has no
         light/dark variants, so it stays fixed at the light-mode ink color
         regardless of this app's own Light/Dark toggle. In Dark mode that's
         near-invisible against the dark sidebar/background (confirmed via
         computed style: #1A2A33 on both a dark and light background). Force
         all three to the active theme's text color instead. */
      [data-testid="stWidgetLabel"] p, [data-testid="stWidgetLabel"] span,
      [data-testid="stCaptionContainer"] p, [data-testid="stCaptionContainer"] span,
      [data-testid="stExpander"] summary, [data-testid="stExpander"] summary span,
      [data-testid="stExpander"] summary p,
      [data-testid="stRadioOption"] p, [data-testid="stRadioOption"] span,
      [data-testid="stMetricLabel"] [data-testid="stMarkdownContainer"] p,
      [data-testid="stMetricLabel"] p, [data-testid="stMetricLabel"] span,
      [data-testid="stMetricValue"]
        {{ color: {t['text']} !important; }}
      /* the sidebar specifically hardcodes its own heading color regardless of
         this app's theme (confirmed: the sidebar's "SCRAP-AI" h3 stayed at the
         static ink color even though an identical h3 in the main content area
         correctly inherited the theme's text color) -- override headings
         scoped to the sidebar only, since main-content headings already work. */
      section[data-testid="stSidebar"] h1, section[data-testid="stSidebar"] h2,
      section[data-testid="stSidebar"] h3, section[data-testid="stSidebar"] h4
        {{ color: {t['text']} !important; }}
      /* help/tooltip (?) icons paint via SVG stroke, not the CSS `color`
         property, so the general text-color fix above doesn't reach them --
         same static light-mode ink color otherwise, invisible on dark. */
      [data-testid="stTooltipIcon"] svg, [data-testid="stTooltipIcon"] path
        {{ stroke: {t['text']} !important; opacity: .7; }}
      /* an OPEN expander's summary gets a hardcoded near-white background from
         Streamlit regardless of theme (collapsed ones are already transparent,
         which is correct) -- force the same transparent look so light text
         stays legible instead of light-on-near-white. */
      [data-testid="stExpander"] summary {{ background: transparent !important; }}
      /* sidebar collapse toggle sits right at the top edge by default -- nudge it
         down a bit for visibility/click comfort (reported cramped near the very
         top-left corner). Streamlit also hides this button until hover by
         default; forcing it permanently visible per the user's explicit ask. */
      [data-testid="stSidebarHeader"] {{ padding-top: 14px; }}
      [data-testid="stSidebarCollapseButton"] {{ visibility: visible !important; }}
      /* toolbarMode=minimal (see .streamlit/config.toml) leaves this header bar
         completely empty -- confirmed no status widget or sidebar-collapse
         control lives inside it -- so collapse it instead of showing blank
         white space above the banner. */
      header[data-testid="stHeader"] {{ height: 0rem; min-height: 0; }}
      .block-container {{ padding-top: 1.5rem; }}
      /* header banner -- IBM Plex Sans is scoped to the banner only, the rest of the app stays Inter. */
      .kady-header {{
        background: linear-gradient(120deg, {t['header_from']} 0%, {t['header_to']} 100%);
        border-radius: 14px; padding: 18px 24px; margin: 2px 0 14px 0;
        box-shadow: 0 6px 22px rgba(11,46,79,0.18);
        font-family: 'IBM Plex Sans', 'Inter', system-ui, sans-serif;
      }}
      .kady-header h1 {{ color: {t['header_text']}; font-size: 24px; font-weight: 700;
        font-family: 'IBM Plex Sans', 'Inter', system-ui, sans-serif;
        margin: 0; letter-spacing: .2px; line-height: 1.25; }}
      .kady-header .kady-brand {{ color: {t['header_accent']}; opacity: .85; font-weight: 600;
        font-size: 15px; letter-spacing: .8px; text-transform: uppercase; margin: 0 0 3px 0;
        font-family: 'IBM Plex Sans', 'Inter', system-ui, sans-serif; }}
      .kady-header p {{ color: {t['header_sub']}; font-size: 14px; margin: 5px 0 0 0;
        font-family: 'IBM Plex Sans', 'Inter', system-ui, sans-serif; }}
      /* metric cards */
      .kady-cards {{ display: flex; gap: 12px; margin: 2px 0 8px 0; flex-wrap: wrap; }}
      .kady-card {{ flex: 1; min-width: 150px; background: {t['surface']};
        border: 1px solid {t['border']}; border-radius: 12px; padding: 12px 16px;
        box-shadow: 0 2px 8px rgba(11,46,79,0.06); }}
      .kady-card .v {{ font-size: 26px; font-weight: 700; color: {SEA}; line-height: 1.1; }}
      .kady-card .l {{ font-size: 11.5px; color: {t['muted']}; margin-top: 2px;
        text-transform: uppercase; letter-spacing: .4px; }}
      /* chips */
      .kady-chips {{ margin: 2px 0 10px 0; }}
      .kady-chip {{ display: inline-block; padding: 3px 11px; border-radius: 999px;
        font-size: 12px; font-weight: 600; margin-right: 8px; margin-bottom: 4px;
        border: 1px solid {t['border']}; background: {t['surface2']}; color: {t['text']}; }}
      .chip-ok {{ background: rgba(0,168,150,0.14); color: {SEA}; border-color: rgba(0,168,150,0.4); }}
      .chip-warn {{ background: rgba(226,149,120,0.16); color: #C46B4E; border-color: rgba(226,149,120,0.5); }}
      .chip-info {{ background: rgba(28,93,153,0.12); color: {BLUE}; border-color: rgba(28,93,153,0.35); }}
      .badge {{ display:inline-block; padding:2px 10px; border-radius:8px; font-size:12px;
        font-weight:600; }}
      .badge-drug {{ background: rgba(0,168,150,0.16); color:{SEA}; }}
      .badge-none {{ background: {t['surface2']}; color:{t['muted']}; }}
      /* LLM narrative card (Insights tab) -- a key-scoped st.container(border=True),
         see ui.card_container(). Left accent bar signals "model-generated", distinct
         from the plain deterministic narrative markdown above it. */
      .st-key-llm-narrative-card {{
        background: {t['surface']} !important; border: 1px solid {t['border']} !important;
        border-left: 4px solid {TEAL} !important; border-radius: 12px !important;
        box-shadow: 0 2px 10px rgba(11,46,79,0.08);
      }}
      .kady-card-title {{ font-size: 12.5px; font-weight: 700; letter-spacing: .5px;
        text-transform: uppercase; color: {TEAL}; margin: 0 0 8px 2px; }}
      /* tabs */
      /* Tab bar: minimal ghost-button style. Verified selectors for this Streamlit
         version (the old button[data-baseweb="tab"] rule was dead -- this version
         renders tabs as div[data-testid="stTab"] with role="tab", not a <button>). */
      div[role="tablist"] {{ gap: 4px; border-bottom: none; }}
      [data-testid="stTab"] {{
        border-radius: 6px; padding: 5px 12px;
        background: transparent; border: none;
        transition: background .15s;
      }}
      [data-testid="stTab"] p {{ color: {t['muted']}; font-weight: 500; font-size: 13px; }}
      [data-testid="stTab"][aria-selected="true"] {{ background: {t['surface2']}; }}
      [data-testid="stTab"][aria-selected="true"] p {{ color: {t['text']}; font-weight: 700; }}
      [data-testid="stTab"] .react-aria-SelectionIndicator {{ display: none; }}
      /* sidebar section label */
      .kady-side {{ font-size: 12px; font-weight: 700; color: {TEAL};
        text-transform: uppercase; letter-spacing: .6px; margin: 6px 0 2px 0; }}
      /* primary button */
      .stButton>button[kind="primary"] {{ background: {TEAL}; border: 0; font-weight: 600; }}
      .stButton>button[kind="primary"]:hover {{ background: {NAVY}; }}
    </style>
    """, unsafe_allow_html=True)


def header(title: str, subtitle: str, attribution: str | None = None) -> None:
    icon = _b64(_ASSETS / "scrap_ai_logo_v27.png")
    org_logo = _b64(_ASSETS / "sjcrh_logo_white.png")
    img = (f'<img src="data:image/png;base64,{icon}" alt="SCRAP-AI logo" '
           'style="height:168px;width:auto;flex:0 0 auto;"/>') if icon else ""
    org_img = (f'<img src="data:image/png;base64,{org_logo}" alt="St. Jude Children\'s '
              'Research Hospital" style="height:84px;width:auto;flex:0 0 auto;'
              'margin-left:auto;padding-left:18px"/>') if org_logo else ""
    attribution_html = (f'<p style="opacity:.72;font-size:14px;margin-top:2px">{attribution}</p>'
                        if attribution else "")
    st.markdown(f"""<div class="kady-header" style="display:flex;align-items:center;gap:18px">
        {img}
        <div style="min-width:0">
          <p class="kady-brand">SCRAP-AI</p>
          <h1>{title}</h1>
          <p>{subtitle}</p>
          {attribution_html}
        </div>
        {org_img}
        </div>""", unsafe_allow_html=True)


def metric_cards(cards: list[tuple[str, str]]) -> None:
    """cards = [(value, label), ...]"""
    html = '<div class="kady-cards">'
    for v, l in cards:
        html += f'<div class="kady-card"><div class="v">{v}</div><div class="l">{l}</div></div>'
    html += "</div>"
    st.markdown(html, unsafe_allow_html=True)


def card_container(key: str, title: str, icon: str = ""):
    """A bordered, branded card (a `key`-scoped `st.container(border=True)`,
    styled via the `.st-key-{key}` hook — see apply_theme()'s `.kady-llm-card`
    rule) with a small title bar. Returns the container; write further content
    into it with a `with` block, e.g.:

        with ui.card_container("llm-narrative-card", "LLM narrative", icon="\U0001f916"):
            st.markdown(text)
    """
    c = st.container(border=True, key=key)
    with c:
        label = f"{icon}&nbsp;&nbsp;{title}" if icon else title
        st.markdown(f'<div class="kady-card-title">{label}</div>', unsafe_allow_html=True)
    return c


_HIGHLIGHT_RE = re.compile(r"==(.+?)==")


def _escape_html(text: str) -> str:
    """Minimal manual HTML escaper (not the stdlib `html` module -- this file
    already uses `html` as a local variable name in a couple of places, so
    importing it module-wide would be a real shadowing footgun)."""
    return (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def render_highlighted_html(text: str) -> str:
    """Pure transform: escape any stray HTML in `text`, then turn this app's
    bounded `==highlighted text==` syntax into a `<mark>` span. Split out from
    markdown_with_highlights() so the substitution logic is directly testable
    without a live Streamlit script-run context. The model/deterministic text
    never emits raw HTML: escaping runs FIRST, so the `<mark>` this function
    inserts is the only HTML that ever reaches the page -- not an
    arbitrary-HTML-injection surface.
    """
    escaped = _escape_html(text)
    return _HIGHLIGHT_RE.sub(
        r'<mark style="background:#FBE9EE;color:#17324F;padding:0 3px;'
        r'border-radius:3px;font-weight:600;">\1</mark>',
        escaped,
    )


def markdown_with_highlights(text: str) -> None:
    """Render narrative markdown, treating a bounded `==highlighted text==`
    syntax as a soft highlight -- this app's ONLY inline-emphasis channel
    beyond standard **bold**/*italic*, defined and parsed entirely by this
    app's own code (see gene_narrative_llm()'s system prompt and
    gene_narrative_deterministic()).
    """
    st.markdown(render_highlighted_html(text), unsafe_allow_html=True)


def chips(items: list[tuple[str, str]]) -> None:
    """items = [(text, kind), ...] kind in {ok,warn,info,''}"""
    html = '<div class="kady-chips">'
    for text, kind in items:
        cls = f"kady-chip chip-{kind}" if kind else "kady-chip"
        html += f'<span class="{cls}">{text}</span>'
    html += "</div>"
    st.markdown(html, unsafe_allow_html=True)


def sidebar_section(label: str) -> None:
    st.sidebar.markdown(f'<div class="kady-side">{label}</div>', unsafe_allow_html=True)


def sidebar_credits(team_label: str, github_url: str | None = None) -> None:
    """Small credits card at the bottom of the sidebar: team label + optional
    GitHub link. GitHub link renders only if `github_url` is non-empty, so a
    blank/unset value never shows a fake or dead link."""
    t = current_theme()
    link_html = (f'<a href="{github_url}" target="_blank" style="color:{TEAL};'
                'font-weight:600;text-decoration:none;font-size:12.5px;">'
                '\U0001f517 View source on GitHub</a>') if github_url else ""
    st.sidebar.markdown(f"""
    <div style="background:{t['surface']};border:1px solid {t['border']};border-radius:10px;
                padding:12px 14px;margin-top:14px;font-size:13px;">
      <div style="font-weight:700;font-size:11px;color:{t['muted']};text-transform:uppercase;
                  letter-spacing:.4px;margin-bottom:6px;">{team_label}</div>
      {link_html}
    </div>""", unsafe_allow_html=True)


def druggability_badge(n_drugs: int, n_pdb: int) -> str:
    if n_drugs > 0:
        return (f'<span class="badge badge-drug">💊 {n_drugs} drug mechanism(s)'
                f'{" · " + str(n_pdb) + " PDB" if n_pdb else ""}</span>')
    if n_pdb > 0:
        return f'<span class="badge badge-none">🧬 {n_pdb} PDB · no approved drug</span>'
    return '<span class="badge badge-none">no drug / structure</span>'

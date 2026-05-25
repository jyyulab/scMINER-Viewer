"""Shared page chrome (header + footer) for the single-study viewer and
multi-study browser.

The header carries the scMINER logo plus the "scMINER Viewer" wordmark;
the footer mirrors the scMINER portal's footer (St. Jude logo +
copyright + Yu Lab link + GitHub icon). Both render as plain HTML so
they work uniformly across the bare `run_app` page and the
`run_browser` shell. Logos live next to this module under `assets/`
and are served via Shiny's `static_assets` route mounted at
`/scviewer-assets` (see `static_assets_mount()`).
"""

from __future__ import annotations

import datetime
from pathlib import Path

from shiny import ui


# Mount point for the package's static assets. Picked to avoid clashing
# with any user-provided static dir or shiny's own asset routes.
_STATIC_PREFIX = "/scviewer-assets"

_ASSETS_DIR = Path(__file__).parent / "assets"


def static_assets_mount() -> dict[str, Path]:
    """Return the {prefix: dir} mapping to pass to `App(static_assets=...)`.

    Both the single-study viewer (`build_app`) and the multi-study
    browser (`build_browser`) need this so the header/footer logos
    resolve over HTTP.
    """
    return {_STATIC_PREFIX: _ASSETS_DIR}


_CHROME_CSS = f"""
  :root {{
    --scv-color-header-bg: #e6e9ed;
    --scv-color-footer-bg: #f2f5f9;
    --scv-color-text-primary: #2c3e50;
    --scv-color-brand-primary: #0f75bc;
    --scv-color-border: #d0d7de;
  }}

  /* Header */
  .scv-header {{
      background-color: var(--scv-color-header-bg);
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.08);
      padding: 10px 24px;
      display: flex;
      align-items: center;
      gap: 14px;
      border-bottom: 1px solid var(--scv-color-border);
  }}
  .scv-header-logo {{ height: 44px; }}
  .scv-header-title {{
      font-size: 20px;
      font-weight: 600;
      color: var(--scv-color-text-primary);
      letter-spacing: 0.3px;
  }}

  /* Body wrapper -- keep footer at viewport bottom on short pages. */
  .scv-page {{
      display: flex;
      flex-direction: column;
      min-height: 100vh;
  }}
  .scv-content {{ flex: 1 0 auto; padding: 12px 16px; }}

  /* Footer */
  .scv-footer {{
      background-color: var(--scv-color-footer-bg);
      border-top: 1px solid var(--scv-color-border);
      padding: 12px 20px;
      flex-shrink: 0;
  }}
  .scv-footer-inner {{
      max-width: 1400px;
      margin: 0 auto;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      min-height: 50px;
      flex-wrap: wrap;
  }}
  .scv-footer-logo {{ height: 38px; flex-shrink: 0; }}
  .scv-footer-github {{ height: 24px; flex-shrink: 0; }}
  .scv-footer-text {{
      flex: 1 1 auto;
      color: var(--scv-color-text-primary);
      font-size: clamp(11px, 1vw, 14px);
      text-align: center;
  }}
  .scv-footer-text a {{
      color: var(--scv-color-brand-primary);
      text-decoration: none;
  }}
  .scv-footer-text a:hover {{ text-decoration: underline; }}

  @media (max-width: 576px) {{
      .scv-header {{ padding: 8px 12px; }}
      .scv-header-title {{ font-size: 16px; }}
      .scv-header-logo {{ height: 32px; }}
  }}
"""


def page_chrome_css() -> ui.Tag:
    """Return a `<style>` tag with the header/footer CSS.

    Embed once per page (the `:root` vars are idempotent if duplicated,
    but the rules themselves are wasted bytes when repeated).
    """
    return ui.tags.style(_CHROME_CSS)


def page_header() -> ui.Tag:
    """scMINER logo + 'scMINER Viewer' wordmark, sitting at the page top."""
    return ui.tags.header(
        {"class": "scv-header"},
        ui.tags.img(
            {"class": "scv-header-logo",
             "src": f"{_STATIC_PREFIX}/scMINER_logo.png",
             "alt": "scMINER"},
        ),
        ui.tags.span({"class": "scv-header-title"}, "scMINER Viewer"),
    )


def page_footer() -> ui.Tag:
    """St. Jude logo + copyright + Yu Lab link + GitHub icon.

    Mirrors `scminer-portal-frontend/src/components/pagefooter.vue`.
    """
    year = datetime.date.today().year
    return ui.tags.footer(
        {"class": "scv-footer"},
        ui.div(
            {"class": "scv-footer-inner"},
            ui.tags.img(
                {"class": "scv-footer-logo",
                 "src": f"{_STATIC_PREFIX}/StJude.png",
                 "alt": "St. Jude"},
            ),
            ui.div(
                {"class": "scv-footer-text"},
                f"Copyright © {year} | Powered by ",
                ui.tags.a(
                    {"href": "https://www.stjude.org/research/labs/yu-lab.html",
                     "target": "_blank", "rel": "noopener"},
                    "the Yu Lab",
                ),
                " at St. Jude Children's Research Hospital",
            ),
            ui.tags.a(
                {"href": "https://github.com/jyyulab",
                 "target": "_blank", "rel": "noopener"},
                ui.tags.img(
                    {"class": "scv-footer-github",
                     "src": f"{_STATIC_PREFIX}/github-mark.png",
                     "alt": "GitHub"},
                ),
            ),
        ),
    )

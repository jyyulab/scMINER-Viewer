#' Shared page chrome (header + footer) for the single-study viewer and
#' multi-study browser.
#'
#' The header carries the scMINER logo plus the "scMINER Viewer"
#' wordmark; the footer mirrors the scMINER portal's footer (St. Jude
#' logo + copyright + Yu Lab link + GitHub icon). Logos ship under
#' `inst/assets/` and are served via `shiny::addResourcePath("scviewer-
#' assets", ...)` from `.register_assets()`, which is called from
#' `.onLoad` (zzz.R) so the route is mounted as soon as the package is
#' loaded.

#' Mount point for the package's static assets. Matches the Python
#' package's `_STATIC_PREFIX` so `<img src="scviewer-assets/...">` works
#' identically on both sides.
#' @noRd
.assets_prefix <- "scviewer-assets"

#' Register the package's `inst/assets/` directory with Shiny so the
#' header/footer logos resolve over HTTP. Idempotent.
#' @noRd
.register_assets <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) return(invisible(NULL))
  path <- system.file("assets", package = "scminerViewer")
  if (!nzchar(path) || !dir.exists(path)) return(invisible(NULL))
  shiny::addResourcePath(.assets_prefix, path)
  invisible(NULL)
}

#' Inline CSS for the header/footer chrome. Returns a `shiny::tags$style`
#' tag — embed once per page (the `:root` vars are idempotent if
#' duplicated, but the rules are wasted bytes when repeated).
#' @noRd
.page_chrome_css <- function() {
  shiny::tags$style(shiny::HTML("
    :root {
      --scv-color-header-bg: #e6e9ed;
      --scv-color-footer-bg: #f2f5f9;
      --scv-color-text-primary: #2c3e50;
      --scv-color-brand-primary: #0f75bc;
      --scv-color-border: #d0d7de;
    }

    /* Header */
    .scv-header {
        background-color: var(--scv-color-header-bg);
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.08);
        padding: 10px 24px;
        display: flex;
        align-items: center;
        gap: 14px;
        border-bottom: 1px solid var(--scv-color-border);
    }
    .scv-header-logo { height: 44px; }
    .scv-header-title {
        font-size: 20px;
        font-weight: 600;
        color: var(--scv-color-text-primary);
        letter-spacing: 0.3px;
    }

    .scv-content { padding: 12px 16px; }

    /* Footer */
    .scv-footer {
        background-color: var(--scv-color-footer-bg);
        border-top: 1px solid var(--scv-color-border);
        padding: 12px 20px;
    }
    .scv-footer-inner {
        max-width: 1400px;
        margin: 0 auto;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 16px;
        min-height: 50px;
        flex-wrap: wrap;
    }
    .scv-footer-logo { height: 38px; flex-shrink: 0; }
    .scv-footer-github { height: 24px; flex-shrink: 0; }
    .scv-footer-text {
        flex: 1 1 auto;
        color: var(--scv-color-text-primary);
        font-size: clamp(11px, 1vw, 14px);
        text-align: center;
    }
    .scv-footer-text a {
        color: var(--scv-color-brand-primary);
        text-decoration: none;
    }
    .scv-footer-text a:hover { text-decoration: underline; }

    @media (max-width: 576px) {
        .scv-header { padding: 8px 12px; }
        .scv-header-title { font-size: 16px; }
        .scv-header-logo { height: 32px; }
    }
  "))
}

#' scMINER logo + "scMINER Viewer" wordmark, sitting at the page top.
#' @noRd
.page_header <- function() {
  shiny::tags$header(
    class = "scv-header",
    shiny::tags$img(
      class = "scv-header-logo",
      src   = paste0(.assets_prefix, "/scMINER_logo.png"),
      alt   = "scMINER"
    ),
    shiny::tags$span(class = "scv-header-title", "scMINER Viewer")
  )
}

#' St. Jude logo + copyright + Yu Lab link + GitHub icon. Mirrors
#' scminer-portal-frontend/src/components/pagefooter.vue.
#' @noRd
.page_footer <- function() {
  year <- format(Sys.Date(), "%Y")
  shiny::tags$footer(
    class = "scv-footer",
    shiny::div(
      class = "scv-footer-inner",
      shiny::tags$img(
        class = "scv-footer-logo",
        src   = paste0(.assets_prefix, "/StJude.png"),
        alt   = "St. Jude"
      ),
      shiny::div(
        class = "scv-footer-text",
        sprintf("Copyright © %s | Powered by ", year),
        shiny::tags$a(
          href   = "https://www.stjude.org/research/labs/yu-lab.html",
          target = "_blank", rel = "noopener",
          "the Yu Lab"
        ),
        " at St. Jude Children's Research Hospital"
      ),
      shiny::tags$a(
        href   = "https://github.com/jyyulab",
        target = "_blank", rel = "noopener",
        shiny::tags$img(
          class = "scv-footer-github",
          src   = paste0(.assets_prefix, "/github-mark.png"),
          alt   = "GitHub"
        )
      )
    )
  )
}

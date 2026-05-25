.onLoad <- function(libname, pkgname) {
  # Mount the package's header/footer logo assets at /scviewer-assets/
  # so the viewer + browser HTML can resolve `<img src="scviewer-assets/
  # scMINER_logo.png">` etc. Safe to call before any Shiny session
  # starts (shiny::addResourcePath is process-global).
  .register_assets()
}

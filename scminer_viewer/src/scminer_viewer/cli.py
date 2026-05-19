"""CLI entry point — `scminer-viewer {run, browse, info, list}`."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scminer-viewer",
        description="Serve scMINER study bundles in the browser.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # ---- run: single study ------------------------------------------------
    run = sub.add_parser("run",
                          help="Launch the single-study Shiny webui.")
    run.add_argument("bundle_path",
                     help="Path to a .scminer.h5 bundle.")
    run.add_argument("--shard-dir", default=None,
                     help="Directory containing the shard tree if it is "
                          "not co-located with the bundle. Defaults to "
                          "the bundle's parent directory.")
    run.add_argument("--host", default="127.0.0.1")
    run.add_argument("--port", type=int, default=8000)
    run.add_argument("--no-browser", action="store_true",
                     help="Don't open a browser window.")
    run.add_argument("--allow-remote", action="store_true",
                     help="Shortcut for --host 0.0.0.0 --no-browser. "
                          "Use on HPC / remote nodes when SSH-tunneling "
                          "isn't an option. No authentication is added; "
                          "see README for the recommended SSH-tunnel "
                          "pattern instead.")

    # ---- browse: multi-study card-grid index ------------------------------
    browse = sub.add_parser("browse",
                             help="Launch the multi-study browser.")
    browse.add_argument("root_dir",
                        help="Directory containing <studyID>/<studyID>"
                             ".scminer.h5 subfolders.")
    browse.add_argument("--shard-dir", default=None,
                        help="Override the per-study shard root for every "
                             "bundle under root_dir.")
    browse.add_argument("--host", default="127.0.0.1")
    browse.add_argument("--port", type=int, default=8000)
    browse.add_argument("--no-browser", action="store_true",
                        help="Don't open a browser window.")
    browse.add_argument("--allow-remote", action="store_true",
                        help="Shortcut for --host 0.0.0.0 --no-browser. "
                             "Same caveats as in `run`.")

    # ---- info: print one-line summary -------------------------------------
    info = sub.add_parser("info", help="Print summary of a bundle.")
    info.add_argument("bundle_path",
                      help="Path to a .scminer.h5 bundle.")
    info.add_argument("--shard-dir", default=None)

    # ---- list: enumerate studies under a root -----------------------------
    lst = sub.add_parser("list",
                          help="List every study bundle under a root dir.")
    lst.add_argument("root_dir",
                     help="Directory containing <studyID>/<studyID>"
                          ".scminer.h5 subfolders.")

    # ---- prepare: build a study from a YAML config ------------------------
    prep = sub.add_parser(
        "prepare",
        help="Run prepare_study() from a YAML config. Writes graph layout "
             "+ .scminer.h5 bundle under config.output/<studyID>/.",
    )
    prep.add_argument("config_path",
                      help="Path to a YAML config (see "
                           "scminer_viewer.prepare.load_study_config).")
    prep.add_argument("--emit", default="graph,bundle",
                      help="Comma-separated subset of {graph,bundle}. "
                           "Default: graph,bundle.")
    prep.add_argument("--quiet", action="store_true",
                      help="Suppress per-shard progress messages.")

    return p


def _apply_allow_remote(args: argparse.Namespace) -> None:
    """Resolve `--allow-remote` into the underlying --host / --no-browser
    settings and print a security warning to stderr. Honors any explicit
    --host the user passed alongside the shortcut."""
    if not getattr(args, "allow_remote", False):
        return
    if args.host == "127.0.0.1":
        args.host = "0.0.0.0"
    args.no_browser = True
    print(
        f"[scminer-viewer] --allow-remote: binding on {args.host}:"
        f"{args.port}. The viewer has no built-in authentication; "
        f"prefer an SSH tunnel on shared / public-routable networks "
        f"(see README \"Exposing to a remote host\").",
        file=sys.stderr,
    )


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    _apply_allow_remote(args)

    if args.cmd == "info":
        from .data import load_study
        study = load_study(args.bundle_path, shard_dir=args.shard_dir)
        print(study)
        return 0

    if args.cmd == "list":
        from .browser import discover_studies
        df = discover_studies(args.root_dir)
        if df.empty:
            print(f"No studies found at: {args.root_dir}")
            return 1
        # Compact table
        cols = ["studyID", "studyAbbr", "shortTitle", "species",
                "n_cells", "n_genes", "n_clusters", "bundle_path"]
        print(df[cols].to_string(index=False))
        return 0

    if args.cmd == "run":
        path = Path(args.bundle_path)
        if not path.exists():
            print(f"Bundle not found: {path}", file=sys.stderr)
            return 1
        from .app import run_app
        run_app(
            bundle_path=str(path),
            shard_dir=args.shard_dir,
            host=args.host,
            port=args.port,
            launch_browser=not args.no_browser,
        )
        return 0

    if args.cmd == "browse":
        root = Path(args.root_dir)
        if not root.exists():
            print(f"Root dir not found: {root}", file=sys.stderr)
            return 1
        from .browser import run_browser
        run_browser(
            root_dir=str(root),
            shard_dir=args.shard_dir,
            host=args.host,
            port=args.port,
            launch_browser=not args.no_browser,
        )
        return 0

    if args.cmd == "prepare":
        cfg_path = Path(args.config_path)
        if not cfg_path.exists():
            print(f"Config not found: {cfg_path}", file=sys.stderr)
            return 1
        emit_tokens = [t.strip() for t in args.emit.split(",") if t.strip()]
        from .prepare import prepare_study
        result = prepare_study(
            str(cfg_path),
            emit=tuple(emit_tokens),
            verbose=not args.quiet,
        )
        print(f"Output dir:  {result['out_dir']}")
        if result.get("bundle_path"):
            print(f"Bundle path: {result['bundle_path']}")
        return 0

    return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

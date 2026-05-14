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

    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

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

    return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

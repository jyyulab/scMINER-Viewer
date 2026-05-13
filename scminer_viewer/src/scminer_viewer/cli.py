"""CLI entry point — `scminer-viewer run path/to/bundle.scminer.h5`."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scminer-viewer",
        description="Serve an scMINER study bundle in a browser.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Launch the Shiny webui.")
    run.add_argument("bundle_path",
                     help="Path to a .scminer.h5 bundle.")
    run.add_argument("--host", default="127.0.0.1")
    run.add_argument("--port", type=int, default=8000)
    run.add_argument("--no-browser", action="store_true",
                     help="Don't open a browser window.")

    info = sub.add_parser("info", help="Print summary of a bundle.")
    info.add_argument("bundle_path",
                      help="Path to a .scminer.h5 bundle.")

    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.cmd == "info":
        from .data import load_study

        study = load_study(args.bundle_path)
        print(study)
        return 0

    if args.cmd == "run":
        path = Path(args.bundle_path)
        if not path.exists():
            print(f"Bundle not found: {path}", file=sys.stderr)
            return 1
        from .app import run_app

        run_app(
            bundle_path=str(path),
            host=args.host,
            port=args.port,
            launch_browser=not args.no_browser,
        )
        return 0

    return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

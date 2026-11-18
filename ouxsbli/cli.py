"""
Command line front end of :class:`ouxsbli.Case`, usable from any directory::

    ouxsbli root                                  # which source tree is used
    ouxsbli list                                  # template cases of that tree
    ouxsbli init 2D_solver/EVC ./evc nx=66 ny=66  # copy a template and patch it
    ouxsbli build ./evc                           # cmake + make in ./evc/build
    ouxsbli run ./evc -n 2                        # build if needed, then mpirun

``key=value`` arguments are the keyword arguments of ``Case`` (aliases such as
``scheme``, ``accuracy``, ``prec`` go to ``config.fypp``; other names are
``mod_globals.f90`` parameters). ``build`` and ``run`` work on any case
directory, including one of your own outside the source tree.
"""
import argparse
import sys

from .case import Case
from .paths import repo_root, root_origin, template_cases


def _value(text: str):
    """'65' -> 65, '1e-3' -> 0.001, 'true' -> True, anything else stays a string."""
    low = text.lower()
    if low in ("true", ".true."):
        return True
    if low in ("false", ".false."):
        return False
    for conv in (int, float):
        try:
            return conv(low.replace("d", "e") if conv is float else text)
        except ValueError:
            pass
    return text


def _params(items):
    out = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"expected key=value, got {item!r}")
        key, val = item.split("=", 1)
        out[key.strip()] = _value(val.strip())
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="ouxsbli", description="Set up, build and run OUxSBLI cases.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("root", help="print the source tree in use")
    sub.add_parser("list", help="list the template cases of the source tree")
    p = sub.add_parser("init", help="copy a template case to a new directory and patch it")
    p.add_argument("template", help="e.g. 2D_solver/EVC, or a path to a case directory")
    p.add_argument("dir", help="new case directory (replaced if it exists)")
    p.add_argument("params", nargs="*", metavar="key=value")
    p = sub.add_parser("build", help="cmake + make in <dir>/build")
    p.add_argument("dir", nargs="?", default=".")
    p = sub.add_parser("run", help="build if needed, then mpirun in <dir>")
    p.add_argument("dir", nargs="?", default=".")
    p.add_argument("-n", "--nranks", type=int, default=2)
    p.add_argument("--rebuild", action="store_true", help="run cmake + make even if build/a.out exists")
    args = ap.parse_args(argv)

    try:
        if args.cmd == "root":
            print(f"{repo_root()}  (from {root_origin()})")
        elif args.cmd == "list":
            print("\n".join(template_cases()))
        elif args.cmd == "init":
            case = Case(args.template, args.dir, **_params(args.params))
            case._setup()
            print(f"{case.workdir}  (from {case.source})")
        elif args.cmd == "build":
            Case.from_dir(args.dir, verbose=True).build()
        elif args.cmd == "run":
            case = Case.from_dir(args.dir, verbose=True)
            if args.rebuild or not case._built:
                case.build()
            case.run(nranks=args.nranks)
    except (RuntimeError, FileNotFoundError, ValueError) as err:
        print(f"ouxsbli: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

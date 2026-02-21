#!/usr/bin/env python3
"""
Root a phylogenetic tree using midpoint or outgroup rooting.

Usage:
  root_tree.py --input tree.nwk --output tree.rooted.nwk --method midpoint
  root_tree.py --input tree.nwk --output tree.rooted.nwk --method outgroup --outgroup "taxon_name"
"""

import argparse
import sys

from ete3 import Tree


def parse_args():
    parser = argparse.ArgumentParser(description="Root a phylogenetic tree")
    parser.add_argument("--input", required=True, help="Input Newick tree")
    parser.add_argument("--output", required=True, help="Output rooted Newick tree")
    parser.add_argument("--method", choices=["midpoint", "outgroup"], default="midpoint")
    parser.add_argument("--outgroup", default=None, help="Outgroup taxon name (if method=outgroup)")
    return parser.parse_args()


def main():
    args = parse_args()

    # Read tree
    t = Tree(args.input)

    if args.method == "midpoint":
        # Midpoint rooting: root at the midpoint of the longest path
        midpoint = t.get_midpoint_outgroup()
        if midpoint:
            t.set_outgroup(midpoint)
            print(f"[root] Midpoint rooting applied")
        else:
            print(f"[root] Could not determine midpoint, tree unchanged", file=sys.stderr)

    elif args.method == "outgroup":
        if not args.outgroup:
            print("[error] --outgroup required when method=outgroup", file=sys.stderr)
            sys.exit(1)
        try:
            outgroup_node = t & args.outgroup
            t.set_outgroup(outgroup_node)
            print(f"[root] Outgroup rooting applied: {args.outgroup}")
        except Exception as e:
            print(f"[error] Outgroup '{args.outgroup}' not found: {e}", file=sys.stderr)
            # Fallback to midpoint
            midpoint = t.get_midpoint_outgroup()
            if midpoint:
                t.set_outgroup(midpoint)
                print(f"[root] Fallback to midpoint rooting")

    # Write rooted tree
    t.write(outfile=args.output, format=1)
    print(f"[root] Saved: {args.output}")


if __name__ == "__main__":
    main()

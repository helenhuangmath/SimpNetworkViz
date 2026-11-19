"""The `snv` argument parser, with help text made safe for reStructuredText.

sphinx-argparse reads each option's help string as reStructuredText, where a pair of
vertical bars (as in "|log2FC|") means a substitution. This wrapper escapes them so the
package source does not have to change.
"""
import argparse

from simpnetworkviz.cli import build_parser as _build_parser


def _escape(parser):
    for action in parser._actions:
        if isinstance(action.help, str):
            action.help = action.help.replace("|", "\\|")
        if isinstance(action, argparse._SubParsersAction):
            for sub in action.choices.values():
                _escape(sub)


def build_parser():
    parser = _build_parser()
    _escape(parser)
    return parser

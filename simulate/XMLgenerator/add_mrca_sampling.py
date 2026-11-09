#!/usr/bin/env python3
"""Add a reflected MRCA-height proposal to a generated fixed-LCA PHYFUM XML.

Python 3.6.8 or newer; standard library only. Analysis filenames, including
MLE reader references, receive --suffix before their final extension.
The optional second positional argument is the exact output XML path.
"""

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET


class CommentPreservingTreeBuilder(ET.TreeBuilder):
    """Preserve in-document comments without Python 3.8's insert_comments option."""

    def __init__(self):
        super().__init__()
        self.depth = 0

    def start(self, tag, attrs):
        self.depth += 1
        return super().start(tag, attrs)

    def end(self, tag):
        element = super().end(tag)
        self.depth -= 1
        return element

    def comment(self, text):
        if self.depth:
            super().start(ET.Comment, {})
            super().data(text)
            super().end(ET.Comment)


def suffixed_name(name, suffix):
    """Keep directories and insert the suffix before the final extension."""
    path = Path(name)
    return str(path.with_name(path.stem + suffix + path.suffix))


def add_sampling(tree, suffix="_MRCASampling"):
    """Validate the generated fixed-LCA configuration, then modify it in place."""
    if not suffix or "/" in suffix or "\\" in suffix:
        raise ValueError("--suffix must be nonempty and cannot contain path separators")
    root = tree.getroot()
    ids = {}
    for element in root.iter():
        identifier = element.get("id")
        if identifier:
            if identifier in ids:
                raise ValueError("Duplicate XML ID: " + identifier)
            ids[identifier] = element

    def resolve(element):
        if element is not None and element.get("idref"):
            return ids.get(element.get("idref"))
        return element

    likelihoods = [e for e in root.iter("cenancestorTreeLikelihood") if not e.get("idref")]
    schedules = [e for e in root.iter("operators") if not e.get("idref")]
    if len(likelihoods) != 1 or len(schedules) != 1:
        raise ValueError("Expected one defined cenancestorTreeLikelihood and one operator schedule")
    likelihood = likelihoods[0]
    if likelihood.get("heightRules", "false").lower() != "true":
        raise ValueError("The LCA must be fixed using heightRules='true'")
    height = resolve(likelihood.find("cenancestorHeight/parameter"))
    if height is None or height.tag != "parameter" or height.get("dimension", "1") != "1":
        raise ValueError("Expected a scalar cenancestorHeight parameter")
    try:
        upper = float(height.get("value", ""))
    except ValueError:
        raise ValueError("cenancestorHeight needs an explicit scalar value") from None
    if not math.isfinite(upper) or upper <= 0:
        raise ValueError("cenancestorHeight must be positive and finite")
    tree_model = resolve(likelihood.find("treeModel"))
    parameter = tree_model.find("rootHeight/parameter") if tree_model is not None else None
    if parameter is None or not parameter.get("id"):
        raise ValueError("Expected an identified treeModel rootHeight parameter")
    root_id = parameter.get("id")
    schedule = schedules[0]
    for entry in schedule:
        operator = resolve(entry)
        if operator is None:
            raise ValueError("Unresolved operator reference")
        references = [e.get("idref") for e in operator.iter("parameter")]
        if height.get("id") and height.get("id") in references:
            raise ValueError("cenancestorHeight is sampled; this script requires a fixed LCA")
        if operator.tag == "randomWalkOperator" and root_id in references:
            raise ValueError("A randomWalkOperator already targets the root height")

    # Omit 'lower': TreeModel supplies the live oldest-child height bound.
    operator = ET.Element("randomWalkOperator", {
        "windowSize": "2.5", "weight": "3.0", "boundaryCondition": "reflecting",
        "upper": height.get("value"), "autoOptimize": "true"})
    ET.SubElement(operator, "parameter", {"idref": root_id})
    operator.text = "\n            "
    operator[0].tail = "\n        "
    operator.tail = "\n    "
    schedule.append(operator)

    # Keep MLE/other analysis readers aligned with the renamed log writers.
    for element in root.iter():
        for attribute in ("fileName", "operatorAnalysis"):
            name = element.get(attribute)
            if name:
                element.set(attribute, suffixed_name(name, suffix))
    return tree


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_xml", type=Path, help="Generated PHYFUM XML with a fixed LCA")
    parser.add_argument("output_xml", nargs="?", type=Path,
                        help="Exact output XML path (default: input path with suffix)")
    parser.add_argument("--suffix", default="_MRCASampling",
                        help="Suffix before filename extensions (default: _MRCASampling)")
    args = parser.parse_args(argv)
    try:
        output = args.output_xml or Path(suffixed_name(str(args.input_xml), args.suffix))
        if output.resolve() == args.input_xml.resolve():
            raise ValueError("Output XML must differ from the input XML")
        tree = ET.parse(str(args.input_xml), parser=ET.XMLParser(
            target=CommentPreservingTreeBuilder()))
        add_sampling(tree, args.suffix)
        # Exclusive creation protects both existing results and the original XML.
        with output.open("xb") as stream:
            tree.write(stream, encoding="utf-8", xml_declaration=True)
        print(output)
    except (OSError, ValueError, ET.ParseError) as error:
        parser.exit(2, "error: %s\n" % error)


if __name__ == "__main__":
    main()

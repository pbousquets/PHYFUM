#!/usr/bin/env python3
"""Convert a generated PHYFUM two-rate XML to four independent rates.

Uses only Python's standard library. Nu inherits mu's initial value, bounds,
proposals and prior; zeta inherits gamma's. Lambda is unchanged. Output paths
(including MLE logs and analyses) receive the same suffix as the XML filename.
"""

import argparse
from copy import deepcopy
import math
from pathlib import Path
import xml.etree.ElementTree as ET


class CommentPreservingTreeBuilder(ET.TreeBuilder):
    """Keep comments on Python versions predating TreeBuilder(insert_comments=...)."""

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
        # ElementTree serializes a single root, so retain comments inside it only.
        if self.depth:
            super().start(ET.Comment, {})
            super().data(text)
            super().end(ET.Comment)


def suffixed_name(name, suffix):
    """Insert the suffix before the final extension, preserving the directory."""
    path = Path(name)
    return str(path.with_name(path.stem + suffix + path.suffix))


def insert_after(parent, original, addition):
    parent.insert(list(parent).index(original) + 1, addition)


def convert_tree(tree, suffix):
    """Convert a parsed tree in place; reject unsupported or incomplete inputs."""
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
    models = [model for model in root.iter("flipflopModel") if model.get("idref") is None]
    if len(models) != 1:
        raise ValueError("Expected exactly one defined flipflopModel")
    model = models[0]
    if model.find("nu") is not None or model.find("zeta") is not None:
        raise ValueError("The model already defines nu or zeta; expected a two-rate XML")

    def clone(element, old_id, new_id, rate):
        """Give copied operators/priors their own IDs and repair internal references."""
        addition = deepcopy(element)
        replacements = {old_id: new_id}
        for child in addition.iter():
            identifier = child.get("id")
            if identifier:
                candidate = identifier + "." + rate
                if candidate in ids:
                    raise ValueError("New XML ID already exists: " + candidate)
                replacements[identifier] = candidate
                ids[candidate] = child
                child.set("id", candidate)
        for child in addition.iter():
            reference = child.get("idref")
            if reference in replacements:
                child.set("idref", replacements[reference])
        return addition

    for old_rate, new_rate in (("mu", "nu"), ("gamma", "zeta")):
        wrapper = model.find(old_rate)
        parameter = wrapper.find("parameter") if wrapper is not None else None
        if parameter is not None and parameter.get("idref"):
            parameter = ids.get(parameter.get("idref"))
        if parameter is None or parameter.tag != "parameter" or not parameter.get("id"):
            raise ValueError("Expected an identified parameter for " + old_rate)
        old_id = parameter.get("id")
        try:
            value = float(parameter.get("value", ""))
        except ValueError:
            raise ValueError(old_rate + " needs an explicit scalar initial value") from None
        if not math.isfinite(value) or value <= 0:
            raise ValueError(old_rate + " needs a positive finite initial value for scale proposals")
        if list(parameter) or parameter.get("dimension", "1") != "1":
            raise ValueError(old_rate + " must be a scalar parameter")
        new_id = (old_id[:-len(old_rate)] + new_rate
                  if old_id.endswith(old_rate) else old_id + "." + new_rate)
        if new_id in ids:
            raise ValueError("New XML ID already exists: " + new_id)
        new_wrapper = ET.Element(new_rate)
        new_wrapper.text = wrapper.text
        new_wrapper.tail = wrapper.tail
        new_parameter = deepcopy(parameter)
        new_parameter.set("id", new_id)
        new_parameter.attrib.pop("idref", None)
        new_wrapper.append(new_parameter)
        insert_after(model, wrapper, new_wrapper)
        ids[new_id] = new_parameter

        # Copy actual configured settings, including any customized prior or tuning.
        parents = {child: parent for parent in root.iter() for child in parent}
        operators = []
        priors = []
        for element in list(root.iter()):
            if element.find("parameter[@idref='%s']" % old_id) is None:
                continue
            if element.tag == "scaleOperator":
                if len(element) != 1:
                    raise ValueError("Expected a single-parameter scaleOperator for " + old_id)
                operators.append(element)
            elif isinstance(element.tag, str) and element.tag.endswith("Prior"):
                if len(element.findall("parameter")) != 1:
                    raise ValueError("Expected a single-parameter prior for " + old_id)
                priors.append(element)
        if not operators or not priors:
            raise ValueError("Missing scaleOperator or prior for " + old_id)
        for element in operators + priors:
            insert_after(parents[element], element, clone(element, old_id, new_id, new_rate))

        # Extend joint rate/time proposals, and log the new independently sampled rates.
        for operator in root.iter("upDownOperator"):
            for group in operator:
                for reference in list(group.findall("parameter")):
                    if reference.get("idref") == old_id:
                        insert_after(group, reference, clone(reference, old_id, new_id, new_rate))
        for log in root.iter("log"):
            for entry in list(log):
                if (entry.tag == "parameter" and entry.get("idref") == old_id
                        or entry.tag == "column" and entry.find(
                            "parameter[@idref='%s']" % old_id) is not None):
                    addition = clone(entry, old_id, new_id, new_rate)
                    if addition.tag == "column" and addition.get("label"):
                        addition.set("label", new_rate)
                    insert_after(log, entry, addition)

    # fileName also occurs in MLE analysis readers: update those references together
    # with their log writers, so all stages continue to use the same suffixed files.
    for element in root.iter():
        for attribute in ("fileName", "operatorAnalysis"):
            name = element.get(attribute)
            if name:
                element.set(attribute, suffixed_name(name, suffix))
    return tree


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_xml", type=Path, help="Existing two-rate PHYFUM XML")
    parser.add_argument("output_xml", nargs="?", type=Path,
                        help="Exact output XML path (default: input with suffix before extension)")
    parser.add_argument("--suffix", default="_4param",
                        help="Suffix for XML and analysis output filenames (default: _4param)")
    args = parser.parse_args(argv)
    try:
        output = args.output_xml or Path(suffixed_name(str(args.input_xml), args.suffix))
        if output.resolve() == args.input_xml.resolve():
            raise ValueError("Output XML must differ from the input XML")
        # Preserve comments, including commented alternative model configurations.
        tree = ET.parse(args.input_xml, parser=ET.XMLParser(
            target=CommentPreservingTreeBuilder()))
        convert_tree(tree, args.suffix)
        # Exclusive creation protects existing analyses and earlier conversion results.
        with output.open("xb") as stream:
            tree.write(stream, encoding="utf-8", xml_declaration=True)
        print(output)
    except (OSError, ValueError, ET.ParseError) as error:
        parser.exit(2, "error: %s\n" % error)


if __name__ == "__main__":
    main()

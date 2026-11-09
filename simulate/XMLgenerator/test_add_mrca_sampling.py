"""Dependency-free tests; run with unittest discover -p test_add_mrca_sampling.py."""

from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
import xml.etree.ElementTree as ET

from add_mrca_sampling import add_sampling, CommentPreservingTreeBuilder, suffixed_name


SCRIPT = Path(__file__).with_name("add_mrca_sampling.py")
EXAMPLES = Path(__file__).resolve().parents[2] / "examples/release/flipflop"


class AddMRCASamplingTest(unittest.TestCase):
    def test_examples_preserve_model_and_prior(self):
        for name in ("sim1Absolute.xml", "sim1AbsoluteFourRates.xml"):
            tree = ET.parse(str(EXAMPLES / name))
            root = tree.getroot()
            model = ET.tostring(root.find("flipflopModel"))
            prior = ET.tostring(root.find("mcmc/posterior/prior"))
            original_operators = list(root.find("operators"))
            filenames = [(e, attr, e.get(attr)) for e in root.iter()
                         for attr in ("fileName", "operatorAnalysis") if e.get(attr)]
            add_sampling(tree)
            self.assertEqual(model, ET.tostring(root.find("flipflopModel")))
            self.assertEqual(prior, ET.tostring(root.find("mcmc/posterior/prior")))
            self.assertEqual(list(root.find("operators"))[:-1], original_operators)
            move = root.find("operators")[-1]
            self.assertEqual(move.tag, "randomWalkOperator")
            self.assertEqual(move.get("upper"), "45")
            self.assertEqual(move.get("boundaryCondition"), "reflecting")
            self.assertIsNone(move.get("lower"))
            self.assertEqual(move.find("parameter").get("idref"), "treeModel.rootHeight")
            for element, attribute, name in filenames:
                self.assertEqual(element.get(attribute), suffixed_name(name, "_MRCASampling"))

    def test_references_custom_ids_and_analysis_readers(self):
        tree = ET.parse(str(EXAMPLES / "sim1Absolute.xml"))
        root = tree.getroot()
        likelihood = root.find("cenancestorTreeLikelihood")
        height = likelihood.find("cenancestorHeight/parameter")
        height.set("value", "67.5")
        likelihood.find("cenancestorHeight").remove(height)
        root.append(height)
        ET.SubElement(likelihood.find("cenancestorHeight"), "parameter", idref="luca_height")
        root.find("treeModel/rootHeight/parameter").set("id", "custom.root")
        for tag in ("log", "pathSamplingAnalysis", "steppingStoneSamplingAnalysis"):
            ET.SubElement(root, tag, fileName="results/run.MLE.log")
        add_sampling(tree, "_custom")
        move = root.find("operators")[-1]
        self.assertEqual(move.get("upper"), "67.5")
        self.assertEqual(move[0].get("idref"), "custom.root")
        self.assertEqual(root.find("pathSamplingAnalysis").get("fileName"),
                         "results/run.MLE_custom.log")
        self.assertEqual(root.find("steppingStoneSamplingAnalysis").get("fileName"),
                         "results/run.MLE_custom.log")

    def test_invalid_inputs(self):
        for case in ("duplicate", "branch_rules", "moving_lca", "nan", "negative", "vector",
                     "missing_root", "empty_suffix", "path_suffix"):
            tree = ET.parse(str(EXAMPLES / "sim1Absolute.xml"))
            root = tree.getroot()
            suffix = "_test"
            height = root.find("cenancestorTreeLikelihood/cenancestorHeight/parameter")
            if case == "duplicate":
                add_sampling(tree)
            elif case == "branch_rules":
                root.find("cenancestorTreeLikelihood").attrib.pop("heightRules")
            elif case == "moving_lca":
                op = ET.SubElement(root.find("operators"), "scaleOperator")
                ET.SubElement(op, "parameter", idref="luca_height")
            elif case in ("nan", "negative", "vector"):
                height.set("value", {"nan": "nan", "negative": "-1", "vector": "4 5"}[case])
            elif case == "missing_root":
                root.find("treeModel/rootHeight/parameter").attrib.pop("id")
            elif case == "empty_suffix":
                suffix = ""
            elif case == "path_suffix":
                suffix = "../out"
            with self.subTest(case=case):
                with self.assertRaises(ValueError):
                    add_sampling(tree, suffix)

    def test_cli_defaults_explicit_output_and_overwrite_protection(self):
        def run(*args):
            return subprocess.run([sys.executable, str(SCRIPT)] + list(args),
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  universal_newlines=True)
        with TemporaryDirectory() as directory:
            source = Path(directory) / "analysis.xml"
            source.write_bytes((EXAMPLES / "sim1Absolute.xml").read_bytes())
            before = source.read_bytes()
            result = run(str(source))
            self.assertEqual(result.returncode, 0, result.stderr)
            output = Path(directory) / "analysis_MRCASampling.xml"
            self.assertTrue(output.is_file())
            self.assertIn(b"<!--", output.read_bytes())
            saved = output.read_bytes()
            self.assertNotEqual(run(str(source)).returncode, 0)
            self.assertEqual(output.read_bytes(), saved)
            explicit = Path(directory) / "chosen.xml"
            result = run(str(source), str(explicit), "--suffix", "_trial")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(ET.parse(str(explicit)).find("mcmc").get("operatorAnalysis")
                            .endswith("_trial.ops"))
            self.assertNotEqual(run(str(source), str(source)).returncode, 0)
            self.assertEqual(source.read_bytes(), before)

    def test_comment_preservation_and_final_extension(self):
        tree = ET.fromstring("<beast><!--keep--><x/></beast>", parser=ET.XMLParser(
            target=CommentPreservingTreeBuilder()))
        self.assertIn(b"<!--keep-->", ET.tostring(tree))
        self.assertEqual(suffixed_name("dir/run.log.gz", "_test"), "dir/run.log_test.gz")
        self.assertEqual(suffixed_name("dir/noextension", "_test"), "dir/noextension_test")


if __name__ == "__main__":
    unittest.main()

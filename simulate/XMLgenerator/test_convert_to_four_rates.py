"""Standard-library tests: python -m unittest discover -s simulate/XMLgenerator -p test_convert_to_four_rates.py"""

from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
import xml.etree.ElementTree as ET

from convert_to_four_rates import CommentPreservingTreeBuilder, convert_tree, suffixed_name


SCRIPT = Path(__file__).with_name("convert_to_four_rates.py")
EXAMPLE = Path(__file__).resolve().parents[2] / "examples/release/flipflop/sim1Absolute.xml"


class ConvertFourRatesTest(unittest.TestCase):
    def test_comments_without_python38_treebuilder_option(self):
        xml = "<!--before--><beast>text<!--inside--><child/>tail<!--last--></beast><!--after-->"
        root = ET.fromstring(xml, parser=ET.XMLParser(target=CommentPreservingTreeBuilder()))
        self.assertEqual(ET.tostring(root, encoding="unicode"),
                         "<beast>text<!--inside--><child />tail<!--last--></beast>")

    def test_complete_conversion_and_custom_settings(self):
        tree = ET.parse(EXAMPLE)
        root = tree.getroot()
        original_alignment = ET.tostring(root.find("afalignment"))
        # The converter must inherit the user's actual values, priors and tuning.
        root.find("flipflopModel/mu/parameter").set("value", "0.037")
        root.find("flipflopModel/gamma/parameter").set("value", "0.081")
        root.find("operators/scaleOperator/parameter[@idref='flipflop.mu']/..").set("scaleFactor", "0.7")
        prior = root.find("mcmc/posterior/prior/halfNormalPrior/parameter[@idref='flipflop.gamma']/..")
        prior.set("stdev", "0.12")
        prior.set("id", "gammaPrior")
        convert_tree(tree, "_4param")
        self.assertEqual(ET.tostring(root.find("afalignment")), original_alignment)
        for old, new in (("mu", "nu"), ("gamma", "zeta")):
            old_param = root.find(f"flipflopModel/{old}/parameter")
            new_param = root.find(f"flipflopModel/{new}/parameter")
            self.assertEqual(new_param.get("value"), old_param.get("value"))
            self.assertEqual(new_param.get("lower"), old_param.get("lower"))
            self.assertEqual(new_param.get("id"), f"flipflop.{new}")
            for path in ("operators/scaleOperator", "mcmc/posterior/prior/halfNormalPrior"):
                old_node = root.find(f"{path}/parameter[@idref='flipflop.{old}']/..")
                new_node = root.find(f"{path}/parameter[@idref='flipflop.{new}']/..")
                self.assertEqual({k: v for k, v in old_node.attrib.items() if k != "id"},
                                 {k: v for k, v in new_node.attrib.items() if k != "id"})
            self.assertIsNotNone(root.find(f"operators/upDownOperator/up/parameter[@idref='flipflop.{new}']"))
            self.assertIsNotNone(root.find(f"mcmc/log[@fileName]/parameter[@idref='flipflop.{new}']"))
        ids = [e.get("id") for e in root.iter() if e.get("id")]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertTrue(all(e.get("idref") in ids for e in root.iter() if e.get("idref")))

    def test_parameter_references_and_custom_ids(self):
        tree = ET.parse(EXAMPLE)
        root = tree.getroot()
        wrapper = root.find("flipflopModel/mu")
        parameter = wrapper.find("parameter")
        wrapper.remove(parameter)
        parameter.set("id", "custom.methylation")
        root.insert(0, parameter)
        ET.SubElement(wrapper, "parameter", idref="custom.methylation")
        for e in root.iter():
            if e.get("idref") == "flipflop.mu":
                e.set("idref", "custom.methylation")
        convert_tree(tree, "_new")
        new_parameter = root.find("flipflopModel/nu/parameter")
        self.assertEqual(new_parameter.get("id"), "custom.methylation.nu")
        self.assertEqual(new_parameter.get("value"), "0.1")

    def test_suffixes_include_mle_readers(self):
        tree = ET.parse(EXAMPLE)
        root = tree.getroot()
        for tag in ("log", "pathSamplingAnalysis", "steppingStoneSamplingAnalysis"):
            ET.SubElement(root, tag, fileName="results/run.MLE.log")
        convert_tree(tree, "_independent")
        for element in root:
            if element.get("fileName"):
                self.assertEqual(element.get("fileName"), "results/run.MLE_independent.log")
        self.assertEqual(root.find("mcmc").get("operatorAnalysis"), "sim1Absolute_independent.ops")
        self.assertEqual(root.find("mcmc/log[@fileName]").get("fileName"), "sim1Absolute_independent.log")
        self.assertEqual(root.find("mcmc/logTree").get("fileName"), "sim1Absolute_independent.trees")
        self.assertEqual(suffixed_name("a.b/run.xml", "_4param"), "a.b/run_4param.xml")
        self.assertEqual(suffixed_name("run", "_4param"), "run_4param")

    def test_invalid_inputs(self):
        for failure in ("nu", "prior", "operator", "value", "collision"):
            with self.subTest(failure=failure):
                tree = ET.parse(EXAMPLE)
                root = tree.getroot()
                if failure == "nu":
                    ET.SubElement(root.find("flipflopModel"), "nu")
                elif failure == "prior":
                    prior = root.find("mcmc/posterior/prior")
                    prior.remove(prior.find("halfNormalPrior"))
                elif failure == "operator":
                    operators = root.find("operators")
                    operators.remove(operators.find("scaleOperator/parameter[@idref='flipflop.mu']/.."))
                elif failure == "value":
                    root.find("flipflopModel/mu/parameter").set("value", "nan")
                else:
                    ET.SubElement(root, "parameter", id="flipflop.nu", value="0.2")
                with self.assertRaises(ValueError):
                    convert_tree(tree, "_4param")
        for suffix in ("", "../bad", "a\\b"):
            with self.assertRaises(ValueError):
                convert_tree(ET.parse(EXAMPLE), suffix)

    def test_cli_paths_and_overwrite_protection(self):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "original.xml"
            source.write_bytes(EXAMPLE.read_bytes())
            original = source.read_bytes()

            def run(*args):
                return subprocess.run([sys.executable, str(SCRIPT), str(source), *map(str, args)],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                      universal_newlines=True)

            result = run()
            self.assertEqual(result.returncode, 0, result.stderr)
            default_output = source.with_name("original_4param.xml")
            self.assertTrue(default_output.exists())
            self.assertIn("Manual alternative", default_output.read_text())  # Comments preserved.
            self.assertEqual(source.read_bytes(), original)
            converted = default_output.read_bytes()
            self.assertNotEqual(run().returncode, 0)
            self.assertEqual(default_output.read_bytes(), converted)
            custom_output = Path(directory) / "chosen.xml"
            result = run(custom_output, "--suffix", "_custom")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(ET.parse(custom_output).find("mcmc").get("operatorAnalysis"),
                             "sim1Absolute_custom.ops")
            self.assertNotEqual(run(source).returncode, 0)
            self.assertEqual(source.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()

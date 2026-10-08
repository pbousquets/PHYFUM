"""Check four-rate XML wiring and backwards-compatible generator defaults.

Run with: python -m unittest discover -s simulate/XMLgenerator -p 'test_*.py'
Requires the generator's usual click, pandas, and yattag dependencies.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import xml.etree.ElementTree as ET

from click.testing import CliRunner

from createXML import createXML
from methylationBetas2xml import main


class FourRateXMLTest(unittest.TestCase):
    def generate(self, **kwargs):
        generator = createXML(age=45, mu=0.03, gamma=0.07, **kwargs)
        generator.addSamples(SimpleNamespace(
            sample_names=["A", "B", "C"], samples=["0.1,0.5", "0.2,0.6", "0.3,0.7"]
        ))
        generator.buildDoc("test", iterations=10, sampling=1, screenSampling=1)
        return ET.fromstring(str(generator))

    def assert_rate_wiring(self, root):
        for new, old in (("nu", "mu"), ("zeta", "gamma")):
            parameter = root.find(f"flipflopModel/{new}/parameter")
            self.assertEqual(parameter.get("id"), f"flipflop.{new}")
            self.assertEqual(parameter.get("lower"), "0.0")
            self.assertIsNone(parameter.get("idref"))
            # Each rate must have one independent proposal and one matching proper prior.
            for parent, tag in (("operators", "scaleOperator"), ("mcmc/posterior/prior", "halfNormalPrior")):
                new_nodes = root.findall(f"{parent}/{tag}/parameter[@idref='flipflop.{new}']/..")
                old_nodes = root.findall(f"{parent}/{tag}/parameter[@idref='flipflop.{old}']/..")
                self.assertEqual(len(new_nodes), 1)
                self.assertEqual(new_nodes[0].attrib, old_nodes[0].attrib)
            self.assertIsNotNone(root.find(
                f"operators/upDownOperator/up/parameter[@idref='flipflop.{new}']"
            ))
            self.assertEqual(len(root.findall(
                f"mcmc/log[@fileName]/parameter[@idref='flipflop.{new}']"
            )), 1)
        ids = [node.get("id") for node in root.iter() if node.get("id")]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertTrue(all(node.get("idref") in ids for node in root.iter() if node.get("idref")))

    def test_original_model_remains_default(self):
        root = self.generate()
        self.assertIsNone(root.find("flipflopModel/nu"))
        self.assertIsNone(root.find("flipflopModel/zeta"))
        for rate in ("nu", "zeta"):
            self.assertEqual(root.findall(f".//*[@idref='flipflop.{rate}']"), [])

    def test_four_rate_defaults_and_overrides(self):
        for overrides in ({}, {"nu": 0.04}, {"zeta": 0.08}, {"nu": 0.04, "zeta": 0.08}):
            with self.subTest(overrides=overrides):
                root = self.generate(four_rates=True, **overrides)
                self.assert_rate_wiring(root)
                for rate, default in (("nu", 0.03), ("zeta", 0.07)):
                    self.assertEqual(float(root.find(f"flipflopModel/{rate}/parameter").get("value")),
                                     overrides.get(rate, default))

    def test_cli_and_mode_validation(self):
        runner = CliRunner()
        sample = str(Path(__file__).with_name("test.csv"))
        with TemporaryDirectory() as directory:
            output = str(Path(directory) / "four")
            args = ["--input", sample, "--age", "45", "--output", output,
                    "--four-rates", "--nu", "0.04", "--zeta", "0.08"]
            result = runner.invoke(main, args)
            self.assertEqual(result.exit_code, 0, result.output)
            root = ET.parse(output + ".xml").getroot()
            self.assert_rate_wiring(root)
            self.assertEqual(float(root.find("flipflopModel/nu/parameter").get("value")), 0.04)
            self.assertEqual(float(root.find("flipflopModel/zeta/parameter").get("value")), 0.08)
            for rate in ("nu", "zeta"):
                result = runner.invoke(main, ["--input", sample, "--age", "45", f"--{rate}", "0.04"])
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn("require --four-rates", result.output)
                with self.assertRaises(ValueError):
                    createXML(age=45, **{rate: 0.04})

    def test_example_preserves_original_analysis(self):
        examples = Path(__file__).resolve().parents[2] / "examples/release/flipflop"
        original = ET.parse(examples / "sim1Absolute.xml").getroot()
        modified = ET.parse(examples / "sim1AbsoluteFourRates.xml").getroot()
        self.assert_rate_wiring(modified)
        self.assertEqual(modified.find("mcmc").get("operatorAnalysis"), "sim1AbsoluteFourRates.ops")
        for tag in ("log", "logTree"):
            self.assertTrue(modified.find(f"mcmc/{tag}[@fileName]").get("fileName")
                            .startswith("sim1AbsoluteFourRates."))
        # Removing the new parameters and wiring must recover the original analysis.
        model = modified.find("flipflopModel")
        for rate in ("nu", "zeta"):
            model.remove(model.find(rate))
            for parent in modified.iter():
                for child in list(parent):
                    if child.get("idref") == f"flipflop.{rate}":
                        parent.remove(child)
            for path in ("operators", "mcmc/posterior/prior"):
                parent = modified.find(path)
                for child in list(parent):
                    if child.tag in ("scaleOperator", "halfNormalPrior") and len(child) == 0:
                        parent.remove(child)
        for tree in (original, modified):
            for node in tree.iter():
                node.text = (node.text or "").strip()
                node.tail = ""
                for name in ("fileName", "operatorAnalysis"):
                    if name in node.attrib:
                        node.set(name, node.get(name).replace("sim1AbsoluteFourRates.", "sim1Absolute."))
        self.assertEqual(ET.tostring(original), ET.tostring(modified))


if __name__ == "__main__":
    unittest.main()

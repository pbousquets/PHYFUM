"""Fixed-LCA XML wiring for the bounded, reflected root-height proposal."""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import xml.etree.ElementTree as ET

from createXML import createXML
from click.testing import CliRunner
from methylationBetas2xml import main


class ReflectedRootXMLTest(unittest.TestCase):
    def assert_root_move(self, root):
        height = root.find("cenancestorTreeLikelihood/cenancestorHeight/parameter")
        self.assertEqual(root.find("cenancestorTreeLikelihood").get("heightRules"), "true")
        moves = root.findall("operators/randomWalkOperator")
        self.assertEqual(len(moves), 1)
        move = moves[0]
        self.assertEqual(move.get("boundaryCondition"), "reflecting")
        self.assertEqual(float(move.get("upper")), float(height.get("value")))
        # Do not override the root parameter's live child-height lower bound.
        self.assertIsNone(move.get("lower"))
        self.assertEqual(float(move.get("windowSize")), 2.5)
        self.assertEqual(float(move.get("weight")), 3.0)
        self.assertEqual(move.get("autoOptimize"), "true")
        self.assertEqual(move.find("parameter").get("idref"), "treeModel.rootHeight")
        self.assertIsNotNone(root.find("operators/scaleOperator/parameter[@idref='treeModel.rootHeight']"))
        self.assertIsNotNone(root.find("operators/upDownOperator"))
        self.assertIsNotNone(root.find("operators/subtreeSlide"))
        self.assertEqual(root.findall("operators/.//parameter[@idref='luca_branch']"), [])
        self.assertEqual(root.findall("operators/.//parameter[@idref='luca_height']"), [])

    def test_both_rate_modes_and_different_ages(self):
        for four_rates in (False, True):
            for age in (0.1, 12.5, 45, 90):
                with self.subTest(four_rates=four_rates, age=age):
                    generator = createXML(age=age, four_rates=four_rates, mrca_sampling=True)
                    generator.addSamples(SimpleNamespace(
                        sample_names=["A", "B"], samples=["0.1,0.5", "0.2,0.6"]))
                    generator.buildDoc("test", iterations=10, sampling=1, screenSampling=1)
                    self.assert_root_move(ET.fromstring(str(generator)))

    def test_cli_default_and_opt_in(self):
        sample = Path(__file__).with_name("test.csv")
        with TemporaryDirectory() as directory:
            for four_rates in (False, True):
                for enable in (False, True):
                    with self.subTest(four_rates=four_rates, enable=enable):
                        output = Path(directory) / f"test-{four_rates}-{enable}"
                        args = ["--input", str(sample), "--age", "45", "--output", str(output)]
                        if four_rates:
                            args.append("--four-rates")
                        if enable:
                            args.append("--mrca-sampling")
                        result = CliRunner().invoke(main, args)
                        self.assertEqual(result.exit_code, 0, result.output)
                        root = ET.parse(str(output) + ".xml").getroot()
                        if enable:
                            self.assert_root_move(root)
                        else:
                            self.assertIsNone(root.find("operators/randomWalkOperator"))
                            self.assertIsNotNone(root.find(
                                "operators/scaleOperator/parameter[@idref='treeModel.rootHeight']"))

    def test_four_rate_example(self):
        example = Path(__file__).resolve().parents[2] / "examples/release/flipflop/sim1AbsoluteFourRates.xml"
        self.assert_root_move(ET.parse(example).getroot())


if __name__ == "__main__":
    unittest.main()

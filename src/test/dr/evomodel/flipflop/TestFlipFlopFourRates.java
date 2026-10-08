package test.dr.evomodel.flipflop;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodelxml.substmodel.FlipFlopModelParser;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.XMLObject;
import dr.xml.XMLParser;
import junit.framework.TestCase;

import java.io.StringReader;
import javax.xml.parsers.DocumentBuilderFactory;

/** Regression tests for separate methylation and demethylation rates and XML defaults. */
public class TestFlipFlopFourRates extends TestCase {
    private static final double TOLERANCE = 1e-10;

    private FlipFlopModel model(Parameter mu, Parameter gamma, Parameter nu, Parameter zeta,
                               boolean normalize, boolean useFrequencies) {
        AFsequence dataType = new AFsequence(6);
        FrequencyModel frequencies = new FrequencyModel(dataType,
                new double[]{0.1, 0.2, 0.1, 0.3, 0.2, 0.1});
        Parameter cells = new Parameter.Default(2.0);
        Parameter lambda = new Parameter.Default(0.7);
        if (nu == null && zeta == null) {
            return new FlipFlopModel("legacy", dataType, cells, gamma, lambda, mu,
                    normalize, useFrequencies, frequencies);
        }
        return new FlipFlopModel("fourRates", dataType, cells, gamma, lambda, mu, nu, zeta,
                normalize, useFrequencies, frequencies);
    }

    private double[] probabilities(FlipFlopModel model) {
        double[] result = new double[36];
        model.getTransitionProbabilities(0.8, result);
        return result;
    }

    private void assertProbabilitiesEqual(double[] expected, double[] actual) {
        for (int i = 0; i < expected.length; i++) {
            assertEquals("Transition entry " + i, expected[i], actual[i], TOLERANCE);
        }
    }

    public void testIndependentRateMatrix() {
        FlipFlopModel model = model(new Parameter.Default(0.3), new Parameter.Default(0.2),
                new Parameter.Default(0.4), new Parameter.Default(0.5), false, false);
        model.setupMatrix();
        // Check the explicit diagonal before setupMatrix's makeValid can repair it.
        model.storeIntoAmat();
        // Source states are rows: (0,0), (1,0), (2,0), (0,1), (1,1), (0,2).
        // Hand-calculated from the four-rate-params reference (whose matrix is transposed).
        double[][] expected = {
                {-1.2, 1.2, 0, 0, 0, 0},
                {1.2, -2.9, 1.3, 0.4, 0, 0},
                {0, 1.0, -1.8, 0, 0.8, 0},
                {0.7, 0.4, 0, -2.4, 0.6, 0.7},
                {0, 0, 1.1, 0.5, -2.7, 1.1},
                {0, 0, 0, 0, 0.8, -0.8}
        };
        for (int i = 0; i < expected.length; i++) {
            double sum = 0;
            for (int j = 0; j < expected.length; j++) {
                double rate = model.getAmat()[i][j];
                assertEquals("Rate " + i + "," + j, expected[i][j], rate, TOLERANCE);
                sum += rate;
            }
            assertEquals(0.0, sum, TOLERANCE);
        }
        double[] transition = probabilities(model);
        for (int i = 0; i < 6; i++) {
            double sum = 0;
            for (int j = 0; j < 6; j++) {
                assertTrue(transition[6*i+j] >= 0.0);
                sum += transition[6*i+j];
            }
            assertEquals(1.0, sum, TOLERANCE);
        }
    }

    public void testLegacyConstructorKeepsRatesTied() {
        for (boolean normalize : new boolean[]{false, true}) {
            for (boolean useFrequencies : new boolean[]{false, true}) {
                Parameter mu = new Parameter.Default(0.3);
                Parameter gamma = new Parameter.Default(0.2);
                FlipFlopModel legacy = model(mu, gamma, null, null, normalize, useFrequencies);
                probabilities(legacy); // Populate the cache before changing the shared parameters.
                mu.setParameterValue(0, 0.6);
                gamma.setParameterValue(0, 0.8);
                FlipFlopModel explicit = model(mu, gamma, mu, gamma, normalize, useFrequencies);
                assertEquals(4, legacy.getVariableCount());
                assertProbabilitiesEqual(probabilities(explicit), probabilities(legacy));
            }
        }
    }

    public void testIndependentUpdatesAndRestore() {
        Parameter mu = new Parameter.Default(0.3);
        Parameter gamma = new Parameter.Default(0.2);
        Parameter nu = new Parameter.Default(0.4);
        Parameter zeta = new Parameter.Default(0.5);
        FlipFlopModel model = model(mu, gamma, nu, zeta, false, false);
        double[] original = probabilities(model);
        model.storeModelState();
        nu.setParameterValue(0, 0.9);
        double[] changed = probabilities(model);
        assertTrue(Math.abs(original[0] - changed[0]) > 1e-6);
        assertEquals(0.9, model.getAmat()[1][3], TOLERANCE);
        assertEquals(1.2, model.getAmat()[1][0], TOLERANCE);
        zeta.setParameterValue(0, 0.8);
        probabilities(model);
        assertEquals(1.5, model.getAmat()[1][0], TOLERANCE);
        assertEquals(0.9, model.getAmat()[1][3], TOLERANCE);
        assertEquals(1.2, model.getAmat()[0][1], TOLERANCE);
        assertEquals(0.8, model.getAmat()[5][4], TOLERANCE);
        model.restoreModelState();
        assertEquals(0.4, nu.getParameterValue(0), TOLERANCE);
        assertEquals(0.5, zeta.getParameterValue(0), TOLERANCE);
        assertProbabilitiesEqual(original, probabilities(model));
    }

    private XMLParser parseModel(String extraRates) throws Exception {
        XMLParser parser = new XMLParser(false, true);
        parser.addXMLObjectParser(new ParameterParser());
        parser.addXMLObjectParser(new FlipFlopModelParser());
        // Supply an already parsed pattern list; the model and parameters use real XML parsing.
        XMLObject patterns = new XMLObject(DocumentBuilderFactory.newInstance()
                .newDocumentBuilder().newDocument().createElement("patterns"));
        patterns.setNativeObject(new Patterns(new AFsequence(6)));
        parser.getObjectStore().put("patterns", patterns);
        String xml = "<beast><flipflopModel id='model' normalize='false'>"
                + "<patterns idref='patterns'/>"
                + "<stemCells><parameter value='2'/></stemCells>"
                + "<gamma><parameter id='gamma' value='0.2'/></gamma>"
                + "<lambda><parameter value='0.7'/></lambda>"
                + "<mu><parameter id='mu' value='0.3'/></mu>"
                + extraRates + "</flipflopModel></beast>";
        parser.parse(new StringReader(xml), false);
        return parser;
    }

    public void testXmlDefaultsAndIndependentOverrides() throws Exception {
        // Each optional rate may be supplied independently.
        for (boolean hasNu : new boolean[]{false, true}) {
            for (boolean hasZeta : new boolean[]{false, true}) {
                XMLParser parser = parseModel(
                        (hasNu ? "<nu><parameter id='nu' value='0.4'/></nu>" : "")
                        + (hasZeta ? "<zeta><parameter id='zeta' value='0.5'/></zeta>" : ""));
                FlipFlopModel parsed = (FlipFlopModel) parser.getObjectStore().get("model").getNativeObject();
                assertEquals(4 + (hasNu ? 1 : 0) + (hasZeta ? 1 : 0), parsed.getVariableCount());
                probabilities(parsed);
                assertEquals(hasNu ? 0.4 : 0.3, parsed.getAmat()[1][3], TOLERANCE);
                assertEquals(0.7 + (hasZeta ? 0.5 : 0.2), parsed.getAmat()[1][0], TOLERANCE);
                Parameter mu = (Parameter) parser.getObjectStore().get("mu").getNativeObject();
                Parameter gamma = (Parameter) parser.getObjectStore().get("gamma").getNativeObject();
                mu.setParameterValue(0, 0.6);
                gamma.setParameterValue(0, 0.8);
                probabilities(parsed);
                assertEquals(hasNu ? 0.4 : 0.6, parsed.getAmat()[1][3], TOLERANCE);
                assertEquals(0.7 + (hasZeta ? 0.5 : 0.8), parsed.getAmat()[1][0], TOLERANCE);
            }
        }
    }

    public void testXmlExplicitSharedReferences() throws Exception {
        XMLParser parser = parseModel("<nu><parameter idref='mu'/></nu>"
                + "<zeta><parameter idref='gamma'/></zeta>");
        FlipFlopModel parsed = (FlipFlopModel) parser.getObjectStore().get("model").getNativeObject();
        assertEquals(4, parsed.getVariableCount());
        assertProbabilitiesEqual(probabilities(model(new Parameter.Default(0.3),
                new Parameter.Default(0.2), null, null, false, false)), probabilities(parsed));
    }
}

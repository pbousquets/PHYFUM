package test.dr.evomodel.flipflop;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.AFsequence;
import dr.evolution.io.NewickImporter;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockCenancestorBranchRates;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.*;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.CenancestorTreeLikelihood;
import dr.evomodel.treelikelihood.FlipFlopErrorModel;
import dr.evomodelxml.substmodel.FlipFlopNeoplasmCenancestorFrequencyParser;
import dr.evomodelxml.treelikelihood.CenancestorTreeLikelihoodParser;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import junit.framework.TestCase;

import java.io.StringReader;
import javax.xml.parsers.DocumentBuilderFactory;

/** Strict XML parsing and integration with the actual cenancestor likelihood. */
public class TestFlipFlopNeoplasmParser extends TestCase {
    private static final String HEIGHTS = "<originHeight><parameter id='origin' value='5'/></originHeight>"
            + "<cenancestorHeight><parameter id='height' value='2'/></cenancestorHeight>";
    private static final String ROOT = "<frequencies><parameter id='root' value='0.5 0 0.5'/></frequencies>";

    private void store(XMLParser parser, String id, Object object) throws Exception {
        XMLObject xo = new XMLObject(DocumentBuilderFactory.newInstance().newDocumentBuilder()
                .newDocument().createElement(id));
        xo.setNativeObject(object);
        parser.getObjectStore().put(id, xo);
    }

    private String modelXml(String attributes, String origin, String time, String extraRates) {
        return "<flipflopNeoplasmCenancestorFrequency id='neoplasm' " + attributes + ">"
                + "<patterns idref='patterns'/>" + origin + time
                + "<mu><parameter id='mu' value='0.03'/></mu>"
                + "<gamma><parameter id='gamma' value='0.02'/></gamma>"
                + extraRates + "</flipflopNeoplasmCenancestorFrequency>";
    }

    private XMLParser parse(String attributes, String origin, String time, String extraRates) throws Exception {
        return parse(attributes, origin, time, extraRates, new Patterns(new AFsequence(6)));
    }

    private XMLParser parse(String attributes, String origin, String time, String extraRates,
            Patterns patterns) throws Exception {
        XMLParser parser = new XMLParser(false, true);
        parser.addXMLObjectParser(new ParameterParser());
        parser.addXMLObjectParser(new FlipFlopNeoplasmCenancestorFrequencyParser());
        store(parser, "patterns", patterns);
        parser.parse(new StringReader("<beast>" + modelXml(attributes, origin, time, extraRates) + "</beast>"), false);
        return parser;
    }

    private Patterns datedPatterns(Date... dates) {
        Taxa taxa = new Taxa();
        for (int i = 0; i < dates.length; i++) {
            Taxon taxon = new Taxon("sample" + i);
            // Supply metadata without changing BEAST's process-wide date cache.
            if (dates[i] != null) {
                taxon.setAttribute("date", dates[i]);
            }
            taxa.addTaxon(taxon);
        }
        return new Patterns(new AFsequence(6), taxa);
    }

    public void testOriginInferredFromLargestAlignmentAge() throws Exception {
        for (Date[] dates : new Date[][]{
                {new Date(40, Units.Type.YEARS, false), new Date(45, Units.Type.YEARS, false),
                        new Date(42, Units.Type.YEARS, false)},
                {new Date(45, Units.Type.YEARS, false), new Date(45, Units.Type.YEARS, false)}
        }) {
            XMLParser inferred = parse("", ROOT,
                    "<cenancestorHeight><parameter id='height' value='40'/></cenancestorHeight>", "", datedPatterns(dates));
            XMLParser expected = parse("", ROOT, HEIGHTS.replace("value='5'", "value='7'"), "");
            assertSameFrequencies(frequencies(expected).getFrequencies(), frequencies(inferred).getFrequencies());
            // The inferred origin is fixed, but the shared cenancestor height stays live.
            parameter(inferred, "height").setParameterValue(0, 44);
            expected = parse("", ROOT, HEIGHTS.replace("value='5'", "value='3'"), "");
            assertSameFrequencies(frequencies(expected).getFrequencies(), frequencies(inferred).getFrequencies());
        }
    }

    public void testZeroSampleAgeAndInferredOriginBeforeCenancestor() throws Exception {
        Patterns zeroAge = datedPatterns(new Date(0, Units.Type.YEARS, false));
        XMLParser parser = parse("", ROOT,
                "<cenancestorHeight><parameter value='0'/></cenancestorHeight>", "", zeroAge);
        assertSameFrequencies(new double[]{0.5, 0, 0, 0, 0, 0.5}, frequencies(parser).getFrequencies());
        try {
            parse("", ROOT, "<cenancestorHeight><parameter value='1'/></cenancestorHeight>", "", zeroAge);
            fail("Negative inferred duration accepted");
        } catch (XMLParseException expected) {
            assertTrue(expected.getMessage().contains("originHeight must be greater than or equal to cenancestorHeight"));
        }
    }

    public void testExplicitOriginOverridesSampleAges() throws Exception {
        double[] expected = frequencies(parse("", ROOT, HEIGHTS, "")).getFrequencies();
        for (Patterns patterns : new Patterns[]{
                datedPatterns(new Date(45, Units.Type.YEARS, false)),
                datedPatterns((Date) null),
                datedPatterns(new Date(45, Units.Type.YEARS, true))
        }) {
            assertSameFrequencies(expected, frequencies(parse("", ROOT, HEIGHTS, "", patterns)).getFrequencies());
        }
    }

    public void testIncompatibleDatesRequireExplicitOrigin() throws Exception {
        Date valid = new Date(45, Units.Type.YEARS, false);
        Patterns[] invalid = {
                new Patterns(new AFsequence(6)), datedPatterns(), datedPatterns(valid, null),
                datedPatterns(new Date(45, Units.Type.YEARS, true)),
                datedPatterns(new Date(-1, Units.Type.YEARS, false)),
                datedPatterns(new Date(Double.NaN, Units.Type.YEARS, false)),
                datedPatterns(new Date(Double.POSITIVE_INFINITY, Units.Type.YEARS, false)),
                datedPatterns(valid, new Date(46, Units.Type.MONTHS, false)),
                datedPatterns(valid, Date.createTimeSinceOrigin(46, Units.Type.YEARS, 1.0))
        };
        for (Patterns patterns : invalid) {
            try {
                parse("", ROOT, "<cenancestorHeight><parameter value='2'/></cenancestorHeight>", "", patterns);
                fail("Incompatible sample dates accepted");
            } catch (XMLParseException expected) {
                assertTrue(expected.getMessage().contains("explicit originHeight"));
            }
        }
    }

    private Object object(XMLParser parser, String id) {
        return parser.getObjectStore().get(id).getNativeObject();
    }

    private Parameter parameter(XMLParser parser, String id) {
        return (Parameter) object(parser, id);
    }

    private FlipFlopNeoplasmCenancestorFrequencyModel frequencies(XMLParser parser) {
        return (FlipFlopNeoplasmCenancestorFrequencyModel) object(parser, "neoplasm");
    }

    private void assertSameFrequencies(double[] expected, double[] actual) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], actual[i], 1e-12);
        }
    }

    public void testHeightsDetermineDuration() throws Exception {
        XMLParser heights = parse("", ROOT, HEIGHTS, "");
        Parameter mu = new Parameter.Default(0.03);
        Parameter gamma = new Parameter.Default(0.02);
        FlipFlopNeoplasmCenancestorFrequencyModel expected = new FlipFlopNeoplasmCenancestorFrequencyModel(
                new AFsequence(6), new Parameter.Default(new double[]{0.5, 0, 0.5}),
                new Parameter.Default(3.0), mu, gamma, mu, gamma);
        assertSameFrequencies(expected.getFrequencies(), frequencies(heights).getFrequencies());
        parameter(heights, "height").setParameterValueQuietly(0, 5);
        assertSameFrequencies(new double[]{0.5, 0, 0, 0, 0, 0.5}, frequencies(heights).getFrequencies());
    }

    public void testRemovedElapsedTimeRejectedWithoutStrictChecking() throws Exception {
        for (String time : new String[]{"<elapsedTime><parameter value='3'/></elapsedTime>",
                HEIGHTS + "<elapsedTime><parameter value='3'/></elapsedTime>"}) {
            XMLParser parser = new XMLParser(false, false);
            parser.addXMLObjectParser(new ParameterParser());
            parser.addXMLObjectParser(new FlipFlopNeoplasmCenancestorFrequencyParser());
            store(parser, "patterns", new Patterns(new AFsequence(6)));
            try {
                parser.parse(new StringReader("<beast>" + modelXml("", ROOT, time, "") + "</beast>"), false);
                fail("Removed elapsedTime option accepted");
            } catch (XMLParseException expected) {
                // No legacy elapsedTime input should silently change the intended duration.
                assertTrue(expected.getMessage().contains("cenancestorHeight"));
            }
        }
    }

    public void testTiedAndIndependentRates() throws Exception {
        for (String extra : new String[]{"", "<nu><parameter value='0.07'/></nu>",
                "<zeta><parameter value='0.04'/></zeta>",
                "<nu><parameter value='0.07'/></nu><zeta><parameter value='0.04'/></zeta>"}) {
            XMLParser parser = parse("", ROOT, HEIGHTS, extra);
            parameter(parser, "mu").setParameterValue(0, 0.06);
            parameter(parser, "gamma").setParameterValue(0, 0.08);
            FlipFlopNeoplasmCenancestorFrequencyModel expected = new FlipFlopNeoplasmCenancestorFrequencyModel(
                    new AFsequence(6), new Parameter.Default(new double[]{0.5, 0, 0.5}),
                    new Parameter.Default(3.0), new Parameter.Default(0.06), new Parameter.Default(0.08),
                    new Parameter.Default(extra.contains("<nu>") ? 0.07 : 0.06),
                    new Parameter.Default(extra.contains("<zeta>") ? 0.04 : 0.08));
            assertSameFrequencies(expected.getFrequencies(), frequencies(parser).getFrequencies());
        }
    }

    public void testProportionSchemesAndNormalization() throws Exception {
        double[] expected = frequencies(parse("", ROOT, HEIGHTS, "")).getFrequencies();
        assertSameFrequencies(expected, frequencies(parse("methylatedProportion='0.5'", "", HEIGHTS, "")).getFrequencies());
        XMLParser attribute = parse("methylatedProportion='0.5'",
                "<frequencies><parameter id='root' value='1'/></frequencies>", HEIGHTS, "");
        assertEquals(3, parameter(attribute, "root").getDimension());
        assertSameFrequencies(expected, frequencies(attribute).getFrequencies());
        XMLParser live = parse("", "<methylatedProportion><parameter id='p' value='0.5'/></methylatedProportion>", HEIGHTS, "");
        assertSameFrequencies(expected, frequencies(live).getFrequencies());
        parameter(live, "p").setParameterValue(0, 0.25);
        assertSameFrequencies(frequencies(parse("methylatedProportion='0.25'", "", HEIGHTS, "")).getFrequencies(),
                frequencies(live).getFrequencies());
        assertSameFrequencies(expected, frequencies(parse("normalize='true'",
                "<frequencies><parameter value='2 0 2'/></frequencies>", HEIGHTS, "")).getFrequencies());
    }

    public void testInvalidXmlRejected() throws Exception {
        String[][] invalid = {
                {"", "", HEIGHTS, ""},
                {"", "<frequencies><parameter value='0.5 0.5'/></frequencies>", HEIGHTS, ""},
                {"methylatedProportion='1.1'", "", HEIGHTS, ""},
                {"normalize='true'", "<frequencies><parameter value='0 0 0'/></frequencies>", HEIGHTS, ""},
                {"", ROOT, "<elapsedTime><parameter value='-1'/></elapsedTime>", ""},
                {"", ROOT, HEIGHTS + "<elapsedTime><parameter value='3'/></elapsedTime>", ""},
                {"", ROOT, "<originHeight><parameter value='5'/></originHeight>", ""},
                {"", ROOT, HEIGHTS, "<nu><parameter value='-0.1'/></nu>"},
                {"", ROOT + "<methylatedProportion><parameter value='0.5'/></methylatedProportion>", HEIGHTS, ""}
        };
        for (String[] input : invalid) {
            try {
                parse(input[0], input[1], input[2], input[3]);
                fail("Invalid neoplasm XML accepted");
            } catch (XMLParseException expected) {
                // Fail at parsing, rather than later in the MCMC.
            }
        }
    }

    public void testTimingXmlValidationAndInference() throws Exception {
        for (String coordinate : new String[]{"height", "branch", "fraction", "invalid"}) {
            XMLParser parser = new XMLParser(false, true);
            parser.addXMLObjectParser(new ParameterParser());
            parser.addXMLObjectParser(new dr.evomodelxml.tree.NeoplasmTimingParser());
            parser.addXMLObjectParser(new dr.evomodelxml.tree.NeoplasmTimingParameterParser());
            parser.addXMLObjectParser(new dr.evomodelxml.operators.NeoplasmTreeScaleOperatorParser());
            store(parser, "patterns", datedPatterns(new Date(5, Units.Type.YEARS, false)));
            store(parser, "treeModel", new TreeModel(new NewickImporter("(A:1,B:1);").importTree(null)));
            String xml = "<beast><neoplasmTiming id='timing' priorCoordinate='" + coordinate + "'>"
                    + "<treeModel idref='treeModel'/><patterns idref='patterns'/>"
                    + "<lcaFraction><parameter id='q' value='0.25'/></lcaFraction></neoplasmTiming>"
                    + "<neoplasmTimingParameter id='derived' quantity='height'><neoplasmTiming idref='timing'/></neoplasmTimingParameter>"
                    + "<neoplasmTreeScale id='scale' scaleFactor='0.8' weight='3'><neoplasmTiming idref='timing'/></neoplasmTreeScale></beast>";
            if (coordinate.equals("invalid")) {
                try { parser.parse(new StringReader(xml), false); fail(); }
                catch (XMLParseException expected) { }
            } else {
                parser.parse(new StringReader(xml), false);
                dr.evomodel.tree.NeoplasmTimingModel timing = (dr.evomodel.tree.NeoplasmTimingModel) object(parser, "timing");
                assertEquals(5.0, timing.getOriginHeight(), 0);
                assertSame(timing.getHeightParameter(), parameter(parser, "derived"));
                assertEquals(coordinate.equals("fraction") ? 0 : Math.log(4), timing.getLogLikelihood(), 1e-12);
                assertEquals(2.0, parameter(parser, "derived").getParameterValue(0), 1e-12);
            }
        }
    }

    public void testLikelihoodUpdatesInAllTimingModes() throws Exception {
        for (String mode : new String[]{"height", "branch", "fraction"}) {
            boolean fractionMode = mode.equals("fraction");
            boolean heightRules = !mode.equals("branch");
            XMLParser parser = new XMLParser(false, true);
            parser.addXMLObjectParser(new ParameterParser());
            parser.addXMLObjectParser(new FlipFlopNeoplasmCenancestorFrequencyParser());
            parser.addXMLObjectParser(new CenancestorTreeLikelihoodParser());
            parser.addXMLObjectParser(new dr.evomodelxml.tree.NeoplasmTimingParser());
            AFsequence type = new AFsequence(6);
            Taxa taxa = new Taxa();
            taxa.addTaxon(new Taxon("A"));
            taxa.addTaxon(new Taxon("B"));
            Patterns patterns = new Patterns(type, taxa);
            patterns.addPattern(new int[]{20, 40});
            TreeModel tree = new TreeModel(new NewickImporter("(A:1,B:1);").importTree(taxa));
            Parameter cells = new Parameter.Default(2.0);
            FlipFlopModel substitution = new FlipFlopModel("flipflop", type, cells,
                    new Parameter.Default(0.02), new Parameter.Default(0.7), new Parameter.Default(0.03),
                    false, false, new FrequencyModel(type, new double[]{0.1, 0.2, 0.2, 0.1, 0.2, 0.2}));
            store(parser, "patterns", patterns);
            store(parser, "treeModel", tree);
            store(parser, "siteModel", new GammaSiteModel(substitution));
            store(parser, "branchRates", new StrictClockCenancestorBranchRates(new Parameter.Default(1.0)));
            store(parser, "errorModel", new FlipFlopErrorModel(taxa, new Taxa(), cells,
                    new Parameter.Default(0.05), new Parameter.Default(0.95), new Parameter.Default(30.0)));
            String timingXml = "<neoplasmTiming id='timing'><treeModel idref='treeModel'/>"
                    + "<originHeight><parameter id='origin' value='5'/></originHeight>"
                    + "<lcaFraction><parameter id='q' value='0.25'/></lcaFraction></neoplasmTiming>";
            String xml = "<beast>" + (fractionMode ? timingXml : "")
                    + modelXml("", ROOT, fractionMode ? "<neoplasmTiming idref='timing'/>" : HEIGHTS,
                    "<nu><parameter id='nu' value='0.07'/></nu>")
                    + "<cenancestorTreeLikelihood id='likelihood'" + (fractionMode ? "" : " heightRules='" + heightRules + "'") + ">"
                    + "<patterns idref='patterns'/><treeModel idref='treeModel'/><siteModel idref='siteModel'/>"
                    + "<branchRates idref='branchRates'/><errorModel idref='errorModel'/>"
                    + "<flipflopNeoplasmCenancestorFrequency idref='neoplasm'/>"
                    + (fractionMode ? "" : "<cenancestorHeight><parameter idref='height'/></cenancestorHeight>"
                    + "<cenancestorBranch><parameter id='branch' value='1'/></cenancestorBranch>")
                    + "</cenancestorTreeLikelihood></beast>";
            parser.parse(new StringReader(xml), false);
            CenancestorTreeLikelihood likelihood = (CenancestorTreeLikelihood) object(parser, "likelihood");
            double original = likelihood.getLogLikelihood();
            assertTrue(Double.isFinite(original));
            // Follow the MCMC lifecycle: store, one evaluated proposal, then restore.
            // Only founder frequencies depend on nu here, testing listener propagation.
            Runnable[] proposals = {
                    () -> parameter(parser, "nu").setParameterValue(0, 0.3),
                    () -> parameter(parser, fractionMode ? "q" : heightRules ? "height" : "branch")
                            .setParameterValue(0, fractionMode ? 0.6 : heightRules ? 3 : 2),
                    () -> tree.setNodeHeight(tree.getRoot(), 1.5),
                    () -> parameter(parser, fractionMode ? "q" : "origin").setParameterValue(0, fractionMode ? 0.8 : 6),
                    () -> {
                        parameter(parser, "root").setParameterValue(0, 0.7);
                        parameter(parser, "root").setParameterValue(2, 0.3);
                    }
            };
            for (int i = 0; i < proposals.length; i++) {
                likelihood.storeModelState();
                proposals[i].run();
                double changed = likelihood.getLogLikelihood();
                assertTrue("Proposal " + i, Math.abs(changed - original) > 1e-6);
                assertLikelihoodMatchesPartials("Proposal " + i, likelihood, frequencies(parser));
                likelihood.restoreModelState();
                assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
                assertLikelihoodMatchesPartials("Restored proposal " + i, likelihood, frequencies(parser));
                likelihood.makeDirty();
                assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
            }
            if (fractionMode) {
                dr.evomodel.tree.NeoplasmTimingModel timing = frequencies(parser).getTimingModel();
                likelihood.storeModelState();
                new dr.evomodel.operators.NeoplasmTreeScaleOperator(timing, 0.8, 1,
                        dr.inference.operators.CoercionMode.COERCION_OFF).propose(1.1);
                double scaled = likelihood.getLogLikelihood();
                likelihood.makeDirty();
                assertEquals(scaled, likelihood.getLogLikelihood(), 1e-12);
                likelihood.restoreModelState();
                likelihood.makeDirty();
                assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
                likelihood.setUsed();
                try {
                    likelihood.getLogLikelihood();
                    fail("Omitted timing prior must not silently change the target");
                } catch (IllegalStateException expected) { }
                timing.setUsed();
                assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
            }
            // A defective founder generator must give a rejected proposal (-infinity),
            // not NaNs or an exception, and restoration must recover the likelihood.
            likelihood.storeModelState();
            parameter(parser, "gamma").setParameterValue(0, 0.0); // zeta shares gamma.
            parameter(parser, "nu").setParameterValue(0, 0.06); // equals 2*mu.
            assertEquals(Double.NEGATIVE_INFINITY, likelihood.getLogLikelihood(), 0.0);
            likelihood.restoreModelState();
            assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
            likelihood.makeDirty();
            assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
            likelihood.storeModelState();
            if (fractionMode) tree.setNodeHeight(tree.getRoot(), 6);
            else parameter(parser, heightRules ? "height" : "branch").setParameterValue(0, 6);
            assertEquals(Double.NEGATIVE_INFINITY, likelihood.getLogLikelihood(), 0.0);
            likelihood.restoreModelState();
            likelihood.makeDirty();
            assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
            if (!fractionMode && heightRules) {
                likelihood.storeModelState();
                tree.setNodeHeight(tree.getRoot(), 3); // Valid intree, but MRCA older than fixed C=2.
                assertEquals(Double.NEGATIVE_INFINITY, likelihood.getLogLikelihood(), 0.0);
                likelihood.restoreModelState();
                likelihood.makeDirty();
                assertEquals(original, likelihood.getLogLikelihood(), 1e-12);
            }
        }
    }

    private void assertLikelihoodMatchesPartials(String stage, CenancestorTreeLikelihood likelihood,
            FlipFlopNeoplasmCenancestorFrequencyModel model) {
        double actual = likelihood.getLogLikelihood();
        double[] partials = likelihood.getCenancestorPartials();
        double[] frequencies = model.getFrequencies();
        double probability = 0;
        for (int i = 0; i < frequencies.length; i++) {
            probability += frequencies[i] * partials[i];
        }
        assertEquals(stage, Math.log(probability), actual, 1e-12);
    }
}

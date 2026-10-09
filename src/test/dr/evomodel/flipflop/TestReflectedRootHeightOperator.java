package test.dr.evomodel.flipflop;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.AFsequence;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.StrictClockCenancestorBranchRates;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralFrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.CenancestorTreeLikelihood;
import dr.evomodel.treelikelihood.FlipFlopErrorModel;
import dr.inference.model.Parameter;
import dr.inference.operators.RandomWalkOperator;
import dr.inferencexml.operators.RandomWalkOperatorParser;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import dr.xml.XMLParser;
import junit.framework.TestCase;

import java.io.StringReader;
import javax.xml.parsers.DocumentBuilderFactory;

/** Exercise the existing BEAST operator as configured for PHYFUM's fixed LCA. */
public class TestReflectedRootHeightOperator extends TestCase {
    private static final double LCA_HEIGHT = 45.0;

    private Taxa taxa() {
        Taxa taxa = new Taxa();
        for (String name : new String[]{"A", "B", "C"}) {
            taxa.addTaxon(new Taxon(name));
        }
        return taxa;
    }

    private TreeModel tree() throws Exception {
        return new TreeModel(new NewickImporter("((A:30,B:30):14,C:44);").importTree(taxa()));
    }

    private RandomWalkOperator operator(TreeModel tree, double window) throws Exception {
        XMLParser parser = new XMLParser(false, true);
        parser.addXMLObjectParser(new RandomWalkOperatorParser());
        XMLObject root = new XMLObject(DocumentBuilderFactory.newInstance().newDocumentBuilder()
                .newDocument().createElement("parameter"));
        root.setNativeObject(tree.getRootHeightParameter());
        parser.getObjectStore().put("treeModel.rootHeight", root);
        parser.parse(new StringReader("<beast><randomWalkOperator id='rootWalk' windowSize='"
                + window + "' weight='3' upper='45' boundaryCondition='reflecting' autoOptimize='true'>"
                + "<parameter idref='treeModel.rootHeight'/></randomWalkOperator></beast>"), false);
        return (RandomWalkOperator) parser.getObjectStore().get("rootWalk").getNativeObject();
    }

    private double oldestChild(TreeModel tree) {
        NodeRef root = tree.getRoot();
        return Math.max(tree.getNodeHeight(tree.getChild(root, 0)),
                tree.getNodeHeight(tree.getChild(root, 1)));
    }

    private void assertValidRoot(TreeModel tree) throws Exception {
        double height = tree.getNodeHeight(tree.getRoot());
        assertTrue(Double.isFinite(height));
        assertTrue(height >= oldestChild(tree));
        assertTrue(height <= LCA_HEIGHT);
        tree.checkTreeIsValid();
    }

    /** Independent triangular-wave expression, including multiple reflections. */
    private double reflected(double value, double lower) {
        double width = LCA_HEIGHT - lower;
        double offset = (value - lower) % (2.0 * width);
        if (offset < 0) {
            offset += 2.0 * width;
        }
        return lower + (offset <= width ? offset : 2.0 * width - offset);
    }

    public void testParsedOperatorReflectsBothBoundariesAndLargeSteps() throws Exception {
        TreeModel tree = tree();
        RandomWalkOperator operator = operator(tree, 1000.0);
        assertSame(tree.getRootHeightParameter(), operator.getParameter());
        assertEquals(3.0, operator.getWeight(), 0.0);
        boolean crossedLower = false;
        boolean crossedUpper = false;
        for (int i = 0; i < 500; i++) {
            // Include proposals starting exactly at both valid endpoints.
            if (i < 2) {
                tree.setNodeHeight(tree.getRoot(), i == 0 ? 30.0 : LCA_HEIGHT);
            }
            long seed = 12345L + i;
            MathUtils.setSeed(seed);
            MathUtils.nextInt(1); // The operator first chooses the parameter dimension.
            double candidate = tree.getNodeHeight(tree.getRoot())
                    + (2.0 * MathUtils.nextDouble() - 1.0) * operator.getWindowSize();
            crossedLower |= candidate < oldestChild(tree);
            crossedUpper |= candidate > LCA_HEIGHT;
            double expected = reflected(candidate, oldestChild(tree));
            MathUtils.setSeed(seed);
            assertEquals(0.0, operator.doOperation(), 0.0);
            assertEquals(expected, tree.getNodeHeight(tree.getRoot()), 1e-10);
            assertValidRoot(tree);
        }
        assertTrue(crossedLower);
        assertTrue(crossedUpper);
    }

    public void testBoundsFollowChildrenAndRootChangingTopology() throws Exception {
        TreeModel tree = tree();
        Parameter originalRootParameter = tree.getRootHeightParameter();
        RandomWalkOperator operator = operator(tree, 2.5);
        NodeRef oldRoot = tree.getRoot();
        NodeRef inner = tree.getChild(oldRoot, 0);
        tree.setNodeHeight(inner, 43.0);
        MathUtils.setSeed(56789);
        for (int i = 0; i < 100; i++) {
            assertEquals(0.0, operator.doOperation(), 0.0);
            assertValidRoot(tree);
        }

        // Rotate ((A,B),C) into (A,(B,C)), making a different node the root.
        NodeRef b = tree.getChild(inner, 1);
        NodeRef c = tree.getChild(oldRoot, 1);
        tree.beginTreeEdit();
        tree.removeChild(oldRoot, inner);
        tree.removeChild(oldRoot, c);
        tree.removeChild(inner, b);
        tree.addChild(oldRoot, b);
        tree.addChild(oldRoot, c);
        tree.addChild(inner, oldRoot);
        tree.setNodeHeight(oldRoot, 35.0);
        tree.setNodeHeight(inner, 44.0);
        tree.setRoot(inner);
        tree.endTreeEdit();
        assertSame(inner, tree.getRoot());
        assertSame(originalRootParameter, tree.getRootHeightParameter());
        assertEquals(35.0, oldestChild(tree), 0.0);
        for (int i = 0; i < 100; i++) {
            operator.doOperation();
            assertValidRoot(tree);
            assertEquals(35.0, tree.getNodeHeight(oldRoot), 0.0);
        }
    }

    public void testNarrowIntervalAndDatedTips() throws Exception {
        TreeModel tree = new TreeModel(new NewickImporter("((A:20,B:30):14,C:40);")
                .importTree(taxa()));
        NodeRef inner = tree.getChild(tree.getRoot(), 0);
        tree.setNodeHeight(tree.getRoot(), LCA_HEIGHT);
        tree.setNodeHeight(inner, LCA_HEIGHT - 1e-6);
        double[] tips = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tips.length; i++) {
            tips[i] = tree.getNodeHeight(tree.getExternalNode(i));
        }
        RandomWalkOperator operator = operator(tree, 2.5);
        MathUtils.setSeed(98765);
        for (int i = 0; i < 100; i++) {
            operator.doOperation();
            assertValidRoot(tree);
            for (int j = 0; j < tips.length; j++) {
                assertEquals(tips[j], tree.getNodeHeight(tree.getExternalNode(j)), 0.0);
            }
        }
    }

    public void testDerivedStemLikelihoodAndStateRestoration() throws Exception {
        TreeModel tree = tree();
        Taxa taxa = taxa();
        AFsequence type = new AFsequence(6);
        Patterns patterns = new Patterns(type, taxa);
        patterns.addPattern(new int[]{20, 40, 60});
        Parameter cells = new Parameter.Default(2.0);
        FlipFlopModel substitution = new FlipFlopModel("flipflop", type, cells,
                new Parameter.Default(0.02), new Parameter.Default(0.7), new Parameter.Default(0.03),
                false, false, new FrequencyModel(type, new double[]{0.1, 0.2, 0.2, 0.1, 0.2, 0.2}));
        Parameter height = new Parameter.Default(LCA_HEIGHT);
        Parameter stem = new Parameter.Default(1.0);
        height.addBounds(new Parameter.DefaultBounds(LCA_HEIGHT, 0.0, 1));
        stem.addBounds(new Parameter.DefaultBounds(LCA_HEIGHT, 0.0, 1));
        CenancestorTreeLikelihood likelihood = new CenancestorTreeLikelihood(patterns, tree,
                new GammaSiteModel(substitution),
                new StrictClockCenancestorBranchRates(new Parameter.Default(1.0)),
                new FlipFlopErrorModel(taxa, new Taxa(), cells, new Parameter.Default(0.05),
                        new Parameter.Default(0.95), new Parameter.Default(30.0)),
                height, stem, new GeneralFrequencyModel(type, new double[]{0.5, 0, 0, 0, 0, 0.5}),
                false, false, true, true, false, true);
        RandomWalkOperator operator = operator(tree, 2.5);
        MathUtils.setSeed(45678);
        for (int i = 0; i < 100; i++) {
            double beforeRoot = tree.getNodeHeight(tree.getRoot());
            double beforeStem = stem.getParameterValue(0);
            double beforeLikelihood = likelihood.getLogLikelihood();
            likelihood.storeModelState();
            assertEquals(0.0, operator.doOperation(), 0.0);
            assertValidRoot(tree);
            assertEquals(LCA_HEIGHT - tree.getNodeHeight(tree.getRoot()), stem.getParameterValue(0), 0.0);
            assertEquals(LCA_HEIGHT, height.getParameterValue(0), 0.0);
            double proposed = likelihood.getLogLikelihood();
            assertTrue(Double.isFinite(proposed));
            likelihood.makeDirty();
            assertEquals(proposed, likelihood.getLogLikelihood(), 1e-10);
            if (i % 2 == 0) {
                likelihood.restoreModelState();
                assertEquals(beforeRoot, tree.getNodeHeight(tree.getRoot()), 0.0);
                assertEquals(beforeStem, stem.getParameterValue(0), 0.0);
                assertEquals(beforeLikelihood, likelihood.getLogLikelihood(), 1e-10);
                likelihood.makeDirty();
                assertEquals(beforeLikelihood, likelihood.getLogLikelihood(), 1e-10);
            } else {
                likelihood.acceptModelState();
            }
        }
    }
}

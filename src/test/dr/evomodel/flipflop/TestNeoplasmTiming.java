package test.dr.evomodel.flipflop;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.*;
import dr.evomodel.operators.*;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import junit.framework.TestCase;

/** Geometry, reversibility, topology changes and prior-only sampling checks. */
public class TestNeoplasmTiming extends TestCase {
    private NeoplasmTimingModel timing(String newick, double origin, double fraction, boolean jacobian) throws Exception {
        return new NeoplasmTimingModel(new TreeModel(new NewickImporter(newick).importTree(null)),
                new Parameter.Default(origin), new Parameter.Default(fraction), jacobian);
    }
    private NeoplasmTreeScaleOperator scale(NeoplasmTimingModel t) {
        return new NeoplasmTreeScaleOperator(t, 0.8, 1, CoercionMode.COERCION_OFF);
    }
    private void checkGeometry(NeoplasmTimingModel t) {
        assertTrue(t.isValid());
        assertEquals(t.getOriginHeight(), t.getRootHeight() + t.getBranchLength() + t.getFounderDuration(), 1e-12);
        assertTrue(t.getBranchLength() >= 0 && t.getFounderDuration() >= 0);
    }

    public void testDerivedTimesAndStateRestoration() throws Exception {
        NeoplasmTimingModel t = timing("(A:20,B:20);", 45, 0.4, true);
        assertEquals(30.0, t.getCenancestorHeight(), 1e-12);
        assertEquals(10.0, t.getBranchLength(), 1e-12);
        t.storeModelState();
        t.getTree().setNodeHeight(t.getTree().getRoot(), 25);
        assertEquals(33.0, t.getHeightParameter().getParameterValue(0), 1e-12);
        assertEquals(8.0, t.getBranchParameter().getParameterValue(0), 1e-12);
        t.getFractionParameter().setParameterValue(0, 0);
        assertEquals(25.0, t.getCenancestorHeight(), 0);
        t.getFractionParameter().setParameterValue(0, 1);
        assertEquals(0.0, t.getFounderDuration(), 0);
        t.restoreModelState();
        assertEquals(30.0, t.getCenancestorHeight(), 1e-12);
        assertEquals(Math.log(25), t.getLogLikelihood(), 1e-12);
        try {
            t.getHeightParameter().setParameterValueQuietly(0, 32);
            fail("Derived height must be read-only");
        } catch (UnsupportedOperationException expected) { }
    }

    public void testInvalidSupportAndInitialStates() throws Exception {
        for (double q : new double[]{-0.1, 1.1, Double.NaN}) {
            try { timing("(A:20,B:20);", 45, q, true); fail(); }
            catch (IllegalArgumentException expected) { }
        }
        for (double h : new double[]{0, 20, 19, Double.NaN}) {
            try { timing("(A:20,B:20);", h, 0.4, true); fail(); }
            catch (IllegalArgumentException expected) { }
        }
        NeoplasmTimingModel t = timing("(A:20,B:20);", 45, 0.4, true);
        for (double root : new double[]{45, 46, Double.NaN, -1}) {
            t.storeModelState();
            t.getTree().setNodeHeight(t.getTree().getRoot(), root);
            assertEquals(Double.NEGATIVE_INFINITY, t.getLogLikelihood(), 0);
            t.restoreModelState();
            checkGeometry(t);
        }
        t.storeModelState();
        t.getOriginParameter().setParameterValue(0, 50);
        assertFalse(t.isValid());
        t.restoreModelState();
        checkGeometry(t);
    }

    public void testJointScaleInverseAndHastings() throws Exception {
        NeoplasmTimingModel t = timing("((A:10,B:10):10,C:20);", 45, 0.4, true);
        NeoplasmTreeScaleOperator op = scale(t);
        double forward = op.propose(1.1);
        assertEquals(22.0, t.getRootHeight(), 1e-12);
        assertEquals(11.0, t.getBranchLength(), 1e-12);
        assertEquals(33.0, t.getCenancestorHeight(), 1e-12);
        assertEquals(11.0 / 23, t.getFractionParameter().getParameterValue(0), 1e-12);
        assertEquals(Math.log(1.1) + Math.log(25.0 / 23), forward, 1e-12);
        assertEquals(0.0, forward + op.propose(1 / 1.1), 1e-12);
        assertEquals(20.0, t.getRootHeight(), 1e-12);
        assertEquals(0.4, t.getFractionParameter().getParameterValue(0), 1e-12);
        try { op.propose(2); fail("Origin crossing accepted"); }
        catch (OperatorFailedException expected) { }
        assertEquals(20.0, t.getRootHeight(), 1e-12); // Early rejection did not mutate.
    }

    public void testDatedTipsAndEndpointScales() throws Exception {
        NeoplasmTimingModel t = timing("(A:10,B:20);", 45, 0.5, true);
        double[] tips = new double[t.getTree().getExternalNodeCount()];
        for (int i = 0; i < tips.length; i++) tips[i] = t.getTree().getNodeHeight(t.getTree().getExternalNode(i));
        try { scale(t).propose(0.25); fail("Parent below dated tip accepted"); }
        catch (OperatorFailedException expected) { }
        scale(t).propose(0.8);
        for (int i = 0; i < tips.length; i++) assertEquals(tips[i], t.getTree().getNodeHeight(t.getTree().getExternalNode(i)), 0);
        for (double q : new double[]{0, 1}) {
            NeoplasmTimingModel endpoint = timing("(A:20,B:20);", 45, q, true);
            double hr = scale(endpoint).propose(0.8);
            assertEquals(0, hr + scale(endpoint).propose(1.25), 1e-12);
            assertEquals(q, endpoint.getFractionParameter().getParameterValue(0), 1e-12);
        }
    }

    public void testRootChangingTopologyMoves() throws Exception {
        MathUtils.setSeed(92341);
        NeoplasmTimingModel t = timing("((A:10,B:10):10,(C:5,D:5):15);", 45, 0.4, true);
        TreeModel tree = t.getTree();
        SubtreeSlideOperator op = new SubtreeSlideOperator(tree, 1, 15, true, false, false, false, CoercionMode.COERCION_OFF);
        int replacements = 0;
        for (int i = 0; i < 500; i++) {
            t.storeModelState();
            NodeRef root = tree.getRoot();
            double oldHeight = t.getCenancestorHeight();
            try {
                double hr = op.doOperation();
                if (Double.isFinite(hr) && t.isValid()) {
                    if (root.getNumber() != tree.getRoot().getNumber()) replacements++;
                    checkGeometry(t);
                    assertEquals(tree.getNodeHeight(tree.getRoot()) + 0.4 * t.getGap(), t.getCenancestorHeight(), 1e-12);
                }
            } catch (OperatorFailedException expected) { }
            t.restoreModelState();
            assertEquals(root.getNumber(), tree.getRoot().getNumber());
            assertEquals(oldHeight, t.getCenancestorHeight(), 1e-12);
        }
        assertTrue("Exercise actual root replacements", replacements > 0);
    }

    public void testPriorOnlySamplerPreservesHeightDensity() throws Exception {
        // Uniform joint density in (R,C), 0<R<C<H, implies E[R]=H/3, E[C]=2H/3.
        // A uniform q target instead gives E[R]=H/2, E[C]=3H/4. Detect missing Jacobians.
        for (boolean jacobian : new boolean[]{true, false}) {
            MathUtils.setSeed(jacobian ? 9123 : 9133);
            NeoplasmTimingModel t = timing("(A:2,B:2);", 10, 0.5, jacobian);
            Parameter root = t.getTree().getRootHeightParameter();
            root.addBounds(new Parameter.DefaultBounds(10, 0, 1));
            RandomWalkOperator rootMove = new RandomWalkOperator(root, 3, RandomWalkOperator.BoundaryCondition.reflecting, 1, CoercionMode.COERCION_OFF);
            RandomWalkOperator local = new RandomWalkOperator(t.getFractionParameter(), 0.4, RandomWalkOperator.BoundaryCondition.reflecting, 1, CoercionMode.COERCION_OFF);
            NeoplasmTreeScaleOperator joint = scale(t);
            double sumR = 0, sumC = 0, sumQ = 0;
            int count = 0;
            for (int i = 0; i < 120000; i++) {
                double oldLog = t.getLogLikelihood();
                t.storeModelState();
                double hr;
                try {
                    hr = i % 3 == 0 ? rootMove.doOperation() : i % 3 == 1 ? local.doOperation() : joint.doOperation();
                } catch (OperatorFailedException rejected) { hr = Double.NEGATIVE_INFINITY; }
                if (Math.log(MathUtils.nextDouble()) < t.getLogLikelihood() - oldLog + hr) t.acceptModelState();
                else t.restoreModelState();
                if (i >= 20000) {
                    sumR += t.getRootHeight(); sumC += t.getCenancestorHeight();
                    sumQ += t.getFractionParameter().getParameterValue(0); count++;
                }
            }
            assertEquals(jacobian ? 10.0 / 3 : 5, sumR / count, 0.12);
            assertEquals(jacobian ? 20.0 / 3 : 7.5, sumC / count, 0.12);
            assertEquals(0.5, sumQ / count, 0.015);
        }
    }
}

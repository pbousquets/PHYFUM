package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.NeoplasmTimingModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;

/** Scale all internal heights and the LCA-MRCA stem, leaving origin and tips fixed. */
public class NeoplasmTreeScaleOperator extends AbstractCoercableOperator {
    private final NeoplasmTimingModel timing;
    private double scaleFactor;

    public NeoplasmTreeScaleOperator(NeoplasmTimingModel timing, double scaleFactor, double weight, CoercionMode mode) {
        super(mode);
        if (!Double.isFinite(scaleFactor) || scaleFactor <= 0 || scaleFactor >= 1
                || !Double.isFinite(weight) || weight <= 0) {
            throw new IllegalArgumentException("Require 0 < scaleFactor < 1 and positive finite weight");
        }
        this.timing = timing;
        this.scaleFactor = scaleFactor;
        setWeight(weight);
    }

    public double doOperation() throws OperatorFailedException {
        return propose(scaleFactor + MathUtils.nextDouble() * (1 / scaleFactor - scaleFactor));
    }

    /** Deterministic part of the proposal, exposed for inverse/Jacobian tests. */
    public double propose(double scale) throws OperatorFailedException {
        if (!timing.isValid() || !Double.isFinite(scale) || scale <= 0) throw new OperatorFailedException("Invalid timing/scale");
        TreeModel tree = timing.getTree();
        Parameter q = timing.getFractionParameter();
        double oldGap = timing.getGap();
        double newGap = timing.getOriginHeight() - scale * timing.getRootHeight();
        double newFraction = scale * timing.getBranchLength() / newGap;
        if (!(newGap > 0) || !Double.isFinite(newFraction) || newFraction < 0 || newFraction > 1
                || newFraction < q.getBounds().getLowerLimit(0) || newFraction > q.getBounds().getUpperLimit(0)) {
            throw new OperatorFailedException("Scaled LCA exceeds the fixed origin or fraction bounds");
        }
        // Validate final geometry before changing anything, including fixed dated tips.
        double[] heights = new double[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            heights[node.getNumber()] = tree.getNodeHeight(node) * (tree.isExternal(node) ? 1 : scale);
            if (!Double.isFinite(heights[node.getNumber()])) throw new OperatorFailedException("Nonfinite scaled height");
        }
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            NodeRef parent = tree.getParent(node);
            if (parent != null && heights[node.getNumber()] > heights[parent.getNumber()]) {
                throw new OperatorFailedException("Scaled parent is younger than its child");
            }
        }
        tree.beginTreeEdit();
        try {
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                tree.setNodeHeight(node, heights[node.getNumber()]);
            }
            q.setParameterValue(0, newFraction);
        } finally {
            tree.endTreeEdit();
        }
        // Account for any extra height bounds installed by the user's model.
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            double value = tree.getNodeHeight(node);
            if (value < tree.getNodeHeightLower(node) || value > tree.getNodeHeightUpper(node)) {
                throw new OperatorFailedException("Scaled tree violates height bounds");
            }
        }
        // State determinant s^(K+1)*D/D'; inversion of the uniform scale adds s^-2.
        return (tree.getInternalNodeCount() - 1) * Math.log(scale) + Math.log(oldGap / newGap);
    }

    public String getOperatorName() { return "neoplasmTreeScale(" + timing.getId() + ")"; }
    public double getCoercableParameter() { return Math.log(1 / scaleFactor - 1); }
    public void setCoercableParameter(double value) { scaleFactor = 1 / (Math.exp(value) + 1); }
    public double getRawParameter() { return scaleFactor; }
    public double getTargetAcceptanceProbability() { return 0.234; }
    public String getPerformanceSuggestion() { return ""; }
}

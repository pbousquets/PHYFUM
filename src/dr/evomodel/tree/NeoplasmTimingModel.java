package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.inference.model.*;

/** Shared origin/LCA/MRCA geometry and the height (or branch)-to-fraction Jacobian. */
public class NeoplasmTimingModel extends AbstractModelLikelihood {
    private final TreeModel tree;
    private final Parameter fraction, origin;
    private final double fixedOrigin;
    private final boolean timeCoordinatePrior;
    private final DerivedParameter height = new DerivedParameter("cenancestorHeight", 0);
    private final DerivedParameter branch = new DerivedParameter("cenancestorBranch", 1);
    private final DerivedParameter duration = new DerivedParameter("founderDuration", 2);

    public NeoplasmTimingModel(TreeModel tree, Parameter origin, Parameter fraction, boolean timeCoordinatePrior) {
        super("neoplasmTiming");
        this.tree = tree;
        this.origin = origin;
        this.fraction = fraction;
        this.timeCoordinatePrior = timeCoordinatePrior;
        if (origin.getDimension() != 1 || fraction.getDimension() != 1 || origin == fraction) {
            throw new IllegalArgumentException("Origin and LCA fraction must be distinct scalar parameters");
        }
        fixedOrigin = origin.getParameterValue(0);
        if (!isValid()) throw new IllegalArgumentException("Require finite heights, valid tree, 0 <= LCA fraction <= 1 and MRCA < origin");
        origin.addBounds(new Parameter.DefaultBounds(fixedOrigin, fixedOrigin, 1));
        fraction.addBounds(new Parameter.DefaultBounds(1, 0, 1));
        addModel(tree);
        addVariable(origin);
        addVariable(fraction);
    }

    public TreeModel getTree() { return tree; }
    public Parameter getOriginParameter() { return origin; }
    public Parameter getFractionParameter() { return fraction; }
    public Parameter getHeightParameter() { return height; }
    public Parameter getBranchParameter() { return branch; }
    public Parameter getDurationParameter() { return duration; }
    public double getOriginHeight() { return fixedOrigin; }
    public double getRootHeight() { return tree.getNodeHeight(tree.getRoot()); }
    public double getGap() { return fixedOrigin - getRootHeight(); }
    public double getBranchLength() { return fraction.getParameterValue(0) * getGap(); }
    public double getCenancestorHeight() { return getRootHeight() + getBranchLength(); }
    public double getFounderDuration() { return (1 - fraction.getParameterValue(0)) * getGap(); }

    /** Also catches root replacements and tree proposals that do not consult bounds. */
    public boolean isValid() {
        double q = fraction.getParameterValue(0);
        if (!Double.isFinite(fixedOrigin) || fixedOrigin <= 0 || origin.getParameterValue(0) != fixedOrigin
                || !Double.isFinite(q) || q < 0 || q > 1 || !(getRootHeight() < fixedOrigin)) return false;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            double h = tree.getNodeHeight(node);
            if (!Double.isFinite(h) || h < 0) return false;
            if (h < tree.getNodeHeightLower(node) || h > tree.getNodeHeightUpper(node)) return false;
            NodeRef parent = tree.getParent(node);
            if (parent != null && h > tree.getNodeHeight(parent)) return false;
        }
        return true;
    }

    /** Include this object ONCE in the prior, independently of the number of partitions. */
    public double getLogLikelihood() {
        if (!isValid()) return Double.NEGATIVE_INFINITY;
        return timeCoordinatePrior ? Math.log(getGap()) : 0;
    }
    public boolean evaluateEarly() { return true; }
    public Model getModel() { return this; }
    public void makeDirty() { changed(); }
    protected void handleModelChangedEvent(Model model, Object object, int index) { changed(); }
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) { changed(); }
    private void changed() {
        // Read-only outputs always read the current root, including after topology changes.
        height.fireParameterChangedEvent();
        branch.fireParameterChangedEvent();
        duration.fireParameterChangedEvent();
        fireModelChanged();
    }
    protected void storeState() { }
    protected void restoreState() { }
    protected void acceptState() { }

    private class DerivedParameter extends Parameter.Abstract {
        private final String name;
        private final int kind;
        private Bounds<Double> bounds = new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0, 1);
        DerivedParameter(String name, int kind) { this.name = name; this.kind = kind; }
        public int getDimension() { return 1; }
        public String getParameterName() { return getId() == null ? name : getId(); }
        public double getParameterValue(int dim) {
            if (dim != 0) throw new IndexOutOfBoundsException();
            return kind == 0 ? getCenancestorHeight() : kind == 1 ? getBranchLength() : getFounderDuration();
        }
        private UnsupportedOperationException readOnly() {
            return new UnsupportedOperationException("Derived neoplasm timing: operate on the LCA fraction or tree instead");
        }
        public void setParameterValue(int dim, double value) { throw readOnly(); }
        public void setParameterValueQuietly(int dim, double value) { throw readOnly(); }
        public void setParameterValueNotifyChangedAll(int dim, double value) { throw readOnly(); }
        public void addDimension(int index, double value) { throw readOnly(); }
        public double removeDimension(int index) { throw readOnly(); }
        protected void adoptValues(Parameter source) { throw readOnly(); }
        protected void storeValues() { }
        protected void restoreValues() { }
        protected void acceptValues() { }
        public Bounds<Double> getBounds() { return bounds; }
        public void addBounds(Bounds newBounds) {
            final Bounds<Double> previous = bounds;
            bounds = new Bounds<Double>() {
                public int getBoundsDimension() { return 1; }
                public Double getLowerLimit(int dim) { return Math.max(previous.getLowerLimit(dim), (Double) newBounds.getLowerLimit(dim)); }
                public Double getUpperLimit(int dim) { return Math.min(previous.getUpperLimit(dim), (Double) newBounds.getUpperLimit(dim)); }
            };
        }
    }
}

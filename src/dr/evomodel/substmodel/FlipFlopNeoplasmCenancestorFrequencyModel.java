package dr.evomodel.substmodel;

import dr.evolution.datatype.AFsequence;
import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.evomodel.tree.NeoplasmTimingModel;

/**
 * Evolves a single founder cell before error-free expansion to a stem-cell niche.
 * Initial frequencies are ordered unmethylated, hemimethylated, fully methylated.
 * Rates are absolute per unit of elapsed time; no replacement or flux normalization
 * is used during the single-cell phase.
 */
public class FlipFlopNeoplasmCenancestorFrequencyModel extends AbstractGeneralFrequencyModel {
    private final Parameter initialFrequencies;
    private final boolean useMethylatedProportion;
    private final Parameter elapsedTime;
    private final Parameter mu, gamma, nu, zeta;
    private final FounderSubstitutionModel founderModel;
    private final int stateCount;
    private final int stemCells;
    private NeoplasmTimingModel timing;

    /** Bind the shared geometry; duration must be its live read-only output. */
    public void setTimingModel(NeoplasmTimingModel timing) {
        if (this.timing != null || elapsedTime != timing.getDurationParameter()) {
            throw new IllegalArgumentException("Founder frequencies must use the shared timing duration");
        }
        this.timing = timing;
        addModel(timing);
    }
    public NeoplasmTimingModel getTimingModel() { return timing; }
    public boolean hasValidDuration() {
        double t = elapsedTime.getParameterValue(0);
        return Double.isFinite(t) && t >= 0 && (timing == null || timing.isValid());
    }

    public FlipFlopNeoplasmCenancestorFrequencyModel(DataType dataType, Parameter initialFrequencies,
            Parameter elapsedTime, Parameter mu, Parameter gamma, Parameter nu, Parameter zeta) {
        this(dataType, initialFrequencies, false, elapsedTime, mu, gamma, nu, zeta);
    }

    /** A live proportion p represents the origin distribution (1-p, 0, p). */
    public static FlipFlopNeoplasmCenancestorFrequencyModel fromMethylatedProportion(
            DataType dataType, Parameter proportion, Parameter elapsedTime,
            Parameter mu, Parameter gamma, Parameter nu, Parameter zeta) {
        return new FlipFlopNeoplasmCenancestorFrequencyModel(dataType, proportion, true,
                elapsedTime, mu, gamma, nu, zeta);
    }

    private FlipFlopNeoplasmCenancestorFrequencyModel(DataType dataType, Parameter initialFrequencies,
            boolean useMethylatedProportion, Parameter elapsedTime,
            Parameter mu, Parameter gamma, Parameter nu, Parameter zeta) {
        super(dataType);
        if (!(dataType instanceof AFsequence)) {
            throw new IllegalArgumentException("Neoplasm frequencies require an AFsequence data type");
        }
        stateCount = dataType.getStateCount();
        stemCells = (int) Math.round((Math.sqrt(8.0 * stateCount + 1) - 3) / 2);
        if (stemCells < 1 || (stemCells + 1L) * (stemCells + 2L) / 2 != stateCount) {
            throw new IllegalArgumentException("State count must equal (S+1)*(S+2)/2 for S >= 1");
        }
        this.initialFrequencies = initialFrequencies;
        this.useMethylatedProportion = useMethylatedProportion;
        this.elapsedTime = elapsedTime;
        this.mu = mu;
        this.gamma = gamma;
        this.nu = nu;
        this.zeta = zeta;
        // Validate inputs before constructing the shared eigendecomposition model.
        getOriginFrequencies();
        scalar(elapsedTime, "elapsedTime");
        validateRates();
        founderModel = new FounderSubstitutionModel(mu, gamma, nu, zeta);
        addModel(founderModel); // Includes the eigen cache in MCMC store/restore.
        addVariable(initialFrequencies);
        initialFrequencies.addBounds(new Parameter.DefaultBounds(1.0, 0.0,
                initialFrequencies.getDimension()));
        for (Parameter parameter : new Parameter[]{elapsedTime, mu, gamma, nu, zeta}) {
            addVariable(parameter); // Shared rates are registered only once.
            parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }
    }

    public int getFrequencyCount() {
        return stateCount;
    }

    public double getFrequency(int index) {
        return getFrequencies()[index];
    }

    public void setFrequency(int index, double value) {
        throw new UnsupportedOperationException("Cenancestor frequencies are derived; change the origin frequencies instead");
    }

    public double[] getFrequencies() {
        double[] initial = getOriginFrequencies();
        double time = scalar(elapsedTime, "elapsedTime");
        validateRates();
        double[] transition = new double[9];
        if (time == 0) {
            // At zero duration even a defective generator has the identity transition.
            transition[0] = transition[4] = transition[8] = 1.0;
        } else {
            founderModel.getTransitionProbabilities(time, transition);
        }
        double[] frequencies = new double[stateCount];
        // PHYFUM enumerates m first, then k: clonal states (0,0), (S,0), (0,S).
        int[] clonalStates = {0, stemCells, stateCount - 1};
        for (int destination = 0; destination < 3; destination++) {
            for (int source = 0; source < 3; source++) {
                frequencies[clonalStates[destination]] += initial[source] * transition[3 * source + destination];
            }
        }
        // Only the rate decomposition is cached. Read time and origin frequencies live:
        // branchRules can update cenancestorHeight quietly without a change event.
        return frequencies;
    }

    private void validateRates() {
        double methylation = scalar(mu, "mu");
        double demethylation = scalar(gamma, "gamma");
        double hemiMethylation = scalar(nu, "nu");
        double hemiDemethylation = scalar(zeta, "zeta");
        requireNonnegativeFinite(Math.max(Math.max(2 * methylation, 2 * demethylation),
                hemiMethylation + hemiDemethylation), "Total exit rate");
    }

    private double[] getOriginFrequencies() {
        if (useMethylatedProportion) {
            double proportion = scalar(initialFrequencies, "methylatedProportion");
            if (proportion > 1.0) {
                throw new IllegalArgumentException("methylatedProportion must be between 0 and 1");
            }
            return new double[]{1 - proportion, 0, proportion};
        }
        if (initialFrequencies.getDimension() != 3) {
            throw new IllegalArgumentException("Origin frequencies must have exactly three elements (U, H, M)");
        }
        double[] frequencies = initialFrequencies.getParameterValues();
        double sum = 0;
        for (double value : frequencies) {
            requireNonnegativeFinite(value, "Origin frequency");
            sum += value;
        }
        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Origin frequencies must sum to 1");
        }
        return frequencies;
    }

    private static double scalar(Parameter parameter, String name) {
        if (parameter.getDimension() != 1) {
            throw new IllegalArgumentException(name + " must be a scalar parameter");
        }
        double value = parameter.getParameterValue(0);
        requireNonnegativeFinite(value, name);
        return value;
    }

    private static void requireNonnegativeFinite(double value, String name) {
        if (!Double.isFinite(value) || value < 0) {
            throw new IllegalArgumentException(name + " must be finite and nonnegative");
        }
    }

    /**
     * Row-source generator in (U,H,M) order:
     * [-2mu, 2mu, 0; zeta, -(nu+zeta), nu; 0, 2gamma, -2gamma].
     * Supply only the single-cell generator; all matrix exponentiation, conditioning
     * checks and eigen-cache management are inherited from FlipFlopModel's backend.
     * There is no S-1 replacement term or stationary-frequency rate multiplier.
     */
    private static class FounderSubstitutionModel extends ComplexSubstitutionModel {
        private final Parameter mu, gamma, nu, zeta;

        FounderSubstitutionModel(Parameter mu, Parameter gamma, Parameter nu, Parameter zeta) {
            this(new AFsequence(3), mu, gamma, nu, zeta);
        }

        private FounderSubstitutionModel(DataType type, Parameter mu, Parameter gamma,
                Parameter nu, Parameter zeta) {
            // The backend requires a FrequencyModel, but these placeholders do not
            // enter Q or the origin distribution, and flux normalization is disabled.
            super("neoplasmFounder", type,
                    new FrequencyModel(type, new double[]{1.0 / 3, 1.0 / 3, 1.0 / 3}), null);
            this.mu = mu;
            this.gamma = gamma;
            this.nu = nu;
            this.zeta = zeta;
            for (Parameter rate : new Parameter[]{mu, gamma, nu, zeta}) {
                addVariable(rate); // Invalidates the inherited cache, deduplicating tied rates.
            }
            setNormalization(false);
            // Initialize the eigen storage even if the initial duration is zero,
            // so the first MCMC storeState never encounters uninitialized arrays.
            setupMatrix();
        }

        @Override
        public void storeIntoAmat() {
            // setupMatrix supplies the diagonals from these off-diagonal rates.
            amat[0][1] = 2 * mu.getParameterValue(0);
            amat[0][2] = 0;
            amat[1][0] = zeta.getParameterValue(0);
            amat[1][2] = nu.getParameterValue(0);
            amat[2][0] = 0;
            amat[2][1] = 2 * gamma.getParameterValue(0);
        }
    }
}

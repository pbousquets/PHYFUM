package dr.evomodel.substmodel;

import dr.evolution.datatype.AFsequence;
import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;

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
    private final int stateCount;
    private final int stemCells;

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
        // Validate all inputs before registering them with the model.
        getFrequencies();
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
        double methylation = scalar(mu, "mu");
        double demethylation = scalar(gamma, "gamma");
        double hemiMethylation = scalar(nu, "nu");
        double hemiDemethylation = scalar(zeta, "zeta");
        double[][] transition = transitionMatrix(time, methylation, demethylation,
                hemiMethylation, hemiDemethylation);
        double[] frequencies = new double[stateCount];
        // PHYFUM enumerates m first, then k: clonal states (0,0), (S,0), (0,S).
        int[] clonalStates = {0, stemCells, stateCount - 1};
        for (int destination = 0; destination < 3; destination++) {
            for (int source = 0; source < 3; source++) {
                frequencies[clonalStates[destination]] += initial[source] * transition[source][destination];
            }
        }
        // Deliberately uncached: branchRules can update cenancestorHeight quietly.
        // Reading live parameters also makes restoration require no additional model state.
        return frequencies;
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
     * Uniformization computes exp(Q*t) without a singular S-1 replacement term
     * or an eigenvector inverse, including zero rates and absorbing states.
     */
    private static double[][] transitionMatrix(double time, double mu, double gamma, double nu, double zeta) {
        double rate = Math.max(Math.max(2 * mu, 2 * gamma), nu + zeta);
        requireNonnegativeFinite(rate, "Total exit rate");
        double[][] transition = identity();
        if (time == 0 || rate == 0) {
            return transition;
        }
        double[][] jump = {
                {1 - 2 * mu / rate, 2 * mu / rate, 0},
                {zeta / rate, 1 - (nu + zeta) / rate, nu / rate},
                {0, 2 * gamma / rate, 1 - 2 * gamma / rate}
        };
        // Scale until rate*t <= 1, then square back. Logarithms avoid overflow
        // when both the finite rate and finite duration are large.
        int squarings = Math.max(0, (int) Math.ceil((Math.log(rate) + Math.log(time)) / Math.log(2)));
        double scaledTime = rate * Math.scalb(time, -squarings);
        double weight = Math.exp(-scaledTime);
        double[][] power = identity();
        for (int i = 0; i < 3; i++) {
            transition[i][i] = weight;
        }
        // exp(Q*t) = exp(-r*t) * sum_n (r*t)^n / n! * (I + Q/r)^n.
        for (int n = 1; n <= 32; n++) {
            power = multiply(power, jump);
            weight *= scaledTime / n;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    transition[i][j] += weight * power[i][j];
                }
            }
        }
        normalizeRows(transition);
        for (int n = 0; n < squarings; n++) {
            transition = multiply(transition, transition);
            normalizeRows(transition); // Prevent rounding drift under repeated squaring.
        }
        return transition;
    }

    private static double[][] identity() {
        return new double[][]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    }

    private static double[][] multiply(double[][] left, double[][] right) {
        double[][] product = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    product[i][j] += left[i][k] * right[k][j];
                }
            }
        }
        return product;
    }

    private static void normalizeRows(double[][] matrix) {
        for (double[] row : matrix) {
            double sum = row[0] + row[1] + row[2];
            for (int j = 0; j < 3; j++) {
                row[j] /= sum;
            }
        }
    }
}

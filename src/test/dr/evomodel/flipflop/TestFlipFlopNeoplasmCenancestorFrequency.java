package test.dr.evomodel.flipflop;

import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.FlipFlopNeoplasmCenancestorFrequencyModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/** Analytic checks of founder-cell evolution, clonal expansion, and live parameters. */
public class TestFlipFlopNeoplasmCenancestorFrequency extends TestCase {
    private static final double TOLERANCE = 1e-11;

    private FlipFlopNeoplasmCenancestorFrequencyModel model(int cells, double[] initial,
            double time, double mu, double gamma, double nu, double zeta) {
        return new FlipFlopNeoplasmCenancestorFrequencyModel(
                new AFsequence((cells + 1) * (cells + 2) / 2), new Parameter.Default(initial),
                new Parameter.Default(time), new Parameter.Default(mu), new Parameter.Default(gamma),
                new Parameter.Default(nu), new Parameter.Default(zeta));
    }

    private void assertDistribution(double[] expected, double[] actual) {
        assertEquals(expected.length, actual.length);
        double sum = 0;
        for (int i = 0; i < expected.length; i++) {
            assertEquals("State " + i, expected[i], actual[i], TOLERANCE);
            assertTrue(Double.isFinite(actual[i]) && actual[i] >= 0);
            sum += actual[i];
        }
        assertEquals(1.0, sum, TOLERANCE);
    }

    /** Independent assembly through the backend's generic off-diagonal rate vector. */
    private double[] backendFrequencies(double[] initial, double time,
            double mu, double gamma, double nu, double zeta) {
        AFsequence type = new AFsequence(3);
        // Upper triangle in row order, then lower triangle in column order.
        Parameter rates = new Parameter.Default(new double[]{2 * mu, 0, nu, zeta, 0, 2 * gamma});
        ComplexSubstitutionModel backend = new ComplexSubstitutionModel("reference", type,
                new FrequencyModel(type, new double[]{1.0 / 3, 1.0 / 3, 1.0 / 3}), rates) {
            @Override
            protected double[] getPi() {
                // Absolute rates, as in FlipFlopModel with useFrequencyModel=false.
                return new double[]{1, 1, 1};
            }
        };
        backend.setNormalization(false);
        double[] transition = new double[9];
        backend.getTransitionProbabilities(time, transition);
        double[] result = new double[3];
        for (int source = 0; source < 3; source++) {
            for (int destination = 0; destination < 3; destination++) {
                result[destination] += initial[source] * transition[3 * source + destination];
            }
        }
        return result;
    }

    private void assertZeroFrequencies(double[] frequencies) {
        // This is the backend's invalid-decomposition sentinel, not a distribution.
        for (double frequency : frequencies) {
            assertEquals(0.0, frequency, 0.0);
        }
    }

    public void testZeroTimeAndClonalExpansion() {
        for (int cells : new int[]{1, 2, 5, 10}) {
            double[] expected = new double[(cells + 1) * (cells + 2) / 2];
            expected[0] = 0.2;
            expected[cells] = 0.3;
            expected[expected.length - 1] = 0.5;
            assertDistribution(expected, model(cells, new double[]{0.2, 0.3, 0.5},
                    0, 0.3, 0.2, 0.7, 0.4).getFrequencies());
            assertDistribution(expected, model(cells, new double[]{0.2, 0.3, 0.5},
                    100, 0, 0, 0, 0).getFrequencies());
        }
    }

    public void testRateDirectionAtShortTimes() {
        double[][] generator = {{-0.6, 0.6, 0}, {0.4, -1.1, 0.7}, {0, 0.4, -0.4}};
        double time = 1e-7;
        for (int source = 0; source < 3; source++) {
            double[] initial = new double[3];
            initial[source] = 1;
            double[] result = model(1, initial, time, 0.3, 0.2, 0.7, 0.4).getFrequencies();
            for (int destination = 0; destination < 3; destination++) {
                assertEquals(generator[source][destination],
                        (result[destination] - initial[destination]) / time, 1e-7);
            }
        }
    }

    public void testTiedRatesMatchTwoIndependentAlleles() {
        double mu = 0.3;
        double gamma = 0.2;
        for (double time : new double[]{1e-6, 0.5, 20, 1e6}) {
            double decay = Math.exp(-(mu + gamma) * time);
            double fromU = mu / (mu + gamma) * (1 - decay);
            double fromM = mu / (mu + gamma) + gamma / (mu + gamma) * decay;
            for (int source = 0; source < 3; source++) {
                double a = source == 2 ? fromM : fromU;
                double b = source == 0 ? fromU : fromM;
                double[] expected = {(1-a)*(1-b), a*(1-b) + b*(1-a), a*b};
                double[] initial = new double[3];
                initial[source] = 1;
                assertDistribution(expected, model(1, initial, time, mu, gamma, mu, gamma).getFrequencies());
            }
        }
    }

    public void testDiagonalizableAbsorbingStates() {
        // Irreversible U -> H -> M with distinct exit rates 0.8 and 0.5.
        double time = 2;
        double survival = Math.exp(-0.8 * time);
        double hemi = 0.8 / (0.8 - 0.5) * (Math.exp(-0.5 * time) - survival);
        assertDistribution(new double[]{survival, hemi, 1 - survival - hemi},
                model(1, new double[]{1, 0, 0}, time, 0.4, 0, 0.5, 0).getFrequencies());
        assertDistribution(new double[]{0, 0, 1},
                model(1, new double[]{0, 0, 1}, time, 0.4, 0, 0.5, 0).getFrequencies());
    }

    public void testDefectiveGeneratorUsesBackendFailure() {
        // Equal irreversible exit rates form a defective matrix. Do not fall back
        // to a separate exponential: match the shared backend's rejection behavior.
        double[] initial = {1, 0, 0};
        assertZeroFrequencies(backendFrequencies(initial, 2, 0.4, 0, 0.8, 0));
        assertZeroFrequencies(model(1, initial, 2, 0.4, 0, 0.8, 0).getFrequencies());
        assertDistribution(initial, model(1, initial, 0, 0.4, 0, 0.8, 0).getFrequencies());
    }

    public void testAsymmetricStationaryLimitAndLargeTime() {
        // Detailed balance: pi proportional to (zeta*2gamma, 2mu*2gamma, 2mu*nu).
        double sum = 0.16 + 0.24 + 0.42;
        double[] expected = {0.16/sum, 0.24/sum, 0.42/sum};
        assertDistribution(expected, model(1, new double[]{1, 0, 0}, 1e6,
                0.3, 0.2, 0.7, 0.4).getFrequencies());
    }

    public void testExtremeRatesUseBackendFailure() {
        // The shared decomposition does not converge at these extreme rates.
        // This explicitly records the loss of the former uniformization guarantee.
        double[] initial = {1, 0, 0};
        assertZeroFrequencies(backendFrequencies(initial, 1e200,
                0.3e200, 0.2e200, 0.7e200, 0.4e200));
        assertZeroFrequencies(model(1, initial, 1e200,
                0.3e200, 0.2e200, 0.7e200, 0.4e200).getFrequencies());
    }

    public void testMatchesSharedBackend() {
        for (double[] rates : new double[][]{{0.3, 0.2, 0.7, 0.4}, {0.03, 0.02, 0.03, 0.02},
                {1.5, 0.05, 0.4, 0.6}, {0, 0, 0, 0}, {3e-13, 2e-13, 7e-13, 4e-13}}) {
            for (double time : new double[]{1e-5, 0.5, 20, 1e6}) {
                for (int source = 0; source < 3; source++) {
                    double[] initial = new double[3];
                    initial[source] = 1;
                    double[] expected = backendFrequencies(initial, time,
                            rates[0], rates[1], rates[2], rates[3]);
                    double[] actual = model(1, initial, time,
                            rates[0], rates[1], rates[2], rates[3]).getFrequencies();
                    // Compare entries rather than imposing row renormalization:
                    // tiny-rate eigenvalue cutoffs/clipping are inherited unchanged.
                    for (int state = 0; state < 3; state++) {
                        assertEquals(expected[state], actual[state], TOLERANCE);
                    }
                }
            }
        }
    }

    public void testIndependentRateUpdatesCacheAndAcceptRestore() {
        Parameter initial = new Parameter.Default(new double[]{0.2, 0.3, 0.5});
        Parameter time = new Parameter.Default(0.0);
        Parameter[] rates = {new Parameter.Default(0.3), new Parameter.Default(0.2),
                new Parameter.Default(0.7), new Parameter.Default(0.4)};
        FlipFlopNeoplasmCenancestorFrequencyModel model = new FlipFlopNeoplasmCenancestorFrequencyModel(
                new AFsequence(3), initial, time, rates[0], rates[1], rates[2], rates[3]);
        ComplexSubstitutionModel backend = (ComplexSubstitutionModel) model.getModel(0);
        double[][] eigenvectors = backend.getEigenVectors();
        // Store before the first nonzero duration: eigen state must already exist.
        model.storeModelState();
        time.setParameterValueQuietly(0, 2.0);
        assertDistribution(backendFrequencies(initial.getParameterValues(), 2, 0.3, 0.2, 0.7, 0.4),
                model.getFrequencies());
        assertSame(eigenvectors, backend.getEigenVectors());
        model.acceptModelState();
        for (int index = 0; index < rates.length; index++) {
            double[] before = model.getFrequencies();
            model.storeModelState();
            rates[index].setParameterValue(0, rates[index].getParameterValue(0) * 1.5);
            double[] expected = backendFrequencies(initial.getParameterValues(), 2,
                    rates[0].getParameterValue(0), rates[1].getParameterValue(0),
                    rates[2].getParameterValue(0), rates[3].getParameterValue(0));
            assertDistribution(expected, model.getFrequencies());
            assertNotSame(eigenvectors, backend.getEigenVectors());
            model.restoreModelState();
            assertDistribution(before, model.getFrequencies());
            // Accept the same proposal next, then test restoration from that state.
            model.storeModelState();
            rates[index].setParameterValue(0, rates[index].getParameterValue(0) * 1.5);
            assertDistribution(expected, model.getFrequencies());
            model.acceptModelState();
            eigenvectors = backend.getEigenVectors();
        }
    }

    public void testRestoreAfterFailedDecomposition() {
        Parameter nu = new Parameter.Default(0.5);
        FlipFlopNeoplasmCenancestorFrequencyModel model = new FlipFlopNeoplasmCenancestorFrequencyModel(
                new AFsequence(3), new Parameter.Default(new double[]{1, 0, 0}),
                new Parameter.Default(2.0), new Parameter.Default(0.4), new Parameter.Default(0.0),
                nu, new Parameter.Default(0.0));
        double[] original = model.getFrequencies();
        model.storeModelState();
        nu.setParameterValue(0, 0.8);
        assertZeroFrequencies(model.getFrequencies());
        model.restoreModelState();
        assertDistribution(original, model.getFrequencies());
    }

    public void testParameterChangesQuietTimeChangesAndRestore() {
        Parameter initial = new Parameter.Default(new double[]{0.5, 0, 0.5});
        Parameter time = new Parameter.Default(2.0);
        Parameter mu = new Parameter.Default(0.3);
        Parameter gamma = new Parameter.Default(0.2);
        FlipFlopNeoplasmCenancestorFrequencyModel model = new FlipFlopNeoplasmCenancestorFrequencyModel(
                new AFsequence(3), initial, time, mu, gamma, mu, gamma);
        double[] original = model.getFrequencies();
        model.storeModelState();
        mu.setParameterValue(0, 0.6);
        time.setParameterValueQuietly(0, 4.0);
        initial.setParameterValue(0, 0.2);
        initial.setParameterValue(2, 0.8);
        assertDistribution(model(1, new double[]{0.2, 0, 0.8}, 4, 0.6, 0.2, 0.6, 0.2).getFrequencies(),
                model.getFrequencies());
        model.restoreModelState();
        assertDistribution(original, model.getFrequencies());
    }

    public void testLiveMethylatedProportion() {
        Parameter proportion = new Parameter.Default(0.5);
        Parameter mu = new Parameter.Default(0.3);
        Parameter gamma = new Parameter.Default(0.2);
        FlipFlopNeoplasmCenancestorFrequencyModel model =
                FlipFlopNeoplasmCenancestorFrequencyModel.fromMethylatedProportion(new AFsequence(6),
                        proportion, new Parameter.Default(2.0), mu, gamma, mu, gamma);
        for (double p : new double[]{0, 0.25, 1}) {
            proportion.setParameterValue(0, p);
            assertDistribution(model(2, new double[]{1-p, 0, p}, 2, 0.3, 0.2, 0.3, 0.2).getFrequencies(),
                    model.getFrequencies());
        }
    }

    public void testInvalidInputsRejected() {
        double[][] invalidFrequencies = {{0.5, 0.5}, {0.2, 0.2, 0.2}, {-0.1, 0.5, 0.6},
                {Double.NaN, 0, 1}, {Double.POSITIVE_INFINITY, 0, 0}};
        for (double[] frequencies : invalidFrequencies) {
            try {
                model(1, frequencies, 2, 0.3, 0.2, 0.7, 0.4);
                fail("Invalid origin frequencies accepted");
            } catch (IllegalArgumentException expected) {
                // Invalid inputs should fail before a likelihood is evaluated.
            }
        }
        try {
            model(1, new double[]{1, 0, 0}, -1, 0.3, 0.2, 0.7, 0.4);
            fail("Negative elapsed time accepted");
        } catch (IllegalArgumentException expected) {
            // The origin must precede the cenancestor.
        }
    }
}

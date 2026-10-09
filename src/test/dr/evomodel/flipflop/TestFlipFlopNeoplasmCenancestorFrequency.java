package test.dr.evomodel.flipflop;

import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.FlipFlopNeoplasmCenancestorFrequencyModel;
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

    public void testAbsorbingStatesAndRepeatedEigenvalue() {
        // Irreversible U -> H -> M, both steps at rate 0.8: an Erlang waiting time.
        double time = 2;
        double survival = Math.exp(-0.8 * time);
        assertDistribution(new double[]{survival, 0.8*time*survival, 1-(1+0.8*time)*survival},
                model(1, new double[]{1, 0, 0}, time, 0.4, 0, 0.8, 0).getFrequencies());
        assertDistribution(new double[]{0, 0, 1},
                model(1, new double[]{0, 0, 1}, time, 0.4, 0, 0.8, 0).getFrequencies());
    }

    public void testAsymmetricStationaryLimitAndLargeTime() {
        // Detailed balance: pi proportional to (zeta*2gamma, 2mu*2gamma, 2mu*nu).
        double sum = 0.16 + 0.24 + 0.42;
        double[] expected = {0.16/sum, 0.24/sum, 0.42/sum};
        assertDistribution(expected, model(1, new double[]{1, 0, 0}, 1e6,
                0.3, 0.2, 0.7, 0.4).getFrequencies());
        // The product of finite rates and time can overflow; the result must not.
        assertDistribution(expected, model(1, new double[]{1, 0, 0}, 1e200,
                0.3e200, 0.2e200, 0.7e200, 0.4e200).getFrequencies());
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

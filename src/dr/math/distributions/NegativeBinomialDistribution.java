/*
 * NegativeBinomialDistribution.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.math.distributions;
import dr.math.GammaFunction;
import dr.math.UnivariateFunction;
import org.apache.commons.math.MathException;
import org.apache.commons.math.special.Beta;

/**
 * @author Trevor Bedford
 * @version $Id$
 */
public class NegativeBinomialDistribution implements Distribution {

    double mean;
    double stdev;

    public NegativeBinomialDistribution(double mean, double stdev) {
        this.mean = mean;
        this.stdev = stdev;
    }

    public double pdf(double x) {
        if (x < 0)  return 0;
        return Math.exp(logPdf(x));
    }

    public double logPdf(double x) {
        if (x < 0)  return Double.NEGATIVE_INFINITY;
        double r = -1 * (mean*mean) / (mean - stdev*stdev);
        double p = mean / (stdev*stdev);
        return Math.log(Math.pow(1-p,x)) + Math.log(Math.pow(p, r)) + GammaFunction.lnGamma(r+x) - GammaFunction.lnGamma(r) - GammaFunction.lnGamma(x+1);
    }

    public double cdf(double x) {
        double r = -1 * (mean*mean) / (mean - stdev*stdev);
        double p = mean / (stdev*stdev);
        try {
            return Beta.regularizedBeta(p, r, x+1);
        } catch (MathException e) {
            // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
            // exceptions is bad. Better to return NaN and let the calling code deal with it.
            return Double.NaN;
//                throw MathRuntimeException.createIllegalArgumentException(
//                "Couldn't calculate beta cdf for alpha = " + alpha + ", beta = " + beta + ": " +e.getMessage());
        }
    }

    public double quantile(double y) {
        // TB - I'm having trouble implementing this
        return Double.NaN;
    }

    public double mean() {
        return mean;
    }

    public double variance() {
        return stdev*stdev;
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException();
    }

    public static void main(String[] args) {
        System.out.println("Test negative binomial");
        System.out.println("Mean 5, sd 5, x 5, pdf 0.074487, logPdf -2.59713");
        NegativeBinomialDistribution dist = new NegativeBinomialDistribution(5, 5);
        System.out.println("pdf = " + dist.pdf(5));
        System.out.println("logPdf = " + dist.logPdf(5));
    }

}

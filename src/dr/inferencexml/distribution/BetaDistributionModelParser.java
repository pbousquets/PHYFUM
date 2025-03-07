/*
 * BetaDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.BetaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;

/**
 */
public class BetaDistributionModelParser extends DistributionModelParser {

    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";

    public String getParserName() {
        return BetaDistributionModel.BETA_DISTRIBUTION_MODEL;
    }

    ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
        return new BetaDistributionModel(parameters[0], parameters[1]);
    }

    public String[] getParameterNames() {
        return new String[]{ALPHA, BETA};
    }

    public String getParserDescription() {
        return "A model of a beta distribution.";
    }

    public boolean allowOffset() {
        return false;
    }

    public Class getReturnType() {
        return BetaDistributionModel.class;
    }
}

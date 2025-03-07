/*
 * StrictClockCenancestorBranchRatesParser.java
 *
 * Modified by Diego Mallo from StrictClockBranchRatesParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.xml.*;
import dr.evomodel.branchratemodel.StrictClockCenancestorBranchRates;

import java.util.logging.Logger;

/**
 */
public class StrictClockCenancestorBranchRatesParser extends AbstractXMLObjectParser {

    public static final String STRICT_CLOCK_CENANCESTOR_BRANCH_RATES = "strictClockCenancestorBranchRates";
    public static final String RATE = "rate";

    public String getParserName() {
        return STRICT_CLOCK_CENANCESTOR_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);
        Bounds<Double> rateBounds=rateParameter.getBounds();

        if (rateBounds==null|| rateBounds!=null && rateBounds.getLowerLimit(0)<0.00000001){
            Logger.getLogger("dr.evomodel").warning("WARNING:\nIf used with ascertainment bias correction, a lower-unbounded clockrate may generate numerical instability resulting in positive data likelihoods and impeding proper mixing.");
        }

        Logger.getLogger("dr.evomodel").info("Using strict molecular clock model.");

        return new StrictClockCenancestorBranchRates(rateParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a strict clock model. " +
                        "All branches have the same rate of molecular evolution.";
    }

    public Class getReturnType() {
        return StrictClockCenancestorBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(RATE, Parameter.class, "The molecular evolutionary rate parameter", false),
    };
}


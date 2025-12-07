/*
 * CenancestorTreeLikelihoodParser.java
 *
 * Modified from TreeLikelihoodParser.java by DM
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

package dr.evomodelxml.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractGeneralFrequencyModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.Parameter;
import dr.xml.*;
import dr.evomodel.branchratemodel.CenancestorBranchRateModel;
import dr.evomodel.treelikelihood.CenancestorTreeLikelihood;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 */
public class CenancestorTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TREE_LIKELIHOOD = "cenancestorTreeLikelihood";
    public static final String ANCESTRAL_TREE_LIKELIHOOD = "cenancestorAncestralTreeLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String ALLOW_MISSING_TAXA = "allowMissingTaxa";
    public static final String STORE_PARTIALS = "storePartials";
    public static final String SCALING_FACTOR = "scalingFactor";
    public static final String SCALING_THRESHOLD = "scalingThreshold";
    public static final String FORCE_JAVA_CORE = "forceJavaCore";
    public static final String FORCE_RESCALING = "forceRescaling";
    public static final String CENANCESTOR_HEIGHT = "cenancestorHeight";
    public static final String CENANCESTOR_BRANCH = "cenancestorBranch";
    public static final String HEIGHT_RULES = "heightRules";
    public static final String DIVISION_MODEL = "divisionModel";
    //Division models
    public static final String CLONING = "cloning";
    public static final String BUDDING = "budding";

    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
        boolean allowMissingTaxa = xo.getAttribute(ALLOW_MISSING_TAXA, false);
        boolean storePartials = xo.getAttribute(STORE_PARTIALS, true);
        boolean forceJavaCore = xo.getAttribute(FORCE_JAVA_CORE, false);
        boolean heightRules = xo.getAttribute(HEIGHT_RULES,false);

        if (Boolean.valueOf(System.getProperty("java.only"))) {
            forceJavaCore = true;
        }

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
        AbstractGeneralFrequencyModel cenancestorFrequencyModel = (AbstractGeneralFrequencyModel) xo.getChild(AbstractGeneralFrequencyModel.class);

        CenancestorBranchRateModel branchRateModel = (CenancestorBranchRateModel) xo.getChild(CenancestorBranchRateModel.class);

        Parameter cenancestor = (Parameter) xo.getElementFirstChild(CENANCESTOR_HEIGHT);

        Parameter cenancestorBranch = (Parameter) xo.getElementFirstChild(CENANCESTOR_BRANCH);

        //Parameter AsStatistic=(Parameter) xo.getChild(USE_AS_STATISTIC);

        TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);
        if (tipStatesModel != null && tipStatesModel.getPatternList() != null) {
            throw new XMLParseException("The same sequence error model cannot be used for multiple partitions");
        }
        if (tipStatesModel != null && tipStatesModel.getModelType() == TipStatesModel.Type.STATES) {
            throw new XMLParseException("The state emitting TipStateModel requires BEAGLE");
        }

        boolean forceRescaling = xo.getAttribute(FORCE_RESCALING, false);

        String division_model = "default";
        if (xo.hasAttribute(DIVISION_MODEL)) {
            division_model = (String) xo.getAttribute(DIVISION_MODEL);
            Set<String> allowed = new HashSet<>(Arrays.asList(CLONING, BUDDING));

            if (!allowed.contains(division_model)) {
                Logger.getLogger("dr.evolution").info("\nWARNING: cenancestorTreeLikelihood division model " +
                        division_model +
                        " was not recognized, using the default model");
                division_model = "default";
            }
        }

        return new CenancestorTreeLikelihood(
                patternList,
                treeModel,
                siteModel,
                branchRateModel,
                tipStatesModel,
                cenancestor,
                cenancestorBranch,
                cenancestorFrequencyModel,
                division_model,
                //AsStatistic,
                useAmbiguities, allowMissingTaxa, storePartials, forceJavaCore, forceRescaling,
                heightRules);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree with an extra branch to the cenancestor given the site model.";
    }

    public Class getReturnType() {
        return CenancestorTreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            AttributeRule.newBooleanRule(ALLOW_MISSING_TAXA, true),
            AttributeRule.newBooleanRule(STORE_PARTIALS, true),
            AttributeRule.newBooleanRule(FORCE_JAVA_CORE, true),
            AttributeRule.newBooleanRule(FORCE_RESCALING, true),
            AttributeRule.newBooleanRule(HEIGHT_RULES, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(SiteModel.class),
            new ElementRule(CenancestorBranchRateModel.class, true),
            new ElementRule(TipStatesModel.class, true),
            new ElementRule(FrequencyModel.class, true),
            new ElementRule(CENANCESTOR_HEIGHT,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(CENANCESTOR_BRANCH,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new StringAttributeRule(DIVISION_MODEL, "Crypt division model", new String[]{CLONING,BUDDING}, true),
            //new ElementRule(USE_AS_STATISTIC,
            //new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
    };
}


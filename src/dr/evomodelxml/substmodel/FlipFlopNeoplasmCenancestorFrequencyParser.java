package dr.evomodelxml.substmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Date;
import dr.evomodel.substmodel.FlipFlopNeoplasmCenancestorFrequencyModel;
import dr.inference.model.DifferenceParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Arrays;

/** XML interface for single-cell evolution before clonal niche expansion. */
public class FlipFlopNeoplasmCenancestorFrequencyParser extends AbstractXMLObjectParser {
    public static final String FREQUENCY_MODEL = "flipflopNeoplasmCenancestorFrequency";
    public static final String FREQUENCIES = "frequencies";
    public static final String METHYLATED_PROPORTION = "methylatedProportion";
    public static final String ELAPSED_TIME = "elapsedTime";
    public static final String ORIGIN_HEIGHT = "originHeight";
    public static final String CENANCESTOR_HEIGHT = "cenancestorHeight";

    public String getParserName() {
        return FREQUENCY_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        PatternList patterns = (PatternList) xo.getChild(PatternList.class);
        DataType dataType = patterns.getDataType();
        Parameter mu = parameter(xo, "mu");
        Parameter gamma = parameter(xo, "gamma");
        Parameter nu = xo.hasChildNamed("nu") ? parameter(xo, "nu") : mu;
        Parameter zeta = xo.hasChildNamed("zeta") ? parameter(xo, "zeta") : gamma;

        boolean hasElapsed = xo.hasChildNamed(ELAPSED_TIME);
        boolean hasOrigin = xo.hasChildNamed(ORIGIN_HEIGHT);
        boolean hasCenancestor = xo.hasChildNamed(CENANCESTOR_HEIGHT);
        if (hasElapsed ? (hasOrigin || hasCenancestor) : !hasCenancestor) {
            throw new XMLParseException("Specify elapsedTime OR cenancestorHeight with an optional originHeight");
        }
        Parameter elapsedTime;
        if (hasElapsed) {
            elapsedTime = parameter(xo, ELAPSED_TIME);
        } else {
            Parameter origin = hasOrigin ? parameter(xo, ORIGIN_HEIGHT) : inferOriginHeight(patterns);
            Parameter cenancestor = parameter(xo, CENANCESTOR_HEIGHT);
            if (origin.getDimension() != 1 || cenancestor.getDimension() != 1
                    || !Double.isFinite(origin.getParameterValue(0))
                    || !Double.isFinite(cenancestor.getParameterValue(0))
                    || origin.getParameterValue(0) < 0 || cenancestor.getParameterValue(0) < 0) {
                throw new XMLParseException("originHeight and cenancestorHeight must be finite nonnegative scalars");
            }
            // BEAST heights increase backwards in time. For contemporaneous samples,
            // originHeight is their age measured from origin time zero.
            elapsedTime = new DifferenceParameter(Arrays.asList(origin, cenancestor));
        }

        boolean hasFrequencies = xo.hasChildNamed(FREQUENCIES);
        boolean hasProportion = xo.hasChildNamed(METHYLATED_PROPORTION);
        boolean hasProportionAttribute = xo.hasAttribute(METHYLATED_PROPORTION);
        if (hasProportion && (hasFrequencies || hasProportionAttribute)) {
            throw new XMLParseException("A methylatedProportion parameter cannot be combined with origin frequencies or the proportion attribute");
        }
        if (!hasFrequencies && !hasProportion && !hasProportionAttribute) {
            throw new XMLParseException("Specify three origin frequencies (U,H,M) or methylatedProportion");
        }
        try {
            if (hasProportion) {
                return FlipFlopNeoplasmCenancestorFrequencyModel.fromMethylatedProportion(dataType,
                        parameter(xo, METHYLATED_PROPORTION), elapsedTime, mu, gamma, nu, zeta);
            }
            Parameter frequencies = hasFrequencies ? parameter(xo, FREQUENCIES) : new Parameter.Default(3);
            if (hasProportionAttribute) {
                double p = xo.getDoubleAttribute(METHYLATED_PROPORTION);
                if (!Double.isFinite(p) || p < 0 || p > 1) {
                    throw new IllegalArgumentException("methylatedProportion must be between 0 and 1");
                }
                // Match the existing parser's attribute scheme: initialize the supplied
                // parameter, but keep just three origin states instead of niche states.
                frequencies.setDimension(3);
                frequencies.setParameterValueQuietly(0, 1 - p);
                frequencies.setParameterValueQuietly(1, 0);
                frequencies.setParameterValue(2, p);
            }
            if (xo.getAttribute("normalize", false)) {
                double sum = 0;
                for (double value : frequencies.getParameterValues()) {
                    if (!Double.isFinite(value) || value < 0) {
                        throw new IllegalArgumentException("Origin frequencies must be finite and nonnegative");
                    }
                    sum += value;
                }
                if (!Double.isFinite(sum) || sum <= 0) {
                    throw new IllegalArgumentException("Origin frequencies must have a positive finite sum");
                }
                for (int i = 0; i < frequencies.getDimension(); i++) {
                    frequencies.setParameterValue(i, frequencies.getParameterValue(i) / sum);
                }
            }
            return new FlipFlopNeoplasmCenancestorFrequencyModel(dataType, frequencies,
                    elapsedTime, mu, gamma, nu, zeta);
        } catch (IllegalArgumentException error) {
            throw new XMLParseException(error.getMessage());
        }
    }

    private static Parameter parameter(XMLObject xo, String name) throws XMLParseException {
        return (Parameter) xo.getElementFirstChild(name);
    }

    /** Infer once from this alignment's ages, avoiding BEAST's global taxon date cache. */
    private static Parameter inferOriginHeight(PatternList patterns) throws XMLParseException {
        final String guidance = "; supply an explicit originHeight instead";
        final int count;
        try {
            count = patterns.getTaxonCount();
        } catch (RuntimeException error) {
            throw new XMLParseException("Cannot infer originHeight without alignment taxa" + guidance);
        }
        if (count == 0) {
            throw new XMLParseException("Cannot infer originHeight from an empty taxon list" + guidance);
        }
        Date reference = null;
        double maximumAge = 0;
        for (int i = 0; i < count; i++) {
            Date date = patterns.getTaxon(i).getDate();
            if (date == null || date.isBackwards() || !Double.isFinite(date.getTimeValue())
                    || date.getTimeValue() < 0 || !Double.isFinite(date.getOrigin())) {
                throw new XMLParseException("Cannot infer originHeight: sample " + patterns.getTaxonId(i)
                        + " must have a finite nonnegative forwards age" + guidance);
            }
            if (reference != null && (date.getUnits() != reference.getUnits()
                    || date.getOrigin() != reference.getOrigin())) {
                throw new XMLParseException("Cannot infer originHeight: sample ages must use the same units and origin"
                        + guidance);
            }
            reference = date;
            maximumAge = Math.max(maximumAge, date.getTimeValue());
        }
        // The dates are fixed input metadata; the cenancestor height remains a live parameter.
        return new Parameter.Default(maximumAge);
    }

    public String getParserDescription() {
        return "Evolve origin frequencies (unmethylated, hemimethylated, methylated) through a single-cell "
                + "phase, then expand clonally to the alignment's niche size. Use elapsedTime or the difference "
                + "originHeight - cenancestorHeight in BEAST height units. Omitted originHeight defaults to the largest "
                + "alignment sample age (forwards dates with common units and origin). Rates are absolute; omitted nu/zeta share mu/gamma.";
    }

    public Class getReturnType() {
        return FlipFlopNeoplasmCenancestorFrequencyModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PatternList.class),
            new ElementRule(FREQUENCIES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(METHYLATED_PROPORTION, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            AttributeRule.newDoubleRule(METHYLATED_PROPORTION, true),
            AttributeRule.newBooleanRule("normalize", true),
            new XORRule(
                    new ElementRule(ELAPSED_TIME, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new AndRule(
                            new ElementRule(ORIGIN_HEIGHT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
                            new ElementRule(CENANCESTOR_HEIGHT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}))),
            new ElementRule("mu", new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule("gamma", new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule("nu", new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule("zeta", new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
    };
}

package dr.evomodelxml.operators;

import dr.evomodel.operators.NeoplasmTreeScaleOperator;
import dr.evomodel.tree.NeoplasmTimingModel;
import dr.inference.operators.CoercionMode;
import dr.xml.*;

public class NeoplasmTreeScaleOperatorParser extends AbstractXMLObjectParser {
    public String getParserName() { return "neoplasmTreeScale"; }
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        try {
            return new NeoplasmTreeScaleOperator((NeoplasmTimingModel) xo.getChild(NeoplasmTimingModel.class),
                    xo.getDoubleAttribute("scaleFactor"), xo.getDoubleAttribute("weight"), CoercionMode.parseMode(xo));
        } catch (IllegalArgumentException e) { throw new XMLParseException(e.getMessage()); }
    }
    public String getParserDescription() { return "Joint scale of internal heights and LCA-MRCA stem at fixed origin and tip dates."; }
    public Class getReturnType() { return NeoplasmTreeScaleOperator.class; }
    public XMLSyntaxRule[] getSyntaxRules() { return new XMLSyntaxRule[]{
        new ElementRule(NeoplasmTimingModel.class), AttributeRule.newDoubleRule("scaleFactor"),
        AttributeRule.newDoubleRule("weight"), AttributeRule.newBooleanRule("autoOptimize", true)
    }; }
}

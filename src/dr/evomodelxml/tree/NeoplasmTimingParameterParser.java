package dr.evomodelxml.tree;

import dr.evomodel.tree.NeoplasmTimingModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/** Read-only parameters for existing prior and logger interfaces. */
public class NeoplasmTimingParameterParser extends AbstractXMLObjectParser {
    public String getParserName() { return "neoplasmTimingParameter"; }
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        NeoplasmTimingModel timing = (NeoplasmTimingModel) xo.getChild(NeoplasmTimingModel.class);
        String quantity = xo.getStringAttribute("quantity");
        if (quantity.equals("height")) return timing.getHeightParameter();
        if (quantity.equals("branch")) return timing.getBranchParameter();
        if (quantity.equals("duration")) return timing.getDurationParameter();
        throw new XMLParseException("quantity must be height, branch or duration");
    }
    public String getParserDescription() { return "Read-only derived neoplasm timing parameter for priors/logging, not operators."; }
    public Class getReturnType() { return Parameter.class; }
    public XMLSyntaxRule[] getSyntaxRules() { return new XMLSyntaxRule[]{
        new ElementRule(NeoplasmTimingModel.class), AttributeRule.newStringRule("quantity")
    }; }
}

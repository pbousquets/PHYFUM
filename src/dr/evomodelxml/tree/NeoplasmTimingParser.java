package dr.evomodelxml.tree;

import dr.evolution.alignment.PatternList;
import dr.evomodel.tree.NeoplasmTimingModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.substmodel.FlipFlopNeoplasmCenancestorFrequencyParser;
import dr.inference.model.Parameter;
import dr.xml.*;

/** A single source of fixed origin and derived LCA times for all partitions. */
public class NeoplasmTimingParser extends AbstractXMLObjectParser {
    public String getParserName() { return "neoplasmTiming"; }
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter origin;
        if (xo.hasChildNamed("originHeight")) {
            origin = (Parameter) xo.getElementFirstChild("originHeight");
        } else {
            PatternList patterns = (PatternList) xo.getChild(PatternList.class);
            if (patterns == null) throw new XMLParseException("neoplasmTiming requires originHeight or an alignment with sample ages");
            origin = FlipFlopNeoplasmCenancestorFrequencyParser.inferOriginHeight(patterns);
        }
        String coordinate = xo.getAttribute("priorCoordinate", "height");
        if (!coordinate.equals("height") && !coordinate.equals("branch") && !coordinate.equals("fraction")) {
            throw new XMLParseException("priorCoordinate must be height, branch or fraction");
        }
        try {
            return new NeoplasmTimingModel((TreeModel) xo.getChild(TreeModel.class), origin,
                    (Parameter) xo.getElementFirstChild("lcaFraction"), !coordinate.equals("fraction"));
        } catch (IllegalArgumentException e) { throw new XMLParseException(e.getMessage()); }
    }
    public String getParserDescription() {
        return "LCA fraction (C-R)/(H-R). Include once in the prior: height (default) or branch coordinates add log(H-R); fraction coordinates do not.";
    }
    public Class getReturnType() { return NeoplasmTimingModel.class; }
    public XMLSyntaxRule[] getSyntaxRules() { return new XMLSyntaxRule[]{
        new ElementRule(TreeModel.class), new ElementRule(PatternList.class, true),
        new ElementRule("originHeight", new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
        new ElementRule("lcaFraction", new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
        AttributeRule.newStringRule("priorCoordinate", true)
    }; }
}

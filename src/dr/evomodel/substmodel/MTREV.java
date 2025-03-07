/*
 * MTREV.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.AminoAcids;
import dr.util.Author;
import dr.util.Citation;

import java.util.*;

/**
 * MTREV24 model of amino acid evolution
 * (complete sequence data of mtDNA from 24 vertebrate species)
 * Adachi, J., and Hasegawa, M. 1996. J. Mol. Evol. 42:459-468.
 *
 * @version $Id: MTREV.java,v 1.4 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class MTREV extends EmpiricalRateMatrix.AbstractAminoAcid {
	
	public static final MTREV INSTANCE = new MTREV();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private MTREV() { super("mtREV24");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];
		
		// Q matrix
		rate[0][1]=1.2199217606346e+01; 	rate[0][2]=1.4182139942122e+01;
		rate[0][3]=9.2985091873208e+00; 	rate[0][4]=3.1542792981957e+01;
		rate[0][5]=1.0025852846688e+00; 	rate[0][6]=5.1418866803338e+00;
		rate[0][7]=6.3531246495131e+01; 	rate[0][8]=7.3137132861715e+00;
		rate[0][9]=5.0782382656186e+01; 	rate[0][10]=1.3399741808481e+01;
		rate[0][11]=4.4021672780560e+00; 	rate[0][12]=7.4673480520104e+01;
		rate[0][13]=3.3513021631978e+00; 	rate[0][14]=2.8582502221773e+01;
		rate[0][15]=2.0413623195312e+02; 	rate[0][16]=2.5301305153906e+02;
		rate[0][17]=1.0000000000000e+00; 	rate[0][18]=3.4084158197615e+00;
		rate[0][19]=1.0266468401249e+02; 
	
		rate[1][2]=6.9661274444534e+00; 	rate[1][3]=1.0000000000000e+00;
		rate[1][4]=5.4384584796568e+01; 	rate[1][5]=1.1631134513343e+02;
		rate[1][6]=1.0000000000000e+00; 	rate[1][7]=1.2122831341194e+01;
		rate[1][8]=8.6961067087353e+01; 	rate[1][9]=1.0000000000000e+00;
		rate[1][10]=8.1976829394538e+00; 	rate[1][11]=7.4423215395318e+01;
		rate[1][12]=1.0000000000000e+00; 	rate[1][13]=2.4659158338099e+00;
		rate[1][14]=1.2439947713615e+01; 	rate[1][15]=3.1791814866372e+00;
		rate[1][16]=1.0935327216119e+00; 	rate[1][17]=1.1550775790126e+01;
		rate[1][18]=1.0000000000000e+00; 	rate[1][19]=4.0211417480338e+00;

		rate[2][3]=4.1809325468160e+02; 	rate[2][4]=3.1020979842967e+01;
		rate[2][5]=9.1349622725361e+01; 	rate[2][6]=3.3185663516310e+01;
		rate[2][7]=2.8052324651124e+01; 	rate[2][8]=2.6112087577885e+02;
		rate[2][9]=1.4261453863336e+01; 	rate[2][10]=7.9775653461977e+00;
		rate[2][11]=3.2036829276162e+02; 	rate[2][12]=3.4424354918739e+01;
		rate[2][13]=7.9996445145608e+00; 	rate[2][14]=3.8586541461044e+01;
		rate[2][15]=2.6020426225852e+02; 	rate[2][16]=1.2550758780474e+02;
		rate[2][17]=5.6207759736659e+00; 	rate[2][18]=1.0071406219571e+02;
		rate[2][19]=1.0000000000000e+00; 
	
		rate[3][4]=1.0000000000000e+00; 	rate[3][5]=2.9097352675564e+01;
		rate[3][6]=3.0713149855302e+02; 	rate[3][7]=2.9877072751897e+01;
		rate[3][8]=5.9995408885817e+01; 	rate[3][9]=2.2827096245105e+00;
		rate[3][10]=1.0000000000000e+00; 	rate[3][11]=1.2183938185384e+00;
		rate[3][12]=1.0000000000000e+00; 	rate[3][13]=2.6221929413096e+00;
		rate[3][14]=7.0708004204733e+00; 	rate[3][15]=3.6327934317139e+01;
		rate[3][16]=1.4743408713748e+01; 	rate[3][17]=1.0453246057102e+01;
		rate[3][18]=1.1165627147496e+01; 	rate[3][19]=1.0000000000000e+00;

		rate[4][5]=3.9599394038972e+01; 	rate[4][6]=1.0000000000000e+00;
		rate[4][7]=1.6163581056674e+01; 	rate[4][8]=7.4467985406234e+01;
		rate[4][9]=3.3018175376623e+01; 	rate[4][10]=1.3500725995091e+01;
		rate[4][11]=1.0000000000000e+00; 	rate[4][12]=3.2504095376923e+00;
		rate[4][13]=3.7264767083096e+01; 	rate[4][14]=1.6454136037822e+01;
		rate[4][15]=1.4581783243113e+02; 	rate[4][16]=9.4720031458442e+01;
		rate[4][17]=1.7684087896962e+01; 	rate[4][18]=1.3409157685926e+02;
		rate[4][19]=1.0000000000000e+00; 
	
		rate[5][6]=1.6503249008836e+02; 	rate[5][7]=3.5530760735494e+00;
		rate[5][8]=3.0652523140859e+02; 	rate[5][9]=4.3905393139325e+00;
		rate[5][10]=2.0895470525345e+01; 	rate[5][11]=2.4504076430724e+02;
		rate[5][12]=2.4931300477797e+01; 	rate[5][13]=1.0059428264289e+01;
		rate[5][14]=7.2256314165467e+01; 	rate[5][15]=2.8480937892158e+01;
		rate[5][16]=4.9962974409828e+01; 	rate[5][17]=1.0000000000000e+00;
		rate[5][18]=2.0430790980529e+01; 	rate[5][19]=9.9986289000676e+00;

		rate[6][7]=1.4884496769963e+01; 	rate[6][8]=2.5853576435567e+01;
		rate[6][9]=1.7418201388328e+00; 	rate[6][10]=1.0000000000000e+00;
		rate[6][11]=1.6519126809071e+02; 	rate[6][12]=1.0000000000000e+00;
		rate[6][13]=1.4067850525292e+00; 	rate[6][14]=6.7547121641947e+00;
		rate[6][15]=2.8794794140840e+01; 	rate[6][16]=7.8001372062558e+00;
		rate[6][17]=1.0000000000000e+00; 	rate[6][18]=6.9067239183061e+00;
		rate[6][19]=1.1127702362585e+01; 
	
		rate[7][8]=1.0000000000000e+00; 	rate[7][9]=3.1466649021550e+00;
		rate[7][10]=1.2699794194865e+00; 	rate[7][11]=1.1962111069278e+01;
		rate[7][12]=1.0000000000000e+00; 	rate[7][13]=1.0000000000000e+00;
		rate[7][14]=1.0000000000000e+00; 	rate[7][15]=6.6277950574411e+01;
		rate[7][16]=5.8800079133028e+00; 	rate[7][17]=5.7494182626674e+00;
		rate[7][18]=1.6887657206208e+00; 	rate[7][19]=1.3320553471351e+00;

		rate[8][9]=6.4536986087271e+00; 	rate[8][10]=6.0472584534958e+00;
		rate[8][11]=6.7197196398961e+01; 	rate[8][12]=6.2977633277779e+00;
		rate[8][13]=2.5347805183364e+01; 	rate[8][14]=3.2089868698728e+01;
		rate[8][15]=4.0766987134407e+01; 	rate[8][16]=2.3570850628539e+01;
		rate[8][17]=3.7286635325194e+00; 	rate[8][18]=3.5270764890474e+02;
		rate[8][19]=1.0000000000000e+00; 
	
		rate[9][10]=1.7320653206333e+02; 	rate[9][11]=1.0298655619743e+01;
		rate[9][12]=2.7262244199514e+02; 	rate[9][13]=4.4561065036310e+01;
		rate[9][14]=1.0856482766156e+01; 	rate[9][15]=2.5107659603898e+01;
		rate[9][16]=1.9391167162525e+02; 	rate[9][17]=1.0000000000000e+00;
		rate[9][18]=1.3161329199391e+01; 	rate[9][19]=6.4365086389428e+02;

		rate[10][11]=7.8314019154706e+00; 	rate[10][12]=2.8290920517725e+02;
		rate[10][13]=1.1371735519833e+02; 	rate[10][14]=2.1105885757279e+01;
		rate[10][15]=3.8741359395934e+01; 	rate[10][16]=6.6524559321657e+01;
		rate[10][17]=1.7071378554833e+01; 	rate[10][18]=2.3234516108847e+01;
		rate[10][19]=4.8247261078055e+01; 
	
		rate[11][12]=4.8092094826036e+01; 	rate[11][13]=3.3887559483420e+00;
		rate[11][14]=2.6368577564199e+01; 	rate[11][15]=5.5679895711418e+01;
		rate[11][16]=7.1750284708933e+01; 	rate[11][17]=1.2631893872825e+01;
		rate[11][18]=2.6932728996777e+01; 	rate[11][19]=1.0000000000000e+00;

		rate[12][13]=4.7798798034572e+01; 	rate[12][14]=9.9165053447429e+00;
		rate[12][15]=5.8505442466161e+01; 	rate[12][16]=2.7798190504760e+02;
		rate[12][17]=1.1427000119701e+01; 	rate[12][18]=2.1029990530586e+01;
		rate[12][19]=2.0397078683768e+02; 
	
		rate[13][14]=9.1089574817139e+00; 	rate[13][15]=3.3835737720574e+01;
		rate[13][16]=1.7815549567056e+01; 	rate[13][17]=4.1272404968214e+00;
		rate[13][18]=2.4504156395152e+02; 	rate[13][19]=3.3435675442163e+00;

		rate[14][15]=8.9421193040709e+01; 	rate[14][16]=6.7485067008375e+01;
		rate[14][17]=2.2161693733113e+00; 	rate[14][18]=8.5338209390745e+00;
		rate[14][19]=4.3342126659660e+00; 
	
		rate[15][16]=3.1432036618746e+02; 	rate[15][17]=2.0305343047059e+01;
		rate[15][18]=3.4167877957799e+01; 	rate[15][19]=1.0000000000000e+00;

		rate[16][17]=5.2559565123081e+00; 	rate[16][18]=2.0382362288681e+01;
		rate[16][19]=1.0765527137500e+02; 
	
		rate[17][18]=1.3814733274637e+01; 	rate[17][19]=2.8259139240676e+00;

		rate[18][19]=1.0000000000000e+00;
		
		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		double[] f = new double[n];
		f[0]=0.072; f[1]=0.019; f[2]=0.039; f[3]=0.019; f[4]=0.006;
		f[5]=0.025; f[6]=0.024; f[7]=0.056; f[8]=0.028; f[9]=0.088;
		f[10]=0.168; f[11]=0.023; f[12]=0.054; f[13]=0.061; f[14]=0.054;
		f[15]=0.072; f[16]=0.086; f[17]=0.029; f[18]=0.033; f[19]=0.043;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "MTREV amino acid substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("J", "Adachi"),
                    new Author("M", "Hasegawa")
            },
            "Model of amino acid substitution in proteins encoded by mitochondrial DNA",
            1996, "J Mol Evol", 42, 459, 468
    );
}

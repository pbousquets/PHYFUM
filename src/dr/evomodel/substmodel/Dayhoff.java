/*
 * Dayhoff.java
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
 * Dayhoff model for amino acid evolution
 * Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1978)
 * A model of evolutionary change in proteins.
 * Dayhoff, M.O. (ed.) Atlas of Protein Sequence Structur., Vol5, Suppl. 3,
 * National Biomedical Research Foundation, Washington DC, pp. 345-352.
 *
 * @version $Id: Dayhoff.java,v 1.3 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Dayhoff extends EmpiricalRateMatrix.AbstractAminoAcid {
	
	public static final Dayhoff INSTANCE = new Dayhoff();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private Dayhoff() { super("dayhoff");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];
		
		// Q matrix
		rate[0][1]=9.6472567159749e-01; 	rate[0][2]=3.5927991886410e+00;
		rate[0][3]=4.3200552414656e+00; 	rate[0][4]=1.3184584178499e+00;
		rate[0][5]=3.2267534963169e+00; 	rate[0][6]=7.0141987829615e+00;
		rate[0][7]=8.5773867857875e+00; 	rate[0][8]=8.1434196396611e-01;
		rate[0][9]=2.3518447453539e+00; 	rate[0][10]=1.4735711728911e+00;
		rate[0][11]=9.3940162271805e-01; 	rate[0][12]=2.5490196078431e+00;
		rate[0][13]=6.5922920892495e-01; 	rate[0][14]=8.9189834148670e+00;
		rate[0][15]=1.4540712836859e+01; 	rate[0][16]=1.3411904595370e+01;
		rate[0][17]=3.8517964118027e-02; 	rate[0][18]=8.7897227856660e-01;
		rate[0][19]=7.4036511156187e+00; 
	
		rate[1][2]=1.1890243902439e+00; 	rate[1][3]=5.9525626545377e-02;
		rate[1][4]=8.4778922655537e-01; 	rate[1][5]=8.8348561504191e+00;
		rate[1][6]=5.5954088952654e-02; 	rate[1][7]=3.1434881434075e-01;
		rate[1][8]=8.4753987678285e+00; 	rate[1][9]=2.2684090115941e+00;
		rate[1][10]=5.5954088952654e-01; 	rate[1][11]=1.6681312769010e+01;
		rate[1][12]=3.1707317073171e+00; 	rate[1][13]=4.8959827833572e-01;
		rate[1][14]=3.6754156468900e+00; 	rate[1][15]=5.4755072760812e+00;
		rate[1][16]=9.6472567159749e-01; 	rate[1][17]=7.5538020086083e+00;
		rate[1][18]=2.7977044476327e-01; 	rate[1][19]=8.6083213773314e-01;

		rate[2][3]=3.2459324155194e+01; 	rate[2][4]=7.3852625416383e-02;
		rate[2][5]=3.7732198142415e+00; 	rate[2][6]=5.3911764705882e+00;
		rate[2][7]=5.0264375413087e+00; 	rate[2][8]=1.9061418685121e+01;
		rate[2][9]=2.7901430842607e+00; 	rate[2][10]=1.2482698961938e+00;
		rate[2][11]=1.1542279411765e+01; 	rate[2][12]=1.9117647058824e-01;
		rate[2][13]=5.0183823529412e-01; 	rate[2][14]=1.5181660899654e+00;
		rate[2][15]=1.7697478991597e+01; 	rate[2][16]=8.3557302231237e+00;
		rate[2][17]=8.6029411764706e-01; 	rate[2][18]=3.4411764705882e+00;
		rate[2][19]=5.7352941176471e-01; 
	
		rate[3][4]=2.5534152404601e-02; 	rate[3][5]=4.8811013767209e+00;
		rate[3][6]=4.0561952440551e+01; 	rate[3][7]=4.4423506911730e+00;
		rate[3][8]=3.0865788117500e+00; 	rate[3][9]=8.5749078239692e-01;
		rate[3][10]=2.5926985518518e-02; 	rate[3][11]=2.5930851063830e+00;
		rate[3][12]=1.1667143483333e-01; 	rate[3][13]=1.2963492759259e-02;
		rate[3][14]=4.7853935065891e-01; 	rate[3][15]=3.4167709637046e+00;
		rate[3][16]=2.3984722282163e+00; 	rate[3][17]=3.2408731898147e-02;
		rate[3][18]=8.1351689612015e-02; 	rate[3][19]=6.3829787234043e-01;

		rate[4][5]=2.1864264103535e-02; 	rate[4][6]=1.4770525083277e-02;
		rate[4][7]=3.9055458751427e-01; 	rate[4][8]=1.0223340673168e+00;
		rate[4][9]=1.5970515970516e+00; 	rate[4][10]=3.9098448749850e-02;
		rate[4][11]=8.0776309049169e-03; 	rate[4][12]=1.4155086538140e-01;
		rate[4][13]=8.6898395721925e-02; 	rate[4][14]=6.8155604487784e-01;
		rate[4][15]=5.8097784568373e+00; 	rate[4][16]=5.9929928084086e-01;
		rate[4][17]=3.4759358288770e-01; 	rate[4][18]=3.4759358288770e+00;
		rate[4][19]=1.7647058823529e+00; 
	
		rate[5][6]=2.5476780185759e+01; 	rate[5][7]=1.0174974779977e+00;
		rate[5][8]=2.1573939173192e+01; 	rate[5][9]=6.5266504894988e-01;
		rate[5][10]=2.6634492806410e+00; 	rate[5][11]=5.5466331269350e+00;
		rate[5][12]=4.0247678018576e+00; 	rate[5][13]=1.8038017885416e-02;
		rate[5][14]=5.5044618466582e+00; 	rate[5][15]=2.0267580716497e+00;
		rate[5][16]=1.9256432155439e+00; 	rate[5][17]=9.6202762055552e-02;
		rate[5][18]=1.0061919504644e-01; 	rate[5][19]=1.2538699690402e+00;

		rate[6][7]=2.8869795109055e+00; 	rate[6][8]=1.5519031141869e+00;
		rate[6][9]=2.1701112877583e+00; 	rate[6][10]=4.0484429065744e-01;
		rate[6][11]=2.9823529411765e+00; 	rate[6][12]=1.0705882352941e+00;
		rate[6][13]=1.9801735189768e-02; 	rate[6][14]=1.7993079584775e+00;
		rate[6][15]=2.8184873949580e+00; 	rate[6][16]=1.2261663286004e+00;
		rate[6][17]=7.3114099162219e-02; 	rate[6][18]=7.6470588235294e-01;
		rate[6][19]=1.3058823529412e+00; 
	
		rate[7][8]=3.7906768788150e-01; 	rate[7][9]=2.3128004846840e-02;
		rate[7][10]=2.5776602775942e-01; 	rate[7][11]=9.6662260409782e-01;
		rate[7][12]=6.0145406477198e-01; 	rate[7][13]=5.4775280898876e-01;
		rate[7][14]=1.2382877804129e+00; 	rate[7][15]=8.2853366065527e+00;
		rate[7][16]=1.1110604644803e+00; 	rate[7][17]=1.2888301387971e-01;
		rate[7][18]=1.7114723586662e-02; 	rate[7][19]=1.9233311302049e+00;

		rate[8][9]=2.7354343963341e-01; 	rate[8][10]=1.5876246692449e+00;
		rate[8][11]=9.6993944636678e-01; 	rate[8][12]=1.2544085640577e-01;
		rate[8][13]=1.6868512110727e+00; 	rate[8][14]=3.3075513942601e+00;
		rate[8][15]=1.2530894710826e+00; 	rate[8][16]=8.1434196396611e-01;
		rate[8][17]=1.0121107266436e+00; 	rate[8][18]=4.4982698961938e+00;
		rate[8][19]=1.5570934256055e+00; 
	
		rate[9][10]=9.2275320303002e+00; 	rate[9][11]=1.6663354531002e+00;
		rate[9][12]=1.1780604133545e+01; 	rate[9][13]=6.9753577106518e+00;
		rate[9][14]=4.2551201720752e-01; 	rate[9][15]=8.8575970928912e-01;
		rate[9][16]=6.8951811852420e+00; 	rate[9][17]=9.8802836705702e-02;
		rate[9][18]=1.3434022257552e+00; 	rate[9][19]=3.1526232114467e+01;

		rate[10][11]=6.5787197231834e-01; 	rate[10][12]=1.8622837370242e+01;
		rate[10][13]=5.6340830449827e+00; 	rate[10][14]=1.1377976796255e+00;
		rate[10][15]=6.1690558576372e-01; 	rate[10][16]=1.2098794893211e+00;
		rate[10][17]=1.7543252595156e+00; 	rate[10][18]=1.0346020761246e+00;
		rate[10][19]=6.2906574394464e+00; 
	
		rate[11][12]=8.6029411764706e+00; 	rate[11][13]=6.6640454965565e-03;
		rate[11][14]=1.2089100346021e+00; 	rate[11][15]=3.4411764705882e+00;
		rate[11][16]=4.9442190669371e+00; 	rate[11][17]=3.4272233982290e-02;
		rate[11][18]=4.7794117647059e-01; 	rate[11][19]=3.7500000000000e-01;

		rate[12][13]=3.2500000000000e+00; 	rate[12][14]=5.9976931949250e-01;
		rate[12][15]=2.1848739495798e+00; 	rate[12][16]=3.6916835699797e+00;
		rate[12][17]=1.6247577591604e-01; 	rate[12][18]=1.1508700794053e-01;
		rate[12][19]=9.0588235294118e+00; 
	
		rate[13][14]=3.9359861591695e-01; 	rate[13][15]=1.6386554621849e+00;
		rate[13][16]=4.9442190669371e-01; 	rate[13][17]=2.8676470588235e+00;
		rate[13][18]=2.4852941176471e+01; 	rate[13][19]=4.4117647058824e-01;

		rate[14][15]=8.6431043005437e+00; 	rate[14][16]=2.8308077795013e+00;
		rate[14][17]=3.5840244687362e-02; 	rate[14][18]=4.3804743506776e-02;
		rate[14][19]=1.7301038062284e+00; 
	
		rate[15][16]=1.9663865546218e+01; 	rate[15][17]=2.7857142857143e+00;
		rate[15][18]=1.2016806722689e+00; 	rate[15][19]=1.0840336134454e+00;

		rate[16][17]=4.2019597219666e-02; 	rate[16][18]=1.5162271805274e+00;
		rate[16][19]=5.6592292089249e+00; 
	
		rate[17][18]=2.2941176470588e+00; 	rate[17][19]=1.2654363316538e-01;

		rate[18][19]=1.0000000000000e+00; 
		
		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		double[] f = new double[n];
		f[0] = 0.087; f[1] = 0.041; f[2] = 0.040; f[3] = 0.047;
		f[4] = 0.033; f[5] = 0.038; f[6] = 0.05; f[7] = 0.089;
		f[8] = 0.034; f[9] = 0.037; f[10] = 0.085; f[11] = 0.08;
		f[12] = 0.015; f[13] = 0.04; f[14] = 0.051; f[15] = 0.07;
		f[16] = 0.058; f[17] = 0.01; f[18] = 0.03; f[19] = 0.065;

		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "Dayhoff amino acid substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("MO", "Dayhoff"),
                    new Author("RM", "Schwartz"),
                    new Author("BC", "Orcutt")
            },
            "A model of evolutionary change in proteins",
            1972,
            "in Dayhoff, M.O. (ed.) Atlas of Protein Sequence Structur., Vol 5, Suppl. 3",
            5,
            345, 352,
            Citation.Status.PUBLISHED
    );

}

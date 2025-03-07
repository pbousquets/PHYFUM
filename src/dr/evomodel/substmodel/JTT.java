/*
 * JTT.java
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
 * JTT model for amino acid evolution
 * D.T. Jones, W.R. Taylor, and J.M. Thornton
 * The rapid generation of mutation data matrices from protein sequences
 * CABIOS  vol. 8 no. 3 1992 pp. 275-282
 *
 * @version $Id: JTT.java,v 1.3 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class JTT extends EmpiricalRateMatrix.AbstractAminoAcid {
	
	public static final JTT INSTANCE = new JTT();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private JTT() { super("JTT");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];
		
		// Q matrix
		rate[0][1]=3.1628651460584e+00; 	rate[0][2]=3.2804935927860e+00;
		rate[0][3]=4.8477237048666e+00; 	rate[0][4]=3.4612244897959e+00;
		rate[0][5]=3.3130910900946e+00; 	rate[0][6]=6.3199473337722e+00;
		rate[0][7]=1.0440154440154e+01; 	rate[0][8]=1.3061224489796e+00;
		rate[0][9]=2.1726844583987e+00; 	rate[0][10]=1.8443597219107e+00;
		rate[0][11]=2.2137668626773e+00; 	rate[0][12]=2.7210884353741e+00;
		rate[0][13]=8.3265306122449e-01; 	rate[0][14]=1.1537414965986e+01;
		rate[0][15]=2.2838213546288e+01; 	rate[0][16]=2.7007955724663e+01;
		rate[0][17]=5.1311953352770e-01; 	rate[0][18]=8.3673469387755e-01;
		rate[0][19]=1.7474335188621e+01; 
	
		rate[1][2]=2.6598918637222e+00; 	rate[1][3]=9.1014867485456e-01;
		rate[1][4]=6.1624649859944e+00; 	rate[1][5]=1.8036482885837e+01;
		rate[1][6]=1.8924731182796e+00; 	rate[1][7]=8.1810886516769e+00;
		rate[1][8]=1.9119717452198e+01; 	rate[1][9]=1.4410687351864e+00;
		rate[1][10]=2.2211961707760e+00; 	rate[1][11]=3.9239234676922e+01;
		rate[1][12]=2.5060690943044e+00; 	rate[1][13]=3.9439775910364e-01;
		rate[1][14]=4.1953094963476e+00; 	rate[1][15]=5.9016766126741e+00;
		rate[1][16]=3.8437069743152e+00; 	rate[1][17]=7.6766706682673e+00;
		rate[1][18]=1.4173669467787e+00; 	rate[1][19]=1.0308123249300e+00;

		rate[2][3]=3.2226935854843e+01; 	rate[2][4]=1.8710963455150e+00;
		rate[2][5]=4.5351268130622e+00; 	rate[2][6]=3.3951344979102e+00;
		rate[2][7]=4.5987249708180e+00; 	rate[2][8]=2.3693774375271e+01;
		rate[2][9]=2.9235880398671e+00; 	rate[2][10]=8.0960899565551e-01;
		rate[2][11]=1.5024269384537e+01; 	rate[2][12]=1.9003322259136e+00;
		rate[2][13]=4.3853820598007e-01; 	rate[2][14]=7.1083317047749e-01;
		rate[2][15]=2.9456208772690e+01; 	rate[2][16]=1.3735908553410e+01;
		rate[2][17]=1.6706217370669e-01; 	rate[2][18]=4.1661129568106e+00;
		rate[2][19]=9.7452934662237e-01; 
	
		rate[3][4]=6.2857142857143e-01; 	rate[3][5]=3.0662020905923e+00;
		rate[3][6]=4.5450549450549e+01; 	rate[3][7]=7.5402435402435e+00;
		rate[3][8]=6.0544672718586e+00; 	rate[3][9]=6.8808114961961e-01;
		rate[3][10]=3.6130902064968e-01; 	rate[3][11]=1.6718197057180e+00;
		rate[3][12]=1.0879120879121e+00; 	rate[3][13]=1.9340659340659e-01;
		rate[3][14]=7.3949579831933e-01; 	rate[3][15]=3.4196528109572e+00;
		rate[3][16]=2.4749487800335e+00; 	rate[3][17]=3.4536891679749e-01;
		rate[3][18]=2.6895604395604e+00; 	rate[3][19]=1.8608058608059e+00;

		rate[4][5]=5.5191637630662e-01; 	rate[4][6]=3.2442396313364e-01;
		rate[4][7]=3.3297297297297e+00; 	rate[4][8]=4.3726708074534e+00;
		rate[4][9]=9.1868131868132e-01; 	rate[4][10]=9.9466248037677e-01;
		rate[4][11]=2.9830508474576e-01; 	rate[4][12]=2.4095238095238e+00;
		rate[4][13]=4.1485714285714e+00; 	rate[4][14]=7.3949579831933e-01;
		rate[4][15]=1.2862939958592e+01; 	rate[4][16]=2.8125907990315e+00;
		rate[4][17]=6.8244897959184e+00; 	rate[4][18]=1.2885714285714e+01;
		rate[4][19]=3.7714285714286e+00; 
	
		rate[5][6]=2.0316061593796e+01; 	rate[5][7]=1.3922214897825e+00;
		rate[5][8]=3.3861536130889e+01; 	rate[5][9]=4.7172339855267e-01;
		rate[5][10]=4.2320327755868e+00; 	rate[5][11]=1.7835941652395e+01;
		rate[5][12]=2.6573751451800e+00; 	rate[5][13]=2.7595818815331e-01;
		rate[5][14]=9.4992143198743e+00; 	rate[5][15]=3.2350653941322e+00;
		rate[5][16]=3.0973838067678e+00; 	rate[5][17]=1.0512692882031e+00;
		rate[5][18]=1.5331010452962e+00; 	rate[5][19]=1.0778164924506e+00;

		rate[6][7]=6.6857641051189e+00; 	rate[6][8]=1.4458024443999e+00;
		rate[6][9]=6.7068415455512e-01; 	rate[6][10]=5.7932850559579e-01;
		rate[6][11]=1.0365070686558e+01; 	rate[6][12]=1.0138248847926e+00;
		rate[6][13]=2.6359447004608e-01; 	rate[6][14]=1.1291226167887e+00;
		rate[6][15]=1.8337006611901e+00; 	rate[6][16]=1.9520424900414e+00;
		rate[6][17]=6.9519420671494e-01; 	rate[6][18]=3.8018433179723e-01;
		rate[6][19]=2.7772657450077e+00; 
	
		rate[7][8]=1.2113479939567e+00; 	rate[7][9]=3.2670032670033e-01;
		rate[7][10]=4.1817641817642e-01; 	rate[7][11]=1.6354950592239e+00;
		rate[7][12]=7.6447876447876e-01; 	rate[7][13]=3.0579150579151e-01;
		rate[7][14]=1.2391551215081e+00; 	rate[7][15]=1.1138492529797e+01;
		rate[7][16]=1.8888816176952e+00; 	rate[7][17]=3.3491450634308e+00;
		rate[7][18]=3.1853281853282e-01; 	rate[7][19]=2.8416988416988e+00;

		rate[8][9]=1.0931677018634e+00; 	rate[8][10]=3.2194389461470e+00;
		rate[8][11]=3.1498052426571e+00; 	rate[8][12]=1.9130434782609e+00;
		rate[8][13]=2.7329192546584e+00; 	rate[8][14]=6.7304834977469e+00;
		rate[8][15]=4.3726708074534e+00; 	rate[8][16]=2.8162964522581e+00;
		rate[8][17]=7.8083407275954e-01; 	rate[8][18]=3.5118012422360e+01;
		rate[8][19]=7.2877846790890e-01; 
	
		rate[9][10]=1.4069798333535e+01; 	rate[9][11]=1.2292791953809e+00;
		rate[9][12]=2.8366300366300e+01; 	rate[9][13]=4.7384615384615e+00;
		rate[9][14]=5.8780435251023e-01; 	rate[9][15]=2.4105749323141e+00;
		rate[9][16]=1.5243062022723e+01; 	rate[9][17]=8.2888540031397e-01;
		rate[9][18]=1.8434065934066e+00; 	rate[9][19]=5.7699633699634e+01;

		rate[10][11]=8.8039805231089e-01; 	rate[10][12]=2.2425954997384e+01;
		rate[10][13]=1.5099529042386e+01; 	rate[10][14]=6.2626896912611e+00;
		rate[10][15]=3.4917298022888e+00; 	rate[10][16]=1.6109411169944e+00;
		rate[10][17]=3.2366001345593e+00; 	rate[10][18]=1.4505494505495e+00;
		rate[10][19]=1.0557823129252e+01; 
	
		rate[11][12]=3.6577885391445e+00; 	rate[11][13]=1.4915254237288e-01;
		rate[11][14]=1.2868062479229e+00; 	rate[11][15]=2.8162964522581e+00;
		rate[11][16]=5.7494151926786e+00; 	rate[11][17]=5.4790729851263e-01;
		rate[11][18]=5.3268765133172e-01; 	rate[11][19]=7.4899112187248e-01;

		rate[12][13]=2.5666666666667e+00; 	rate[12][14]=9.4491129785247e-01;
		rate[12][15]=1.6397515527950e+00; 	rate[12][16]=1.2180790960452e+01;
		rate[12][17]=1.1972789115646e+00; 	rate[12][18]=1.1130952380952e+00;
		rate[12][19]=1.7746031746032e+01; 
	
		rate[13][14]=8.8739495798319e-01; 	rate[13][15]=5.6298136645963e+00;
		rate[13][16]=8.3099273607748e-01; 	rate[13][17]=3.3224489795918e+00;
		rate[13][18]=3.3392857142857e+01; 	rate[13][19]=3.6000000000000e+00;

		rate[14][15]=1.6261762676085e+01; 	rate[14][16]=6.8852490148602e+00;
		rate[14][17]=4.2256902761104e-01; 	rate[14][18]=6.7787114845938e-01;
		rate[14][19]=1.2549019607843e+00; 
	
		rate[15][16]=2.7891216619293e+01; 	rate[15][17]=1.8740017746229e+00;
		rate[15][18]=3.7349896480331e+00; 	rate[15][19]=2.4182194616977e+00;

		rate[16][17]=4.8702870978900e-01; 	rate[16][18]=1.1985472154964e+00;
		rate[16][19]=6.7925746569814e+00; 
	
		rate[17][18]=4.6020408163265e+00; 	rate[17][19]=1.4693877551020e+00;

		rate[18][19]=1.0000000000000e+00; 
		
		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		double[] f = new double[n];
		f[0] = 0.077; f[1] = 0.051; f[2] = 0.043; f[3] = 0.052;
		f[4] = 0.02; f[5] = 0.041; f[6] = 0.062; f[7] = 0.074;
		f[8] = 0.023; f[9] = 0.052; f[10] = 0.091; f[11] = 0.059;
		f[12] = 0.024; f[13] = 0.04; f[14] = 0.051; f[15] = 0.069;
		f[16] = 0.059; f[17] = 0.014; f[18] = 0.032; f[19] = 0.066;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}


	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "JTT amino acid substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("DT", "Jones"),
                    new Author("WR", "Taylor"),
                    new Author("JM", "Thornton")
            },
            "The rapid generation of mutation data matrices from protein sequences",
            1992,
            "CABIOS",
            8,
            275, 282,
            Citation.Status.PUBLISHED
    );
}

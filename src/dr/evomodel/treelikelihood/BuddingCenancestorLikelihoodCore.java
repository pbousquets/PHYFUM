/*
 * BuddingCenancestorLikelihoodCore.java
 *
 * Diego Mallo
 *
 * This file is part of PHYFUM.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * PHYFUM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  PHYFUM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with PHYFUM; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treelikelihood;

import java.util.Arrays;

/**
 * BuddingCenancestorLikelihoodCore - An implementation of LikelihoodCore for small cell populations that reproduce by budding
 *
 * @author Diego Mallo
 */

public class BuddingCenancestorLikelihoodCore extends GeneralCenancestorLikelihoodCore {
    protected int S;
    protected double [][] pBud;
    protected final int i_state_o = 0;
    protected int i_state_k;
    protected int i_state_m;

    //Vars for calculations so that we are not allocating memory all the time
    protected double [] pLA,pLB;

    /**
     * Constructor
     *
     * @param stateCount number of states
     */
    public BuddingCenancestorLikelihoodCore(int stateCount) {
        super(stateCount);
        this.S = (int) (Math.sqrt(2 * stateCount + 1.0/4.0)-3.0/2.0);
        generatepBud();
        this.i_state_k = S; //For clarity
        this.i_state_m = stateCount - 1;
        pLA = new double[stateCount];
        pLB = new double[stateCount];
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials2, double[] matrices2,
                                                    double[] partials3)
    {
        double sumAbud , sumBbud; //Likelihood sums of the valid budding options on the [A,B] branch
        int u = 0;
        int v = 0;

        for (int l = 0; l < matrixCount; l++) {

            for (int k = 0; k < patternCount; k++) {

                int w = l * matrixSize;

                Arrays.fill(pLA,0.0);
                Arrays.fill(pLB, 0.0);

                //Calculate the partials at the end of the branch of each daughter branch A,B
                for (int i = 0; i < stateCount; i++) { //i = from state
                    for (int j = 0; j < stateCount; j++) { // j = to state
                        pLA[i] += matrices1[w] * partials1[v + j];
                        pLB[i] += matrices2[w] * partials2[v + j];
                        w++;
                    }
                }

                //Calculate the partials for the paternal node, taking into account budding combinations
                //For other implementations this will be a double loop, but for budding it is always 3 so I am doing it manually
                for (int i = 0; i < stateCount; i++) {
                    sumAbud = pLA[i_state_o] * pBud[i][0] + pLA[i_state_k] * pBud[i][1] + pLA[i_state_m] * pBud[i][2]; // pL of each possible budding state * probability of that bud
                    sumBbud = pLB[i_state_o] * pBud[i][0] + pLB[i_state_k] * pBud[i][1] + pLB[i_state_m] * pBud[i][2];
                    partials3[u] = pLA[i] * sumBbud + pLB[i] * sumAbud;
                    u++;
                }

                v += stateCount;
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials2, double[] matrices2,
                                                    double[] partials3, int[] matrixMap)
    {
        double sumAbud , sumBbud; //Likelihood sums of the valid budding options on the [A,B] branch

        int u = 0;
        int v = 0;

        for (int k = 0; k < patternCount; k++) {

            int w = matrixMap[k] * matrixSize;

            Arrays.fill(pLA,0.0);
            Arrays.fill(pLB, 0.0);

            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    pLA[i] += matrices1[w] * partials1[v + j];
                    pLB[i] += matrices2[w] * partials2[v + j];
                    w++;
                }
            }

            //Calculate the partials for the paternal node, taking into account budding combinations
            //For other implementations this will be a double loop, but for budding it is always 3 so I am doing it manually
            for (int i = 0; i < stateCount; i++) {
                sumAbud = pLA[i_state_o] * pBud[i][0] + pLA[i_state_k] * pBud[i][1] + pLA[i_state_m] * pBud[i][2]; // pL of each possible budding state * the number of cells that can generate such bud in the current from state i
                sumBbud = pLB[i_state_o] * pBud[i][0] + pLB[i_state_k] * pBud[i][1] + pLB[i_state_m] * pBud[i][2];
                partials3[u] = pLA[i] * sumBbud + pLB[i] * sumAbud;
                u++;
            }
            v += stateCount;
        }
    }

    //This would assume budding also happens at origing, I do not think it is a good idea, but it could also be selectable by the user
 /*   *//**
     * Calculates partial likelihoods at a node when both children have partials.
     *//*
    protected void calculatePartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials3)
    {
        double sumAbud; //Likelihood sums of the valid budding options on the [A,B] branch
        int u = 0;
        int v = 0;

        for (int l = 0; l < matrixCount; l++) {

            for (int k = 0; k < patternCount; k++) {

                int w = l * matrixSize;
                Arrays.fill(pLA,0.0);

                //Calculate the partials at the end of the branch of each daughter branch A,B
                for (int i = 0; i < stateCount; i++) { //i = from state
                    for (int j = 0; j < stateCount; j++) { // j = to state
                        pLA[i] += matrices1[w] * partials1[v + j];
                        w++;
                    }
                }

                //Calculate the partials for the paternal node, taking into account budding combinations
                //For other implementations this will be a double loop, but for budding it is always 3 so I am doing it manually
                for (int i = 0; i < stateCount; i++) {
                    sumAbud = pLA[i_state_o] * pBud[i][0] + pLA[i_state_k] * pBud[i][1] + pLA[i_state_m] * pBud[i][2]; // pL of each possible budding state * the number of cells that can generate such bud in the current from state i
                    partials3[u] = pLA[i] * 0.5 + sumAbud; //The 0.5 is the probability of not budding. It is in pBud, but since with only one daughter pLA is not multiplied by sumXbud, we need to do it manually
                    u++;
                }

                v += stateCount;
            }
        }
    }

    *//**
     * Calculates partial likelihoods at a node when both children have partials.
     *//*
    protected void calculatePartialsPruning(double[] partials1, double[] matrices1,
                                            double[] partials3, int[] matrixMap) {
        double sumAbud;
        int u = 0;
        int v = 0;

        for (int k = 0; k < patternCount; k++) {

            int w = matrixMap[k] * matrixSize;
            Arrays.fill(pLA,0.0);

            //Calculate the partials at the end of the branch of each daughter branch A,B
            for (int i = 0; i < stateCount; i++) { //i = from state
                for (int j = 0; j < stateCount; j++) { // j = to state
                    pLA[i] += matrices1[w] * partials1[v + j];
                    w++;
                }
            }

            //Calculate the partials for the paternal node, taking into account budding combinations
            //For other implementations this will be a double loop, but for budding it is always 3 so I am doing it manually
            for (int i = 0; i < stateCount; i++) {
                sumAbud = pLA[i_state_o] * pBud[i][0] + pLA[i_state_k] * pBud[i][1] + pLA[i_state_m] * pBud[i][2]; // pL of each possible budding state * the number of cells that can generate such bud in the current from state i
                partials3[u] = pLA[i] * 0.5 + sumAbud; //The 0.5 is the probability of not budding. It is in pBud, but since with only one daughter pLA is not multiplied by sumXbud, we need to do it manually
                u++;
            }

            v += stateCount;
        }

    }*/

    /**
     * Pending to implement. These are not currently in use in PHYFUM but I should still implement them, just not right now
     */

    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                int[] states2, double[] matrices2,
                                                double[] partials3)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }
    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */

    protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
                                                      double[] partials2, double[] matrices2,
                                                      double[] partials3)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }
    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                int[] states2, double[] matrices2,
                                                double[] partials3, int[] matrixMap)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
                                                      double[] partials2, double[] matrices2,
                                                      double[] partials3, int[] matrixMap)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }
    /**
     * Calculates partial likelihoods at a node when the only child has states.
     */
    protected void calculateStatesPruning(int[] states1, double[] matrices1,
                                          double[] partials3)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    /**
     * Calculates partial likelihoods at a node when the only child has states.
     */
    protected void calculateStatesPruning(int[] states1, double[] matrices1,
                                          double[] partials3, int[] matrixMap)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }


    /**
     * This method will create a matrix with the probability of each bud.
     * For each state, we count the number of fully demethylated cells S-k-m, half-methylated k, and fully-methylated m
     * We divide that by the number of total cells * 2 since both sides can bud
     * For S=2, possible states are for (k,m): 2,0,0 1,0,1 0,0,2 1,1,0, 0,2,0 0,1,1
     **/
    private void generatepBud(){
        this.pBud = new double[stateCount][];
        double combs = S*2; //Any cell can bud, of each A or B daughter
        int ii = 0;
        for (int m = 0; m <= S; m++){
            for (int k = 0; k+m <= S; k++){
                    this.pBud[ii] = new double[]{(S-k-m)/combs, k/combs, m/combs};
                    ii++;
            }
        }
    }
}

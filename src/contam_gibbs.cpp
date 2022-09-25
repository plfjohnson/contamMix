/* $Id
  -----------------------------------------------------------------------------
  Copyright 2012,2013 Philip Johnson.

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published
     by the Free Software Foundation, either version 3 of the License,
     or (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License at <http://www.gnu.org/licenses/> for
     more details.

  Gibbs sampler routine for contamination estimation. (This version
  assumes contamination is first genome)

  R CMD SHLIB contam_gibbs.cpp
  -----------------------------------------------------------------------------
*/

#define R_NO_REMAP
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include "Rwrappers.h"
#include <vector>
#include <iostream>

using namespace std;

// based off gsl_ran_dirichlet -- this way can directly access R's
// gamma rng, etc
void ran_dirichlet (double priorAlpha, const int K, const int alpha[],
                    double theta[])
{
  for (int i = 0; i < K; ++i) {
      theta[i] = rgamma(alpha[i] + priorAlpha, 1.0);
  }

  double norm = 0.0;
  for (int i = 0; i < K; ++i) {
      norm += theta[i];
  }
  for (int i = 0; i < K; ++i) {
      theta[i] /= norm;
  }
}


// Exported C entrypoints for calling from R
extern "C" {
    // 2-stage gibbs sampler that ASSUMES the authentic proportion is
    // first in the list
    SEXP gibbs_step(SEXP s_data, SEXP s_proportions, SEXP s_dirichAlpha) {
        GetRNGstate();
        CRMatrix<double> data(s_data);
        CRVector<double> prop(s_proportions);
        if (data.ncol() != prop.size()) {
            Rf_error("mismatch in dimensions for data parameter");
        }
        CRVector<double> dirichAlpha(s_dirichAlpha);
        if (dirichAlpha.size() != 1) {
            Rf_error("prior parameter (alpha) should be a scalar value");
        }
        double priorAlpha = dirichAlpha[0];

        // STEP ONE: augment data with random assignment of READs
        // (R_i) to GENOME (j) accordingly to proportion (p):

        // Z_i ~ k(j|R_i,p) \propto Pr(R_i|j) Pr(j|p_j)

        // STEP TWO: summarize Z_i as mtCounts_j = \sum_i(Z_i==j)
        vector<int> mtCounts(prop.size(), 0);//initialize to 0
        for (unsigned int readI = 0;  readI < data.nrow();  ++readI) {
            double tot = 0;
            for (unsigned int mtI = 0;  mtI < data.ncol();  ++mtI) {
                tot += prop[mtI] * data(readI, mtI);
            }

            double r = unif_rand()*tot;
            tot = 0;
            for (unsigned int mtI = 0;  mtI < data.ncol();  ++mtI) {
                tot += prop[mtI] * data(readI, mtI);
                if (tot > r) {
                    //cout << mtI << " ";
                    ++mtCounts[mtI];
                    break;
                }
            }
        }
        //cout << endl;
        /*
        for (unsigned int i = 0;  i < mtCounts.size();  ++i) {
            cout << mtCounts[i] << " ";
        }
        cout << endl;
        */

        // STEP THREE: draw proportions (p) according to their
        // posterior probability.  Use Dirichelet prior for everything
        // *except* authentic, which is unif(0,1). (actually should be
        // 0.5 - 1 given assumptions of the model, but getting the
        // maths to work out seems awkward)

        // p_a ~ Beta
        // p ~ Dirichlet(mtCounts + alpha) = Pr(p | Z, alpha)
        CRVector<double> newProp(prop.size(), true);

        // prior is beta(1,1) to be uniform [0,1] --> posterior is
        // beta(1+a, 1+(n-a))
        int numAuthentic = mtCounts.data()[0];
        newProp[0] =
            rbeta(1 + numAuthentic, 1 + data.nrow() - numAuthentic);

        // prior is dirichlet for potential contaminate genomes -->
        // posterior is alpha=priorAlpha+counts
        ran_dirichlet(priorAlpha, mtCounts.size()-1, mtCounts.data()+1,
                      newProp.data()+1);
        // scale down dirichlet to be in 1-authentic
        for (unsigned int i = 1;  i < newProp.size();  ++i) {
            //cout << newProp.data()[i] << " ";
            newProp.data()[i] *= 1-newProp[0];
        }
        //cout << endl;
        PutRNGstate();
        return newProp;
    }
}

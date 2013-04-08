#include "assert.h"
#include "math.h"
#include "mex.h"

#define ARG_UNARY_COST_SUPPORT 0
#define ARG_UNARY_COST_METACLASS 1
#define ARG_TC_V 2
#define ARG_TC_H 3
#define ARG_MIN_HEIGHTS 4
#define ARG_MAX_HEIGHTS 5
#define ARG_MIN_H_DISTS 6
#define ARG_CONSTS 7

#define RET_COSTS 0
#define RET_A 1
#define RET_B 2
#define RET_AEQ 3
#define RET_BEQ 4

// Notes on the implementation:
//
// W[V|H] is organized in the following way: Each matrix is a 4x4xRpxR matrix.
// So each block of 16 represents the cost of t2xt1. Next, each 16 block is
// layed out such that we navigate the Rp's, then navigates the Rs.
// 



// Adds the costs from the data.
void AddCostsDataTerm(double* costs,
                      const unsigned int num_unknowns_S,
                      const unsigned int num_unknowns_M,
                      const double* unary_cost_support,
                      const double* unary_cost_metaclass,
                      const double alpha,
                      const double beta,
                      const unsigned int r,
                      const unsigned int g) {

  int rp = r + 1;
  
  // First term, local support costs.
  for (int ii = 0; ii < num_unknowns_S; ++ii) {
    costs[ii] += alpha * unary_cost_support[ii];
  }
  
  // Move the location of consts to the beginning of the M matrix.
  costs = costs + num_unknowns_S;
  
  // Second term, local metaclass costs.
  for (int ii = 0; ii < num_unknowns_M; ++ii) {
    costs[ii] += beta * unary_cost_metaclass[ii];
  }
}

// Adds the costs from the global consistency term.
void AddCostsGlobalGroundConsistency(double* costs,
                                     const unsigned int r,
                                     const double* min_heights,
                                     const unsigned int offset_M,
                                     const double kappa) {
  
  const unsigned int rp = r + 1;
  
  // First, find the min of the min heights.
  double min_all = 10e10;
  for (int ii = 0; ii < r; ++ii) {
    if (min_all > min_heights[ii]) {
      min_all = min_heights[ii];
    }
  }
  
  // Next add cost proportional to the distance from the min height.
  for (int ii = 0; ii < r; ++ii) {
    if ((min_all + .03) < min_heights[ii]) {
      costs[offset_M + ii * 4] += kappa;
    }
  }
}

void AddCostsTransition(double* costs,
                        const unsigned int r,
                        const unsigned int g,
                        const double* TC_V,
                        const double* TC_H,
                        const double lambda,
                        const unsigned int offset_WV,
                        const unsigned int offset_WH) {

  const unsigned int row_size = r * 16 + 4;
  double cost;
  
  for (int ii = 0; ii < r; ++ii) {
    
    // Add costs for each WV observed region.    
    for (int jj = 0; jj < r; ++jj) {
      for (int t1 = 0; t1 < 4; ++t1) {
        for (int t2 = 0; t2 < 4; ++t2) {
          
          if (ii == jj) {
            cost = 10e10;
          } else {
            cost = lambda * TC_V[t1 * 4 + t2];
          }
          
          costs[offset_WV + ii * row_size + jj * 16 + t1 * 4 + t2] = cost;
        }
      }
    }
    
    // Add costs for the unobserved regions.
    for (int t1 = 0; t1 < 4; ++t1) {
      costs[offset_WV + ii * row_size + r * 16 + t1] = lambda * TC_V[t1 * 4];
    }

    // Add costs for each WH observed region.    
    for (int jj = 0; jj < r; ++jj) {
      for (int t1 = 0; t1 < 4; ++t1) {
        for (int t2 = 0; t2 < 4; ++t2) {
          
          if (ii == jj) {
            cost = lambda;
          } else {
            cost = lambda * TC_H[t1 * 4 + t2];
          }
          
          costs[offset_WH + ii * row_size + jj * 16 + t1 * 4 + t2] = cost;
        }
      }
    }
    
    // Add costs for the unobserved regions.
    for (int t1 = 0; t1 < 4; ++t1) {
      costs[offset_WH + ii * row_size + r * 16 + t1] = lambda * TC_H[t1 * 4];
    }
  }
}

// Adds the costs from the Consistency cost term.
void AddCostsConsistency(double* costs,
                         const unsigned int r, 
                         const unsigned int g,
                         const double* min_heights,
                         const double* max_heights,
                         const double* min_horz_dists,
                         const double tau,
                         const double upsilon) {
 
  const unsigned int rp = r + 1;
  
  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < rp; ++jj) {
      costs[ii * g + jj] += tau * pow(min_heights[ii] - max_heights[jj], 2);
    }
    
    for (int jj = rp; jj < 2 * rp; ++jj) {
      int sup_ind = jj-r-1;
      double min_horz_dist = min_horz_dists[ii * rp + sup_ind];
      bool bad_height = min_heights[ii] < min_heights[sup_ind]
          || max_heights[ii] > max_heights[sup_ind];
      costs[ii * g + jj] += upsilon * (pow(min_horz_dist,2) + 
           (bad_height ? 1 : 0));
    }
  }
}

// Sets up several groups of Inequalities:
//   1. All hidden variables must be greater than or equal to 0
//   2. All hidden variables must be less than or equal to 1.
//   3. S(i,r') + M(k,1) <= 1 for all i and k
//   4. S(i,2r') + M(k,1) <= 1 for all i and k
//
// Args:
//   plhs - the return arguments.
//   num_unknowns - the number of unknowns in the LP problem.
//   r - the number of regions in the current image.
//   g - equal to 2*rp + 1
void SetInequalities(mxArray* plhs[],
    const unsigned int num_unknowns,
    const unsigned int r,
    const unsigned int g,
    const unsigned int offset_M,
    const unsigned int offset_WV,
    const unsigned int offset_WH) {
  
  const unsigned int rp = r + 1;
  
  
  int num_equations = 
            num_unknowns // geq 0
          + num_unknowns // leq 1
          + 4 * r  // Each pair (i,t1) <= to M(i,t1) over all j, t2.
          + r * (2 * (4 * r)) // Limit to columns of J inside each block.
//           + 4 * r  // Each pair (j,t2) sums to M(j,t2) over all i, t1.
          + r * r * 2 // Each pair (i,rp) + M(i,1) < 2
          + r * 2 * rp // Each block in W <= S(i,j)
          + 4 * r * r + 4 * r * r + 4 * r + 4 * r
          + 1; 
  
  int nnz_ineq = 
            num_unknowns // geq 0
          + num_unknowns // leq 1
          + 4 * r * (1 + 2 * (4 * r + 1))  // Each pair (i,t1) <= to M(i,t1) over all j, t2.
          + 5 * r * (2 * (4 * r)) // Limit to columns of J inside each block.
//           + 4 * r * (4 * r * 2 + 1)        // Each pair (j,t2) sums to M(j,t2) over all i, t1
          + r * r * 2 * 2 // Each pair (i,rp) + M(i,1) < 2
          + r * 2 * (17 * r + 5)    // Each block in W <= S(i,j)
          + (4 * r * r) * 5 + (4 * r * r) * 5 + (4 * r) * 2 + (4 * r) * 2
          + r * 2;
  

          
  plhs[RET_A] = mxCreateSparse(num_unknowns, num_equations, nnz_ineq, mxREAL);
  
  double* A = mxGetPr(plhs[RET_A]);
  mwIndex* A_r = mxGetIr(plhs[RET_A]);
  mwIndex* A_c = mxGetJc(plhs[RET_A]);

  plhs[RET_B] = mxCreateDoubleMatrix(num_equations, 1, mxREAL);
  double* b = (double*) mxGetData(plhs[RET_B]);
  
  // Make sure all unknowns are geq to zero.
  int p = 0;
  int equation_ind = 0;
  for (int ii = 0; ii < num_unknowns; ++ii, ++p, ++equation_ind) {
    A_c[equation_ind] = p;
    
    A[p] = -1;
    A_r[p] = ii;
  }
    
  // Make sure all unknowns are leq to one.
  for (int ii = 0; ii < num_unknowns; ++ii, ++p, ++equation_ind) {
    A_c[equation_ind] = p;
  
    A[p] = 1;
    A_r[p] = ii;    
    
    b[p] = 1;
  }

  // *********************************************************************************
  //   Add an inequality constraint such that sum_j sum_t2 W_{i,j,t1,t2} = M_{i,t1}
  // *********************************************************************************
  unsigned int row_size_w = r * 16 + 4;
  
  for (int ii = 0; ii < r; ++ii) {
    for (int t1 = 0; t1 < 4; ++t1, ++equation_ind) {
      A_c[equation_ind] = p;
      b[equation_ind] = 0;
      
      // First, set the negative M element.
      A[p] = -1;
      A_r[p] = offset_M + ii * 4 + t1;
      ++p;
      
      // Next, grab the observed verticals.
      for (int jj = 0; jj < r; ++jj) {
        for (int t2 = 0; t2 < 4; ++t2, ++p) {
          A[p] = 1;
          A_r[p] = offset_WV + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }      
      }
      
      // Next, grab the hidden floor vertical.
      A[p] = 1;
      A_r[p] = offset_WV + ii * row_size_w + r * 16 + t1;
      ++p;
      
      // Next, grab the observed horizontals.
      for (int jj = 0; jj < r; ++jj) {
        for (int t2 = 0; t2 < 4; ++t2, ++p) {
          A[p] = 1;
          A_r[p] = offset_WH + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }      
      }
              
      // Finally, grab the hidden floor horizontal.
      A[p] = 1;
      A_r[p] = offset_WH + ii * row_size_w + r * 16 + t1;
      ++p;
    }
  }

  // Limits the value of each cell inside each 4-block.
  for (int ii = 0; ii < r; ++ii) {
    // First the observed verticals.
    for (int jj = 0; jj < r; ++jj) {
      for (int t1 = 0; t1 < 4; ++t1, ++equation_ind) {
        A_c[equation_ind] = p;
        b[equation_ind] = 0;
        
        // First set the negative M element.
        A[p] = -1;
        A_r[p] = offset_M + ii * 4 + t1;
        ++p;
        
        // Next sum over each of the jj.
        for (int t2 = 0; t2 < 4; ++t2, ++p) {
          A[p] = 1;
          A_r[p] = offset_WV + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }
      }
    }
    
    // Next the hidden verticals
    for (int t1 = 0; t1 < 4; ++t1, ++equation_ind) {
      A_c[equation_ind] = p;
      b[equation_ind] = 0;

      // First set the negative M element.
      A[p] = -1;
      A_r[p] = offset_M + ii * 4 + t1;
      ++p;

      // Next the hidden t1.
      A[p] = 1;
      A_r[p] = offset_WV + ii * row_size_w + r * 16 + t1;
      ++p;
    }
    
    // Next the observed horizontals.
    for (int jj = 0; jj < r; ++jj) {
      for (int t1 = 0; t1 < 4; ++t1, ++equation_ind) {
        A_c[equation_ind] = p;
        b[equation_ind] = 0;
        
        // First set the negative M element.
        A[p] = -1;
        A_r[p] = offset_M + ii * 4 + t1;
        ++p;
        
        // Next sum over each of the jj.
        for (int t2 = 0; t2 < 4; ++t2, ++p) {
          A[p] = 1;
          A_r[p] = offset_WH + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }
      }
    }
    
    
    // Next the hidden horizontal
    for (int t1 = 0; t1 < 4; ++t1, ++equation_ind) {
      A_c[equation_ind] = p;
      b[equation_ind] = 0;

      // First set the negative M element.
      A[p] = -1;
      A_r[p] = offset_M + ii * 4 + t1;
      ++p;

      // Next the hidden t1.
      A[p] = 1;
      A_r[p] = offset_WV + ii * row_size_w + r * 16 + t1;
      ++p;
    }
  }
  
  // Limits the value of each column of J inside each 4-block.
  for (int ii = 0; ii < r; ++ii) {
    // First, the observed verticals.
    for (int jj = 0; jj < r; ++jj) {
      for (int t2 = 0; t2 < 4; ++t2, ++equation_ind) {
        A_c[equation_ind] = p;
        b[equation_ind] = 0;
        
        // First, set the negative M element.
        A[p] = -1;
        A_r[p] = offset_M + jj * 4 + t2;
        ++p;
      
        // Next set the values for each t1.
        for (int t1 = 0; t1 < 4; ++t1, ++p) {
          A[p] = 1;
          A_r[p] = offset_WV + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }
      }
    }
    
    // Next the observed hiddens.
    for (int jj = 0; jj < r; ++jj) {
      for (int t2 = 0; t2 < 4; ++t2, ++equation_ind) {
        A_c[equation_ind] = p;
        b[equation_ind] = 0;
        
        // First, set the negative M element.
        A[p] = -1;
        A_r[p] = offset_M + jj * 4 + t2;
        ++p;
      
        // Next set the values for each t1.
        for (int t1 = 0; t1 < 4; ++t1, ++p) {
          A[p] = 1;
          A_r[p] = offset_WH + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }
      }
    }
  }
  
  // Second sum  (k, t2)
//   
//   for (int jj = 0; jj < r; ++jj) {
//     for (int t2 = 0; t2 < 4; ++t2, ++equation_ind) {
//       A_c[equation_ind] = p;
//       b[equation_ind] = 0;
//       
//       // First, set the negative M element.
//       A[p] = -1;
//       A_r[p] = offset_M + jj * 4 + t2;
//       ++p;
//       
//       // Next grab the observed verticals.
//       for (int ii = 0; ii < r; ++ii) {
//         for (int t1 = 0; t1 < 4; ++t1, ++p) {
//           A[p] = 1;
//           A_r[p] = offset_WV + ii * row_size_w + jj * 16 + t1 * 4 + t2;
//         }
//       }
//       
//       // Next, grab the observed horizontals.
//       for (int ii = 0; ii < r; ++ii) {
//         for (int t1 = 0; t1 < 4; ++t1, ++p) {
//           A[p] = 1;
//           A_r[p] = offset_WH + ii * row_size_w + jj * 16 + t1 * 4 + t2;
//         }      
//       }
//     }
//   }

  // Each block of W <= S(i,j)
  for (int ii = 0; ii < r; ++ii) {
    
    // First, the WV observed blocks.
    for (int jj = 0; jj < r; ++jj, ++equation_ind) {
      A_c[equation_ind] = p;
      b[equation_ind] = 0;
      
      // First, set the negative S element.
      A[p] = -1;
      A_r[p] = ii * g + jj;
      ++p;
      
      // Next, turn on all of the elements in the current block.
      for (int t1 = 0; t1 < 4; ++t1) {
        for (int t2 = 0; t2 < 4; ++t2, ++p) {
          A[p] = 1;
          A_r[p] = offset_WV + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }
      }
    }
    
    // Next, the WV hidden block for the current ii.
    A_c[equation_ind] = p;
    b[equation_ind] = 0;
    ++equation_ind;
    
    // First, set the negative S element.
    A[p] = -1;
    A_r[p] = ii * g + r;
    ++p;
    
    for (int t1 = 0; t1 < 4; ++t1, ++p) {
      A[p] = 1;
      A_r[p] = offset_WV + ii * row_size_w + r * 16 + t1;
    }
    
    
    // Next, the WH observed blocks.
    for (int jj = 0; jj < r; ++jj, ++equation_ind) {
      A_c[equation_ind] = p;
      b[equation_ind] = 0;
      
      // First, set the negative S element.
      A[p] = -1;
      A_r[p] = ii * g + rp + jj;
      ++p;
      
      // Next, turn on all of the elements in the current block.
      for (int t1 = 0; t1 < 4; ++t1) {
        for (int t2 = 0; t2 < 4; ++t2, ++p) {
          A[p] = 1;
          A_r[p] = offset_WH + ii * row_size_w + jj * 16 + t1 * 4 + t2;
        }
      }
    }
    
    // Finally, the WH hidden block for the current ii.
    A_c[equation_ind] = p;
    b[equation_ind] = 0;
    ++equation_ind;
    
    // First, set the negative S element.
    A[p] = -1;
    A_r[p] = ii * g + 2 * rp - 1;
    ++p;
    
    for (int t1 = 0; t1 < 4; ++t1, ++p) {
      A[p] = 1;
      A_r[p] = offset_WH + ii * row_size_w + r * 16 + t1;
    }
  }
  
  // Add constraints for support from below by hidden region.
  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < r; ++jj, ++equation_ind) {
      A_c[equation_ind] = p;
      b[equation_ind] = 1;

      A[p] = 1;
      A_r[p] = ii * g + rp - 1;
      ++p;

      A[p] = 1;
      A_r[p] = offset_M + jj * 4;
      ++p;
    }
  }
  
  // Add constraints for support from behind by hidden region.
  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < r; ++jj, ++equation_ind) {
      A_c[equation_ind] = p;
      b[equation_ind] = 1;

      A[p] = 1;
      A_r[p] = ii * g + 2 * rp - 1;
      ++p;

      A[p] = 1;
      A_r[p] = offset_M + jj * 4;
      ++p;
    }
  }
  
  A_c[equation_ind] = p;
  b[equation_ind] = -1;
  ++equation_ind;
  for (int ii = 0; ii < r; ++ii) {
    A[p] = -1;
    A_r[p] = ii * g + rp - 1;
    ++p;
    
    A[p] = -1;
    A_r[p] = offset_M + ii * 4;
    ++p;
  }
  
  assert(num_equations == equation_ind);
  assert(nnz_ineq == p);
  A_c[num_equations] = p; 
}


// Sets up 4 groups of equalities:
//   1. All rows of S must sum to one.
//   2. All rows of M must sum to one.
//   3. All rows of WV must sum to one.
//   4. All rows of WH must sum to one.
//   5. The ground terms for S must be equal to the floor terms of M.
// 
//
// Args:
//   plhs - the return arguments.
//   num_unknowns - the number of unknowns in the LP problem.
//   r - the number of regions in the current image.
//   g - equal to 2*rp + 1
void SetEqualities(mxArray* plhs[],
    const int num_unknowns,
    const int r,
    const int g,
    const unsigned int offset_M,
    const unsigned int offset_WV,
    const unsigned int offset_WH) {

  const unsigned int rp = r + 1;
  
  // Calculate the number of indices per row of W that are turned on.
  const unsigned int num_per_row_w = 2 * (16 * r + 4) + 1;
  
  // Calculate the number of equations:
  const unsigned int num_equations =
       r      // Rows of S sum to 1.
     + r      // Rows of M sum to 1.
     + r      // Last column of S equals last column of M.
     + r      // Number of elements in WV+WH+S(:,g)
     + r * rp;// Blocks of WV sum to S(i,j)
  
  // Calculate the number of non-zero entries in the matrix.
  const int nnz_eq =
        r * g  // Rows of S sum to 1.
      + r * 4  // Rows of M sum to 1.
      + 2 * r  // Last column of S equals last column of M.
      + r * num_per_row_w // Number of elements in WV+WH+S(:,g)
      + r * rp + r * r * 16 + r * 4; // Blocks of WV sum to S(i,j)
          
  plhs[RET_AEQ] = mxCreateSparse(num_unknowns, num_equations, nnz_eq, mxREAL);
  
  double* Aeq = mxGetPr(plhs[RET_AEQ]);
  mwIndex* Aeq_r = mxGetIr(plhs[RET_AEQ]);
  mwIndex* Aeq_c = mxGetJc(plhs[RET_AEQ]);

  plhs[RET_BEQ] = mxCreateDoubleMatrix(num_equations, 1, mxREAL);
  double* beq = (double*) mxGetData(plhs[RET_BEQ]);
  
  int equation_ind = 0;
  
  // **************************************
  //   Ensure that the rows of S sum to 1.
  // **************************************
  int p = 0;
  for (int ii = 0; ii < r; ++ii, ++equation_ind) {
    for (int jj = 0; jj < g; ++jj, ++p) {
      Aeq[p] = 1;
      Aeq_r[p] = ii * g + jj;
    }
    
    assert(ii * g == p);
    beq[equation_ind] = 1;
    Aeq_c[equation_ind] = ii * g;
  }
  
  // **************************************
  //   Ensure that the rows of M sum to 1.
  // **************************************
  for (int ii = 0; ii < r; ++ii, ++equation_ind) {
    beq[equation_ind] = 1;
    Aeq_c[equation_ind] = p;
    
    for (int jj = 0; jj < 4; ++jj, ++p) {
      Aeq[p] = 1;
      Aeq_r[p] = r * g + ii * 4 + jj;
    }
  }
  
  // *******************************************************************
  //   Ensure that the last column of S matches the first column of M.
  // *******************************************************************
  for (int ii = 0; ii < r; ++ii, ++equation_ind) {
    Aeq_c[equation_ind] = p;
    beq[equation_ind] = 0;
    
    // Ground in S
    Aeq[p] = 1;
    Aeq_r[p] = (ii + 1) * g - 1;
    ++p;
    
    // Floor in M
    Aeq[p] = -1;
    Aeq_r[p] = offset_M + ii * 4;
    ++p;  
  }
  
  const unsigned int row_size_W = 16 * r + 4;
    
  // **************************************
  //   Ensure that the rows of W sum to 1.
  // **************************************
  for (int ii = 0; ii < r; ++ii, ++equation_ind) {
    beq[equation_ind] = 1;
    Aeq_c[equation_ind] = p;
    
    // First the Ground.
    Aeq[p] = 1;
    Aeq_r[p] = (ii + 1) * g - 1;
    ++p;
    
    // Then the current row of WV.
    for (int jj = 0; jj < row_size_W; ++jj, ++p) {
      Aeq[p] = 1;
      Aeq_r[p] = offset_WV + ii * row_size_W + jj;
    }
    
    for (int jj = 0; jj < row_size_W; ++jj, ++p) {
      Aeq[p] = 1;
      Aeq_r[p] = offset_WH + ii * row_size_W + jj;
    }
  }
  
  // **********************************************
  //   Ensure that the blocks of WV sum to S(i,j)
  // **********************************************
  for (int ii = 0; ii < r; ++ii) {
    for (int jj = 0; jj < r; ++jj, ++equation_ind) {
      beq[equation_ind] = 0;
      Aeq_c[equation_ind] = p;
      
      // First, turn the S(i,j) term on.
      Aeq[p] = -1;
      Aeq_r[p] = ii * g + jj;
      ++p;
      
      // Next, turn on the WV-block.
      for (int t1 = 0; t1 < 4; ++t1) {
        for (int t2 = 0; t2 < 4; ++t2, ++p) {
          Aeq[p] = 1;
          Aeq_r[p] = offset_WV + ii * row_size_W + jj * 16 + t1 * 4 + t2;
        }
      }
    }
    
    // Finally, add the equation for the final 4-vector column
    beq[equation_ind] = 0;
    Aeq_c[equation_ind] = p;
    ++equation_ind;
    
    // First, turn the S(i,j) term on.
    Aeq[p] = -1;
    Aeq_r[p] = ii * g + r;
    ++p;
    
    for (int t1 = 0; t1 < 4; ++t1, ++p) {
      Aeq[p] = 1;
      Aeq_r[p] = offset_WV + ii * row_size_W + r * 16 + t1;
    }
  }
  
  
  Aeq_c[num_equations] = p;
}

void SetValidS(double* valid_s, const unsigned int r, const unsigned int g,
    const double* min_heights, const double* max_heights,
    const double* min_horz_dists) {
  
  const unsigned int rp = r + 1;
  
  for (int ii = 0; ii < r; ++ii) {
    
    for (int jj = 0; jj < r; ++jj) {
      if (ii == jj) {
        continue;
      }
      
      // Is the current support lower than the current prop?
      double min_horz_dist = min_horz_dists[ii * rp + jj];
      
      if (min_heights[jj] <= min_heights[ii] && min_horz_dist < 3) {
        valid_s[ii * g + jj] = 1;
        
//         valid_s[ii * g + jj] = true;
//         valid_s[ii * g + rp + jj] = true;
      }
    }
  }
}

// Returns:
//   f - the costs for each variable.
//   A - NxI where N is the number of unknowns are I is the number of inequalities.
//   b - Ix1 where I is the number of inequalities.
//   Aeq - NxIeq matrix where N is the number of unknowns and Ieq is the number
//         of equalities.
//   beq - Ieq x 1 matrix where Ieq is the number of equality constraints.
void mexFunction(int nlhs, mxArray* plhs[],
                 const int nrhs, const mxArray* prhs[]) {

  if (nrhs != 8) {
    mexErrMsgTxt("Exactly 8 arguments required");
  }

  // ********************************
  //   Grab the unary support costs
  // ********************************
  if (mxGetNumberOfDimensions(prhs[ARG_UNARY_COST_SUPPORT]) != 2) {
    mexErrMsgTxt("unaryCostSupport must be of size GxR");
  } else if (mxGetClassID(prhs[ARG_UNARY_COST_SUPPORT]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("unaryCostSupport must be of type 'double'");
  }
  
  const unsigned int g = mxGetM(prhs[ARG_UNARY_COST_SUPPORT]);
  const unsigned int r = mxGetN(prhs[ARG_UNARY_COST_SUPPORT]);
  const unsigned int rp = r + 1;
  const double* unary_cost_support = 
          (double*) mxGetData(prhs[ARG_UNARY_COST_SUPPORT]);
  
  // ***********************************
  //   Grab the unary meta class costs
  // ***********************************
  if (mxGetNumberOfDimensions(prhs[ARG_UNARY_COST_METACLASS]) != 2) {
    mexErrMsgTxt("unaryCostMetaclass must be of size GxR");
  } else if (mxGetClassID(prhs[ARG_UNARY_COST_METACLASS]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("unaryCostMetaclass must be of type 'double'");
  } else if (mxGetN(prhs[ARG_UNARY_COST_METACLASS]) != r) {
    mexErrMsgTxt("unaryCostMetaclass must be of size 4xR");
  }
  
  const double* unary_cost_metaclass = 
          (double*) mxGetData(prhs[ARG_UNARY_COST_METACLASS]);
  
  // ********************************************
  //   Get the vertical transition cost matrix.
  // ********************************************
  if (mxGetNumberOfDimensions(prhs[ARG_TC_V]) != 2
      || mxGetM(prhs[ARG_TC_V]) != 4
      || mxGetN(prhs[ARG_TC_V]) != 4) {
    mexErrMsgTxt("TC_V must be of size 4x4");
  } else if (mxGetClassID(prhs[ARG_TC_V]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("TC_V must be of type 'double'");
  }
  const double* TC_V = (double*) mxGetData(prhs[ARG_TC_V]);
  
  // **********************************************
  //   Get the horizontal transition cost matrix.
  // **********************************************
  if (mxGetNumberOfDimensions(prhs[ARG_TC_H]) != 2
      || mxGetM(prhs[ARG_TC_H]) != 4
      || mxGetN(prhs[ARG_TC_H]) != 4) {
    mexErrMsgTxt("TC_H must be of size 4x4");
  } else if (mxGetClassID(prhs[ARG_TC_H]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("TC_H must be of type 'double'");
  }
  const double* TC_H = (double*) mxGetData(prhs[ARG_TC_H]);
  
  // ************************
  //   Get the min heights.
  // ************************
  if (mxGetNumberOfDimensions(prhs[ARG_MIN_HEIGHTS]) != 2
      || mxGetNumberOfElements(prhs[ARG_MIN_HEIGHTS]) != r + 1) {
    mexErrMsgTxt("minHeights must be of size (r+1)x1");
  } else if (mxGetClassID(prhs[ARG_MIN_HEIGHTS]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("minHeights must be of class 'double'");
  }
  const double* min_heights = (double*) mxGetData(prhs[ARG_MIN_HEIGHTS]);
  
  // ************************
  //   Get the max heights.
  // ************************
  if (mxGetNumberOfDimensions(prhs[ARG_MAX_HEIGHTS]) != 2
      || mxGetNumberOfElements(prhs[ARG_MAX_HEIGHTS]) != r + 1) {
    mexErrMsgTxt("maxHeights must be of size (r+1)x1");
  } else if (mxGetClassID(prhs[ARG_MAX_HEIGHTS]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("maxHeights must be of class 'double'");
  }
  const double* max_heights = (double*) mxGetData(prhs[ARG_MAX_HEIGHTS]);
  
  // *********************************
  //   Get the min horizontal dists.
  // *********************************
  if (mxGetNumberOfDimensions(prhs[ARG_MIN_H_DISTS]) != 2
      || mxGetM(prhs[ARG_MIN_H_DISTS]) != r + 1
      || mxGetN(prhs[ARG_MIN_H_DISTS]) != r ) {
    mexErrMsgTxt("minHorzDists must be of size (r+1)xr");
  } else if (mxGetClassID(prhs[ARG_MIN_H_DISTS]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("minHorzDists must be of class 'double'");
  }
  const double* min_horz_dists = (double*) mxGetData(prhs[ARG_MIN_H_DISTS]);
  
  // ***********************
  //   Get the parameters.
  // ***********************
  if (mxGetClassID(prhs[ARG_CONSTS]) != mxSTRUCT_CLASS) {
    mexErrMsgTxt("Argument 9 must be a struct.");
  } else if (mxGetFieldNumber(prhs[ARG_CONSTS], "alpha") == -1) {
    mexErrMsgTxt("Missing field optConsts.alpha");
  } else if (mxGetFieldNumber(prhs[ARG_CONSTS], "beta") == -1) {
    mexErrMsgTxt("Missing field optConsts.beta");
  } else if (mxGetFieldNumber(prhs[ARG_CONSTS], "lambda") == -1) {
    mexErrMsgTxt("Missing field optConsts.lambda");
  } else if (mxGetFieldNumber(prhs[ARG_CONSTS], "tau") == -1) {
    mexErrMsgTxt("Missing field optConsts.tau");
  } else if (mxGetFieldNumber(prhs[ARG_CONSTS], "upsilon") == -1) {
    mexErrMsgTxt("Missing field optConsts.upsilon");
  } else if (mxGetFieldNumber(prhs[ARG_CONSTS], "kappa") == -1) {
    mexErrMsgTxt("Missing field optConsts.kappa");
  }
  
  const double alpha = (double) mxGetScalar(mxGetField(prhs[ARG_CONSTS], 0, "alpha"));
  const double beta = (double) mxGetScalar(mxGetField(prhs[ARG_CONSTS], 0, "beta"));
  const double kappa = (double) mxGetScalar(mxGetField(prhs[ARG_CONSTS], 0, "kappa"));
  const double lambda = (double) mxGetScalar(mxGetField(prhs[ARG_CONSTS], 0, "lambda"));
  const double tau = (double) mxGetScalar(mxGetField(prhs[ARG_CONSTS], 0, "tau"));
  const double upsilon = (double) mxGetScalar(mxGetField(prhs[ARG_CONSTS], 0, "upsilon"));
  
  // Calculate the total number of unknowns.
  const unsigned int num_unknowns_S = g * r;
  const unsigned int num_unknowns_M = r * 4;
  const unsigned int num_unknowns_WV = r * r * 16 + 4 * r;
  const unsigned int num_unknowns_WH = r * r * 16 + 4 * r;
  const unsigned int num_unknowns = num_unknowns_S + num_unknowns_M
      + num_unknowns_WV + num_unknowns_WH;
  
  const unsigned int offset_S = 0;
  const unsigned int offset_M = offset_S + num_unknowns_S;
  const unsigned int offset_WV = offset_M + num_unknowns_M;
  const unsigned int offset_WH = offset_WV + num_unknowns_WV;
  
  // Setup the output costs
  plhs[RET_COSTS] = mxCreateDoubleMatrix(num_unknowns, 1, mxREAL);
  double* costs = (double*) mxGetData(plhs[RET_COSTS]);
  
  for (int ii = 0; ii < num_unknowns; ++ii) {
    costs[ii] = 0; 
  }
  
  AddCostsDataTerm(costs, num_unknowns_S, num_unknowns_M, unary_cost_support,
      unary_cost_metaclass, alpha, beta, r, g);
  
  AddCostsGlobalGroundConsistency(costs, r, min_heights, offset_M, kappa);
  
  AddCostsTransition(costs, r, g, TC_V, TC_H, lambda, offset_WV, offset_WH);
  
  AddCostsConsistency(costs, r, g, min_heights, max_heights, min_horz_dists,
          tau, upsilon);
  
  SetInequalities(plhs, num_unknowns, r, g, offset_M, offset_WV, offset_WH);
  SetEqualities(plhs, num_unknowns, r, g, offset_M, offset_WV, offset_WH);
  
}

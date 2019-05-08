//  g++ glpk_simple_test2.cpp -L /usr/local/lib -l glpk -o glpk_simple_test2
//  ./glpk_simple_test2
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

int main(void)
{
    glp_prob *lp;
    int row_index[1 + 1000];                    //Row indices of each element
    int col_index[1 + 1000];                    //column indices of each element
    double value[1 + 1000];                     //numerical values of corresponding elements
    double z, x1, x2, x3;
    lp = glp_create_prob();                     //creates a problem object
    glp_set_prob_name(lp, "sample");            //assigns a symbolic name to the problem object
    glp_set_obj_dir(lp, GLP_MIN);               //calls the routine glp_set_obj_dir to set the
                                                //omptimization direction flag,
                                                //where GLP_MAX means maximization
 
    //ROWS
    glp_add_rows(lp, 2);                        //adds three rows to the problem object
    //row 1
    //glp_set_row_name(lp, 1, "p");               //assigns name p to first row
    glp_set_row_bnds(lp, 1, GLP_LO, 1.0, 1.0);   //sets the type and bounds of the first row,
                                                //where GLP_LO means that the row has an lower bound.
                                                //300 <= p <= +inf    
    //row 2
    //glp_set_row_name(lp, 2, "q");               //assigns name q to second row
    glp_set_row_bnds(lp, 2, GLP_LO, 1.0, 1.0);//48 <= q <= +inf 
 
    //COLUMNS
    glp_add_cols(lp, 2);                        //adds three columns to the problem object
    //column 1
    //glp_set_col_name(lp, 1, "x1");              //assigns name x1 to first column
    glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);  //sets the type and bounds to the first column,
                                                //where GLP_LO means that the column has an lower bound
    glp_set_obj_coef(lp, 1, 0.0);               //sets the objective coefficient for thr first column
                                                //P = 5 * c + 2.5 * t
    //column 2
    //glp_set_col_name(lp, 2, "x2");              //assigns name x2 to first column
    glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);  //sets the type and bounds to the second column
    glp_set_obj_coef(lp, 2, 1.0);               //sets the objective coefficient for thr second column
 
    /*
    p = 560 * c + 320 * t
    q = 136 * c + 40 * t
    */
 
    row_index[1] = 1, col_index[1] = 1, value[1] = 1.0;   // a[1,1] = 560.0
    row_index[2] = 1, col_index[2] = 2, value[2] = 1.0;   // a[1,2] = 320.0
    row_index[3] = 2, col_index[3] = 1, value[3] = -1.0;   // a[2,1] = 136.0
    row_index[4] = 2, col_index[4] = 2, value[4] = 1.0;   // a[2,2] = 40.0
     
    for (int i = 1; i < 5; i++) {
        std::cout << value[i];
        std::cout << ((i % 2 == 0) ? "\n" : "\t");
    }
 
    glp_load_matrix(lp, 4, row_index, col_index, value);    //calls the routine glp_load_matrix
                                                            //loads information from three arrays
                                                            //into the problem object
    glp_simplex(lp, NULL);                      //calls the routine glp_simplex
                                                //to solve LP problem
    z = glp_get_obj_val(lp);                    //obtains a computed value of the objective function
    x1 = glp_get_col_prim(lp, 1);               //obtain computed values of structural variables (columns)
    x2 = glp_get_col_prim(lp, 2);               //obtain computed values of structural variables (columns)
 
    printf("\nPrize(z) = %g; chicken(x1) = %g; tofu(x2) = %g;\n", z, x1, x2); //writes out the optimal solution
    glp_delete_prob(lp);                        //calls the routine glp_delete_prob, which frees all the memory

}

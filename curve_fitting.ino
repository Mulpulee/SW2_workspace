// Arduino pin assignment
#define PIN_IR A0

void setup()
{
  Serial.begin(1000000);  // initialize serial port
}

void loop()
{
  unsigned int filtered; // Voltage values from the IR sensor (0 ~ 1023)

  double coef[8];
  double y[7] = { 0, 50, 100, 150, 200, 250, 300 };
  double x[7] = { };

  Serial.println();
  Serial.println();
  Serial.println("This code was written using : https://github.com/nasa/polyfit/blob/main/PolyFit.cpp");
  Serial.println("Details in code");
  Serial.println();

  Serial.print("Polynomial order(under 7) : ");
  while (Serial.available() == 0)
    ;
  int k = Serial.parseInt();
  Serial.print(k);
  Serial.println();
  Serial.println();

  Serial.println("Start measuring.");
  Serial.read();
  
  for (int i = 0; i < 7; i++)
  {
    Serial.println();
    Serial.print("Put the ball on ");
    Serial.print(y[i]);
    Serial.println("mm and press ENTER.");
    
    while (Serial.available() == 0)
    ;
    Serial.read();
    
    filtered = ir_sensor_filtered(20, 0.5, 0); // Replace n with your desired value
    Serial.print("FLT:"); Serial.println(filtered);
    x[i] = filtered;
  }
    
  get_polyfit(k, 7, y, x, coef);

  Serial.print("Curve-fit equation for ");
  Serial.print(k);
  Serial.print(" orders is : ");
  
  for (int i = 0; i < k + 1; i++)
  {
    Serial.print(coef[i], 8);

    switch (i)
    {
      case 0: break;
      case 1: Serial.print(" * x"); break;
      default: Serial.print(" * pow(x, "); Serial.print(i); Serial.print(")"); break;
    }

    if (i < k) Serial.print(" + ");
  }
  
  while (1) ;
}

int compare(const void *a, const void *b) {
  return (*(unsigned int *)a - *(unsigned int *)b);
}

unsigned int ir_sensor_filtered(unsigned int n, float position, int verbose)
{
  // Eliminate spiky noise of an IR distance sensor by repeating measurement and taking a middle value
  // n: number of measurement repetition
  // position: the percentile of the sample to be taken (0.0 <= position <= 1.0)
  // verbose: 0 - normal operation, 1 - observing the internal procedures, and 2 - calculating elapsed time.
  // Example 1: ir_sensor_filtered(n, 0.5, 0) => return the median value among the n measured samples.
  // Example 2: ir_sensor_filtered(n, 0.0, 0) => return the smallest value among the n measured samples.
  // Example 3: ir_sensor_filtered(n, 1.0, 0) => return the largest value among the n measured samples.

  // The output of Sharp infrared sensor includes lots of spiky noise.
  // To eliminate such a spike, ir_sensor_filtered() performs the following two steps:
  // Step 1. Repeat measurement n times and collect n * position smallest samples, where 0 <= postion <= 1.
  // Step 2. Return the position'th sample after sorting the collected samples.

  // returns 0, if any error occurs

  unsigned int *ir_val, ret_val;
  unsigned int start_time;
 
  if (verbose >= 2)
    start_time = millis(); 

  if ((n == 0) || (n > 100) || (position < 0.0) || (position > 1))
    return 0;
    
  if (position == 1.0)
    position = 0.999;

  if (verbose == 1) {
    Serial.print("n: "); Serial.print(n);
    Serial.print(", position: "); Serial.print(position); 
    Serial.print(", ret_idx: ");  Serial.println((unsigned int)(n * position)); 
  }

  ir_val = (unsigned int *)malloc(sizeof(unsigned int) * n);
  if (ir_val == NULL)
    return 0;

  if (verbose == 1)
    Serial.print("IR:");
  
  for (int i = 0; i < n; i++) {
    ir_val[i] = analogRead(PIN_IR);
    if (verbose == 1) {
        Serial.print(" ");
        Serial.print(ir_val[i]);
    }
  }

  if (verbose == 1)
    Serial.print  ("  => ");

  qsort(ir_val, n, sizeof(unsigned int), compare);
  ret_val = ir_val[(unsigned int)(n * position)];

  if (verbose == 1) {
    for (int i = 0; i < n; i++) {
        Serial.print(" ");
        Serial.print(ir_val[i]);
    }
    Serial.print(" :: ");
    Serial.println(ret_val);
  }
  free(ir_val);

  if (verbose >= 2) {
    Serial.print("Elapsed time: "); Serial.print(millis() - start_time); Serial.println("ms");
  }
  
  return ret_val;
}


// ********************************************************************
// * Code PolyFit                                                     *
// * Written by Ianik Plante                                          *
// *                                                                  *
// * KBR                                                              *
// * 2400 NASA Parkway, Houston, TX 77058                             *
// * Ianik.Plante-1@nasa.gov                                          *
// *                                                                  *
// * This code is used to fit a series of n points with a polynomial  *
// * of degree k, and calculation of error bars on the coefficients.  *
// * If error is provided on the y values, it is possible to use a    *
// * weighted fit as an option. Another option provided is to fix the *
// * intercept value, i.e. the first parameter.                       *
// *                                                                  *
// * This code has been written partially using data from publicly    *
// * available sources.                                               *
// *                                                                  *  
// * The code works to the best of the author's knowledge, but some   *   
// * bugs may be present. This code is provided as-is, with no        *
// * warranty of any kind. By using this code you agree that the      * 
// * author, the company KBR or NASA are not responsible for possible *
// * problems related to the usage of this code.                      * 
// *                                                                  *   
// * The program has been reviewed and approved by export control for * 
// * public release. However some export restriction exists. Please   *    
// * respect applicable laws.                                         *
// *                                                                  *   
// ********************************************************************


/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */


// Initialize a 2D array
// **************************************************************
double** Make2DArray(const size_t rows, const size_t cols) {

    double** array;

    array = new double* [rows];
    for (size_t i = 0; i < rows; i++) {
        array[i] = new double[cols];
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            array[i][j] = 0.;
        }
    }

    return array;

}

// Free a 2D array
// **************************************************************
void Free2DArray(double** array, const size_t rows) {
    for (size_t i = 0; i < rows; i++) {
        delete[] array[i];
    }
    delete[] array;
}

// Transpose a 2D array
// **************************************************************
double** MatTrans(double** array, const size_t rows, const size_t cols) {
    double** arrayT = Make2DArray(cols, rows);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            arrayT[j][i] = array[i][j];
        }
    }

    return arrayT;

}

// Perform the multiplication of matrix A[m1,m2] by B[m2,m3]
// **************************************************************
double** MatMul(const size_t m1, const size_t m2, const size_t m3, double** A, double** B) {
    double** array = Make2DArray(m1, m3);

    for (size_t i = 0; i < m1; i++) {
        for (size_t j = 0; j < m3; j++) {
            array[i][j] = 0.;
            for (size_t m = 0; m < m2; m++) {
                array[i][j] += A[i][m] * B[m][j];
            }
        }
    }
    return array;

}

// Perform the multiplication of matrix A[m1,m2] by vector v[m2,1]
// **************************************************************
void MatVectMul(const size_t m1, const size_t m2, double** A, double* v, double* Av) {


    for (size_t i = 0; i < m1; i++) {
        Av[i] = 0.;
        for (size_t j = 0; j < m2; j++) {
            Av[i] += A[i][j] * v[j];
        }
    }


}

double determinant(double** a, const size_t n) {
    double det = 1.0;
    for (size_t i = 0; i < n; i++) {
        size_t pivot = i;
        for (size_t j = i + 1; j < n; j++) {
            if (abs(a[j][i]) > abs(a[pivot][i])) {
                pivot = j;
            }
        }
        if (pivot != i) {
            double* temp = a[i];
            a[i] = a[pivot];
            a[pivot] = temp;
            det *= -1;
        }
        if (a[i][i] == 0) {
            return 0;
        }
        det *= a[i][i];
        for (size_t j = i + 1; j < n; j++) {
            double factor = a[j][i] / a[i][i];
            for (size_t k = i + 1; k < n; k++) {
                a[j][k] -= factor * a[i][k];
            }
        }
    }
    return det;
}

// Perform the 
// **************************************************************
void transpose(double** num, double** fac, double** inverse, const size_t r) {

    double** b = Make2DArray(r, r);
    double deter;

    for (size_t i = 0; i < r; i++) {
        for (size_t j = 0; j < r; j++) {
            b[i][j] = fac[j][i];
        }
    }

    deter = determinant(num, r);

    for (size_t i = 0; i < r; i++) {
        for (size_t j = 0; j < r; j++) {
            inverse[i][j] = b[i][j] / deter;
        }
    }

    Free2DArray(b, r);

}

// Calculates the cofactors 
// **************************************************************
void cofactor(double** num, double** inverse, const size_t f)
{

    double** b = Make2DArray(f, f);
    double** fac = Make2DArray(f, f);

    size_t m;
    size_t n;

    for (size_t q = 0; q < f; q++) {

        for (size_t p = 0; p < f; p++) {

            m = 0;
            n = 0;

            for (size_t i = 0; i < f; i++) {

                for (size_t j = 0; j < f; j++) {

                    if (i != q && j != p) {

                        b[m][n] = num[i][j];

                        if (n < (f - 2)) {
                            n++;
                        }
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
        }
    }

    transpose(num, fac, inverse, f);

    Free2DArray(b, f);
    Free2DArray(fac, f);
}


// Perform the fit of data n data points (x,y) with a polynomial of order k
// **************************************************************
void PolyFit(const double* x, double* y, const size_t n, const size_t k, const bool fixedinter,
    const double fixedinterval, double* beta, double** Weights, double** XTWXInv) {

    // Definition of variables
    // **************************************************************
    double** X = Make2DArray(n, k + 1);           // [n,k+1]
    double** XT;                               // [k+1,n]
    double** XTW;                              // [k+1,n]
    double** XTWX;                             // [k+1,k+1]

    double* XTWY = new double[k + 1];
    double* Y = new double[n];

    size_t begin = 0;
    if (fixedinter) begin = 1;

    // Initialize X
    // **************************************************************
    for (size_t i = 0; i < n; i++) {
        for (size_t j = begin; j < (k + 1); j++) {  // begin
            X[i][j] = pow(x[i], j);
        }
    }

    // Matrix calculations
    // **************************************************************
    XT = MatTrans(X, n, k + 1);                 // Calculate XT
    XTW = MatMul(k + 1, n, n, XT, Weights);         // Calculate XT*W
    XTWX = MatMul(k + 1, n, k + 1, XTW, X);           // Calculate (XTW)*X

    if (fixedinter) XTWX[0][0] = 1.;

    cofactor(XTWX, XTWXInv, k + 1);             // Calculate (XTWX)^-1

    for (size_t m = 0; m < n; m++) {
        if (fixedinter) {
            Y[m] = y[m] - fixedinterval;
        }
        else {
            Y[m] = y[m];
        }
    }
    MatVectMul(k + 1, n, XTW, Y, XTWY);             // Calculate (XTW)*Y
    MatVectMul(k + 1, k + 1, XTWXInv, XTWY, beta);    // Calculate beta = (XTWXInv)*XTWY

    if (fixedinter) beta[0] = fixedinterval;

    Free2DArray(X, n);
    delete[] XTWY;
    delete[] Y;
    Free2DArray(XT, k + 1);
    Free2DArray(XTW, k + 1);
    Free2DArray(XTWX, k + 1);

}


// The main program
// **************************************************************
void get_polyfit(int pk, int pn, double y[7], double x[7], double* coefbeta) {

    const size_t k = pk;                                    // Polynomial order
    const size_t n = pn;                                    // Number of data points
    bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
    double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
    double alphaval = 0.05;                          // Critical apha value

    double** XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
    double** Weights;                                // Matrix Weights [n,n]


    XTWXInv = Make2DArray(k + 1, k + 1);
    Weights = Make2DArray(n, n);

    // Build the weight matrix
    // **************************************************************
    for (size_t i = 0; i < n; i++) Weights[i][i] = 1.;

    if (determinant(Weights, n) == 0. || k > n) {
        return -1;
    }

    // Calculate the coefficients of the fit
    // **************************************************************
    PolyFit(x, y, n, k, fixedinter, fixedinterval, coefbeta, Weights, XTWXInv);


    Free2DArray(XTWXInv, k + 1);
    Free2DArray(Weights, n);

}

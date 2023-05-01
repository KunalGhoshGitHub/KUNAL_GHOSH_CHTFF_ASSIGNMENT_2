#include <iostream>
#include <string>
#include <cassert>
#include <cmath>
#include<algorithm>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>
#include<cassert>

using namespace std;

typedef double real;

double const pi = 4.0*atan(1.0);

// Generate 1D grid
void generate_1D_grid(vector<real> &x, real X0, real Xl, real N)
{
    for (int i = 0; X0 + ((i/N)*(Xl-X0)) <= Xl;i++)
    {
        x.push_back(X0 + ((i/N)*(Xl-X0)));
        //cout<<X0 + ((i/N)*(Xl-X0))<<endl;
    }
}

// This function calculates the transpose of a given matrix
void matrix_transpose(vector<vector<real>>& A, vector<vector<double>> &A_Transpose)
{
    int rows = A.size();
    int cols = A[0].size();

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            A_Transpose[col][row] = A[row][col];
        }
    }
}

// This function calculates the eigen values , eigen vector matrix and inverse of eigen vector matrix
void compute_eigen_value_vector_inverse(vector<vector<real>> &P, vector<vector<real>> &PINV, vector<real> &eigen_values, vector<real> &x,real Nx,real hx)
{
    // Set up eigenvalues and matrices
    int i, j;
	for(i = 1 ; i < Nx ; i++) {
		eigen_values.push_back(-4.0*sin(pi*0.5*x[i])*sin(pi*0.5*x[i]));
		for(j = 1 ; j < Nx ; j++) {
			P[i-1][j-1] = sin(i*pi*x[j]) ;
			PINV[i-1][j-1] = 2.0*hx*sin(j*pi*x[i]) ;
		}
	}
}

// This function calculate the RHS of the equation
void compute_rhs(vector<vector<real>> &RHS, vector<real> &x, vector<real> &y, real h)
{
    int rows = RHS.size();
    int cols = RHS[0].size();
    real pi = 4.0*atan(1);

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            RHS[row][col] = h*h*(1000.0 - (200.0*pi*pi))*sin(10.0*pi*x[row+1])*cos(10.0*pi*y[col+1]);
        }
    }
}

// This function computes the product of two given matrix
void matrix_multiply(vector<vector<real>> &A, vector<vector<double>> &B, vector<vector<double>> &C)
{
    int m = A.size();
    int n = B.size();
    int p = B[0].size();

    // vector<vector<double>> C(m, vector<double>(p, 0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    //return C;
}

// This function computes the V matrix
void compute_V(vector<vector <real>> &F_tilde, vector<real> &eigen_values, vector<vector<real>> &V, real h)
{
    int rows = F_tilde.size();
    int cols = F_tilde[0].size();

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            V[row][col] = (F_tilde[row][col])/((eigen_values[row] + eigen_values[col])+(1000.0*h*h));
        }
    }
}

// Function to save a vector<double> to a file of given name
void write_to_file(vector<real> &u, string str)
{
    ofstream file;
    // Open the file
    file.open(str);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    cout<<"Output file "<< str <<" is opened."<<endl;

    // Writing the vector values to the file in scientific notation
    for(int i=0;i<u.size();i++)
    {
        file<<u[i]<<scientific<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

// Function to save a vector<vector <double>> (MATRIX) to a file of given name
void write_to_file(vector<vector <real> > &u, string str)
{
    ofstream file;
    // Open the file
    file.open(str);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    cout<<"Output file "<< str <<" is opened."<<endl;

    int rows = u.size();
    int cols = u[0].size();

    // Writing the vector values to the file in scientific notation

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            file<<u[row][col]<<scientific<<",";
        }
        file<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

// This function read the inputs from a file
vector<real> input_parameters(string file_name)
{
    // Vector to read the input parameters from a input file
    vector<real> input;
    string item_name;
    int i = 0;

    // Open the file
    ifstream file;
    file.open(file_name);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: Input file could not be opened" << endl;
        exit(1);
    }

    cout<<"Input file "<<file_name<<" is opened."<<endl;

    string line;
    while (getline(file, line))
    {
        // Classifying a string as double if the first character is a numeric
        if(line[0]!= '/')
        {
            // To ensure that we are not reading white space characters
            if(isdigit(line[0]))
            {
                input.push_back(stod(line));
            }
        }
    }

    // Closing the input file
    file.close();

    cout<<"Input file "<<file_name<<" is closed."<<endl;
    cout<<endl;

    return input;
}

int main()
{
    real X0,Xl,Y0,Yl,h;
    int N;

    // Given the domain of the solution
    X0 = 0;
    Xl = 1;
    Y0 = 0;
    Yl = 1;
    vector<real> L2_error_Vector, max_error_vector,Values_of_N;

    // Reading the values of N from a file
    Values_of_N = input_parameters("Values_of_N.txt");
    for (int g= 0;g < Values_of_N.size() ;g++)
    {
        N = Values_of_N[g];
        h = 1.0/N;

        int rows, cols;
        rows = N-1;
        cols = N-1;

        vector<vector <real> > P(rows,vector<real> (cols,0.0));
        vector<real> x;
        vector<real> y;
        vector<vector <real> > RHS(rows,vector<real> (cols,0.0));
        vector<vector <real> > P_Transpose(rows,vector<real> (cols,0.0));
        vector<vector <real> > P_INV_Transpose(rows,vector<real> (cols,0.0));
        vector<vector <real> > P_INV(rows,vector<real> (cols,0.0));
        vector<vector <real> > temp(rows,vector<real> (cols,0.0));
        vector<vector <real> > F_tilde(rows,vector<real> (cols,0.0));
        vector<vector <real> > temp1(rows,vector<real> (cols,0.0));
        vector<vector <real> > V(rows,vector<real> (cols,0.0));
        vector<vector <real> > U(rows,vector<real> (cols,0.0));
        vector<vector <real> > U_analytical(rows,vector<real> (cols,0.0));
        vector<real> eigen_values;
        vector<vector <real>> U_extended(N,vector<real> (N,0.0));

        // Generating the Grids
        generate_1D_grid(x,X0,Xl,N);
        generate_1D_grid(y,Y0,Yl,N);

        // Computing the eigen values, eigen vectors and the inverse of the matrix of eigen vectors
        compute_eigen_value_vector_inverse(P,P_INV,eigen_values,x,N,h);

        // Calculating the Transpose of P such that A = P Lambda P^T
        matrix_transpose(P,P_Transpose);

        // Calculating the Transpose of P^(-1)
        matrix_transpose(P_INV,P_INV_Transpose);

        // Calculating the RHS of the given equation
        compute_rhs(RHS,x,y,h);

        // Applying the boundary conditions
        for (int row = 0;row< rows; row++)
        {
            RHS[row][0] = RHS[row][0] - sin(10.0*pi*x[row+1]);
            RHS[row][rows-1] = RHS[row][rows-1] - sin(10.0*pi*x[row+1]);
        }

        // Multiplying: P^(-1) RHS = temp
        matrix_multiply(P_INV,RHS,temp);

        // Multiplying: F_tilde = temp (P^(-1))^T
        matrix_multiply(temp,P_INV_Transpose,F_tilde);

        // Calculating V
        compute_V(F_tilde, eigen_values,V,h);

        // Multiplying: P V = temp1
        matrix_multiply(P,V,temp1);

        // Multiplying: temp1 P^T = U
        matrix_multiply(temp1,P_Transpose,U);

        vector<vector <real> > U_Numerical(rows+2,vector<real> (cols+2,0.0));
        vector<vector <real> > U_Analytical(rows+2,vector<real> (cols+2,0.0));

        // Saving the matrix U (Numerical Solution)
        for (int row = 0; row < rows+2; row++)
        {
            for (int col = 0; col < cols+2; col++)
            {
                // Using the boundary condition
                if (col == 0)
                {
                    U_Numerical[row][col] = sin(10.0*pi*x[row]);
                }
                // Using the boundary condition
                if (col == rows +1)
                {
                    U_Numerical[row][col] = sin(10.0*pi*x[row]);
                }
                // Copying the U previously computed
                if (row < rows && col < cols)
                {
                    U_Numerical[row+1][col+1] = U[row][col];
                }
            }
        }

        // Computing the analytical solution
        // Saving the matrix U (Analytical Solution)
        for (int row = 0; row < rows+2; row++)
        {
            for (int col = 0; col < cols+2; col++)
            {
                U_Analytical[row][col] = sin(10.0*pi*x[row])*cos(10.0*pi*y[col]);
            }
        }


        // Writing the Numerical Solution to the file
        write_to_file(U_Numerical,"Q_3_Numerical_Solution_"+to_string(N)+"_.csv");

        // Writing the analytical solution
        write_to_file(U_Analytical,"Q_3_Analytical_Solution_"+to_string(N)+"_.csv");

        // Write x to a file
        write_to_file(x,"Q_3_x_grid_"+to_string(N)+"_.csv");

        // Write y to a file
        write_to_file(y,"Q_3_y_grid_"+to_string(N)+"_.csv");


        // Write a Output file showing values of x,y, Numerical Solution and Analytical Solution
        ofstream File("Q_3_Output_"+to_string(N)+"_.dat", ios::out) ;
        File.flags( ios::dec | ios::scientific );
        File.precision(16) ;
        if(!File) {cerr<< "Error: Output file couldnot be opened.\n";}

        real Max_Error, L2_Error;
        Max_Error = L2_Error = 0.0 ;

        File << "TITLE = Flow" << endl << "VARIABLES = X, Y, u, Exact " << endl;
        File << "Zone T = psi I = " << N+1 << " J = " << N+1 << endl ;
        int i,j;

        for(i = 0 ; i < N-1 ; i++) {
            for(j = 0 ; j < N-1 ; j++) {
                if( fabs(U_Numerical[i][j] - U_Analytical[i][j] ) > Max_Error)
                {
                    Max_Error = fabs( U_Numerical[i][j] - U_Analytical[i][j] );
                }
                L2_Error += (U_Numerical[i][j] - U_Analytical[i][j])*(U_Numerical[i][j] - U_Analytical[i][j] )/ ( ( N+1.0 )*(N+1.0) ) ;

                File << x[i] << "\t" << y[j] << "\t" << U_Numerical[i][j] << "\t" << U_Analytical[i][j] << endl ;
            }
        }


        L2_Error = sqrt(L2_Error) ;
        File.close() ;
        // Printing the value of N
        cout<<"For N = "<<N<<endl;

        // Printing the value of L2 Error and Maximum Error
        cout << "\n L2 : " << L2_Error << "\t Max : " << Max_Error <<  endl ;
        L2_error_Vector.push_back(L2_Error);
        max_error_vector.push_back(Max_Error);

	}

	// Writing the L2 Errors in a file
    write_to_file(L2_error_Vector,"Q_3_L2_Error.csv");

    // Writing the Max Errors in a file
	write_to_file(max_error_vector,"Q_3_max_error.csv");
}

#include <iostream>
#include <string>
#include <cassert>
#include <cmath>
#include<vector>
#include<fstream>
#include<cassert>

using namespace std;

typedef double real;

double const pi = 4.0*atan(1);

// This function generate 1-dimensional the grid of x between X0 and Xl. (X0 < Xl)
// N = Number of steps
// So, N+1 is the number of the grid points
void generate_1D_grid(vector<real> &x, real X0, real Xl, real N)
{
    for (int i = 0; X0 + ((i/N)*(Xl-X0)) <= Xl;i++)
    {
        x.push_back(X0 + ((i/N)*(Xl-X0)));
    }
}

// This function will be computing the RHS of the equtaion
// u" = − sin(beta*pi*x), where, x = (0, 1)
void compute_rhs(vector<real> &RHS, vector<real> &x, real beta, real h)
{
    // This loop will iterate over different rows of the RHS
    for(int i = 1;i < x.size()-1;i++)
    {
        // Populating the RHS vector with the value of RHS at the respetive grid points
        // WITHOUT implementing the boundary condition
        RHS.push_back(-h*h*sin(beta*pi*x[i]));
    }
}

// This function will compute the exact solution for the ODE: u" = − sin(beta*pi*x)
// At the specified grid points
void compute_exact_sol(vector<real> &u_exact, vector<real> &x, real beta)
{
    // This loop will iterate over the all the grid points for some N
    for (int i = 0; i < x.size();i++)
    {
        // Calculating and storing the value of the exact solution at the respective grid points
        u_exact.push_back((sin(beta*pi*x[i])/(beta*pi*beta*pi)) - ((sin(beta*pi)/(beta*pi*beta*pi))*x[i]));
    }
}

// This function will compute the L2 error
real compute_L2_error(vector<real> &u_exact, vector<real> &sol)
{
    // Verifying the size of the numerical solution and exact solution vectors
    assert(u_exact.size() == sol.size());

    real sum,L2_error;
    sum = 0;

    // This loop will iterate over each element of the numerical solution and exact solution vectors
    // and adding up the square of the error
    for (int i = 0; i < sol.size();i++)
    {
        sum  = sum + ((u_exact[i]-sol[i])*(u_exact[i]-sol[i]));
    }

    // Taking the average of the sum of the squares
    sum = sum/u_exact.size();

    // Taking the square root of the average of the sum of the squares
    L2_error = pow(sum,0.5);

    return L2_error;
}

// This function will compute the maximum error
real compute_max_error(vector<real> &u_exact, vector<real> &sol)
{
    // Verifying the size of the numerical solution and exact solution vectors
    assert(u_exact.size() == sol.size());

    real max_error;
    max_error = 0;

    // This loop will iterate over each element of the numerical solution and exact solution vectors
    // and calculating the maximum of the error
    for (int i = 0; i < sol.size();i++)
    {
        if (abs(u_exact[i] - sol[i])> max_error)
        {
            max_error = abs(u_exact[i] - sol[i]);
        }
    }
    return max_error;
}

// Thomas_Algorithm_Matix_Solver(vector<vector <real> > &A, vector<real> &B, vector<real> &X)
// solves a system of equations: AX = B, where A is a Tridiagonal matrix
void Thomas_Algorithm_Matix_Solver(vector<vector <real> > &A, vector<real> &B, vector<real> &X)
{
    // Forward Sweep
    for (int row = 1;row < B.size(); row++)
    {
        A[row][row] = A[row][row] - (A[row-1][row]*(A[row][row-1]/A[row-1][row-1]));
        B[row] = B[row] - (B[row-1]*(A[row][row-1]/A[row-1][row-1]));
        A[row][row-1] = 0.0;
    }

    // Backward Sweep
    for (int row = B.size()-1; row >= 1; row--)
    {
        B[row-1] = B[row-1] - (B[row]*(A[row-1][row]/A[row][row]));
        A[row-1][row] = 0.0;
    }

    for (int row = 0;row < B.size(); row++)
    {
        X[row] = B[row]/A[row][row];
    }
}

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

// Function to save a vector<double> to a file of given name
void write_to_file(vector<double> &u, string str)
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

int main()
{
    real X0,Xl,u0,ul,h, beta;
    int N;
    vector<real> Value_of_N, Boundary_Conditions, Interval, L2_error, max_error, Values_of_Beta;

    // Opening the input file to read the value of beta
    Values_of_Beta = input_parameters("Values_of_Beta.txt");

    // Opening the input file to read the values of N
    Value_of_N = input_parameters("Value_of_N.txt");

    // Opening the input file to read the values of u(X0) and  u(Xl) (Boundary Conditions)
    Boundary_Conditions = input_parameters("Boundary_Conditions.txt");

    // Opening the input file to read the Interval
    Interval = input_parameters("Interval.txt");

    // Interval starts at X0 and ends at Xl
    X0 = Interval[0];
    Xl = Interval[1];

    // Assigning the value of boundary conditions to separate variables
    u0 = Boundary_Conditions[0];
    ul = Boundary_Conditions[1];

    // Assigning the value of N
    N = Value_of_N[0];

    // This loop will iterate over different values of beta
    for (int i=0;i < Values_of_Beta.size();i++)
    {
        // Assigning the value of beta to a variable
        beta = Values_of_Beta[i];

        // Calculating the value of h (step size)
        h = (Xl -X0)/N;

        // Declaring vectors for the RHS, grid location (x) and exact solution of the equation
        vector<real> RHS, x, u_exact;

        // Generating the 1 Dimensional grid
        generate_1D_grid(x,X0,Xl,N);

        // Computing the RHS of the equation (WITHOUT implementing the boundary conditions)
        compute_rhs(RHS,x,beta,h);

        // Implementing the boundary conditions x = X0 and x = Xl
        RHS[0] = RHS[0] - (u0);
        RHS[N -2] = RHS[N-2] - (ul);

        // Declaring the rows and columns for the matix in the LHS
        int rows, cols;
        rows = N-1;
        cols = N-1;

        // Declaring the matrix
        vector<vector <real> > A(rows,vector<real> (cols,0.0));

        // Assigning the values of first row of A
        A[0][0] = -2.0;
        A[0][1] = 1.0;

        // Assigning the values of last row of A
        A[N-2][N-2] = -2.0;
        A[N-2][N-3] = 1.0;

        // This loop will iterate over all the rows of A except first and last row
        for(int i = 1;i < N-2;i++)
        {
            A[i][i-1] = 1.0;
            A[i][i] = -2.0;
            A[i][i+1] = 1.0;
        }

        // Declaring a vector (sol) to store the value of the numerical solution
        vector<real> sol(N-1);

        // Solving the equation using Thomas Algorithm
        Thomas_Algorithm_Matix_Solver(A,RHS,sol);

        // Assigning the value of u at x = X0, using the boundary condition
        sol.insert(sol.begin(),u0);

        // Assigning the value of u at x = Xl, using the boundary condition
        sol.push_back(ul);

        // Computing the exact solution
        compute_exact_sol(u_exact,x,beta);

        // Computing the L2 error
        L2_error.push_back(compute_L2_error(u_exact,sol));

        // Computing the maximum error
        max_error.push_back(compute_max_error(u_exact,sol));

        // Write the grid to a csv file
        write_to_file(x,"Q_1_b_Grid_Points_beta_is_"+to_string(beta)+".csv");

        // Write the numerical solution to a csv file
        write_to_file(sol,"Q_1_b_Numerical_Solution_beta_is_"+to_string(beta)+".csv");

        // Write the exact solution to a csv file
        write_to_file(u_exact,"Q_1_b_Exact_Solution_beta_is_"+to_string(beta)+".csv");
    }

    // The order of the L2 error will be same as that of the N in the input file "Values_of_N.txt"
    // Write the L2 to a csv file
    write_to_file(L2_error,"Q_1_b_L2_Error.csv");

    // The order of the maximum error will be same as that of the N in the input file "Values_of_N.txt"
    // Write the maximum error to a csv file
    write_to_file(max_error,"Q_1_b_Max_Error.csv");

    return 0;
}

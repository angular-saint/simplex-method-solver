/*
Author: Abtin Fooladkhah
Student of Computer Science at University of Tabriz
Project: Simplex Method Solver (C++)

Acknowledgements:
- Prof. Jaber Karimpoor
- Dr. Kheyri
for their guidance and support.

Description:
This program implements the Simplex algorithm for solving linear
programming problems. The user inputs the objective function and
constraints via console, and the program computes the optimal solution.
*/

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
using namespace std;

const double EPS = 1e-9;

// Print simplex tableau with Z row at the top

void printTableau(const vector<vector<double>> &T, const vector<string> &varNames, 
                  const vector<string> &basicVar, int iter){
    int rows = T.size();
    int cols = T[0].size();
    cout << "\n===== Iteration " << iter << " =====\n";
    cout << setw(6) << "BV";
    for(int j=0;j<cols-1;j++)
        cout << setw(8) << varNames[j];
    cout << setw(8) << "RHS\n";

    // Print Z row first
    cout << setw(6) << basicVar[rows-1];
    for(int j=0;j<cols;j++)
        cout << setw(8) << fixed << setprecision(3) << T[rows-1][j];
    cout << endl;

    // Print remaining rows
    for(int i=0;i<rows-1;i++){
        cout << setw(6) << basicVar[i];
        for(int j=0;j<cols;j++)
            cout << setw(8) << fixed << setprecision(3) << T[i][j];
        cout << endl;
    }
}

// Select pivot column
int selectPivotColumn(const vector<vector<double>> &T){
    int cols = T[0].size();
    int rows = T.size();
    int pivotCol = -1;
    double minVal = -EPS;
    for(int j=0;j<cols-1;j++){
        if(T[rows-1][j] < minVal){
            minVal = T[rows-1][j];
            pivotCol = j;
        }
    }
    return pivotCol;
}

// Select pivot row
int selectPivotRow(const vector<vector<double>> &T, int pivotCol){
    int rows = T.size();
    int cols = T[0].size();
    int pivotRow = -1;
    double minRatio = 1e18;
    for(int i=0;i<rows-1;i++){
        if(T[i][pivotCol] > EPS){
            double ratio = T[i][cols-1]/T[i][pivotCol];
            if(ratio < minRatio){
                minRatio = ratio;
                pivotRow = i;
            }
        }
    }
    return pivotRow;
}

// Simplex algorithm
void simplex(vector<vector<double>> &T, vector<string> &varNames, vector<string> &basicVar){
    int iter=0;
    int rows = T.size();
    int cols = T[0].size();
    while(true){
        iter++;
        printTableau(T,varNames,basicVar,iter);

        int pivotCol = selectPivotColumn(T);
        if(pivotCol==-1){
            cout << "\n✔ Optimal solution reached\n";
            cout << "Optimal Z = " << fixed << setprecision(3) << T[rows-1][cols-1] << endl;
            return;
        }

        int pivotRow = selectPivotRow(T,pivotCol);
        if(pivotRow==-1){
            cout << "\n❌ Unbounded solution\n";
            return;
        }

        cout << "\nEntering Variable: " << varNames[pivotCol];
        cout << "\nLeaving Variable: " << basicVar[pivotRow] << endl;

        basicVar[pivotRow] = varNames[pivotCol];

        double pivot = T[pivotRow][pivotCol];
        for(int j=0;j<cols;j++)
            T[pivotRow][j]/=pivot;
        for(int i=0;i<rows;i++){
            if(i==pivotRow) continue;
            double factor = T[i][pivotCol];
            for(int j=0;j<cols;j++)
                T[i][j]-=factor*T[pivotRow][j];
        }
    }
}

// Parse compact objective/constraint string into coefficients
vector<double> parseCoeffs(const string &line, int &n){
    vector<double> coeffs;
    int i=0;
    while(i<line.size()){
        double num=0;
        int sign=1;
        if(line[i]=='+'){ sign=1;i++; }
        else if(line[i]=='-'){ sign=-1;i++; }

        size_t start=i;
        while(i<line.size() && (isdigit(line[i])||line[i]=='.')) i++;
        if(start<i) num=stod(line.substr(start,i-start))*sign;
        else num=sign;

        if(i<line.size() && line[i]=='x'){
            i++;
            size_t startIdx=i;
            while(i<line.size() && isdigit(line[i])) i++;
            int idx=stoi(line.substr(startIdx,i-startIdx))-1;
            if(coeffs.size()<=idx) coeffs.resize(idx+1,0);
            coeffs[idx]=num;
        }
    }
    n=coeffs.size();
    return coeffs;
}

int main(){
    string line;
    vector<vector<double>> A;
    vector<double> b,c;
    vector<string> varNames,basicVar;

    // Objective function
    cout << "Enter objective function ( All input without spaces): ";
    getline(cin,line);
    int n;
    c=parseCoeffs(line,n);

    // Number of constraints
    int m;
    cout << "Enter number of constraints: ";
    cin >> m; cin.ignore();

    for(int i=0;i<m;i++){
        cout << "Constraint " << i+1 << " constraints: ";
        getline(cin,line);
        size_t pos = line.find("<=");
        double rhs = 0;
        string lhs=line;
        if(pos!=string::npos){
            rhs = stod(line.substr(pos+2));
            lhs=line.substr(0,pos);
        }
        int tmp;
        vector<double> coeffs = parseCoeffs(lhs,tmp);
        coeffs.resize(n,0);
        A.push_back(coeffs);
        b.push_back(rhs);
    }

    // Build initial simplex tableau
    int rows = m+1;
    int cols = n+m+1;
    vector<vector<double>> T(rows, vector<double>(cols,0));

    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            T[i][j]=A[i][j];
    for(int i=0;i<m;i++)
        T[i][n+i]=1;
    for(int i=0;i<m;i++)
        T[i][cols-1]=b[i];
    for(int j=0;j<n;j++)
        T[m][j]=-c[j];

    varNames.resize(n);
    for(int j=0;j<n;j++)
        varNames[j]="x"+to_string(j+1);
    for(int j=0;j<m;j++)
        varNames.push_back("s"+to_string(j+1));

    for(int i=0;i<m;i++)
        basicVar.push_back("s"+to_string(i+1));
    basicVar.push_back("Z");

    simplex(T,varNames,basicVar);
    return 0;
}


﻿#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip> 

using namespace std;


bool transformToDiagonallyDominant(vector<vector<double>>& matrix) {
    int num_Rows = matrix.size(); //Get the number of rows in the matrix

    for (int i = 0; i < num_Rows; i++) {
        // Find row that has the maximum absolute value in current column
        int max_row = i;
        double max_val = abs(matrix[i][i]);

        for (int j = i + 1; j < num_Rows; j++) {
            if (abs(matrix[j][i]) > max_val) {
                max_row = j;
                max_val = abs(matrix[j][i]);
            }
        }

        if (max_row != i) {
            // Swap rows i and max_row
            swap(matrix[i], matrix[max_row]);
        }

        // See if diagonal element is not dominant in its row
        double diagonal = abs(matrix[i][i]);
        double row_Sum = 0.0;

        for (int j = 0; j < num_Rows; j++) {
            if (i != j) {
                row_Sum += abs(matrix[i][j]);
            }
        }

        if (diagonal <= row_Sum) {
            return false; // Transformation was not successful
        }
    }

    return true; // Transformation was successful
}


// Function to print a matrix
void printMatrix(const vector<vector<double>>& matrix) {
    int num_Rows = matrix.size();
    int num_Cols = matrix[0].size();

    for (int i = 0; i < num_Rows; i++) {
        for (int j = 0; j < num_Cols; j++) {
            cout << matrix[i][j] << ' ';
        }
        cout << '\n';
    }
}

void gaussSeidel(vector<vector<double>>& A, vector<double>& b, vector<double>& start_approx_for_unknowns, double tolerance, string stopping_criterion = "MAE") {
    int num_Rows = A.size();

    // Check if Matrix is a diagonally dominant 
    bool isDiagonallyDominant = true;
    for (int i = 0; i < num_Rows; i++) {
        double diagonal = abs(A[i][i]);
        double row_Sum = 0.0;
        for (int j = 0; j < num_Rows; j++) {
            if (i != j) {
                row_Sum += abs(A[i][j]);
            }
        }
        if (diagonal <= row_Sum) {
            isDiagonallyDominant = false;
            break;
        }
    }

    if (!isDiagonallyDominant) {  // See if the matrix is not diagonally dominant
        cout << "Non-diagonally dominant Matrix Please Transform" << endl;
        cout << "Matrix Before Transformation:" << endl;
        printMatrix(A);  // Print original matrix
        cout << endl;

        if (!transformToDiagonallyDominant(A)) {  // Try to transform the matrix
            cout << "Matrix cannot be transformed into a diagonally dominant form" << endl;
            cout << "----------------------------------------------------------------------------------" << endl;
            return;  // Exit function b/c transformation is not successful
        }
        else {
            cout << "Diagonally dominant Matrix" << endl;
            cout << "Matrix After Transformation:" << endl;
            transformToDiagonallyDominant(A);  // Perform transformation matrix
            printMatrix(A); // Print transformed matix
        }

        cout << endl;
    }
    else {
        cout << "Diagonally dominant Matrix" << endl;  // else matrix is already diagonally dominant
    }


    vector<double> x(num_Rows, 0.0); // store the current approximations of unknowns 
    vector<double> old_x(num_Rows, 0.0); // store the previous approximations of unknowns

    for (int i = 0; i < num_Rows; i++) {
        b[i] = b[i] / A[i][i]; //here we transform biases b dividing them by a diagonal element from a corresponding row
        x[i] = start_approx_for_unknowns[i]; // starting value contains initial approximations; x contains current approximations
        old_x[i] = x[i]; // we store current approximations from x in old_x to be able to use them for finding an abs. approximate error
        for (int j = 0; j < num_Rows; j++) { // here we transform matrix A of our system dividing all elements in each row by a diagonal element from this row
            if (i != j) { //  we do this to avoid redundant operations in the main loop 
                A[i][j] = A[i][j] / A[i][i];
            }
        }
    }

    int iter = 0; // Initialize iteration counter
    double error = 10.0; // this is just to start the while loop below  error is intentionally assigned the large value

    while (error > tolerance) { // This is a main loop where our iterative process is being performed as long as error exceeds tolerance
        error = 0.0; // initialization of error to accumulate the error there 
        for (int i = 0; i < num_Rows; i++) {
            x[i] = b[i]; // we start updating current approximations new_x according to the main formula
            for (int j = 0; j < num_Rows; j++) { //  in this for loop, we finally calculate new approximations and put them in x according to the main formula
                if (i != j) {
                    x[i] = x[i] - A[i][j] * x[j]; //  In Gauss-Seidel method the latest available approximations of unknowns are used
                }
            }

            error += abs(x[i] - old_x[i]); // accumulation of approximate absolute errors calculated for each unknown
            old_x[i] = x[i]; // we store current approximations from x in old_x to be able to use them for finding an abs.approximate error
        }

        // Calculate the stopping error based on the chosen criterion
        double stopping_error;
        cout << "Iteration " << iter + 1 << ": ";
        if (stopping_criterion == "MAE") {
            stopping_error = error / num_Rows; // Mean Absolute Error
            cout << "Mean Absolute Error (MAE): " << fixed << setprecision(4) << stopping_error << endl;
        }
        else if (stopping_criterion == "RMSE") {
            stopping_error = sqrt(pow(error * error, 2) / num_Rows); // Root Mean Square Error
            cout << "Root Mean Square Error (RMSE): " << fixed << setprecision(4) << stopping_error << endl;
        }

        // Check if the stopping error is below the tolerance
        if (stopping_error <= tolerance) {
            break; // If the error is small enough, exit the loop
        }

        error = error / num_Rows; //  approximate mean absolute error is a mean of sum of approximate absolute errors calculated for each unknown

        iter++; // Increment iteration 
    }

    cout << endl;
    cout << "Solution Obtained:" << endl;
    for (int i = 0; i < num_Rows; i++) {
        cout << "x" << (i + 1) << " = " << fixed << setprecision(3) << x[i] << endl;
    }
    cout << "----------------------------------------------------------------------------------" << endl;
}


int main() {

    vector<vector<double>> A = { {3, 1, -4}, {-2, 3, 1}, {2, 0, 5} }; // Matrix A
    vector<double> b = { 7, -5, 10 }; // Part of Matrix A

    vector<vector<double>> B = { {1, -2, 4}, {8, -3, 2}, {-1, 10, 2} }; // Matrix B
    vector<double> c = { 6, 2, 4 }; // Part of Matrix B

    vector<double> startApproxForUnknowns = { 0, 0, 0 };
    double tolerance = 0.001;
    string stoppingCriterion1 = "MAE";
    string stoppingCriterion2 = "RMSE";

    cout << "------------------------------- Gauss-Seidel method ---------------------------------------" << endl;
    vector<vector<double>> A_copy = A;  // Create a copy of the original matrix
    vector<double> b_copy = b;

    // Apply the Gauss-Seidel method to matrix A and vector b with MAE and RMSE stopping criteria
    cout << "(STOPPING CRITERION: MAE for A)" << endl;
    gaussSeidel(A_copy, b_copy, startApproxForUnknowns, tolerance, stoppingCriterion1);
    cout << endl;

    A_copy = A;  // Reset the copy
    b_copy = b;  // Reset the copy
    cout << "(STOPPING CRITERION: RMSE for A)" << endl;
    gaussSeidel(A_copy, b_copy, startApproxForUnknowns, tolerance, stoppingCriterion2);

    cout << endl;
    // Apply the Gauss-Seidel method to matrix B and vector c with MAE and RMSE stopping criteria
    A_copy = B;  // Reset the copy
    b_copy = c;  // Reset the copy
    cout << "(STOPPING CRITERION: MAE for B)" << endl;
    gaussSeidel(A_copy, b_copy, startApproxForUnknowns, tolerance, stoppingCriterion1);
    cout << endl;

    A_copy = B;  // Reset the copy
    b_copy = c;  // Reset the copy
    cout << "(STOPPING CRITERION: RMSE for B)" << endl;
    gaussSeidel(A_copy, b_copy, startApproxForUnknowns, tolerance, stoppingCriterion2);

}



//------------------------------------------------------------------------------MATLAB -------------------------------------------------------------------

ADDED 12/11/23 BY OMAR

clc;
clear all;
A = [3, 1, -4; -2, 3, 1; 2, 0, 5]; % Matrix A
b = [7; -5; 10]; % Part of Matrix A

B = [1, -2, 4; 8, -3, 2; -1, 10, 2]; % Matrix B
c = [6; 2; 4]; % Part of Matrix B

startApproxForUnknowns = [0; 0; 0];
tolerance = 0.001;
stoppingCriterion1 = 'MAE';
stoppingCriterion2 = 'RMSE';

disp('------------------------------- Gauss-Seidel method ---------------------------------------');

disp('(STOPPING CRITERION: MAE for A)');
gaussSeidel(A, b, startApproxForUnknowns, tolerance, stoppingCriterion1);

disp('(STOPPING CRITERION: RMSE for A)');
gaussSeidel(A, b, startApproxForUnknowns, tolerance, stoppingCriterion2);

disp('(STOPPING CRITERION: MAE for B)');
gaussSeidel(B, c, startApproxForUnknowns, tolerance, stoppingCriterion1);

disp('(STOPPING CRITERION: RMSE for B)');
gaussSeidel(B, c, startApproxForUnknowns, tolerance, stoppingCriterion2);



function [success, matrix] = transformToDiagonallyDominant(matrix)
    num_Rows = size(matrix, 1);

    for i = 1:num_Rows
        % Find row that has the maximum absolute value in current column
        [max_val, max_row] = max(abs(matrix(i:num_Rows, i)));
        max_row = max_row + i - 1;
        max_val = max_val(1);

        if max_row ~= i
            % Swap rows i and max_row
            temp = matrix(i, :);
            matrix(i, :) = matrix(max_row, :);
            matrix(max_row, :) = temp;
        end

        % See if diagonal element is not dominant in its row
        diagonal = abs(matrix(i, i));
        row_Sum = sum(abs(matrix(i, :))) - diagonal;

        if diagonal <= row_Sum
            success = false; % Transformation was not successful
            return;
        end
    end

    success = true; % Transformation was successful
end

function printMatrix(matrix)
    [num_Rows, num_Cols] = size(matrix);

    for i = 1:num_Rows
        for j = 1:num_Cols
            fprintf('%f ', matrix(i, j));
        end
        fprintf('\n');
    end
end

function gaussSeidel(A, b, start_approx_for_unknowns, tolerance, stopping_criterion)
    num_Rows = size(A, 1);

    % Check if Matrix is diagonally dominant
    isDiagonallyDominant = true;
    for i = 1:num_Rows
        diagonal = abs(A(i, i));
        row_Sum = sum(abs(A(i, :))) - diagonal;
        if diagonal <= row_Sum
            isDiagonallyDominant = false;
            break;
        end
    end

    if ~isDiagonallyDominant
        disp('Non-diagonally dominant Matrix Please Transform');
        disp('Matrix Before Transformation:');
        printMatrix(A);

        [success, A] = transformToDiagonallyDominant(A);

        if ~success
            disp('Matrix cannot be transformed into a diagonally dominant form');
            disp('----------------------------------------------------------------------------------');
            return;
        else
            disp('Diagonally dominant Matrix');
            disp('Matrix After Transformation:');
            printMatrix(A);
        end
    else
        disp('Diagonally dominant Matrix');
    end

    x = zeros(num_Rows, 1); % Store the current approximations of unknowns
    old_x = zeros(num_Rows, 1); % Store the previous approximations of unknowns

    for i = 1:num_Rows
        b(i) = b(i) / A(i, i);
        x(i) = start_approx_for_unknowns(i);
        old_x(i) = x(i);
        A(i, :) = A(i, :) / A(i, i);
    end

    iter = 0;
    error = 10.0;

    while error > tolerance
        error = 0.0;

        for i = 1:num_Rows
            x(i) = b(i);
            for j = 1:num_Rows
                if i ~= j
                    x(i) = x(i) - A(i, j) * x(j);
                end
            end

            error = error + abs(x(i) - old_x(i));
            old_x(i) = x(i);
        end

        fprintf('Iteration %d: ', iter + 1);
        if strcmp(stopping_criterion, 'MAE')
            stopping_error = error / num_Rows;
            fprintf('Mean Absolute Error (MAE): %.4f\n', stopping_error);
        elseif strcmp(stopping_criterion, 'RMSE')
            stopping_error = sqrt((error * error) / num_Rows);
            fprintf('Root Mean Square Error (RMSE): %.4f\n', stopping_error);
        end

        if stopping_error <= tolerance
            break;
        end

        error = error / num_Rows;
        iter = iter + 1;
    end

    disp(' ');
    disp('Solution Obtained:');
    for i = 1:num_Rows
        fprintf('x%d = %.3f\n', i, x(i));
    end
    disp('----------------------------------------------------------------------------------');
end


//------------------------------------------------------------------------------MATLAB -------------------------------------------------------------------
ADDED 12/11/23
    BY Omar 

clc;
clear all;


A = [3, 1, -4, 7; -2, 3, 1, -5; 2, 0, 5, 10]; % Matrix A

B = [1, -2, 4, 6; 8, -3, 2, 2; -1, 10, 2, 4]; % Matrix B

startApproxForUnknowns = [0; 0; 0];
tolerance = 0.001;
stoppingCriterion1 = 'MAE';
stoppingCriterion2 = 'RMSE';

disp('------------------------------- Gauss-Seidel method ---------------------------------------');

disp('(STOPPING CRITERION: MAE for A)');
gaussSeidel(A, startApproxForUnknowns, tolerance, stoppingCriterion1);

disp('(STOPPING CRITERION: RMSE for A)');
gaussSeidel(A, startApproxForUnknowns, tolerance, stoppingCriterion2);

disp('(STOPPING CRITERION: MAE for B)');
gaussSeidel(B, startApproxForUnknowns, tolerance, stoppingCriterion1);

disp('(STOPPING CRITERION: RMSE for B)');
gaussSeidel(B, startApproxForUnknowns, tolerance, stoppingCriterion2);


function [success, matrix] = transformToDiagonallyDominant(matrix)
    num_Rows = size(matrix, 1);

    for i = 1:num_Rows
        % Find row that has the maximum absolute value in current column
        [max_val, max_row] = max(abs(matrix(i:num_Rows, i)));
        max_row = max_row + i - 1;
        max_val = max_val(1);

        if max_row ~= i
            % Swap rows i and max_row
            temp = matrix(i, :);
            matrix(i, :) = matrix(max_row, :);
            matrix(max_row, :) = temp;
        end

        % See if diagonal element is not dominant in its row
        diagonal = abs(matrix(i, i));
        row_Sum = sum(abs(matrix(i, :))) - diagonal;

        if diagonal <= row_Sum
            success = false; % Transformation was not successful
            return;
        end
    end

    success = true; % Transformation was successful
end

function printMatrix(matrix)
    [num_Rows, num_Cols] = size(matrix);

    for i = 1:num_Rows
        for j = 1:num_Cols
            fprintf('%f ', matrix(i, j));
        end
        fprintf('\n');
    end
end


function gaussSeidel(matrix, start_approx_for_unknowns, tolerance, stopping_criterion)
    num_Rows = size(matrix, 1);

    % Split the matrix into A and b
    A = matrix(:, 1:end-1);
    b = matrix(:, end);

    % Make a copy of A for transformations
    A_copy = A;

    % Check if Matrix A is diagonally dominant
    isDiagonallyDominant = true;
    for i = 1:num_Rows
        diagonal = abs(A_copy(i, i));
        row_Sum = sum(abs(A_copy(i, :))) - diagonal;
        if diagonal <= row_Sum
            isDiagonallyDominant = false;
            break;
        end
    end

    if ~isDiagonallyDominant
        disp('Non-diagonally dominant Matrix Please Transform');
        disp('Matrix Before Transformation:');
        printMatrix(matrix);

        [success, A_copy] = transformToDiagonallyDominant(A_copy);

        if ~success
            disp('Matrix cannot be transformed into a diagonally dominant form');
            disp('----------------------------------------------------------------------------------');
            return;
        else
            disp('Diagonally dominant Matrix');
            disp('Matrix After Transformation:');
            printMatrix([A_copy, b]);
        end
    else
        disp('Diagonally dominant Matrix');
    end

    x = zeros(num_Rows, 1); % Store the current approximations of unknowns
    old_x = zeros(num_Rows, 1); % Store the previous approximations of unknowns

    for i = 1:num_Rows
        b(i) = b(i) / A_copy(i, i);
        x(i) = start_approx_for_unknowns(i);
        old_x(i) = x(i);
        A_copy(i, :) = A_copy(i, :) / A_copy(i, i);
    end

    iter = 0;
    error = 10.0;

    while error > tolerance
        error = 0.0;

        for i = 1:num_Rows
            x(i) = b(i);
            for j = 1:num_Rows
                if i ~= j
                    x(i) = x(i) - A_copy(i, j) * x(j);
                end
            end

            error = error + abs(x(i) - old_x(i));
            old_x(i) = x(i);
        end

        fprintf('Iteration %d: ', iter + 1);
        if strcmp(stopping_criterion, 'MAE')
            stopping_error = error / num_Rows;
            fprintf('Mean Absolute Error (MAE): %.4f\n', stopping_error);
        elseif strcmp(stopping_criterion, 'RMSE')
            stopping_error = sqrt((error * error) / num_Rows);
            fprintf('Root Mean Square Error (RMSE): %.4f\n', stopping_error);
        end

        if stopping_error <= tolerance
            break;
        end

        error = error / num_Rows;
        iter = iter + 1;
    end

    disp(' ');
    disp('Solution Obtained:');
    for i = 1:num_Rows
        fprintf('x%d = %.3f\n', i, x(i));
    end
    disp('----------------------------------------------------------------------------------');
end


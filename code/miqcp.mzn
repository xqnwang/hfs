%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIQCP - Mixed Integer Quadratically Constrained Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data
include "miqcp-data.dzn";

QuadrFloat = true; % specify that the solver supports quadratic constraints

% Parameters
int: n; % number of all series
int: m; % number of bottom-level series

float: lambda; % Lagrange multiplier
int: M; % big-M

int: Nn; % number of rows of X
int: mn; % number of columns of X

set of int: SERIES = 1..n;
set of int: BSERIES = 1..m;
set of int: ROW = 1..Nn;
set of int: COLUMN = 1..mn;

array[ROW] of float: y; % input: y
array[ROW, COLUMN] of float: X; % input: X

% Decision variables
array[COLUMN] of var 0..2: g;
array[SERIES] of var 0..1: z;

% Assertion
constraint assert(mn = m * n, "Error: the input dimension is incorrect");

% Constraints
%% \sum_{i=0}^{m-1} (g_{j+in})^2 - M z_j <= 0 for j = 1,...,n
constraint forall(j in SERIES)(
                  sum (i in BSERIES) (g[j + (i - 1) * n] * g[j + (i - 1) * n]) - (M * z[j]) <= 0);

% Miminize loss function
var float: loss = 
  (1/Nn) * sum (r in ROW) ((y[r] - sum (c in COLUMN) (X[r, c] * g[c])) * (y[r] - sum (c in COLUMN) (X[r, c] * g[c])));
var float: penalty = 
  lambda * sum (j in SERIES) (z[j]);

% Solve item
solve minimize loss + penalty;

% Output item
output [ "z = \(z)\n" ] ++
       [ "G = \n"] ++
       [ show_float(6, 3, g[c]) ++
         if c mod n == 0 then "\n" else "" endif |
         c in COLUMN
       ] ++
       [ "loss = \(loss)\n"] ++
       [ "penalty = \(penalty)\n"] ++
       [ "obj = \(loss + penalty)\n"]

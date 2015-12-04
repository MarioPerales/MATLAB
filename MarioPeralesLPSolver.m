function[x, obj] = MarioPeralesLPSolver(A,b,c)

%
% This is an interactive program demonstrating the Simplex Method 
% via the revised simplex tableau.
% This method is the most efficient of all the 3 versions 
% Written by Ming Gu for Math 170
% March 2007

%%% Note, I (Mario Perales; December 3, 2015) have MODIFIED this original
%%% code. Any sections which are modified are accompinied by three %s 


%%% We must find the number of constraints and unknowns from the Matrix.
%%%
%%% For an mxn matrix, we have m constraints and n unknowns. %%%
[m, n] = size(A);

%%% setup LP for phase 1 
%%% 
%%% We must solve the CMP for phase 1: [A I][0] = b; z >= 0; sum(z) = min
%%%                                         [Z]   
A1 = [A eye(m)]; %%% Creates [A I]
c1 = [zeros(n,1);ones(m,1)]; %%% Created new cost vector
B =(n+1:n+m); %%% Our initial Basis set, B, is the index of the Identity matrix
x = zeros(1,n+m); %%% Our first feasible solution to this is x = b
x(B)= b;
T = A1(:,B)\[b eye(m)]; %%% Creating our simplified Tableau (without the criterion row!)
y = T(:,2:end)*c1(B); %%% Our cost
filler = x(B)*c1(B); 
T = [T;[filler,y']]; %%% Our ST with criterion row
[ITER,B,T,flg] = ILPRevisedSimplex(A1,B,c1,T); %%% Do Phase 1 Calculation to find BFS
x = zeros(n,1);
x(B) = T(1:end-1,1); %%% Update our x to have the BFS
i = 1;
while i <= length(B)
    CB(i) = c(B(i)); %%% We want to find our cost vector with respect to the basis of the BFS.
    i = i + 1;
end
j = 1;
while j <= length(B)
    t(j) = x(B(j)); %%% Same as above, but now for x.
    j = j + 1;
end
t = t'; 

u = T(1:end-1,2:end); %%% Getting the u component from our tableau
y = (CB * u); %%% We recalculate y because we want it in the right Basis
T(end,1) = CB * t; %%% Same for x
T(end,2:end) = y';

%%% We are done with phase 1. We plug this into Phase II.  %%%

% Starting Simplex Method

y = y';
f = ['Starting Phase II Simplex Iteration... '];
format short g;
disp(f);
disp('Initial Basis is');
disp(B');
disp('Displaying Initial solution x, c-A^T*y and their componentwise product');
disp([x c-A'*y x.*(c-A'*y)]);
simplex = 1;
ITER = 0;
obj = c'*x;
pause(2);

while (simplex == 1)
%
% determine the next s and r values.
%
   y        = T(end,2:end)';
   [zmin,s] = min(c-A'*y); 

   t        = T(1:end-1,2:end)*A(:,s);
   [flg,r] = Revisedgetr1(n,s,B,T,t);
   if (r < 1) %%% Switched the order here because we want to output "no lower bound" before infeasibility.
       disp('LP has no lower bound');
       simplex = 0;
       continue;
   end
   if (flg == 1)
       disp('LP is degenerate');
       simplex = 0;
       continue;
   end
   x   = zeros(n,1);
   x(B)= T(1:end-1,1);
   ITER = ITER + 1;
   f = ['Iteration ', num2str(ITER), ' Obj ', num2str(c'*x), '. Smallest component in c-A^T*y: ', ... 
         num2str(zmin), ' at s =', num2str(s), '. Component r = ', num2str(r), ' out of basis'];
   disp(f);
   obj1 = c'*x;

%
% check for convergence.
%
%%% Check for convergence at the end to output the right matrix with the 
%%% corresponding tableau.

   if (abs(zmin) < 1e-14)
       disp('Simplex Method has converged');
       simplex = 0;
       disp('Displaying Optimal Basis');
       disp(B');
       disp('Displaying Optimal solution x, c-A^T*y and their componentwise product');
       disp([x c-A'*y x.*(c-A'*y)]);
       obj = c'*x; %%% Last update of obj. %%%
       continue;
   end
   
%
% update the revised simplex tableau.
%
       
   if (flg == 1)
       disp('LP is degenerate');
       simplex = 0;
       continue;
   end
   [T,B1,flg]=RevisedSimplexTableau(B,r,s,t,zmin,T);  
   B   = B1;
   obj = obj1;
   disp('Current Basis is');
   disp(B');
   pause(1);
   opt = 2;
end
function [ITER,B,T,flg] = ILPRevisedSimplex(A,B,c,T)
%
% This is the Simplex Method using the Revised Simplex Tableau.
% We assume the Tableau T has been initialized on input.
% Written by Ming Gu for Math 170
% March 2007
%
% flg = 0: convergence
% flg = 1: LP unbounded below
% flg = 2: LP degenerate
%
simplex = 1;
ITER    = 0;
[m,n]   = size(A);
x       = zeros(n,1);
x(B)    = T(1:m,1);
flg     = 0;
obj     = c'*x;
ctol    = 1e-13;
while (simplex == 1)
%
% determine the next s and r values.
%
   y        = T(end,2:end)';
   [Y, I]   = sort(c-A'*y);
   mask     = find(Y<-ctol);
%
% check if the iteration has converged.
%
   if (length(mask)<1)
%%%       disp('Simplex Method has converged');
       simplex = 0;
       continue;
   end
   sflg = 1;
   for sk = 1:length(mask)
%
% loop over all possible s choices to reduce
% chance of non-convergence due to degeneracy.
%
       s    = I(sk);
       zmin = Y(sk);
       t    = T(1:end-1,2:end)*A(:,s);
       r    = Revisedgetr(n,s,B,T,t);
       if (r < 1)
          disp('LP has no lower bound');
          simplex = 0;
          flg     = 1;
          continue;
       end
       [Tnew,Bnew,flg] = RevisedSimplexTableau(B,r,s,t,zmin,T);      
       if (flg == 0 & Tnew(end,1)-obj < -ctol)
%
% success. update the revised simplex tableau.
%
           sflg = 0;
           B    = Bnew;
           T    = Tnew;
           obj  = T(end,1);
           ITER = ITER + 1;
           f = ['Iteration ', num2str(ITER), ' Obj ', ...
                num2str(obj), '. c-A^T*y componet = ', ... 
                num2str(zmin), ' at s =', num2str(s), '. Kicked out r = ', ...
                num2str(r)];
%%%           disp(f);
           break;
       end
    end
    if (sflg == 1)
%
%  none of the s choices was good enough to avoid degeneracy. Give up.
%
        disp('LP is degenerate');
        simplex = 0;
        flg     = 2;
%%%        disp(T(1:end-1,1));
        continue;
    end
end
function [flg, r] = Revisedgetr1(n,s,B,T,t)
%
% find the index to kick out
%
% On input: 
% B: current basis index
% T: current Revised Tableau
% t: current pivot column
% s: s-th  column is to JOIN the index.
% n: number of unknowns.
%

% On output: 
%
% r: B(r)-th column is to LEAVE the index.
%    r < 0 indicates unbounded LP.
%

x   = zeros(n,1);
flg = 1;
x(B)= T(1:end-1,1);
if (max(t)<n*eps)
    r = -1;
    return;
end
mask        = find(t>0);
[lamda, r]  = min(x(B(mask))./t(mask));
flg = lamda;
r           = mask(r);
function [Tnew,Bnew,flg] = RevisedSimplexTableau(B,r,s,t,zmin,T)
%
% This function updates a RevisedSimplexTableau
% 
% Written by Ming Gu for Math 170
% March 2007
%
% On input: 
% B: current basis index
% T: current RevisedTableau
% r: B(r)-th column is to LEAVE the index.
% t: Pivot column
% s: s-th  column is to JOIN the index.
% zmin: s-th component in c-A'*y
%

% On output: 
%
% flg:  flg == 0 indicates SUCCESS in updating,
%       flg == 1 indicates FAILURE in updating,
% Bnew: New basis index
% Tnew: New Tableau
%

%
% initialize flg.
%
flg     = 0;
%
% find dimensions of T.
%
[mt,nt] = size(T); 
%
% Set up Bnew
%
B     = B(:);
Bnew  = [B(1:r-1);s;B(r+1:mt-1)];
%
% Setup Tnew
%
Tnew         = zeros(mt,nt); 
if (t(r) == 0)
%
% This is indication of degeneracy. Quit.
%
    flg = 1;
    return;
end
%
% This is the normal case. Proceed. 
%
Temp              = T(r,:)/t(r);
Tnew(1:mt-1,:)    = T(1:mt-1,:) - t*Temp;
Tnew(r,:)         = Temp;
theta0            = zmin/t(r);
Tnew(mt,:)        = T(mt,:) + theta0*T(r,:);
function r = Revisedgetr(n,s,B,T,t)
%
% find the index to kick out
%
% On input: 
% B: current basis index
% T: current Revised Tableau
% t: current pivot column
% s: s-th  column is to JOIN the index.
% n: number of unknowns.
%

% On output: 
%
% r: B(r)-th column is to LEAVE the index.
%    r < 0 indicates unbounded LP.
%

x   = zeros(n,1);
x(B)= T(1:end-1,1);
if (max(t)<n*eps)
    r = -1;
    return;
end
mask        = find(t>0);
[lamda, r]  = min(x(B(mask))./t(mask));
r           = mask(r);
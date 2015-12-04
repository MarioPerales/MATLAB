% This script is designed to solve the 1-norm regression problem.
%
% To solve this linear program, this program employs methods discussed
% discussed on the MA170 class website, taught by Ming Gu at UC Berkeley.
% Most notably, we use a variation of the simplified tableau method.
%
% INPUT(S):
%
% NOTE: Since this is a script, you do not need inputs BUT the script
% reads a textfile which MUST be labeled infile.txt. ***
%
% infile.txt: A text file. The text file must be formatted exactly as
% followed. The first line of the matrix is n m, where n and m are the
% dimensions of your A matrix, respectively. It is then followed by n x m
% lines that look like i j aij. i and j indicate the index of the A matrix
% and aij corresponds to the value of the matrix at A(i,j). Finally, the
% last n lines take the form i bi. This outlines the content and length of
% the b matrix in the 1-norm regression problem.
%
% OUTPUT(S):
%
% outfile.txt: A text file. It displays the time it took to complete the
% script, the optimal solution, optimal dual solution, and optimal object
% value (in this case, A*x). It also indicates the basis in x for which 
% these values correspond to.
%
% Written by: Mario Perales
% December 3, 2015


% We first want to be able to read our infile text. %

tic;
fid = fopen('infile.txt');

% Note that fgetl gets a line from the file and discards discards it afterwards. % 

tline = fgetl(fid);
num = str2num(tline);

% First line gives us the dimensions of the matrix. %

n = num(1);
m = num(2);

% Initiatilize our matrix A to fit our dimensions and fill it in later... %

A = zeros(n,m);

for k=1:(n*m)
    iter_tline = fgetl(fid);
    iter_num = str2num(iter_tline);
    [a, c] = size(iter_num);
    
    % Good o' CS61B error throwing. Just checking if everything is fine. % 
    if(c ~= 3)
        err = MException('',horzcat('It appears line ', num2str(k),' does not have exactly three values!'));
        throw(err);
    end;
    
    % Just reading off the data % 
    i = iter_num(1);
    j = iter_num(2);
    
    A(i,j) = iter_num(3);
end

% Now to initialize our b matrix(i.e. vector)! %
b = zeros(n,1); % Should it be zeroes(1,n)? %

for k=1:n;
    b_iter_tline = fgetl(fid);
    b_iter_num = str2num(b_iter_tline); %#ok<*ST2NM>
    [a, c] = size(b_iter_num);
    if(c ~= 2)
        err = MException('',horzcat('It appears line ', num2str(k),' does not have exactly two values!'));
        throw(err);
    end
    
    i = b_iter_num(1);
    b(i) = b_iter_num(2); 
end;

% Now to start the actual program. %

%  ---------------  %
% | START PHASE 1 | %
%  ---------------  %

% We first want to find an initial index set B for Phase I

B = 1:m;
Bc = setdiff(1:n,B);
t = zeros(n,1);

% Check for degeneracy cases first. % 
% Check if matrix is invertible. %
% Have no idea how to check the second case...%

if(rank(A(B,:)) ~= min(size(A(B,:))))
    err = MException('', 'A(B,:) is non-invertible. Problem is degenerate.');
    throw(err);
end;

% Calculating M and M^-1 %

M = A(B,:);
M_inverse = inv(M);

x = M*b(B);
h = A*x - b;

y(Bc) = sign(h(Bc));
y(B) = -(M_inverse)'*A(Bc,:)'*y(Bc);
obj = abs(A(Bc,:)*x - b(Bc));

%  -------------  %
% | END PHASE I | %
%  -------------  %

while any(abs(y(B)) > 1)
    
    %  ----------------  %
    % | START PHASE II | %
    %  ----------------  %

    s = find(abs(y) > 1,1);
    j = find(B == s);
    e = zeros(m,1);
    e(j) = 1;
    
    t(Bc) = -(sign(y(s)))*(y(Bc)).*(A(Bc,:)*A(B,:)^(-1)*e);
    
    % Need to find r...we employ a method similar to Revisedgetr on the
    % class website. %
    
    mask = find(t(Bc) > 0);
    [lambda, r] = min(abs(h(mask))./(t(Bc(mask))));
    r = Bc(mask(r));
    
    % Time to update B %r
    % B = B\{s}(UNION){r} %
    B = [setdiff(B,s), r];
    Bc = setdiff(1:n,B);
    % Degeneracy check. %
    if(rank(A(B,:)) ~= min(size(A(B,:))))
        err = MException('', 'A(B,:) is non-invertible with new index set. Problem is degenerate.');
        throw(err);
    end;

    M = A(B,:);
    M_inverse = inv(M);

    x = M*b(B);
    h = A*x - b;

    y(Bc) = sign(h(Bc));
    y(B) = -(M_inverse)'*A(Bc,:)'*y(Bc);
    obj = abs(A(Bc,:)*x - b(Bc));
    
end;

%  --------------  %
% | END PHASE II | %
%  --------------  %

time = toc;

% Need to create a text file for this data...%
fid = fopen('outfile.txt','w');

% Display the results! %
fprintf(fid, '%s\r', horzcat('SUCCESS! It only took ', num2str(time), ' seconds to compute the answer!'));
fprintf(fid, '%s\r\n\n\n', 'Here are the results!');
fprintf(fid, '%s\r', 'The optimal primal solution is...');
fprintf(fid, '%f\r', x);
fprintf(fid, '%s\r', 'The optimal dual solution is...');
fprintf(fid, '%f\r', y(B));
fprintf(fid, '%s\r', 'The optimal objective value is...');
fprintf(fid, '%f\r', obj);
fprintf(fid, '%s\r', 'Index set is...');
fprintf(fid, '%f\r', B);


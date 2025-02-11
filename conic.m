clear, clc

%Yalmip set
yalmip('clear');
clear;
f=1;
nt=10; % number of expansion coefficients in \theta
nx=5; % number of expansion coefficients in \xi
nw=10; % number of expansion coefficients in w
xs=[];
ys=[];
%kvals = 1:20;

Pe_vals = logspace(log10(1e-1), log10(1), 100);
for i = 1:length(Pe_vals)
   Pe = Pe_vals(i);
   kvals = 1:10;
   output_lr = asd_function(nt, nx, nw, Pe, kvals, f);
   output_hr = asd_function(nt+20, nx+10, nw+20, Pe, kvals, f);
   % Check convergenge
    while true
        if abs(output_hr.LB-output_lr.LB) <= 0.001 * output_lr.LB
            fname = sprintf('results_Pe_%8.6e.mat', Pe);
            save(fname, "output_hr")
            break
        else
            output_lr = output_hr;
            nt=nt+20;
            nx=nx+10;
            nw=nw+20;
            output_hr = asd_function(nt+20, nx+10, nw+20, Pe, kvals, f);
        end
    end
end

Pe_vals = logspace(log10(1), log10(50), 100);
for i = 1:length(Pe_vals)
    Pe = Pe_vals(i);
    kvals = 1:100  ;
    output_lr = asd_function(nt, nx, nw, Pe, kvals, f);
    output_hr = asd_function(nt+20, nx+10, nw+20, Pe, kvals, f);
    % Check convergenge
    while true
        if abs(output_hr.LB-output_lr.LB) <= 0.001 * output_lr.LB
            fname = sprintf('results_Pe_%8.6e.mat', Pe);
            save(fname, "output_hr")
            break
        else
            output_lr = output_hr;
            nt=nt+20;
            nx=nx+10;
            nw=nw+20;
            output_hr = asd_function(nt+20, nx+10, nw+20, Pe, kvals, f);
        end
    end
end

Pe_vals = logspace(log10(50), log10(100), 100);
for i = 1:length(Pe_vals)
    Pe = Pe_vals(i);
    kvals = 1:130;
    output_lr = asd_function(nt, nx, nw, Pe, kvals, f);
    output_hr = asd_function(nt+20, nx+10, nw+20, Pe, kvals, f);
    % Check convergenge
    while true
        if abs(output_hr.LB-output_lr.LB) <= 0.001 * output_lr.LB
            fname = sprintf('results_Pe_%8.6e.mat', Pe);
            save(fname, "output_hr")
            break
        else
            output_lr = output_hr;
            nt=nt+20;
            nx=nx+10;
            nw=nw+20;
            output_hr = asd_function(nt+20, nx+10, nw+20, Pe, kvals, f);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%
function [output] = asd_function(nt, nx, nw, Pe, kvals, f)
% assert(all(size(X)==[1,n+1]))

yalmip clear
a = sdpvar(1,1);
b = sdpvar(1,1);
v = sdpvar(1,1);
X = sdpvar(nx-1,1);

% Define chebyshev polynomials and their derivatives
n = max([nt, nx, nw]);
T = chebpoly(0:n-1);
Tprime = diff(T);

%Define constraints
Constraints = [];


%Three matrixs for boundary conditions
%First one, for the boundary conditions of \theta
Bt = zeros(nt, nt-1);

Bt(1,:) = -T(1,2:nt)./T(1,1);

for i = 2:nt
    if i-1 <= nt-1
        Bt(i,i-1) = 1;
    end
end

%Second one, for the boundary conditions of \xi
Bx = zeros(nx, nx-1);

Bx(1,:) = -T(1,2:nx)./T(1,1);

for i = 2:nx
    if i-1 <= nx-1
        Bx(i,i-1) = 1;
    end
end

%Third one, for the boundary conditions of w
Bw = zeros(nw, nw-2);
Bw_1=T(-1,2).*T(1,3:nw)-T(1,2).*T(-1,3:nw);
Bw_2=T(1,1).*T(-1,3:nw)-T(-1,1).*T(1,3:nw);
Bw_3=T(-1,1).*T(1,2)-T(1,1).*T(-1,2);
Bw(1,:) = Bw_1/Bw_3;
Bw(2,:) = Bw_2/Bw_3;

for i = 3:nw
    if i-2 <= nw-2
        Bw(i,i-2) = 1;
    end
end

%Constraint 1
E11=Tprime(:,1:nx)'*Tprime(:,1:nx);
E12=-sum(chebfun('(x+1)',[-1 1])*Tprime(:,1:nx));
E22= sum(chebfun('(x+1)',[-1 1])^2);
E= [Bx'*E11*Bx, Bx'*E12';E12*Bx, E22];
E = 0.5 * (E + E');
[V, D] = eig(E);
R = sqrt(D) * V';
Constraints = [Constraints, rcone(R * [X; b], b + 1, v)];

%Constraint 2
% Assemble the block that depends on X
% We need to integrate phi_i*phi_j*phi_l' for every l = 1:length(X)
C = zeros(nt,nw);
X_with_0 = Bx * X;
for k = 1:nx
    Q = T(:,1:nt)' * ( Tprime(:,k) .* T(:,1:nw));
    Q(abs(Q)<1e-12) = 0;
    C = C + Q .* X_with_0(k);
end

% Assemble everything for all values of k
%S = cell(length(n),1);
for i = 1:length(kvals)
    k = kvals(i);
    B11 = Tprime(:,1:nt)'*Tprime(:,1:nt) + k^2 .* T(:,1:nt)'*T(:,1:nt);
    B12 = ( 0.5.*Pe.*k ) .* C;
    B22 = Tprime(:,1:nw)'*Tprime(:,1:nw) + k^2 .* T(:,1:nw)'*T(:,1:nw);
    M = [0.5*(1+b)*B11, 0.5*B12; 0.5*B12.', 0.5*a*B22];
    %Extend M because of the boundary condtion
    N = zeros(nt+nw,nt+nw-3);
    N(1:nt, 1:nt-1) = Bt;
    N(nt+1:nt+nw, nt:nt+nw-3)=Bw;
    S{k}=N'*M*N;
    S{k} = 0.5 * (S{k} + S{k}');
    Constraints = [Constraints, S{k}>=0];
end

% Boundary conditions on `xi
% Constraints = [Constraints, Tprime(-1)*X_with_0==0];

% Define the objective
Objective = (1/2)*f*sum(T(:,1:nx))*X_with_0 - a - (1/4)*v;

% Solve the SDP using MOSEK
options = sdpsettings('solver', 'mosek', 'verbose', 1);
sol = optimize(Constraints, -Objective, options);
pres = check(Constraints);

% Return all variables that we need
output.Pe = Pe;
output.LB = value(Objective);
output.a = value(a);
output.b = value(b);
output.xi = value(X_with_0);
output.sol = sol;
output.nx = nx;
output.nt = nt;
output.nw = nw;
output.kvals = kvals;
output.pres = pres;
end
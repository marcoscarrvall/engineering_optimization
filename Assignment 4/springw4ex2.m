% Two variable valve spring problem - Exercise 4.1
% 1. Visualization of UNCONSTRAINED spring stiffnes and frequency 
% optimization problem
% 2. Computation of steepest descent search direction
% 3. Line search using this search direction: hand controlled optimization cycles,
% including visualization.

% Initialization
clf, hold off, clear
format long

fig = figure(1);
set(fig, 'Color', 'white');

% 1. Problem visualization
% Constant parameter values
springparams1;
w=1;
ktarget=10000; 
frtarget=300;

% Matrix of output values for combinations of design variables D and d 
D = [0.020:0.0005:0.040];
d = [0.002:0.00004:0.005];
for j=1:1:length(d)
  for i=1:1:length(D)
%   Analysis of valve spring.
    [svol,smass,bvol,matc,manc,Lmin,L2,k,F1,F2,Tau1,Tau2,freq1]=...
    springanalysis1(D(i),d(j),L0,L1,n,E,G,rho,Dv,h,p1,p2,nm,ncamfac,nne,matp,bldp);
 	 % Scaled objective function
     fobj(j,i) = ((k-ktarget)/ktarget)^2 + w*((freq1-frtarget)/frtarget)^2; 
     stiffness(j,i) = k;
     freq(j,i) = freq1;
  end
end

% Contour plot of scaled spring optimization problem
%contour(D, d, fobj,[0:0.05:0.2 0.2:0.1:0.5 0.5:0.5:2 2:5:100])
cc = [0.01 0.02 0.05];
contour(D, d, fobj,[cc 10*cc 100*cc 1000*cc 10000*cc 100000*cc 1000000*cc])
ax = gca;
set(ax, 'Color', 'white');       % Background of the plot area
set(ax, 'XColor', [0.1 0.1 0.1]); % Dark grey/black X-axis
set(ax, 'YColor', [0.1 0.1 0.1]); % Dark grey/black Y-axis
set(ax, 'GridColor', [0.8 0.8 0.8]); % Soft light-grey grid lines
set(ax, 'GridAlpha', 0.5);       % Transparency of grid

xlabel('Coil diameter D (m)'), ylabel('Wire diameter d (m)'), ...
title('Figure 1     Spring stiffness and frequency optimization problem (w = 1.0)')
hold on
contour(D,d,stiffness,[10000 10000], "r")
contour(D,d,freq,[300 300], "g")
grid

%end problem visualization


% 1. Initialization for Automation
xq = [0.022  0.0035];   % Initial point [cite: 45]
max_cycles = 100;       % Maximum number of optimization cycles [cite: 55]
tolerance = 1.0e-6;     % Convergence criterion 
cycle = 0;
f_prev = 1e10;          % Initialize with a large value
converged = false;


% Loop over optimization cycle:
while (cycle < max_cycles) && ~converged
   cycle = cycle +1
   % Plot marker in current design point:
   plot(xq(1),xq(2),'o');
   
	% Forward finite diffence gradients of objective function and constraints
	hi=1e-8;
	alpha=0.0;
	sq=[0 0];
	% Objective function in point xq
	fx = springobjw4(alpha,xq,sq,ktarget,frtarget,w);
	% Perturbated objective function values:
	fx1plush = springobjw4(alpha,[xq(1)+hi, xq(2)],sq,ktarget,frtarget,w);
	fx2plush = springobjw4(alpha,[xq(1), xq(2)+hi],sq,ktarget,frtarget,w);
	% Objective function derivatives:
	dfdx1 = (fx1plush - fx)/hi;
	dfdx2 = (fx2plush - fx)/hi;
	 % Gradient vector:
  	df = [dfdx1 dfdx2];
	% Steepest descent search direction:
	sq = -df;
    % Normalize to unity length
    sq = sq / norm(sq);
    scale = 0.01; 
    sq = sq * scale;
    sq
   
	% Setting of options:
	options = optimset('tolx',1.0e-8,'MaxFunEvals',50);

   alphamax = 10.0;

	%Line search (note the lower and upper bound of alfhaq):
   [alphaq, f_current] = ...
      fminbnd('springobjw4', 0, alphamax, options, xq, sq, ktarget, frtarget, w);

   % --- Convergence Check --- 
   if abs(f_current - f_prev) < tolerance
      converged = true;
   end

	% Update values for next iteration
   f_prev = f_current;
   xnew = xq + alphaq * sq;   % Update design point [cite: 33]
    
   % Visualization of path
   plot([xq(1) xnew(1)], [xq(2) xnew(2)], xnew(1), xnew(2), 'o');
   xq = xnew;
    
   % Optional: Print progress to command window
   fprintf('Cycle %d: f = %e\n', cycle, f_current);

end;  % end optimization cycle (while-loop)

%end 

fprintf('\nOptimization finished.\n');
fprintf('Number of optimization cycles necessary for convergence: %d\n', cycle);
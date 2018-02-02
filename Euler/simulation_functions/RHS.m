function du_dt = RHS(a,b,func)
%
% Evaluates the RHS of an ODE func for the each relevant index
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%     a  =  set of positive indices [2:N]
%
%     b  =  set of negative indices [2*M:-1:2*M-N+1]
%
%  func  =  the RHS function

% find how many modes there are being retained and construct output array
N = length(a)+1;
du_dt = zeros(N,N,N,3,4);

% fill in requisite entries (the conjugates do not need to be explicitly
% computed, which cuts simulation time in half
du_dt(a,a,a,:,1) = func(a,a,a); % k = [+,+,+] (conjugate is [-,-,-])
du_dt(a,a,a,:,2) = func(a,a,b); % k = [+,+,-] (conjugate is [-,-,+])
du_dt(a,a,a,:,3) = func(a,b,a); % k = [+,-,+] (conjugate is [-,+,-])
du_dt(a,a,a,:,4) = func(b,a,a); % k = [-,+,+] (conjugate is [+,-,-])

du_dt(1,a,a,:,1) = func(1,a,a); % k = [0,+,+] (conjugate is [0,-,-])
du_dt(a,1,a,:,1) = func(a,1,a); % k = [+,0,+] (conjugate is [-,0,-])
du_dt(a,a,1,:,1) = func(a,a,1); % k = [+,+,0] (conjugate is [-,-,0])

du_dt(1,a,a,:,2) = func(1,a,b); % k = [0,+,-] (conjugate is [0,-,+])
du_dt(a,1,a,:,2) = func(a,1,b); % k = [+,0,-] (conjugate is [-,0,+])
du_dt(a,a,1,:,2) = func(a,b,1); % k = [+,-,0] (conjugate is [-,+,0])

du_dt(1,1,a,:,1) = func(1,1,a); % k = [0,0,+] (conjugate is [0,0,-])
du_dt(1,a,1,:,1) = func(1,a,1); % k = [0,+,0] (conjugate is [0,-,0])
du_dt(a,1,1,:,1) = func(a,1,1); % k = [+,0,0] (conjugate is [-,0,0])

du_dt(1,1,1,:,1) = func(1,1,1); % k = [0,0,0]
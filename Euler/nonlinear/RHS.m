function du_dt = RHS(a,b,func)

N = a(end);
du_dt = zeros(N,N,N,3,4);

du_dt(a,a,a,:,1) = func(a,a,a);
du_dt(a,a,a,:,2) = func(a,a,b);
du_dt(a,a,a,:,3) = func(a,b,a);
du_dt(a,a,a,:,4) = func(b,a,a);

du_dt(1,a,a,:,1) = func(1,a,a);
du_dt(a,1,a,:,1) = func(a,1,a);
du_dt(a,a,1,:,1) = func(a,a,1);

du_dt(1,a,a,:,2) = func(1,a,b);
du_dt(a,1,a,:,2) = func(a,1,b);

du_dt(1,a,a,:,3) = func(1,b,a);
du_dt(a,a,1,:,3) = func(a,b,1);

du_dt(a,1,a,:,4) = func(b,1,a);
du_dt(a,a,1,:,4) = func(b,a,1);

du_dt(1,1,a,:,1) = func(1,1,a);
du_dt(1,a,1,:,1) = func(1,a,1);
du_dt(a,1,1,:,1) = func(a,1,1);

du_dt(1,1,1,:,1) = func(1,1,1);
%This program will solve for the transition path of capital
%from the initial level (k0) to the steady state level (kss)
%for the standard growth model. There is no endogenous labor supply decision. 
%The production function is Cobb-Douglas, with a capital share of alpha.
%The utility function is (c^(1-sigma))/(1-sigma).
%The program assumes that the initial level of capital is below the steady 
%state value. Specifically, it is assumed that k0=lambda*kss.
%T is the number of periods for which the path is computed. The program also allows for
%a tax of tau on income from capital that is rebated lumpsum to the consumer

%Parameter values
delta=.08;alpha=.3;beta=.96;sigma=.999;lambda=.1;T=200;
A=1;tauk=0;

%Steady state capital stock and initial condition

kss=((1/beta-(1-delta))/A/alpha/(1-tauk))^(1/(alpha-1));
k0=lambda*kss;


%Solution for path k(t)

k(1)=k0;


for t=2:T
    
%compute hypothetical path kt for various choices of k(t), given k(t-1)
kcon(1)=k(t-1);
%Step 1:  Set boundaries for interval that k(t) must lie in. 
kmin=k(t-1);kmax=kss;
 %Step 2: Pick a point in the interval as a hypothetical choice for k(t)  
   while abs(kmax-kmin)>.00000015*kss
   kn=.5*(kmin+kmax);
   kcon(2)=kn;
   
 %Step 3: Determine the rest of the path (called kcon for k conditional on a choice of k(1)) given the choice for k(t)
 %stop is simply an indicator to tell us if we need to continue to the next
 %element of the kcon series. i is the index of the kcon series, with i=1 corresponding
 %to k(t-1) and i=2 corresponding to k(t).
   stop=0;
   i=2;
   while stop < 1
      i=i+1;
   %implied value for kcon(i)
         kcon(i)=A*kcon(i-1)^alpha+(1-delta)*kcon(i-1)-...
         (beta*(A*(1-tauk)*alpha*kcon(i-1)^(alpha-1)+(1-delta)))^(1/sigma)*...
        (A*kcon(i-2)^alpha+(1-delta)*kcon(i-2)-kcon(i-1));
    %check to see if kcon is going to zero 
         if kcon(i)<=kcon(i-1), kmin=kn;stop=1;else,kcon(i)=kcon(i);end
    %check to see if kcon is going beyond kss
         if kcon(i)>kss, kmax=kn;stop=1;else,kcon(i)=kcon(i);end
  end
  end
  % We have now determined (within the limits of our approximation) the next value of k(t)
   k(t)=kcon(2);
   end
 
   %We now know the whole sequence k(t). We can also determine the 
   %sequences for consumption and utility;
   
   c(1:(T-1))=A*k(1:(T-1)).^alpha+(1-delta)*k(1:(T-1))-k(2:T);c(T)=c(T-1);
   uv=(c.^(1-sigma))/(1-sigma);betav(1)=1;for j=2:T, betav(j)=beta*betav(j-1);end
   utot=sum(betav.*uv);
   srate=1-c./(A*k.^alpha);
   
   
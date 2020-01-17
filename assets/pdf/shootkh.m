%This program will solve for the transition path of capital
%from the initial level (k0) to the steady state level (kss)
%for the standard growth model with endogenous labor supply decision. 
%The production function is Cobb-Douglas, with a capital share of alpha.
%The utility function is (c^(1-sigma)-1)/(1-sigma)+a*[{(1-h)^(1-gamma)}-1]/(1-gamma).
%The program assumes that the initial level of capital is below the steady 
%state value. Specifically, it is assumed that k0=lambda*kss.
%T is the number of periods for which the path is computed.

%Parameter values
delta=.08;alpha=.3;beta=.96;sigma=0.999;lambda=.1;T=500;
A=1;gamma=.9999;a=2; taul=0;

%Solve for steady state and initial capital

%first solve for steady state ratio of k to h:
khss=((1/beta-(1-delta))/A/alpha)^(1/(alpha-1));
%now determine optimal h
hmin=0;hmax=1;
   while (hmax-hmin) > .00001
       hn=.5*(hmin+hmax);
       cn=hn*(A*khss^alpha-delta*khss);
       foc=a*(1-hn)^(-gamma)/cn^(-sigma)-(1-taul)*A*(1-alpha)*(khss)^alpha;
       if foc > 0, hmax=hn; else, hmin=hn;end
   end
   hss=hn;css=cn;kss=khss*hss;
k0=lambda*kss;
%k0=res0(1);
%Solution for path k(t)

k(1)=k0;


for t=2:T
    
%compute hypothetical path kt for various choices of k(t), given k(t-1)
kt(1)=k(t-1);
%Step 1:  Set boundaries for interval that k(t) must lie in. 
kmin=k(t-1);kmax=kss;
 %Step 2: Pick a point in the interval as a hypothetical choice for k(t)  
   while abs(kmax-kmin)>.00000015*kss
   kn=.5*(kmin+kmax);
   kt(2)=kn;
   %solve for ht(1), ct(1)
   hmin=0;hmax=1;
   while (hmax-hmin) > .00001
       hn=.5*(hmin+hmax);
       cn=A*kt(1)^alpha*hn^(1-alpha)+(1-delta)*kt(1)-kt(2);
       foc=a*(1-hn)^(-gamma)/cn^(-sigma)-(1-taul)*A*(1-alpha)*(kt(1)/hn)^alpha;
       if foc > 0, hmax=hn; else, hmin=hn;end
   end
   ht(1)=hn;
   ct(1)=A*kt(1)^alpha*hn^(1-alpha)+(1-delta)*kt(1)-kt(2);
   
 %Step 3: Determine the rest of the path (called kt, ct and ht) given the 
 %choice for k(t)
 %stop is simply an indicator to tell us if we need to continue to the next
 %element of the kt series. i is the index of the kt series, with i=1 corresponding
 %to k(t-1) and i=2 corresponding to k(t).
   stop=0;
   i=2;
   while stop < 1
      i=i+1;
   %implied value for kt(i)
   %guess a value for ht(i-1), solve for ct(i-1) and hence kt(i)
   
   hmin=0;hmax=1;
   while (hmax-hmin) > .00001
       hn=.5*(hmin+hmax);
       cn=(ct(i-2)^(-sigma)/beta/(alpha*(kt(i-1)/hn)^(alpha-1)+1-delta))^(-1/sigma);
       foc=a*(1-hn)^(-gamma)/cn^(-sigma)-(1-taul)*A*(1-alpha)*(kt(i-1)/hn)^alpha;
       if foc > 0, hmax=hn; else, hmin=hn;end
   end
   ht(i-1)=hn;ct(i-1)=cn;kt(i)=A*kt(i-1)^alpha*ht(i-1)^(1-alpha)+(1-delta)*kt(i-1)-cn;
   
    %check to see if kt is going to zero 
         if kt(i)<=kt(i-1), kmin=kn;stop=1;else,kt(i)=kt(i);end
    %check to see if kt is going to infinity
         if kt(i)>kss, kmax=kn;stop=1;else,kt(i)=kt(i);end
  end
  end
  % We have now determined (within the limits of our approximation) the next value of k(t)
   k(t)=kt(2);c(t-1)=ct(1);h(t-1)=ht(1);
   end
 
   %We now know the whole sequence k(t) and h(t). We can also determine the 
   %sequences for consumption and utility;
   kh=k;
   ca(1:(T-1))=A*k(1:(T-1)).^alpha.*h(1:(T-1)).^(1-alpha)+(1-delta)*k(1:(T-1))-k(2:T);ca(T)=c(T-1);c(T)=c(T-1);h(T)=h(T-1);
   uv=(c.^(1-sigma))/(1-sigma)+a*((1-h).^(1-gamma)-1)/(1-gamma);betav(1)=1;for j=2:T, betav(j)=beta*betav(j-1);end
   utot=sum(betav.*uv);
   ucss=css^(1-sigma)/(1-sigma)/(1-beta);uhss=a*(1-hss)^(1-gamma)/(1-gamma)/(1-beta);
   
   
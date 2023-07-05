$TITLE LAKE MENDOTA PMP MODEL

* Lake Mendota catchment agricultural land-management model-policy simulation
* Weizhe Weng & Kelly M. Cobourn, July 2023
* email: wweng@ufl.edu, kellyc13@vt.edu

***************************************************************************
*****************************  PRELIMINARIES  *****************************
***************************************************************************


**** assume a "calibration year" based on the average of year 2003 to 2014

sets
i        crops /cfc,cfs,sfc,c1fa,c2fa,c3fa,a1fc,a2fc,a3fc/
k        commodities /corn,soy,alfa/
cc(k,i)  commodity-crop mapping /corn.(cfc,cfs,c1fa,c2fa,c3fa),
         soy.(sfc),alfa.(a1fc,a2fc,a3fc)/
r        rotations /ccorn,cornsoy,cornalfa3/
cr(r,i)  rotation-crop mapping /ccorn.(cfc),cornsoy.(cfs,sfc),
         cornalfa3.(c1fa,c2fa,c3fa,a1fc,a2fc,a3fc)/
l        inputs /land,N/
yp       estimated yield function parameters /alpha,beta0,beta2,gamma,rho,beta4/
t        simulation year /2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014/
alias (i,j)

**************************  EXOGENOUS PARAMETERS  **************************

parameters
*Commodity prices are averaged for year 2003 to 2014 from USDA NASS survey data for the state of
*Wisconsin. Units are $/bu for corn and soy and $/ton for alfalfa.

* inflation has been considered in our analysis.
* All the prices has been adjusted to the price level of 2010. Source for inflation rate:http://www.usinflationcalculator.com/

p(i)          average price by crop    /cfc 3.90, cfs 3.90, sfc 9.43, c1fa 3.90,
         c2fa 3.90, c3fa 3.90, a1fc 117.03, a2fc 117.03, a3fc 117.03/

*Exogenous supply elasticities are based on the range of supply elasticity
*estimates published in the economic literature and summarized such that those
*observations that are weighted most heavily are those published in AJAE, those
*estimated more recently, and those for a study region closest to Wisconsin.
*95% confidence intervals are biven by [0.021,0.629] for corn, [0.041,0.719]
*for soybeans, and [0.363,0.633] for alfalfa. Unweighted means of published
*estimates are used in the baseline model simulation.

eta(i)   exogenous supply elasticities for calibration /cfc 0.353, cfs 0.353,
        sfc 0.387, c1fa 0.353, c2fa 0.353, c3fa 0.353, a1fc 0.498, a2fc 0.498,
        a3fc 0.498/


*Gross margins are calculated as revenue less variable and fixed costs. Revenue
*calculations are based on yields and prices published in enterprise budgets
*by the University of Wisconsin Extension. The gross margin for establishment
*year alfalfa is based on an expected yield of 2.2 tons/ac and establishment
*costs. This is not used as the minimum gross margin in the PMP calibration
*because the large second and third year margins on alfalfa imply that alfalfa
*is not the marginal land use.
gmar(i)  gross margins by crop /cfc 76.28, cfs 214.66, sfc 192.9, c1fa 113.99,
         c2fa 111.22, c3fa 76.28, a1fc -36.01, a2fc 351.41, a3fc 351.41/

*Land shares are derived from remote sensing analysis of crop rotations by
*Kemanian, White, and Rozum (2017). The dominant rotations in the Mendota
*catchment are continuous corn (67.4%), corn-soy (8.4%), and corn-corn-corn-
*alfalfa-alfalfa-alfalfa (24.2%).

*Since we are considering fallow in the calibration process, the land shares for fallow is 0.0138, thus,
*total land shares for crops are 0.9862, so we need to adjust the land share of each crop rotation accordingly
share(i) land shares /cfc 0.6646, cfs 0.0414, sfc 0.0414, c1fa 0.0398, c2fa 0.0398,
         c3fa 0.0398, a1fc 0.0398, a2fc 0.0398, a3fc 0.0398/

*Per-unit input costs are calculated using enterprise budgets published by the
*University of Wisconsin Extension. Units for land are acres and include all
*operating costs other than fertilizer, which are assumed to be used in fixed
*proportion to land. Units for N are pounds and include all fertilizer inputs
*where phosphorus (P) and potassium (K) are assumed to be used in fixed
*proportion to nitrogen (N).
table   c(i,l)   per unit costs of inputs
             land        N
cfc         278.12     0.695
cfs         273.186    0.842
sfc         169.0239   2.882
c1fa        277.3849   2.808
c2fa        277.4401   0.856
c3fa        278.1215   0.695
a1fc        169.6961   6.404
a2fc        115.663    5.473
a3fc        115.663    5.473

;



*Each year, Land base calculated using number of pixels in corn , alfalfa , and
*soybeans  from the USDA Cropland Data Layers for 2003-2014.  In here, we take the average of land base.
scalar
b1       average land base in 1000 acres /63.3833/
b2       maximum land base in 1000 acres /68.8843/

;
**************************  ESTIMATED PARAMETERS  **************************
*Yield functions are estimated based on Cycles simulation output for each
*crop in each year from 1980-2014 and with N applications ranging from zero to
*200 kg/ha (0 to 214.1643 lb/ac) in 10-unit increments. Simulated yield data
*is fitted to the nonlinear Mitscherlisch-Baule function with year-specific
*constant terms. Yields for soybeans and alfalfa are invariant to fertilizer
*applications but vary by year.

* In here, we take the average of year-specific constant terms for year 2003-2014



*table    yldp(i,yp)  yield function parameters
*         alpha     beta0       gamma      beta2      rho        beta4
*cfc               507.4714    -0.63903    -3900.29   3900.307   85.76747
*cfs               196.424     -0.12019    -4884.41   4884.464   41.20549
*sfc     48.3814
*c1fa              155.9927    0.076774    317.3992   -292.809   33.67279
*c2fa              304.4659   -0.4478      -0.09376   0.205327   37.86927
*c3fa              330.3981   -0.47785     -0.162499  0.231299   41.62916
*a1fc     3.9722
*a2fc     4.5434
*a3fc     4.3878


;

table    yldp(i,yp)  yield function parameters
         alpha     beta0       gamma      beta2      rho    beta4
cfc               126.7652    0.42239     0.013879   0     0.005102
cfs               38.97159    3.426041    0.025922   0     -9.0881
sfc     55.0926
c1fa              2.05739     80.94295    0.040947   0     -37.9432
c2fa              2.99019     56.44927    0.037627   0     -56.4564
c3fa              4.95003     34.0991     0.029086   0     -78.9817
a1fc     4.4702
a2fc     5.8277
a3fc     5.7339


;

**************************  CALCULATED PARAMETERS  **************************



parameters
fallow           land allocation of fallow
xbar1(i)         baseline land allocation
*The baseline model uses N application rates published in University of
*Wisconsin enterprise budgets as the "observed" level of N fertilizer use.
abarN(i)         recommended N application /cfc 166, cfs 115.3, sfc 16.2,
                 c1fa 27.9, c2fa 94.7, c3fa 166, a1fc 9.9, a2fc 19.8, a3fc 19.8/
*The baseline yield by rotation is calculated using the estimated yield
*functions and abar(i).
ybar(i)          baseline yield by rotation using recommended N application
qbar(i)          baseline production by crop
*The yield elasticity is calculated as (dybar/dabar)*(abar/ybar) where the
*derivative is of the Mitscherlische-Baule function and the elasticity is
*evaluated at the baseline level of N applications and crop yield.
ybarN(i)         yield elasticity with respect to N
ncrop(r)         count of crops within each rotation;


xbar1(i) = share(i)*b2;
ybar(i) = yldp(i,'alpha') + yldp(i,'beta0')*(1 + yldp(i,'gamma') -
          exp(-1*(yldp(i,'beta2') + yldp(i,'rho'))*(yldp(i,'beta4')
          + abarN(i))));
qbar(i) = ybar(i)*xbar1(i);
ybarN(i) = (abarN(i)/ybar(i))*(yldp(i,'beta0')*(yldp(i,'beta2') +
           yldp(i,'rho'))*exp(-1*(yldp(i,'beta2') +
           yldp(i,'rho'))*(yldp(i,'beta4') + abarN(i))));
ncrop(r) = sum(i, 1$cr(r,i));



display ybar, ybarN;



***************************************************************************
**************************  CALIBRATION PROCESS  **************************
***************************************************************************

*******************************  STEP I  **********************************
*The first step in the calibration process is to draw the elasticity of
*substitution parameter from a lognormal distribution with mean 1.15 and
*variance 0.5. Random number generation from the lognormal
*distribution is accomplished with the following algorithm:
*(1) draw a standard normal random variable Z
*(2) define a normally distributed random variable X with mean and std dev:
*m = ln(1.15^2/sqrt(0.5+1.15^2)) and s = sqrt(ln(1+0.5/1.15^2))
*(3) define the lognormally distributed random variable Y = e^X.


sets
G        /G1*G1000/


Parameter
sigma(i)  substitution elasticity by crop
Z(G,i)    underlying normal distribution
Y(G,i)    underlying lognormal distribution;



scalar
m mean of lognormal distribution /-0.02058/
var variance of lognormal distribution /0.56629/;



Z(G,i)=NORMAL(m,var);
Y(G,i)=Exp(Z(G,i));



sigma(i)=sum(G,Y(G,i))/CARD(G);

display sigma;




*******************************  STEP II  *********************************
*The second step in the calibration process is to calibrate the modeled
*supply elasticities against the exogenous supply elasticities.

****** Step II: Part A *****
*Ensure that the two calibration criteria hold per Merel et al. (2011),
*such that all CES delta parameters can be calibrated against the set of
*exogenous supply elasticities.

parameters
cc1(i)           calibration criterion 1
ccv1(i)             violation of calibration criterion 1
b(i)             parameter from Merel et al. 2014
bb(i)            1000*b
psi(i)           term inside calibration criterion 2
cc2(i)           calibration criterion 2
ccv2(i)          violation of calibration criterion 2;

*The first criterion requires that cc1 > 0 for all i.
cc1(i) = eta(i) - ybarN(i)/(1-ybarN(i));
display cc1;
ccv1(i)$(cc1(i) gt 0) = 1;
display ccv1;

*The second criterion requires that cc2 < 0 for all i.
b(i) = sqr(xbar1(i))/(p(i)*qbar(i));
bb(i)=1000*b(i);
psi(i) = sigma(i)*ybarN(i)/(eta(i)*(1-ybarN(i)));
cc2(i) = b(i)*eta(i)*(1-psi(i)) - sum(j$(ord(j) ne ord(i)),
         b(j)*eta(j)*sqr(1 + (1/eta(j)))*(1 + psi(j) - ybarN(j)));
ccv2(i)$(cc2(i) lt 0) = 1;

display bb,cc2,ccv2;


***** Step II: Part B *****
*Use the myopic delta parameters to calculate the initial adjustment term for
*the elasticity calibration system. Solve for the deltas that calibrate against
*exogenous supply elasticities. Iterate over the deltas and adjustment term
*until convergence.



scalars
term             indicator equal to card(i) when all deltas converge /0/
toler            tolerance for convergence /0.001/


parameters
delta0(i)        myopic CES production function parameter
adj(i)           adjustment term using myopic delta
error(i)         absolute value of change in delta
converge(i)      indicator equal to one when delta converges;
delta0(i) = eta(i)/(1+eta(i));
adj(i) = 1 - (b(i)/(delta0(i)*(1-delta0(i))))/(sum(j,
         (b(j)/(delta0(j)*(1-delta0(j)))) +
         (sigma(j)*b(j)*ybarN(j)/(delta0(j)*(delta0(j) - ybarN(j))))));

variables
delta(i)         calibrated CES production function parameter
dummy            dummy objective
;

positive variables delta;

equations
etacalib(i)      calibration to exogenous supply elasticity
edummy           dummy objective function;
etacalib(i).. eta(i) =e= (delta(i))/(1-delta(i))*adj(i);
edummy.. dummy =e= 0;

model stage2 /etacalib,edummy/;


while(term lt card(i),
solve stage2 maximizing dummy using nlp;
*Test for convergence in the deltas
         error(i) = abs(delta0(i) - delta.l(i));
         converge(i)$(error(i) lt toler) = 1 ;
         term = sum(i, converge(i));
*Update values for delta in the adjustment term if convergence test fails
         delta0(i) = delta.l(i);
         adj(i) = 1 - (b(i)/(delta0(i)*(1-delta0(i))))/(sum(j,
                  (b(j)/(delta0(j)*(1-delta0(j)))) +
                  (sigma(j)*b(j)*ybarN(j)/(delta0(j)*(delta0(j) - ybarN(j))))));
);


*******************************  STEP III  ********************************
*The third step in the calibration process is to use the agronomic yield
*elasticities to calibrate the production function parameters mu(i) & beta(i).

parameters
xbar2(i)         baseline N input use by crop
rrho(i)           CES production function parameter;
xbar2(i) = abarN(i)*xbar1(i);
rrho(i) = (sigma(i) - 1)/sigma(i);

display rrho;

scalar
baseN            baseline total N input use;
baseN = sum(i, xbar2(i));

variables
beta(i,l)        CES share parameter
mu(i)            CES scale parameter;

positive variables beta,mu;

equations
yresponse(i)     calibration against agronomic yield elasticity for N
betasum(i)       summation constraint
technology(i)    production function;

yresponse(i).. ybarN(i)*(beta(i,'land')*(xbar1(i)**rrho(i))
               + beta(i,'N')*(xbar2(i)**rrho(i))) =e=
               delta.l(i)*beta(i,'N')*(xbar2(i)**rrho(i));
betasum(i).. sum(l, beta(i,l)) =e= 1;
technology(i).. qbar(i) =e= mu(i)*((beta(i,'land')*(xbar1(i)**rrho(i)) +
                beta(i,'N')*(xbar2(i)**rrho(i)))**(delta.l(i)/rrho(i)));

model stage3 /yresponse,betasum,technology,edummy/;

solve stage3 maximizing dummy using nlp;


*******************************  STEP IV  *********************************
*The fourth step in the calibration process is to solve for the behavioral
*parameters lambda to reflect the crop-specific shadow costs of land and N.
*These are the input cost adjustment terms that rationalize observed economic
*behavior given current prices and the calibrated CES production function. The
*behavioral parameters are chosen such that the first-order necessary
*conditions of the regional profit maximization problem hold at the reference
*allocation. The profit maximization problem is given by:
*        max     sum(i, p(i)*q(i) - ((c(i,'land') +
*                lambda(i,'land'))*x(i,'land') + (c(i,'N') +
*                lambda(i,'N'))*x(i,'N'))
*        s.t.    sum(i, x(i,'land')) <= b1
*                q(i) = mu(i)*(sum(l, beta(i,l)*
*                (x(i,l)**rho(i))))**(delta(i)/rho(i))
*These conditions can hold as long as the calibration conditions
*from Step II of the calibration process hold. Negative values for the lambda
*can be interpreted as hidden social benefits from factor use; positive values
*can be interpreted as a hidden social cost.



*The first-stage shadow value for land is calibrated following Garnache et al. (2017), the shadow value calibration become:
* min sum(i, (p(i)*qbar(i)*(delta(i)-ybarN(i))-(c(i,'land')+lbar*xbar(i)))**2)
* in which we choose the reference shadow values lambda so as to minimise
* the sum of squared deviations between the modelled activity- and input-level expenditures and their observed values in the reference allocation
* this algorithm maximises the informational value of the available data


variables
lbar first stage shadow value for land
lobj objective
lambda(i,l)      factor shadow prices
;

equations
slbar        equations of shadow value for land
foc1(i)          first-order necessary condition wrt land
foc2(i)          first-order necessary condition wrt N;

slbar..lobj=e=sum(i,sqr((p(i)*qbar(i)*(delta.l(i)-ybarN(i))-(c(i,'land')+lbar*xbar1(i))))) ;

foc1(i).. p(i)*qbar(i)*(delta.l(i) - ybarN(i)) =e= (c(i,'land')
          + lambda(i,'land') + lbar.l)*xbar1(i);
foc2(i).. p(i)*qbar(i)*ybarN(i) =e= (c(i,'N') + lambda(i,'N'))*xbar2(i);

model shadowland /slbar/;
model stage4 /foc1,foc2,edummy/;
solve shadowland minimizing lobj using nlp;

solve stage4 maximizing dummy using lp;
display lbar.l, lambda.l ;


***************************************************************************
**************************SIMULATION MODEL (Policy Scenarios)  ***************************
***************************************************************************

parameters

baseleach(t)     baseline nitrogen leaching /2003 205.43, 2004 346.86, 2005 122.01, 2006 340.20, 2007 862.13, 2008 674.72, 2009 169.93, 2010 1264.64, 2011 278.11, 2012 456.82, 2013 431.17, 2014 160.73/;

**************************  POLICY SCENARIOS  *****************************

**************************  Command and control  *****************************
equations
*eunleach(i,t)
*efleach(i,t)
conleach(t)  constraint of leaching
constraintsfc    constraint of N applied
constrainta1fc
constrainta2fc
constrainta3fc ;


constraintsfc..x('sfc','N')/x('sfc','land') =l=15.974;
constrainta1fc..x('a1fc','N')/x('a1fc','land') =l=9.786;
constrainta2fc..x('a2fc','N')/x('a2fc','land') =l=19.573;
constrainta3fc..x('a3fc','N')/x('a3fc','land') =l=19.573;

* 5% reduction in nitrogen
conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.95*baseleach(t);


* 10% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.9*baseleach(t);

* 15% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.85*baseleach(t);

* 20% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.8*baseleach(t);

*25% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.75*baseleach(t);

*30% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.7*baseleach(t);

*35% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.65*baseleach(t);

*40% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.6*baseleach(t);

*45% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.55*baseleach(t);

*50% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.5*baseleach(t);

*55% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.45*baseleach(t);

*60% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.4*baseleach(t);

*65% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.35*baseleach(t);

*70% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.3*baseleach(t);

*75% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.25*baseleach(t);

*80% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.2*baseleach(t);

*85% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.15*baseleach(t);

*90% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.1*baseleach(t);

*95% reduction in nitrogen
*conleach(t)..sum(r, sum(i$cr(r,i),(leach(i,'alpha',t) + leach(i,'beta0',t)*x(i,'N')/x(i,'land')+leach(i,'beta1',t)*sqr(x(i,'N')/x(i,'land')))*x(i,'land')))=l=0.05*baseleach(t);



model conprimal /conleach,profitfallow,rescon,rotcon,cesfun,profit,constraintsfc,constrainta1fc,constrainta2fc,constrainta3fc/;


solve conprimal maximizing tprofit using nlp;


display tprofit.l;

parameters

lland(i)          land allocation by crop
lN(i)
factorc(i)         factor cost of land
land(r)          land allocation by rotation
productionr(r,k) commodity production by rotation
production(k)    commodity production
Napplied(r)      N applications by rotation
unitN(i)         N application per acre
unitleach(i,t)   Per acre NO3 leaching by crop
leaching(i,t)      NO3 leaching by crop
tleach(r,t)        NO3 leaching by rotation
totalleach      NO3 leaching by year
unitemit(i,t)     Per acre N2O emissions by crop
emissions(i,t)     N2O emissions by crop
temit(r,t)         N2O emissions by rotation
rtntoN(r,t)         returns to N
rprofit(r,t)        profit
fallow          land allocation of fallow
totalemit       N2O emission by year

;

factorc(i)=c(i,'land')+lbar.l;
lland(i)= x.l(i,'land');
lN(i)= x.l(i,'N');
land(r) = sum(i$cr(r,i), x.l(i,'land'));
productionr(r,k) = sum(i$(cr(r,i) and cc(k,i)), q.l(i));
production(k) = sum(i$cc(k,i), q.l(i));
Napplied(r) = sum(i$cr(r,i), x.l(i,'N'));
unitN(i)=x.l(i,'N')/x.l(i,'land')  ;
*unitleach(i,t)=  leach(i,'alpha',t) + leach(i,'beta0',t)*unitN(i,t) +leach(i,'beta1',t)*sqr(unitN(i,t)) ;
unitleach(i,t)=  leach(i,'alpha',t) + leach(i,'beta0',t)*unitN(i) +leach(i,'beta1',t)*sqr(unitN(i)) ;
*unitleach(i,t)=  exp(leach(i,'alpha',t) + leach(i,'beta0',t)*unitN(i,t) +leach(i,'beta1',t)*sqr(unitN(i,t))) ;
leaching(i,t) = unitleach(i,t)*x.l(i,'land');
tleach(r,t) = sum(i$cr(r,i), leaching(i,t));
totalleach(t)=sum(r,tleach(r,t))  ;
*unitemit(i,t)=  emit(i,'alpha',t) + emit(i,'beta0',t)*unitN(i,t) +emit(i,'beta1',t)*sqr(unitN(i,t));
unitemit(i,t)=  emit(i,'alpha',t) + emit(i,'beta0',t)*unitN(i) +emit(i,'beta1',t)*sqr(unitN(i));
*unitemit(i,t)=  exp(emit(i,'alpha',t) + emit(i,'beta0',t)*unitN(i,t) +emit(i,'beta1',t)*sqr(unitN(i,t)));
emissions(i,t) = unitemit(i,t)*x.l(i,'land')  ;
temit(r,t) = sum(i$cr(r,i), emissions(i,t));
totalemit(t)=sum(r,temit(r,t));
fallow=fx.l;

rtntoN(r,t) = sum(i$cr(r,i), acp(i)*q.l(i) - c(i,'N')*x.l(i,'N'));
rprofit(r,t) = sum(i$cr(r,i), acp(i)*q.l(i) - sum(l, (c(i,l) +pmplambda(i,l))*x.l(i,l)));



display x.l,q.l,fx.l;
display lland,land,fallow,productionr,production,Napplied,unitN,unitleach,leaching,tleach,totalleach,unitemit,emissions,temit,rtntoN,rprofit;

execute_unload "n5red.gdx"  land fallow production Napplied tleach temit rtntoN rprofit
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=land Rng=land!'
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=fallow Rng=fallow!'
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=production Rng=Production!'
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=Napplied Rng=Napplied!'
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=tleach Rng=leaching!'
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=temit Rng=emissions!'
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=rtntoN Rng=ReturnN!'
execute 'gdxxrw.exe n5red.gdx o=avgn5red.xls Par=rprofit Rng=Profit!'


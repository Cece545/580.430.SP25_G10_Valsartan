function dydt = eqns_v2(t,y,p)
% equations describing the pharmacokinetics of valsartan
% in a one-compartment model with absorption from the gut and
% clearance from the central compartment. 

dydt = zeros(10,1);    % use a column vector 
 
%% EQUATIONS
% 1 - concentration of Valsartan in central compartment with bolus injection,
% first-order clearance, second-order complexation with Receptors
% (first-order wrt. Valsartan), and first order dissociation from
% Receptor-Valsartan complex
% 2 - effective concentration of Receptors in central compartment, with
% second-order complexation with Valsartan (first-order wrt. Receptors),
% second-oder complexation with Angiotensin II (first-order wrt.
% Receptors), first-order dissociation from Valsartan-Receptors Complex,
% and first-order dissociation from Angiotensin II-Receptors Complex
% 3 - concentration of Angiotensin II in central compartment, with
% first order conversion from ACE-Angiotensin I Complex, first-order production rate controlled by Angiotensin II-Receptor Complex , second-order
% complexation with Receptors (first-order wrt. Angiotensin II)
% 4 - concentration of Valsartan-Receptor Complex in central compartment,
% with second-order complexation from Valsartan and Receptor, first order
% dissociation 
% 5 - concentration of Angiotensin II-Receptor Complex, with second order
% complexation from Angiotensin II and Receptor, and first-order
% dissociation
% 6 - concentration of Valsartan-Albumin Complex, with pseudo-first-order
% complexation from Valsartan and first-order dissocaiton 
% 7 - amount of cleared Valsartan in the virtual compartment
% 8 - amount of cleared Angiotensin II in the virtual compartment
% 9 - amount of Angiotensin II production in the virtual compartment
% 10 - concentration of Valsartan in gut


p.C_HSA = 0.526*10^(6); %Concentration of HSA, assumed constant
p.k_CL_1 = 2 / p.Vd/0.05; %Valsartan clearance rate coefficient, 1/h
p.k_off_6 = 3.836 * 10^4 * 3600; %conversion rate from Valsartan-Albumin Complex to Valsartan, h^(-1) 
p.k_on_6 = 1.4*3600 * 0.526*10^(6); %complexation rate constant from Valsartan to Valsartan-Albumin Complex, nM^(-1)*h^(-1) 
p.k_on_4 = 2.686 * 10^(-5)*3600; %complexation rate constant form Valsartan and Receptors to Valsartan-Receptors, nM^(-1)*h^(-1) 
p.k_off_4 =  0.00016576666*3600; %dissociation rate constant from Valsartan-Receptors to Valsartan and Receptors, h^(-1)
p.k_on_5 = 3.93 * 10^(-5)*3600; %complexation rate constant from Receptors and Angiotensin II to Angiotensin II-Receptor Complex, nM^(-1)*h^(-1)
p.k_off_5 = 9.625 * 10^(-4)*3600; %dissociation rate constant from Angiotensin II-Receptor Complex to Angiotensin II and Receptor, h^(-1)
p.k_a = 1.409; %absorption rate constant of Valsartan, 1/h
p.k_CL_3 = 0; % clearance rate constant of Angiotensin II
p.k_feedback = 0;

dydt(1) = -p.k_CL_1*y(1) + p.k_off_6*y(6) - p.k_on_6*y(1) - p.k_on_4*y(1)*y(2) + p.k_off_4*y(4)                                       + p.k_a*y(10); %Valsartan
dydt(2) =                                                 - p.k_on_4*y(1)*y(2) + p.k_off_4*y(4) - p.k_on_5*y(2)*y(3) + p.k_off_5*y(5); %Receptors
dydt(3) = -p.k_CL_3*y(3) + p.k_feedback                                                         - p.k_on_5*y(2)*y(3) + p.k_off_5*y(5); %Angiotensin II                            
dydt(4) =                                                 + p.k_on_4*y(1)*y(2) - p.k_off_4*y(4); %Valsartan-Receptors Complex
dydt(5) =                                                                                       + p.k_on_5*y(2)*y(3) - p.k_off_5*y(5); %Angiotensin II-Receptor Complex 
dydt(6) =                - p.k_off_6*y(6) + p.k_on_6*y(1); %Valsartan-Albumin Complex
dydt(7) = +p.k_CL_1*y(1) ; %"concentration" of cleared Valsartan in the virtual compartment
dydt(8) = +p.k_CL_3*y(3); %%"concentration" of cleared Angiotensin II in the virtual compartment
dydt(9) = -p.k_feedback;
dydt(10)=                                                                                                         - p.k_a*y(10); %Absorption through gut

function [C, CU] = controller_fuzzy3(E, DE)

C = [[0,0,-0.266;1,1,0.2690;0,1,-1.324;1,0,1.3260;0.501,0.501,0.0020;0,0.38127,-0.66;0.311,0,0.220;0.9732,0.5819,0.660;0.6421,0.9832,-0.280;0.4816,0,0.50;0.4515,0.993,-0.590;0,0.5418,-0.830;0.9732,0.3812,0.880;0.5618,0,0.620;0.2709,0.5819,-0.450;0,0.6321,-0.930;0.6020,0.2842,0.390;0.6120,0,0.70;0.6822,0,0.820;0.3311,0.9732,-0.760;0,0.1137,-0.380;0.9732,0.2842,0.980;0.9732,0,1.280;0.9832,0.8026,0.450;0,0.6822,-0.980;0,0.9832,-1.3]];

%Normalizing the inputs to 0-1 interval
EN = (E - (-1.1))/2.2 ;
DEN = (DE - (-2))/4;

%mu=0.3;
%mu1= (1-mu)/mu;

mu=[0.358,0.3549,0.3568,0.3568,0.3549,0.1274,0.1883,0.1743,0.2311,0.050,0.05,0.05,0.05,0.05,0.05,0.05,.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05];
mu1= (1-mu)./mu;


lambda = 4;
e = 0.33;
Noise=0.1;       %0.05;


for i= 1: length(C(:,1))
    
% 1st input type2 Membership function
deltal1 = abs(((EN-(C(i,1)))/ e ) * (1 / (1 + exp(100*((EN-(C(i,1))))))))^lambda;
deltar1 = abs(((EN-(C(i,1)))/ e ) * (1 / (1 + exp(-100*((EN-(C(i,1))))))))^lambda;
G11U(i) = 1/(1+ (mu1(i) * deltal1) + (mu1(i) *deltar1));

deltal2 = abs((((EN-(C(i,1)))- Noise)/ (e) ) * (1 / (1 + exp(100*(((EN-(C(i,1)))- Noise))))))^(lambda);
deltar2 = abs((((EN-(C(i,1)))+ Noise)/ (e) ) * (1 / (1 + exp(-100*(((EN-(C(i,1)))+ Noise))))))^(lambda);
G11L(i) = 1/(1+ (mu1(i) * deltal2) + (mu1(i) *deltar2));

% 2nd Input type2 Membership function
deltal1 = abs(((DEN-(C(i,2)))/ e ) * (1 / (1 + exp(100*((DEN-(C(i,2))))))))^lambda;
deltar1 = abs(((DEN-(C(i,2)))/ e ) * (1 / (1 + exp(-100*((DEN-(C(i,2))))))))^lambda;
G12U(i) = 1/(1+ (mu1(i) * deltal1) + (mu1(i) *deltar1));
 
deltal2 = abs((((DEN-(C(i,2)))- Noise)/ (e) ) * (1 / (1 + exp(100*(((DEN-(C(i,2)))- Noise))))))^(lambda);
deltar2 = abs((((DEN-(C(i,2)))+ Noise)/ (e) ) * (1 / (1 + exp(-100*(((DEN-(C(i,2)))+ Noise))))))^(lambda);
G12L(i) = 1/(1+ (mu1(i) * deltal2) + (mu1(i) *deltar2));

end

for i= 1: length(C(:,1))
    
RU(i) = 1/(((1-G11U(i))/G11U(i)) + ((1-G12U(i))/G12U(i)));
RL(i) = 1/(((1-G11L(i))/G11L(i)) + ((1-G12L(i))/G12L(i)));   
   
end

RUN = RU/(sum(RU));
RLN = RL/(sum(RL));

CU = dot(RUN,C(:,3));
CL = dot(RLN,C(:,3));
C = (CU + CL) / 2 ;

end
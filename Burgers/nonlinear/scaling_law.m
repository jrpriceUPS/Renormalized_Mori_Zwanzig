function coeffs = scaling_law(N,type)
%
% coeffs = scaling_law(N,type)
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%     N  =  resolution of ROM
%
%  type  =  type of ROM (format: 'cn*style*' where n is 1,2,3, or 4, and
%           *style* is B for constant of KdV for algebraically decaying)

% computed on [0,10] April 18
coeffs = zeros(4,1);
if strcmp(type,'c1B')
    coeffs(1) = exp(1.181884962799880)*N^-1.136707559532098;
end

if strcmp(type,'c1KdV')
    coeffs(1) = exp(0.705311893784875)*N^-1.201865905462570;
end

if strcmp(type,'c2B')
    coeffs(1) = exp(1.819785538424278)*N^-1.093945832760621;
    coeffs(2) = -exp(2.463587644947313)*N^-2.325237054919805;
end

if strcmp(type,'c3B')
    coeffs(1) = exp(1.541767557743583)*N^-0.981319248449383;
    coeffs(2) = -exp(1.922964860335419)*N^-1.982058125808136;
    coeffs(3) = exp(0.956206908843237)*N^-2.936870830980371;
end

if strcmp(type,'c4B')
    coeffs(1) = exp(1.905282695041228)*N^-1.066687205300492;
    coeffs(2) = -exp(2.796251624527494)*N^-2.162411635298580;
    coeffs(3) = exp(3.152391615308100)*N^-3.341672392051283;
    coeffs(4) = -exp(2.749000286715377)*N^-4.578886683293675;
end

if strcmp(type,'c2KdV')
    coeffs(1) = exp(1.355925631087403)*N^-1.167283767486007;
    coeffs(2) = -exp(1.625454761197602)*N^-2.495257671455809;
end

if strcmp(type,'c3KdV')
    coeffs(1) = exp(1.127866001886322)*N^-1.080944565654628;
    coeffs(2) = -exp(0.840901100313817)*N^-2.125957363772762;
    coeffs(3) = exp(-2.634244475274525)*N^-2.568963699218120;
end

if strcmp(type,'c4KdV')
    coeffs(1) = exp(1.238239895119442)*N^-1.082866452489714;
    coeffs(2) = -exp(1.577850807631149)*N^-2.233389525099544;
    coeffs(3) = exp(1.514640988759219)*N^-3.506872586205894;
    coeffs(4) = -exp(0.498971763306214)*N^-4.772363054617582;
end
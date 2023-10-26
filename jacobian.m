% Jacobian Function %
% Finite difference %
function[J] = jacobian(data,r,t,roa1)
global frequency;
global rhoa;
global lr;
global lt;

par = 0.15;
%s = frequency;

r2 = r;
for i2 = 1:lr                       % length of resistivity
    r2(i2) = (r(i2)*par)+r(i2);
        %for ii = 1:length(frequency)       % length of measurement
            %s = data(ii);
            %[g] = VES1DFWD(r2,t,s); % forward modeling
            %roa2(ii,:) = g;
        %end
        [roa2, phase] = modelagem1DMT(r2, t, frequency);
    J1(:,i2) = [(roa2-roa1)/(r(i2)*par)]*r(i2)./rhoa;
    r2 = r;
end

t2 = t;
for i3 = 1:lt                       % length of thickness
    t2(i3) = (t(i3)*par)+t(i3);
        %for ii = 1:length(frequency)
        %    s = data(ii);
        %    [g] = VES1DFWD(r,t2,s);
        %    roa3(ii,:) = g;
        %end
        [roa3, phase] = modelagem1DMT(r, t2, frequency);
    J2(:,i3) = [(roa3-roa1)/(t(i3)*par)]*t(i3)./rhoa;
    t2 = t;
end

    J = [J1 J2];
return

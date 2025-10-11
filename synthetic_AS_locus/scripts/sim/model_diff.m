function dy=model_diff(t,y,p, kdeg)

dy = zeros(1,1);

if t < p(2)
    a=kdeg;
else
    a=kdeg*exp(-p(1)*(t-p(2)));
end

dy(1) = a-kdeg*y(1); % protein


end

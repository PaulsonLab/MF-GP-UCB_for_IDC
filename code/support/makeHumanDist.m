function [] = makeHumanDist()

OI = zeros(1,365*24);
for i =1:356*24
    if mod(i,24)>=8 && mod(i,24)<=18
        OI(i) = randi([25 35],1,1);
    end
end

humDist = OI;

save('humDist')



end
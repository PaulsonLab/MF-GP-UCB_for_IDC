%% function x0 = getRand(xl, xu)
% draws a Uniform random variable x0 in [xl,xu]
function x0 = getRand(xl, xu)


x0 = (xu-xl).*rand(size(xl)) + xl;

end





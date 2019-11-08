% Find cube root of array 'target' without using NTHROOT and compare speed to
% using FZERO.

options = optimset('TolX', 1e-9);
target = [(-100:.1:100)' (-1000:1:1000)'];

tic;
xfz = zeros(size(target));
for ii = 1:numel(target)
    xfz(ii) = fzero(@(x) x.^3-target(ii), [-20 20], options);
end
fzero_time = toc

tic;
xbis = bisection(@(x) x.^3, -20, 20, target, options);
bisection_time = toc

fprintf('FZERO took %0.0f times longer than bisection.\n',...
    fzero_time/bisection_time)
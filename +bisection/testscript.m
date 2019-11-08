% Test and verify bisection.
% 
% To run: 
% <a href="matlab: runtests bisection.testscript">runtests bisection.testscript</a>

fk = @(x) -(x+3).*(x-1).^2; % This function kisses but doesn't cross zero at 1.
fi = @(x) fk(-x);

counter reset
[x, fval, exitFlag] = bisection(@countedcube,3,6)
counter

%% Basic usage
[x, fval, exitFlag] = bisection(@sin,-1,2,0,optimset('TolX',1e-4));
assert(abs(x)<1e-4);
assert(abs(fval)<1e-3);
assert(exitFlag>0)

%% Basic DimVar usage
options.TolX = u.um;
options.TolFun = u.nm;
[x, fval, exitFlag] = bisection(@(x) x,-1*u.m,2*u.m,[],options);

%% Check function count for instant-converge
% The number of function counts that are 'overhead' is up to 3:
%   One for checking if f returns an array
%   Two for unbounded root check
% Add to that a couple of evaluations in the main loop.
counter reset
x = bisection(@countedcube,-1,1);
assert(x==0);
counts = counter;
fprintf('%g function evaluations for simple instant-converge\n',counts)

%% Check function count for normal case
counter reset
[x, ~, exitFlag] = bisection(@countedcube,-1,2);
assert(abs(x)<1e-5);
assert(exitFlag == 1);
counts = counter;
if counts > 30
    warning('This should be doable in many fewer function evaluations.')
end
fprintf('%g function evaluations for simple normal case\n',counts)

%% Help example 1
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

% Verify correctness.
[xbis, fval] = bisection(@(x) x.^3, -20, 20, target, options);

exact = nthroot(target,3);
xErr = exact-xbis;
assert(all(abs(xErr(:)) <= options.TolX))
assert(all(abs(fval(:)-target(:)) <= 1000*options.TolX))


%% Help example 2
[A, B] = meshgrid(linspace(1,2,6), linspace(4,12,10));
f = @(x) A.*x.^0.2 + B.*x.^0.87 - 15;
xstar = bisection(f,0,5);

assert(all(all(abs(f(xstar)) < 1e-4)))

%% Help example 1 with expensive function
options = optimset('TolX', 1e-6);
target = [(-100:1:100)' (-1000:10:1000)'];

tic;
xfz = zeros(size(target));
for ii = 1:numel(target)
    xfz(ii) = fzero(@(x) delayedcube(x)-target(ii), [-20 20], options);
end
fzero_time = toc

tic;
xbis = bisection(@(x) delayedcube(x), -20, 20, target, options);
bisection_time = toc

fprintf('FZERO took %0.0f times longer than bisection.\n',...
    fzero_time/bisection_time)

% Verify correctness.
[xbis, fval] = bisection(@(x) x.^3, -20, 20, target, options);

exact = nthroot(target,3);
xErr = exact-xbis;
assert(all(abs(xErr(:)) <= options.TolX))
assert(all(abs(fval(:)-target(:)) <= 1000*options.TolX))

%% Test input robustness
% making sure bad inputs return NaN
f = @(x) -x.^3;
options.TolX = 1e-6;
options.TolFun = 1e-5;
target = [(-100:25:100)' (-1000:250:1000)'];
exact = nthroot(-target,3);
target(3) = nan; % Nan-robustness.
% target(5) = 1i; % Complex totally not supported. TODO: document that fact.
ub = 20*ones(size(target));
lb = -ub;
lb(8) = 21; % Solution not in search interval.
lb(9) = 20; ub(9) = -20; % Flipped search interval - should be okay.
ub(15) = nan; % Nan-robustness.

badInd = [3,15,8];
[x, fval, exitFlag] = bisection(f,lb,ub,target,options);

% Checks:
assert(all(isnan(x(badInd))))
assert(all(isnan(fval(badInd))))
assert(all(exitFlag(badInd) < 0))
assert(all(exitFlag(~badInd) > 0))

%% Trick it - roots at endpoints
[x, fval, exitFlag] = bisection(@(x) x,0,5,0);
assert(exitFlag==1)
[x, fval, exitFlag] = bisection(@(x) x,0,5,5);
assert(exitFlag==1)

[x, fval, exitFlag] = bisection(@(x) -x,0,5,0);
assert(exitFlag==1)
[x, fval, exitFlag] = bisection(@(x) -x,0,5,-5);
assert(exitFlag==1)

%% Trick it - roots just outside endpoints
[x, fval, exitFlag] = bisection(@(x)x,eps,5,0);
assert(exitFlag==-2)
[x, fval, exitFlag] = bisection(@(x)x,0,1,1+eps);
assert(exitFlag==-2)

%% Trick it - instant root find
[x, fval, exitFlag] = bisection(fk,-4,6,0); % 1 is evaluated as zero, stopping.
assert(exitFlag == 2)
assert(x == 1)

%% Unsolvable - no roots
[x, fval, exitFlag] = bisection(fk,-3,5,sqrt(eps));
assert(exitFlag < 0 && isnan(x))

%% Multiple roots
[x, fval, exitFlag] = bisection(fk,-3,5,-sqrt(eps)); % Should be solvable.
assert(exitFlag > 0)

%% NaN tolfun
options.TolFun = nan;
options.TolX = 1e-6;
[x, fval, exitFlag] = bisection(@sin,-1,2,0,options);
assert(exitFlag == -1)

%% NaN tolx
options.TolX = nan;
options.TolFun = 1e-6;
[x, fval, exitFlag] = bisection(@sin,-1,2,0,options);
assert(exitFlag == -1)

%% Many roots, straddling bounds
f = @(x) sin(x)+.5;
[x, fval, exitFlag] = bisection(f,-pi/2,128);
assert(exitFlag == 1)

%% Many roots, non-straddling bounds
f = @(x) sin(x)+.5;
[x, fval, exitFlag] = bisection(f,-pi/2,100);
assert(exitFlag == -2)

%% Make sure f is never evaluated outside bounds.
[x, fval, exitFlag] = bisection(@blowupcube,-20,3,27);

%% Disallow TolFun at kissing root
% 1 is evaluated as zero, but never within a negative tolFun.
options.TolFun = -inf;
options.TolX = 1e-6;
[x, fval, exitFlag] = bisection(fk,-4,6,0,options);
assert(exitFlag == 1)
assert(abs(x+3) < 1e-5)

%% Disallow TolFun at kissing root, other direction
% -1 is evaluated as zero, but never within a negative tolFun.
options.TolFun = -inf;
options.TolX = 1e-6;
[x, fval, exitFlag] = bisection(fi,-6,4,0,options);
assert(exitFlag == 1)
assert(abs(x-3) < 1e-5)

%% Disallow TolFun at kissing root, flipped
% 1 is evaluated as zero, but never within a negative tolFun.
options.TolFun = -inf;
options.TolX = 1e-6;
[x, fval, exitFlag] = bisection(@(x) -fk(x),-4,6,0,options);
assert(exitFlag == 1)
assert(abs(x+3) < 1e-5)

%% Disallow TolFun at kissing root, lower bound is root
options.TolFun = -inf;
options.TolX = 1e-6;
[x, fval, exitFlag] = bisection(fk,-3,5,0,options)
assert(exitFlag == 1)
assert(abs(x+3) < 1e-5)

%% Disallow TolFun at kissing root, upper bound is root
options.TolFun = -inf;
options.TolX = 1e-6;
[x, fval, exitFlag] = bisection(fi,-5,3,0,options)
assert(exitFlag == 1)
assert(abs(x-3) < 1e-5)

%% test functions

function output = countedcube(x)
counter plus
output = x.^3;
end

function output = delayedcube(x)
pause(0.001);
output = x.^3;
end

function output = blowupcube(x)
output = x.^3;
if any((x(:) < -20) | (x(:) > 3))
    error('This function cannot be evaluated outside bounds.')
end
end
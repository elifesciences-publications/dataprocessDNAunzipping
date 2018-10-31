function xmid = bisectionsolve(func,xinit,options)
% use a bisection method to find the zero of a monotonically increasing
% function
% xinit contains initial bracket for x values

opt.Tol = 1e-6;
opt.initscl = 2;
opt.display = 'none';

if (nargin>2)
    opt =  copyStruct(options,opt);
end

xa = xinit(1); xb = xinit(2);

% get a lower bound
fa = func(xa);

while fa>0
    if (xa>0)
        xa = xa-exp(opt.initscl);
    else
        xa = xa*opt.initscl;
    end
    fa = func(xa);
    if (opt.display=='iter')
        [xa fa]
    end
end

% get an upper bound
fb = func(xb);
while fb<0
    if (xb<0)
        xb = xb+exp(opt.initscl);
    else
        xb = xb*opt.initscl;
    end
    fb = func(xb);
    if (opt.display=='iter')
        [xb fb]
    end
end

% bisect interval
while 1
    xmid = (xa+xb)/2;
    fmid = func(xmid);
    if (abs(fmid)<opt.Tol)
        break
    end
    if (fmid>=0)
        xb = xmid;
        fb = func(xb);        
    elseif (fmid<0)
        xa = xmid;
        fa = func(xa);
    else
        error('bad fmid value')
    end
    
    if (opt.display=='iter')
        [xa fa xb fb]
    end
end

end
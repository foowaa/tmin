# tmin
Numerical Optimization based on MATLAB.

See [Link](http://www.cltian.site/2016/04/16/%E5%87%A0%E4%B8%AA%E4%BC%98%E5%8C%96%E7%AE%97%E6%B3%95demo/)

```Matlab
function [x,y,varargout]=tmin(f,x0,options,varargin)
options: 
nmax -- iteration loops (default 5000)
epsilon -- default 1e-5
corrections -- default 100
method -- GD, CG, Newton, BFGS, LBFGS (default LBFGS)
varargout -- varargout{1}=t; varargout{2}=n+1;
```

References:

Wright, Stephen, and Jorge Nocedal. "Numerical optimization." Springer Science 35.67-68 (1999): 7.

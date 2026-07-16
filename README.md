# bisection
Very simple to use and robust method that takes array inputs, so it even has advantages over fzero.

[![Install with mip](https://img.shields.io/badge/install_with-mip-1565C0)](https://mip.sh/packages/bisection) [![View Bisection Method Root Finding on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/28150)

[Install with mip](https://mip.sh/packages/bisection):

```matlab
mip install bisection
mip load bisection
```

BISECTION is a fast, simple-to-use, and robust root-finding method that handles n-dimensional arrays.

Additional optional inputs and outputs for more control and capabilities that don't exist in other implementations of the bisection method or other root finding functions like fzero.
This function really shines in cases where fzero would have to be implemented in a loop to solve multiple cases, in which case this will be much faster.

It can find zero or non-zero roots.

This code can be a bit cryptic. This is for the sake of speed and increased capability. See the many acknowledged other submissions for simpler, easier-to-follow implementations to understand the basics of the bisection method.

# HDCA
Hybrid accelerated DCA for solving the Asymmetric Eigenvalue Complementarity Problem (AEiCP).

This project is supported by the National Natural Science Foundation of China (Grant No: 11601327).

## Instruction

Run `install.m` to install the toolbox on MATLAB.

This toolbox is developed based on the DCAM toolbox, see [here](https://github.com/niuyishuai/DCAM).

Three AEiCP models are considered. Seven DCA-type algorithms are implemented, including: 
* BDCA with exact and inexact line search (BDCAe and BDCAa)
* ADCA
* InDCA
* Hybrid DCA with Line search and Inertial force (HDCA-LI)
* Hybrid DCA with Nesterov's extrapolation and Inertial force (HDCA-NI)
* The classical DCA

See the article [here](https://arxiv.org/abs/2305.12076) for more details about these models and algorithms.

## Citation

```
@Misc{niu2023HDCA,
	title = {Accelerated DC Algorithms for the Asymmetric Eigenvalue Complementarity Problem},
	author = {Yi-Shuai Niu},	
	year = {2023},
	eprint={2305.12076},
        archivePrefix={arXiv},
	url = {https://arxiv.org/abs/2305.12076}
}
```

## Dependencies

This toolbox depends on `MOSEK` for solving the linear and quadratic convex subproblems.

The compared optimization solvers are `KNITRO`, `FILTERSD`, `IPOPT` and MATLAB `FMINCON`. Please make sure that you have installed the corresponding solvers on MATLAB for comparison.


## Samples

The Asymmetric Eigenvalue Complementarity Problem (AEiCP) involves finding complementary eigenvectors $x\in \mathbb{R}^n\setminus \\{0\\}$ and complementary eigenvalues $\lambda\in \mathbb{R}$ that satisfy the following conditions: 

$$\begin{cases}
    w = \lambda B   x - A  x, \\
    x^{\top}  w = 0, \\
    0\neq x\geq 0, w\geq 0,
\end{cases}
$$

where $x^{\top}$ represents the transpose of $x$, $A\in \mathbb{R}^{n\times n}$ denotes an asymmetric real matrix, and $B$ is a positive definite matrix that doesn't necessarily have to be symmetric.

Here is a step-by-step example illustrating how to use different types of DCA algorithms to solve this problem:

* Generate a random AEiCP:
``` Matlab
    n=10;
    randbnd=[-1,1];
    T=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
    mu = -min([0,eigs(T+T',1,'smallestreal')]);
    A=T+mu*eye(n);
    B=eye(n);
```

* Set the initial point $x^0$:

We use random initialization 

``` Matlab
    x0.x = rand(n,1);
    x0.x = x0.x/sum(x0.x);
    x0.y = rand(n,1);
    x0.w = B*x0.x - A*x0.y;
    x0.z = sum(x0.y);
```

* Create a DC function object and a DC programming problem object:
``` Matlab
    dcf=dcfunc;
    dcf.f=@(X,n,A,B,opt)fobj_eval_f1(X,n,A,B,opt);
    mydcp=dcp(dcf,[]);
```
In the function `fobj_eval_f1`, we use the first DC formulation (DCP1) proposed in the paper [here](https://arxiv.org/abs/2305.12076), defined as:

$$0 = \min\\{ f_1(x,y,w,z) := \Vert y-z x \Vert^2 + x^{\top}w :  (x,y,w,z)\in \mathcal{C}_1 \\},$$

where 

$$\mathcal{C}_1: = \\{(x,y,w,z)\in  \mathbb{R}^{n} \times \mathbb{R}^{n} \times \mathbb{R}^{n} \times \mathbb{R}: w = Bx -A y, e^{\top}x = 1, e^{\top}y=z, (x,y,w,z)\geq 0 \\},$$

and $f_1(x,y,w,z)$ has a DC-SOS decomposition $g_1(x,y,w,z) - h_1(x,y,w,z)$ where 

$$
\begin{align}
g_1(x,y,w,z) = &\Vert y\Vert^2 + \frac{((z+1)^2 + \Vert y-x\Vert^2)^2 + ((z-1)^2 + \Vert y+x\Vert^2)^2}{16} + \frac{(z^2 	+ \Vert x\Vert^2)^2}{2} + \frac{\Vert x+w\Vert^2}{4}, \\\
h_1(x,y,w,z) = &\frac{((z+1)^2 + \Vert y+x\Vert^2)^2 + ((z-1)^2 + \Vert y-x\Vert^2)^2}{16} + \frac{z^4 + \Vert x\Vert^4}{2} + \frac{\Vert x-w\Vert^2}{4}.
\end{align}
$$

* Create and setup a DCA-type algorithm object
``` Matlab 
	mydca = dca(mydcp,x0);

	mydca.model='DCP1';	
	mydca.A=A;
	mydca.B=B;
	mydca.verbose=true;
	mydca.tolf=0;
	mydca.tolx=0;
	mydca.maxiter=200;

	mydca.linesearch = false;
	mydca.nesterov = true;
	mydca.inertial = true;
```
`mydca.maxiter=200` is the maximum number of iterations for DCA-type algorithm.
`mydca.tolf` and `mydca.tolx` are stopping tolerences such that we will terminate the algorithm if one of the stopping conditions
$$\Vert x^{k+1}-x^k\Vert\leq \text{mydca.tolx} * (1+\Vert x^{k+1}\Vert) \quad \text{ or } \quad |f(x^{k+1})-f(x^k)|\leq \text{mydca.tolf} * (1+|f(x^{k+1})|)$$
is verified. In this example, we set `mydca.tolf=0`, `mydca.tolf=0` and `mydca.xopt=200`, which means the algorithm will terminate after 200 iterations.

Note that three boolean parameters `mydca.linesearch`, `mydca.nesterov`, and `mydca.inertial` are important for choosing different DCA-type algorithms. The settings in `mydca` for each DCA-type algorithm (DCA|ADCA|BDCAe|BDCAa|ADCA|InDCA|HDCA-LI|HDCA-NI) are summarized below:
  - DCA: linesearch = 0; nesterov = 0; inertial = 0;
  - BDCAe: linesearch = 1; nesterov = 0; inertial = 0; linesearch_type='exact';
  - BDCAa: linesearch = 1; nesterov = 0; inertial = 0; linesearch_type='armijo';
  - ADCA: linesearch = 0; nesterov = 1; inertial = 0; adca_q > 0; restartperiod = inf|>0;
  - InDCA: linesearch = 0; nesterov = 0; inertial = 1;
  - HDCA-LI: linesearch = 1; nesterov = 0; inertial = 1; linesearch_type='exact';
  - HDCA-NI: linesearch = 0; nesterov = 1; inertial = 1; adca_q > 0; restartperiod = inf|>0;

* Optimization
``` Matlab 
	status=mydca.optimize();
```

See more examples and optional settings in the folder `tests`.

## Available Dataset
Three datasets for AEiCP: `RAND1`, `RAND2` and `NEP` are available. 

See `GEN_NEP.m` and `GEN_RANDEICP.m` to generate these datasets.

## License

Released under MIT license

## Contact

niuyishuai82@hotmail.com

# HDCA
Hybrid accelerated DCA for solving the Asymmetric Eigenvalue Complementarity Problem (AEiCP)

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
See the article [here](https://arxiv.org/abs/2301.09098) for more details about these models and algorithms.

## Citation

```
@Misc{niu2023BDCASEiCP,
	title = {Accelerated DC Algorithms for the Asymmetric Eigenvalue Complementarity Problem},
	author = {Yi-Shuai Niu},	
	year = {2023}
}
```

## Dependencies

This toolbox depends on `MOSEK` for solving the linear and quadratic convex subproblems.

The compared optimization solvers are `KNITRO`, `FILTERSD`, `IPOPT` and MATLAB `FMINCON`. Please make sure that you have installed the corresponding solvers on MATLAB for comparison.


## Samples

See some test examples in the folder `tests`.

## Available Dataset
Three datasets for AEiCP: `RAND1`, `RAND2` and `NEP`. 

See `GEN_NEP.m` and `GEN_RANDEICP.m` to generate these datasets.

## License

Released under MIT license

## Contact

niuyishuai82@hotmail.com
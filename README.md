# GaussKonrad


A code for calculating Gauss-Kronrod quadrature nodes and weights. The algorithm is based on [Laurie,1997](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.192.3713 )  and the information in [Pavel Holoborodko's  blog](http://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/ ) .

## Requirements

	* [Eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page )
	* [The GNU Multiple Precision Arithmetic Library](https://gmplib.org/ )
	* [The GNU MPFR Library](http://www.mpfr.org/ )
	* [MPFR C++](http://www.holoborodko.com/pavel/mpfr/ )

## Compilation

	$ mkdir build
	$ cd build
	$ cmake -DEIGEN3_INCLUDE_DIR=path_to_Eigen3 -DGMP_ROOT=path_GMP_root_dir -DMPFR_ROOT=path_to_MPFR_root_dir -DMPFRCPP_ROOT=path_to_MPFRC++_root_dir  path_to_GaussKronrad
	$ make

For example,

	$ cmake -DEIGEN3_INCLUDE_DIR=/Users/sbalan/Library/gcc/eigen-3.2.1/include/eigen3/ -DGMP_ROOT=/Users/sbalan/Library/gcc/gmp-6.0.0/ -DMPFR_ROOT=/Users/sbalan/Library/gcc/mpfr-3.1.2/ -DMPFRCPP_ROOT=/Users/sbalan/Library/gcc/mpfrc++-3.5.9/ ../
	$ make

This will prodcue an executable in the directory `bench` called `bench_PH_mpkronrod` . If you run this program you will get 

	$ ./bench/bench_PH_mpkronrod 
	-0.9956571630258080807355273	0.0116946388673718742780644
	-0.9739065285171717200779640	0.0325581623079647274788190
	-0.9301574913557082260012072	0.0547558965743519960313813
	-0.8650633666889845107320967	0.0750396748109199527670431
	-0.7808177265864168970637176	0.0931254545836976055350655
	-0.6794095682990244062343274	0.1093871588022976418992106
	-0.5627571346686046833390001	0.1234919762620658510779581
	-0.4333953941292471907992659	0.1347092173114733259280540
	-0.2943928627014601981311266	0.1427759385770600807970943
	-0.1488743389816312108848260	0.1477391049013384913748415
	0.0000000000000000000000000	0.1494455540029169056649365
	0.1488743389816312108848260	0.1477391049013384913748415
	0.2943928627014601981311266	0.1427759385770600807970943
	0.4333953941292471907992659	0.1347092173114733259280540
	0.5627571346686046833390001	0.1234919762620658510779581
	0.6794095682990244062343274	0.1093871588022976418992106
	0.7808177265864168970637176	0.0931254545836976055350655
	0.8650633666889845107320967	0.0750396748109199527670431
	0.9301574913557082260012072	0.0547558965743519960313813
	0.9739065285171717200779640	0.0325581623079647274788190
	0.9956571630258080807355273	0.0116946388673718742780644

You can compare the results to the ones in [Pavel's blog](http://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/#Tabulated_Gauss-Kronrod_weights_and_abscissae)


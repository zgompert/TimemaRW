import scipy
import numpy
import dadi
import matplotlib
import matplotlib.pyplot as pyplot
import sys
import re
import pylab
from dadi import Numerics, PhiManip, Integration, Spectrum

## get population ids and 2n sample sizes from the command line
#p1 = sys.argv[1];
#p2 = sys.argv[2];
#n1 = sys.argv[3];
#n2 = sys.argv[4];

n1 = 33;
n2 = 25;
p1 = "bce_cc";
p2 = "bce_rw";

ns = numpy.array([int(n1), int(n2)])
## downsample to 70% of minimum
downs = int(ns.min() * 0.7)

## read data from filtered vcf and create downsampled/projected 2d sfs 
dd = dadi.Misc.make_data_dict_vcf("/uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/variants_rw_plus/morefilter_filtered2x_tcr_rw_knulli_perform.recode.vcf","/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/dating_dadi/IDS_plus.txt")
#dd = dadi.Misc.make_data_dict_vcf("/uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/variants_rw_plus/morefilter_filtered2x_tcr_rw_knulli_variants.vcf","/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_fusion/dating_dadi/IDS_plus.txt")
fs = dadi.Spectrum.from_data_dict(dd,[p1, p2],projections=ns, polarized=False)

of = "dd_"+p1+"_"+p2+"_im_old.txt"
pngf = "dd_"+p1+"_"+p2+"_im_old.png"

ofile = open(of,"w")

ns =  fs.sample_sizes

## define IM model, from dadi website, includes exponential growth
def IM(params, ns, pts):
	s, nu1, nu2, T, m12, m21 = params
	xx = Numerics.default_grid(pts)
	phi = PhiManip.phi_1D(xx)
	phi = PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t : s * (nu1/s) ** (t/T)
	nu2_func = lambda t : (1-s) * (nu2/(1-s)) ** (t/T)
	phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12 = m12, m21 = m21)
	fs = Spectrum.from_phi(phi, ns, (xx, xx))
	return fs
    
# These are the grid point settings will use for extrapolation.
# smallest a bit bigger than sample size, this should be close
pts_l = [20,30,40] ## need to change

func = IM

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: (s, nu1, nu2, T, and m12, m21)
upper_bound = [.95, 20, 20, 30, 1, 1]
lower_bound = [.05, 0.1, 0.1, 0, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [0.7, 0.8, 1.2, 8, .1, .1]
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# trying 3 rounds of optimization
ofile = open(of,"w")
print(ofile)
## initial defaults
ll_opt = -999999999999
theta = 1.0
mpopt = p0
for x in range(20):
	# Perturb our parameters before optimization. This does so by taking each
	# parameter a up to a factor of two up or down.
	p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
	# Do the optimization. By default we assume that theta is a free parameter,
	# since it's trivial to find given the other parameters. If you want to fix
	# theta, add a multinom=False to the call.
	# The maxiter argument restricts how long the optimizer will run. For real 
	# runs, you will want to set this value higher (at least 10), to encourage
	# better convergence. You will also want to run optimization several times
	# using multiple sets of intial parameters, to be confident you've actually
	# found the true maximum likelihood parameters.
	popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound,verbose=1, maxiter=20)
	model = func_ex(popt, ns, pts_l)
	ll_model = dadi.Inference.ll_multinom(model,fs)
	if(ll_model > ll_opt):
		ll_opt = ll_model
		mpopt = popt		
		theta = dadi.Inference.optimal_sfs_scaling(model,fs)

## write round 1 max
ofile.write("{0}".format(ll_opt))
ofile.write(" {0}".format(theta))
for a in range(6):
	ofile.write(" {0}".format(mpopt[a]))
ofile.write("\n")

## start round 2 at max
p0 = mpopt
for x in range(10):
	# Perturb our parameters before optimization. This does so by taking each
	# parameter a up to a factor of two up or down.
	p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
	# Do the optimization. By default we assume that theta is a free parameter,
	# since it's trivial to find given the other parameters. If you want to fix
	# theta, add a multinom=False to the call.
	# The maxiter argument restricts how long the optimizer will run. For real 
	# runs, you will want to set this value higher (at least 10), to encourage
	# better convergence. You will also want to run optimization several times
	# using multiple sets of intial parameters, to be confident you've actually
	# found the true maximum likelihood parameters.
	popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound,verbose=1, maxiter=30)
	model = func_ex(popt, ns, pts_l)
	ll_model = dadi.Inference.ll_multinom(model,fs)
	if(ll_model > ll_opt):
		ll_opt = ll_model
		mpopt = popt		
		theta = dadi.Inference.optimal_sfs_scaling(model,fs)

## write round 2 max
ofile.write("{0}".format(ll_opt))
ofile.write(" {0}".format(theta))
for a in range(6):
	ofile.write(" {0}".format(mpopt[a]))
ofile.write("\n")
	
## start round 3 at max
p0 = mpopt
mtheta = theta
for x in range(5):
	# Perturb our parameters before optimization. This does so by taking each
	# parameter a up to a factor of two up or down.
	p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
	# Do the optimization. By default we assume that theta is a free parameter,
	# since it's trivial to find given the other parameters. If you want to fix
	# theta, add a multinom=False to the call.
	# The maxiter argument restricts how long the optimizer will run. For real 
	# runs, you will want to set this value higher (at least 10), to encourage
	# better convergence. You will also want to run optimization several times
	# using multiple sets of intial parameters, to be confident you've actually
	# found the true maximum likelihood parameters.
	popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound,verbose=1, maxiter=50)
	model = func_ex(popt, ns, pts_l)
	ll_model = dadi.Inference.ll_multinom(model,fs)
	if(ll_model > ll_opt):
		ll_opt = ll_model
		mpopt = popt		
		mtheta = dadi.Inference.optimal_sfs_scaling(model,fs)
	## write each round
	theta = dadi.Inference.optimal_sfs_scaling(model,fs)
	ofile.write("{0}".format(ll_model))
	ofile.write(" {0}".format(theta))
	for a in range(6):
		ofile.write(" {0}".format(popt[a]))
	ofile.write("\n")

## write round 3 max... this is the answer
ofile.write("{0}".format(ll_opt))
ofile.write(" {0}".format(mtheta))
for a in range(6):
	ofile.write(" {0}".format(mpopt[a]))
ofile.write("\n")

model = func_ex(mpopt, ns, pts_l)	
pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model,fs,vmin=1,resid_range=10,show=False)
pylab.savefig(pngf, dpi=400)


ofile.close()

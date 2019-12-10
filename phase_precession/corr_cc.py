from pycircstat.descriptive import mean, resultant_vector_length
import numpy as np
from scipy.stats import norm


def corr_cc(alpha, beta):
	'''computes circular correlation coefficient

	input:
	alpha	sample of angles in radians
	beta	sample of angles in radians

	output:
	rho		correlation coefficient
	pval	significance probability

	references:
	Topics in circular statistics, S.R. Jammalamadaka et al., p. 176

	PHB 3/19/2006 2:02PM

	copyright (c) 2006 philipp berens
	berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens
	distributed under GPL with no liability
	http://www.gnu.org/copyleft/gpl.html
	'''

	n = len(alpha)
	alpha_bar = mean(alpha)
	beta_bar = mean(beta)

	num = sum(np.sin(alpha - alpha_bar) * np.sin(beta - beta_bar))
	den = np.sqrt(sum(np.sin(alpha - alpha_bar)**2) * sum(np.sin(beta - beta_bar)**2))

	rho = num / den	# correlation coefficient

	l20 = mean(np.sin(alpha - alpha_bar)**2)
	l02 = mean(np.sin(beta - beta_bar)**2)
	l22 = mean((np.sin(alpha - alpha_bar)**2) * (np.sin(beta - beta_bar)**2))

	ts = np.sqrt((n * l20 * l02) / l22) * rho
	pval = 2 * (1 - norm.cdf(abs(ts)))

	return rho, pval, ts


def corr_cc_uniform(a, b):
	'''computes a  circular correlation coefficient
	according to Jammmalamadaka 2001, page 177, equation 8.2.4, which
	can deal with uniform distributions of a or b.
	This function is equivalent to circCorrJammalamadaka2.m

	 Input
	 a	angles (samples)
	 b	angles (samples)

	 Output
	 rho   corr. coeff.
	 p	significance probability

	 Written by Richard Kempter Feb 13, 2008,
	 to deal with uniform distributions of a or b (book on page 177 , ii)
	 This function is an extension of the matlab function circCorr.m
	 by philipp berens
	 berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens
	'''

	n = len(a)
	a_bar = mean(a)
	b_bar = mean(b)

	aplusb = a + b
	aminusb = a - b

	aplusb_bar =  mean(aplusb)
	aminusb_bar = mean(aminusb)

	R_aplusb =  resultant_vector_length(aplusb)
	R_aminusb = resultant_vector_length(aminusb)


	den = 2 * np.sqrt(sum(np.sin(a - a_bar)**2) * sum(np.sin(b - b_bar)**2))

	rho =  n* (R_aminusb - R_aplusb) / den

	#RK not sure wheter the next equations on the significance are still valid

	l20 = mean(np.sin(a - a_bar)**2)
	l02 = mean(np.sin(b - b_bar)**2)
	l22 = mean((np.sin(a - a_bar)**2) * (np.sin(b - b_bar)**2))

	ts = np.sqrt((n * l20 * l02) / l22) * rho
	p = 2 * (1 - norm.cdf(abs(ts)))
	return rho, p, ts

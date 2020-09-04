import sys
import os
import numpy as np
import matplotlib.pyplot as plt

"""
=================================================

stochCheck.py

Tom Callister
10/23/14

This program is intended as a pedagogical tool for understanding the output of stochastic.m.

-------------------------------
OVERVIEW OF STOCHASTIC.M OUTPUT
-------------------------------

Stochastic.m outputs,
among other things, the files

	(1) /example_ccspectra.jobX.trialY.dat
	(2) /example_sensints.jobX.trialY.dat
	(3) /example_ccstats.jobX.trialY.dat

For each segment analyzed, File (3) contains the broadband cross-correlation (CC) statistic Y and its associated broadband
uncertainty sigma. Given a number of narrowband measurements Y(f) with associated sigma(f)'s, these are combined
to yield a broadband measurement in the following way:

	Y = Sum( Y(f_i) * sigma^(-2)(f_i) )
		/ Sum( sigma^(-2)(f_i) )

	sigma^(-2) = Sum( sigma^(-2)(f_i) )

where the sums are over the set of frequencies {f_i} analyzed. Stochastic.m instead rather confusingly defines two quanties,
a \weighted\ CC spectrum and set of "sensitivity integrands" (in Files (1) and (2) respectively), which are defined
by the relations

	Y = 3*Re(Sum( CC Spectrum(f_i) df ))

	sigma^(-2) = Sum( SensitivityIntegrand(f_i) df )

-------------
SANITY CHECKS
-------------

To help understand these quantities, stochCheck.py reads output from stochastic.m and performs a set of four sanity checks.
Thse are:

	1. Calculate the broadband CC statistic in two ways, first directly using the CC Spectrum data,
	and secondly with the 'true' narrowband statistics, as derived from the CC Spectrum and Sensitivity Integrands
	These values are then printed alongside the value reported in _ccstats. If everything is
	consistent, these three values should all agree (up to some round-off error).

	2. Similary compute the broandband sigma in two ways, and compare to that reported in _ccspectra.

	3. Display the mean and standard deviation of the set of values {Y(f)/sigma(f)}. If the narrowband
	Y(f) and sigma(f) have been computed correctly, the mean and standard devation should be close to
	zero and one, respectively.

	4. For each segment, plot a histogram of {Y(f)/sigma(f)} in order to graphically confirm the results
	from sanity check #3. Additionally, plot sigma(f) just for fun. These plots are stored in the
	directory ./stochCheck_figures/

---------------------
RUNNING STOCHCHECK.PY
---------------------

To run stochCheck.py, use the following command-line syntax:

	>> python stochCheck.py ~/some/directory/path/ 1

The first command-line argument ('~/some/directory/path/' above) specifies the directory in which stochCheck.py will
look for stochastic.m output files. The second argument is some integer which specifies the job number you wish to
analyze. For instance,

	>> python stochCheck.py ./test/dataDir/ 3

will attempt to read data from the files

	./test/dataDir/example_ccspectra.job3.trial1.dat
	./test/dataDir/example_sensints.job3.trial1.dat
	./test/dataDir/example_ccstats.job3.trial1.dat

--------
RESULTS
--------

stochCheck.py prints to screen the results of the first three sanity checks above. Additionally, it creates a new
directory ./stochCheck_figures/ to store the plots discussed in sanity check #4.

===========================================================
"""

#Object to hold segment data
class DataSegment():

	"""
	=========================================
	An object to hold individual segment data produced by stochastic.m.

	Objects in main() below by reading in the  files:

		./example_ccspectra.jobX.trialY.dat 	(contains ccSpecRe, ccSpecIm)
		./example_sensints.jobX.trialY.dat 	(contains sInts)
		./example_ccstats.jobX.trialY.dat	(contins CC, Sigma)

	where X and Y are specified by the user. This object contains the following instance variables:

	N		#Number of frequency bins.
	startTime	#Segment start time.
	CC		#The broadband cross-correlation (CC) statistic, reported in _ccstats.
	theorSigma	#Standard deviation of the broadband CC statistic, as reported in _ccstats.
	fq		#An array containing the frequency bins at which the narrowband statistics are evaluated
	ccSpecRe	#Array containing the real part of weighted, narrowband CC statistic, as reported in _ccspectra.
			#NOTE: This is NOT identically equal to the standard narrowband CC statistic, but is instead
			#weighted such that twice the integral over this statistc yields the broadband CC value
	ccSpecIm	#Array containing the imaginary part of weighted, narrowband CC statistic. See the comment above.
	sInts		#Array containing sensitivity integrand values, as repored in _sensints. These are defined such that
			#the broadband \sigma^2, the variance of the CC statistic, is given by the reciprocal of the
			#integral over senstivity integrands.

	...and the following methods:

	ccstats_CCStatistic()	#Returns the broadband CC statistic reported in _ccstats
	ccstats_CCSigma()	#Returns the broadband uncertainty (theorSigma) on the CC statistic, as reported in _ccstats
	freq()			#Returns an array of frequencies
	specRe()		#Returns the array ccSpecRe, containing the real, weighted narrowband CC values
	specIm()		#Returns the array ccSpecIm, containing the imaginary, weighted narrowband CC values
	sensInts()		#Returns the array sInts, containing the sensitivity integrands
	getCCStatistic1()	#Directly compute broadband CC statistic using weighted narrowband data (ccSpecRe)
	getCCSigma1()		#Directly compute broadband sigma using sensitivity integrands (sInts)
	getNarrowbandSigmas()	#Recover narrowband sigmas \sigma(f) from the sensitivity integrands
	getNarrowbandCCStats()	#Recover true narrowband CC statistics Y(f) from ccSpecRe and the sensitivity integrands
	getCCStatistic2()	#Compute the broadband CC statistic using the 'true' narrowband \sigma(f) and Y(f) values
	getCCSigma2()		#Compute the broadband sigma using the narrowband \sigma(f)

	=========================================
	"""

	#Default values
	N=0
	startTime=0
	CC=0
	theorSigma=0
	fq=[]
	ccSpecRe=[]
	ccSpecIm=[]
	sInts=[]

	#Initialize & populate object
	def __init__(self,t,sigma,cc,f,Re,Im,SI):
		self.N = len(f)
		self.startTime = t
		self.theorSigma = sigma
		self.CC = cc
		self.fq = np.array(f)
		self.ccSpecRe = np.array(Re)
		self.ccSpecIm = np.array(Im)
		self.sInts = np.array(SI)

	#Return CC Statistic from _ccstats
	def ccstats_CCStatistic(self):
		return self.CC

	#Return CC sigma from _ccstats
	def ccstats_CCSigma(self):
		return self.theorSigma

	#Return array of frequencies
	def freq(self):
		return self.fq

	#Return array of REAL weighted narrowband CC values.
	def specRe(self):
		return self.ccSpecRe

	#Return array of IMAGINARY weighted narrowband CC values.
	def specIm(self):
		return self.ccSpecIm

	#Return sensitivity integrands, as reported in _sensints.
	def sensInts(self):
		return self.sInts


	"""
	The following methods directly compute the broadband CC statistic and its standard deviation
	using the weighted CC spectra and sensitivity integrands reported in the files _ccspecta and
	_sensints.
	"""

	#Directly compute the broadband CC statistic using the weighted CC spectrum.
	def getCCStatistic1(self):
		df = self.fq[1]-self.fq[0]
		return np.sum(2.*df*self.ccSpecRe)

	#Directly compute variance
	def getThSigma1(self):
		df = self.fq[1]-self.fq[0]
		return np.power(np.sum(df*self.sInts),-0.5)

	"""
	The following four methods provide an alternatve means of calculating the broadband CC statistic
	and its std. deviation by first recovering the "true" unweighted narrowband CC values and their
	variances. Broadband values are then obtained by

	Y_broadband = SUM(Y_i*\sigma_i^(-2)) / SUM(\sigma_i^(-2))
	\sigma_broadband^(-2) = SUM( \sigma_i^(-2) )
	"""

	#Return array of /true/ narrowband sigmas
	def getNarrowbandSigmas(self):
		df = self.fq[1]-self.fq[0]
		return np.array(np.power(self.sInts*df,-0.5))

	#Return array of /true/ unweighted (real) narrowband CC values
	def getNarrowbandCCStats(self):
		broadSigma = self.getThSigma2()
		return np.array((2.*self.ccSpecRe)/(np.power(broadSigma,2.)*self.sInts))

	#Combine /true/ variances to compute broadband variance
	def getThSigma2(self):
		return np.power(np.sum(np.power(self.getNarrowbandSigmas(),-2.0)),-0.5)

	#Combine /true/ variances and CC values to compute broadband CC statistic
	def getCCStatistic2(self):
		sigmas=self.getNarrowbandSigmas()
		return np.sum(self.getNarrowbandCCStats()*np.power(sigmas,-2.0))/np.sum(np.power(sigmas,-2.0))


def main(dir,job):

	"""
	========================================================
	main() accepts two arguments:

		dir	#Directory in which python should look
		job	#A stochastic.m job number

	Together, these should specify a set of files ./example_ccspectra.job{JOB}.trial1.dat,
	./example_sensints.job{JOB}.trial1.dat, and ./example_ccstats.job{JOB}.trial1.dat which will be read.

	The data read from these files is organized into DataSegment() objects, and a number of sanity checks
	are performed for each segment.

	========================================================
	"""

	print "\nOpening job #"+job+" in "+dir

	#Read in output data from stochastic.m.
	file1 = dir+'/example_ccspectra.job'+job+'.trial1.dat'
	file2 = dir+'/example_sensints.job'+job+'.trial1.dat'
	file3 = dir+'/example_ccstats.job'+job+'.trial1.dat'

	if (os.path.isfile(file1) == False):
		print "Cannot find",file1
		return

	if (os.path.isfile(file2) == False):
		print "Cannot find",file2
		return

	if (os.path.isfile(file3) == False):
		print "Cannot find",file3
		return

	(sec,freq,ccSpecRe,ccSpecIm) = np.loadtxt(file1,unpack=True,skiprows=3,usecols=(0,2,3,4))
	sensInt = np.loadtxt(file2,skiprows=3,usecols=(3,))
	(crossCorr,thSigma) = np.loadtxt(file3,unpack=True,skiprows=3,usecols=(1,2))

	#Assuming equal segment length, get the length and number of segments.
	segmentLength = len(np.where(sec==sec[0])[0])
	segmentNumber = len(sec)/segmentLength

	#For each segment, create a DataSegment object. Store these objects in the array segments[]
	segments=[DataSegment(sec[i*segmentLength],thSigma[i],crossCorr[i],freq[i*segmentLength:(i+1)*segmentLength],
			ccSpecRe[i*segmentLength:(i+1)*segmentLength],ccSpecIm[i*segmentLength:(i+1)*segmentLength],
			sensInt[i*segmentLength:(i+1)*segmentLength]) for i in range(segmentNumber)]

	#SANITY CHECKS

	#CHECK 1: Print the broadband cross-correlation statistics, as reported by _ccstats and computed
	#using both strategies described above
	print '\nBroad-band Cross-Correlation Statistics:'
	print 'Segment #','\t','"Weighted" Computation','\t','Direct Computation','\t','_ccstats'
	for i,segment in enumerate(segments):
		print "{0} \t\t {1:.7g} \t\t {2:.7g} \t\t {3:.7g}".format(i,segment.getCCStatistic1(),segment.getCCStatistic2(), \
			segment.ccstats_CCStatistic())

	#CHECK 2: Similarly print the broadband sigmas
	print '\nBroad-band Sigmas:'
	print 'Segment #','\t','"Weighted" Computation','\t','Direct Computation','\t','_ccstats'
	for i,segment in enumerate(segments):
		print "{0} \t\t {1:.7g} \t\t {2:.7g} \t\t {3:.7g}".format(i,segment.getThSigma1(),segment.getThSigma2(), \
			segment.ccstats_CCSigma())

	#CHECK 3: Print the means and standard deviations of the narrowband values for Y/sigma.
	meanYonSigma=[]
	stdYonSigma=[]
	print '\nStatistics on Y(f)/sigma(f):'
	print 'Segment #','\t','Mean','\t\t','Std. Deviation'
	for i,segment in enumerate(segments):
		meanYonSigma.append(np.mean(segment.getNarrowbandCCStats()/segment.getNarrowbandSigmas()))
		stdYonSigma.append(np.std(segment.getNarrowbandCCStats()/segment.getNarrowbandSigmas()))
		print "{0} \t\t {1:.4g} \t\t {2:.4g}".format(i,meanYonSigma[i],stdYonSigma[i])
	print ""

	#Create directory to hold figures
	directory = "./stochCheck_figures/"
	if not os.path.exists(directory):
		os.mkdir(directory)

	#CHECK 4: For each segment, create a figure showing Y(f)/sigma(f) and the noise spectrum
	for i,segment in enumerate(segments):
		fig=plt.figure(figsize=(8,3))

		#Plot Y(f)/sigma(f)
		ax=fig.add_subplot(1,2,1)
		plt.hist(segment.getNarrowbandCCStats()/segment.getNarrowbandSigmas(),normed=1,bins=15)
		plt.xlabel(r"$Y(f)/\sigma(f)$")
		plt.title(r"$\mathrm{Histogram\,of\,}Y(f)/\sigma(f)$")

		#Plot sigma(f)
		ax=fig.add_subplot(1,2,2)
		plt.plot(segment.freq(),segment.getNarrowbandSigmas())
		plt.xlim(segment.freq()[0],segment.freq()[-1])
		plt.yscale('log')
		plt.xscale('log')
		plt.ylabel(r"$\sigma(f)$")
		plt.xlabel(r"$f\,(\mathrm{Hz})$")
		plt.title(r"$\mathrm{Noise\,Spectrum}$")

		#plt.tight_layout()
		plt.savefig(directory+"job"+job+"_segment"+str(i)+".ps")

if __name__ == "__main__":

	#If the number of command line arguments is less than or greater than 2, complain and quit
	if (len(sys.argv) < 3):
		print "\nPlease pass a directory and job #. For instance:"
		print ">> python stochCheck.py ./test/directory/ 2\n"

	elif (len(sys.argv) > 3):
		print "\nToo many arguments. Please pass only a directory and job #. For instance:"
		print ">> python stochCheck.py ./test/directory/ 2\n"

	#Otherwise, proceed to main!
	else:
		main(sys.argv[1],sys.argv[2])


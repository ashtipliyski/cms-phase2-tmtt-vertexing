------------------------------------------------------------------------------------
* How to setup the framework and how to run our analysis code:

=== NOTE ===
The software described below is for running the TMT tracking in CMSSW8
and mainly for the tilted barrel geometry.  This differs from the software used
for the demonstrator and December review, which runs in CMSSW6 and for the flat barrel tracker
geometry.  You can find that version of the software here : https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/software/cmssw/trunkSimpleCode/
============

P.S. You can browse the software in https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/software/cmssw/branch/ejclemen_81X/ 

- Setup a CMSSW environment
cmsrel CMSSW_9_1_0_pre1
cd CMSSW_9_1_0_pre1/src
cmsenv

- Software checkout
Copy official code to your svn branch, and then check it out, so you can play with it.

svn copy svn+ssh://svn.cern.ch/reps/UK-TrackTrig/software/cmssw/branch/ejclemen_91X/ svn+ssh://svn.cern.ch/reps/UK-TrackTrig/software/cmssw/branch/{YOURBRANCH} -m "Creating xxx branch of tag:xxx"
# Note the "." at the end of the next line. You need it!
svn co svn+ssh://svn.cern.ch/reps/UK-TrackTrig/software/cmssw/branch/{YOURBRANCH} .

If later on, you want to update the code in your working area with any changes made to Ian_SimpleCode since you checked it out, do:
 a) Backup your code using "cp -pr src/ ~/srcBackup" in your top-level CMSSW directory. 
 b) "svn merge svn+ssh://svn.cern.ch/reps/UK-TrackTrig/software/cmssw/branch/ejclemen_81X/" (to update the files on your disk with changes from the main branch).
 c) Check that there are no SVN conflict messages, and that you understand any differences between your code
after this merge and your backed up code.
 d) If all OK, then "svn commit" to put this updated code into your SVN branch.

- MC samples 
  Some samples are available in Samples81X/
  These are RelVal samples generated in CMSSW_8_1_0_pre16 and do not contain the stubs.
  Please see instructions below on how to produce the stubs on the fly, before running our tracking.
  Alternatively, a subset of the RelVal samples (~1000 events) containing stubs are available privately (ask Emyr)

- Change the configuration parameters:

    edit TMTrackTrigger/TMTrackFinder/python/TMTrackProducer_cff.py

- Compile:
    scram b -j8

- Run the code
    cd TMTrackTrigger/TMTrackFinder/test/
    cmsRun   tmtt_tf_analysis_cfg.py 
    
    or with options:

    cmsRun   tmtt_tf_analysis_cfg.py Events=50 inputMC=../../../Samples81X/tilted_RelVal_ttbar_PU200.txt histFile=outputHistFile.root makeStubs=1

- Look at the printout from the job. At the end, it prints the number of track candidates reconstructed
  and the algorithmic tracking efficiency.

- Look at the analysis histograms (which are explained in '(11) Class "Histos"' below).
    root Hist.root
    TBrowser b

-------------

=== Producing stubs ===:

The RelVal samples (at least from CMSSW_9_0_0_pre5) contain stubs.

To reproduce the stubs, you can do:

cmsRun tmtt_tf_analysis_cfg.py makeStubs=1

This will run the CMSSW modules that produce the stubs (and truth association) before running our tracking algorithms.
Note this will take a long time to run (~100 events per hour for ttbar+200PU), so is only suitable for
debugging on a small number of events.

To produce the stubs, add them to the event content, and output the new collections to file
(along with the existing event content), you can run make_stubs.py:

cmsRun make_stubs.py inputMC=../../../Samples81X/tilted_RelVal_ttbar_PU200.txt Events=50

This will produce an output file output.root, which can then be read by tmtt_tf_analysis_cfg without having to remake the stubs.
i.e. you can then run with (or omit the makeStubs argument as it's False by default):

cmsRun tmtt_tf_analysis_cfg.py makeStubs=0

=== Changing configuration options ===:

a) To change the MC dataset or number of events you run on, either edit tmtt_tf_analysis_cfg.py,
or specify them as arguments after the "cmsRun tmtt_tf_analysis_cfg.py" command.

b) Files TMTrackProducer_Demo_DaisyChain_cff.py, TMTrackProducer_Demo_Thomas_cff.py & 
TMTrackProducer_Demo_Systolic_cff.py, which are all in TMTrackTrigger/TMTrackFinder/python/,  set the 
configuration parameters to match those expected for the May 2016 hardware demonstrator for the Thomas Schuh 
"daisy chain", Thomas Schuh "2-c-bin2, & Systolic Array firmware implementations, respectively. 
These all run the r-phi Hough transform, but disable r-z filters, track fitting etc.

c) File TMTrackTrigger/TMTrackFinder/python/TMTrackProducer_cff.py configures things as they
should be for the Nov. 2016 demonstrator. It includes not only the r-phi HT, but also r-z filters & 
track fitting.

d) All these files include and then override where necessary the default config parameters, which are 
specifiied in TMTrackTrigger/TMTrackFinder/python/TMTrackProducer_Defaults_cfi.py. The comments in that file 
explain what all the config parameters do. 

-------------

=== Software structure ===:

1) Class "TMTrackProducer" -- This is the main routine, which uses classes: "InputData" to unpack the useful
data from the MC dataset, and "Sector" & "HTPair" to do the L1 Hough transform track-finding.
It creates matrices of "Sector" and "HTpair", where the matrix elements each correspond to a different 
(phi,eta) sector. It then uses "TrackFitGeneric" to do the track fitting and optionally "KillDupFitTrks"
to remove duplicate tracks after the fit. It employs "Histos" to create the analysis histograms. 
   To allow comparison of our tracks with those of the AM & Tracklet groups, it also converts our tracks
to the agree common TTTrack format, with this conversion done by the "ConverterToTTTrack" class.

2) Class "InputData" -- This unpacks the most useful information from the Stubs and Tracking Particle 
(truth particles) collections in the MC dataset, and it for convenient access in the "Stub" and "TP"
classes. The "Stub" class uses a class called "DigitalStub" to digitize and then undigitize again
the stub data. This process degrades slightly the resolution, as would happen in the real firmware.
The digitisation is optional. It is called from TMTrackProducer after the stubs have been assigned to
sectors.

3) Class "Sector" -- This knows about the division of the Tracker into (phi,eta) sectors that we use
for the L1 tracking, and decides which of these sectors each stub belongs to.

4) Class "HTpair" -- This does L1 tracking in a single (phi,eta) sector using a pair of Hough transforms:
one in the r-phi plane (implemented with class "HTrphi") and one in the r-z plane (implented with class 
"HTrz"). It runs the r-phi one first, then for each individual track candidate found, it passes the 
associated stubs to the r-z Hough transform, to check that they form a track in that plane too. If they do,
then the information from both Hough transforms is combined to produce 3D tracks. These are stored using 
the "L1track3D" class. 

The r-z Hough transform can optionally be disabled via configuration parameter, in which case, 3D tracks 
are still produced, but in this case, the track r-z helix parameters are set to nominal values.

HTpair can optionally also run r-z track filters after the r-phi Hough transform to clean up the tracks
(by calling class TrkRZfilter).

HTpair can also run duplicate track removal algorithms on the 3D tracks it produces (by calling class
KillDupTrks).

5) Classes "HTrphi" and "HTrz" implement the Hough transforms in the r-phi and r-z planes. They both 
inherit from a base class "HTbase". The Hough transform array is implemented as a matrix of "HTcell" 
objects. The two HT classes both store tracks they find using the "L1track2D" class. They optionally 
use class "KillDupTrks" to attempt to eliminate duplicate tracks. And optionally, class "MuxHToutputs"
can be used to multiplex the tracks found in different HT arrays (sectors) onto a single output
optical link pair.

6) Class "HTcell" -- This represents a single cell in an HT array. It provides functions allowing stubs
to be added to this cell, to check if the stubs in a cell give a good track candidate, and to check
if this matches a tracking particle (truth).

7) Class "L1track2D" represents a 2D track, reconstructed in the r-phi or r-z plane by a Hough transform.
Class "L1track3D" represents a 3D tracks, obtained by combining the information in the 2D tracks
from r-phi and r-z Hough transforms. These classes give access to the stubs associated to each track,
to the reconstructed helix parameters, and to the associated truth particle (if any). They represent
the result of the track finding. Both inherit from a pure virtual class L1trackBase, which contains
no functionality but imposes common function names.

8) Class "KillDupTrks" contains algorithms for killing duplicate tracks found within a single
HT array. Class "KillDupFitTrks" contains algorithms for killing duplicate fitted tracks.

9) Class "TrkRZfilter" contains r-z track filters (and any other filters that must be run after the
Hough transform, because they take too many resources to run inside it).

10) Class "TrackFitGeneric" does one (or more) helix fit(s) to the track candidates, using various
other classes that implement linearized chi2, linear regression or Kalman filter fits. These are:

   - TrackFitLinearAlgo (chi2 linear fit, using code from tracklet group)
   - ChiSquared4ParamsTrackletStyle (chi2 linear fit, written by us)
   - ChiSquared4ParamsApprox (chi2 linear fit, rewritten by us, with maths simplified for easier use in FPGA)
   - KF4ParamsCombIV (Kalman filter fit)
   - globalLinearRegression (linear regression track fit: a simple fit which essentially treats all hits as having the same uncertainty)
   - LinearRegression (older version of code that still gives marginally better performance than globalLinearRegression if you are using the seed filter)

The fit also uses a couple of dedicated utility classes (Matrix & kalmanState & StubCluster).

11) Class "L1fittedTrack" contains the result of running a track fitter (the algorithm for which is 
implemented in class "TrackFitAlgo") on the L1track3D track candidate found by the Hough transform. 
It gives access to the fitted track parameters and chi2, and via a link to the L1track3D candidate 
that produced it, also to the stubs on the track and the associated truth particle (if any). 
It inherits from the pure virutal L1trackBase, ensuring it has some common classes with L1track3D and 
L1track2D.

12) Class "L1fittedTrk4and5" contains a pair of L1fittedTrack objects, containg the result of doing
either a 4 or 5 parameter helix fit to the track, where the former assumes d0 = 0.

13) Class "DataCorrection" -- This is used by class "Stub" to degrade the resolution on the stub
bend information to that expected in the electronics, as opposed to that currently in CMSSW.

14) "Utility" -- contains a few useful functions, that are not part of a class.

15) Class "Settings" -- Reads in the configuration parameters.

16) Class "Histos" -- Books and fills all the histograms. There are several categories of histograms,
with each category being booked/filled by its own function inside "Histos", and being placed inside its
own ROOT directory in the output histogram file. The categories are "InputData" = plots made with the 
Stubs & Tracking Particles; "CheckEtaPhiSectors" = plots checking assignment of stubs to (eta,phi) 
sectors; "HT" = plots checking how stubs are stored in the Hough Transform arrays; "TrackCands" = plots 
showing number of track candidates found & investigating why tracking sometimes failed, 
"EffiAndFakeRate" = plots of tracking efficiency. 

Each user of the code will probably want to book their own set of histograms inside "Histos". So 
just consider the version of this class in SVN as a set of examples of how to do things. Don't feel
obliged to understand what every histogram does.

17) Class "DeadModuleDB" is used both to emulate dead modules by killing stubs in certain tracker
regions, and to recover efficiency caused by dead modules indicating in which sectors looser
track-finding cuts are required.

18) SimTracker/TrackTriggerAssociation/ contains a modification to the official L1 track to TrackingParticle
matching software used by Louise Skinnari's official L1 track analysis code. This modification (made by 
Seb Viret) to TTTrackAssociator.h allows one incorrect hit on L1 tracks, whereas the original matching code
allowed none.

------

To update the dOxygen documentation, just type "doxygen". This creates the web page inside html/
A recent version of this documentation is in http://tomalini.web.cern.ch/tomalini/IanSimpleCode/hierarchy.html .





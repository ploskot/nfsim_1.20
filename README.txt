%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%     @@    @  @@@@@                %
%     @ @   @  @                    %
%     @  @  @  @@@@  ___            %
%     @   @ @  @    /__  | |\ /|    %
%     @    @@  @    ___\ | | v |    %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NFsim - the network free stochastic simulator, v1.20

Pavel Loskot
ZJU-UIUC Institute, Haining, China

This is an updated version of the NFsim simulator version 1.11 from October
2012. It is distributed under the same MIT license as released by the
original contributors in 2016.


################################################################################

Release Notes

Planned Todo's

	(a) Remove the warnings remaining.
	(b) Clean-up source code, remove unused code, refactor some parts.

v.1.20  June, 2022

	(a) The source code was updated to compile under more recent C++
	versions (gcc -std=c++20). Most but not all warnings were removed.
	(b) The file bin/Makefile was created to allow incremental
	recompiling of the source code.
	(c) Added options '-eh' and '-ehb' to record the randomly generated
	reaction events into text file '<output>_eh.gdat' or binary file
	'<output>_eh.dat', respectively.
	(d) The sample biochemical models are stored in models/ folder.
	They are bngl text files processed by the BioNetGen Perl scripts,
	which calls NFsim with appropriate parameters. The models were
	tested for BioNetGen version 2.7.0.


/ original contributors and the description below this line /
################################################################################

NFsim - the network free stochastic simulator, v1.11

michael w. sneddon
justin s. hogg
james r. faeder
thierry emonet

Yale University
University of Pittsburgh
funded by the National Science Foundation


################################################################################

NFsim is a free, open-source, biochemical reaction simulator designed to handle systems 
that have a large or even infinite number of possible molecular interactions or states. 
NFsim also has advanced and flexible options for simulating coarse-grained representations 
of complex nonlinear reaction mechanisms.

NFsim is ideal for modeling polymerization, aggregation, and cooperative reactions that 
cannot be handled with traditional stochastic or ODE simulators. Models are specified in 
the BioNetGen Langauge, providing a powerful model building environment.

If you just want to download and use NFsim, you should simply download a preconfigured
packaged release from http://emonet.biology.yale.edu/nfsim.  If you want to hack on the
code or make contributions, please create a fork and submit pull requests to the dev
branch.

If you use NFsim for your research or work, please cite NFsim as:
Sneddon MW, Faeder JR & Emonet T.  Efficient modeling, simulation and 
coarse-graining of biological complexity with NFsim.  Nature Methods,(2011) 8(2):177-83.

################################################################################

NFsim is released under the MIT License.  See LICENSE.txt for more details.  The git 
repository is hosted on github at: https://github.com/msneddon/nfsim
  
For help with running NFsim, see the user manual, NFsim_manual_[version].pdf,
and open the example model "simple_system.bngl".

Executable files for Windows, Mac and Linux are in the "bin" directory.  Source
code and makefiles for NFsim are in the NFcode directory.  Example models are
in the "models" directory, with README files.  BioNetGen and the ODE and SSA
solvers used by BioNetGen are in the BNG and Network2 directories.  Finally,
a suite of helpful analysis and other modeling tools can be found in the
"NFtools" directory.

Enjoy your new network-free world!

/ release notes for the previous versions were cut below this line /
################################################################################

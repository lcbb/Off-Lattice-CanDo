
=====================================
=                                   =
=   OFF-LATTICE CANDO VERSION 1.1   =
=                                   =
=====================================

Version 1.0 - Base version of Off-Lattice CANDO (Updated March 9th, 2018)

Version 1.1 - Added specification of CNDO topology file inputs 
(Updated March 23rd, 2018)

Authorship: See literature below
Other contributors: Hyungmin Jun

License: The off-lattice CanDo program is free to use with adherance to the 
license terms included (LICENSE.txt).

=========================

Description: Off-lattice CanDo version 1.0 solves the 3D solution shape and 
flexibility of non-scaffolded DNA origami objects. 

=========================

Users of this program are requested to cite the following two references:

P1. K Pan, DN Kim, F Zhang, MR Adendorff, H Yan, M Bathe. "Lattice-free 
	prediction of three-dimensional structure of programmed DNA assemblies." 
	Nat. Commun. 5, 5578 (2014).
	
P2. K Pan, WP Bricker, S Ratanalert, M Bathe. "Structure and conformational 
	dynamics of scaffolded DNA origami nanoparticles." Nucleic Acids Res. 45, 
	6284-6298 (2017).
	
=========================

A. INSTALLATION PROCEDURE

=========================

1. External software requirements

	a. MATLAB and MATLAB Toolkit 
	   (https://www.mathworks.com/products/matlab.html)
	
	b. ADINA finite element software (http://www.adina.com/)
	
	c. UCSF Chimera (https://www.cgl.ucsf.edu/chimera/) 
	
2. Modify personal paths to ADINA and Chimera in 'main_offlattice.m':

    default locations:

	On 32-bit Windows (node-locked or floating license):
		path.ADINA_AUI = '"C:\ADINA93\x32\aui.exe"';
		path.ADINA     = '"C:\ADINA93\x32\adina.exe"';
	
	On 64-bit Windows (node-locked license):
		path.ADINA_AUI = '"C:\ADINA93\x64\aui.exe"';
		path.ADINA     = '"C:\ADINA93\x64\adina.exe"';
	
	On 64-bit Windows (floating license):
		path.ADINA_AUI = '"C:\ADINA93\bin\aui.exe"';
		path.ADINA     = '"C:\ADINA93\bin\adina.exe"';

	path.CHIMERA   = '"C:\Program Files\Chimera 1.10.2\bin\chimera.exe"';
	
3. Modify paths to 'input' and 'output' folders as necessary, also in 
   'main_offlattice.m':

	default locations:
	
	inputDIR = 'input';
	outputDIR = 'output';
	
=========================

B. HOW TO RUN OFF-LATTICE CANDO

=========================

1. Input requirements to solve the 3D solution of a DNA nanoparticle:

	a. A topology file (structure.dat) converted from Tiamat or in the CNDO
	   file format (structure.cndo).
	
	b. A sequence file (structure.txt) with basepairing information
	   (not required unless generating an atomic model). CNDO topology format 
	   already includes sequence information, but Tiamat topology format does 
	   not.
	   
	    - Note that an atomic model can be generated without this sequence 
		file. In this case a randomly generated sequence will be used (Tiamat
		topology format only).
		 
	Please visit http://cando-dna-origami.org/tutorial/ for more information 
	regarding the preparation of input files for off-lattice CanDo.

2a. Convert Tiamat input (.dna) to CanDo input (.dat):

		Use the included executable 'Tiamat2dat.exe' to convert the Tiamat 
		input (.dna) file to a CanDo input (.dat) file. 
	
	a. Open a Windows command prompt window
	
	b. Go to the folder containing your structure.dna file
	
	c. Type the following command and hit Enter:
	   
	   echo structure.dna | Tiamat2dat.exe
	
	d. A new file is generated (dna_ascii.dat) which can be renamed
	
	Please visit http://cando-dna-origami.org/tutorial/ for more information 
	regarding the preparation of input files for off-lattice CanDo.

2b. With CNDO topology format (.cndo), no conversion process necessary.
	
3. Place any input files into the off-lattice CanDo 'inputDIR' as specified in
   Section A3. This includes both topology (structure.dat, structure.cndo) and 
   sequence (structure.txt) files.
   
4. Modify any CanDo adjustable parameters:

	In Section C of 'main_offlattice.m', many adjustable parameters for CanDo
	execution can be set. We highly recommend reading the two papers above 
	before modifying these values or running a simulation.
	
	*** PLEASE READ THIS -- IMPORTANT! ***
	a. Of particular importance is the 'angleHJ' parameter, which sets the 
	   equilibrium value for the Holliday junction twisting angles of DX 
	   crossovers. While the default value is set at 60 degrees, this is not
	   necessarily a useful value for all simulations. For example, in paper P1 
	   above, a value of 60 degrees was used for concentric rings and lattices. 
	   On the other hand, in paper P2, a value of 0 degrees was used for the
	   polyhedral nanoparticles. This value should reflect the expected folding
	   behavior of your nanoparticle design.
	
	*** Parameters for DNA Geometry ***
	b. The default parameters for DNA geometry are those of average B-form DNA.
	   Be careful when modifying these default values.
	
	*** FE Model Resolution ***
	c. The parameter 'jobInfo.model' can be set to 'coarse' or 'fine' (default) 
	   depending on the FE resolution required for your submission. Setting 
	   this value to 'fine' is mandatory for atomic model generation.
	
	*** Atomic Model Generation ***
	d. The parameter 'jobInfo.atomic' can be set to 'true' or 'false' (default)
	   depending on whether you would like to generate an atomic model (.pdb) 
	   of the 3D solution shape. Requires the resolution above to be set to 
	   'fine'.
	
	*** Normal Mode Analysis ***
	e. The parameter 'jobInfo.NMA' can be set to 'true' or 'false' (default)
	   depending on whether you would like to calculate the normal modes of 
	   your 3D solution shape.
	
5. Run off-lattice CanDo simulations:

	Any input files in the 'inputDIR' can be run by the MATLAB script
	'main_offlattice.m'. Multiple inputs can be run sequentially.

	Type the following command into the MATLAB window:
	
	If you have one submission:
		
		main_offlattice({structure},model,atomic)
		
		Here, 'model' is the parameter 'jobInfo.model' above. This should be 
		set to 'coarse' or 'fine'. Also, 'atomic' is the parameter 
		'jobInfo.atomic' above. This should be set to 'true' or 'false'.
	
	If you have multiple submissions:
		
		main_offlattice({structure1,structure2,...,structureN},model,atomic)
	
		See above for setting 'model' and 'atomic' parameters.
	
	For each submission, only the prefix of the input file is required. 
	For example, if your input file is named structure.dat, the input to 
	'main_offlattice' is 'structure'. If your input file is named 
	structure.cndo, the input to 'main_offlattice' is also 'structure'. You 
	may run both Tiamat and CNDO inputs in the same submission batch.
		
=========================

C. OUTPUT OF OFF-LATTICE CANDO SIMULATION

=========================
	
Simulation and output files of the off-lattice CanDo simulations are
re-directed to the 'outputDIR' specified in Section A3. Output files 
include:
   
1. ADINA input and output files (.inp,.por,.dat,...,etc.)

2. MATLAB files (.mat)

3. Images of FE solution (.png) generated by Chimera

4. Atomic structures (.pdb) 

5. Input files for Chimera (.bild)
	
=========================

D. IMPORTANT NOTES AND CONTACT INFORMATION

=========================

*** March 9th, 2018 ***
This version of off-lattice CanDo was tested on a Windows 10 64-bit platform 
using a full MATLAB 2015b installation, ADINA version 9.3.4, and Chimera 
version 10.02.
	
*** Questions ***

Please direct any questions regarding use of this software to Mark Bathe 
(http://lcbb.mit.edu/people/mark-bathe/) at mbathe@mit.edu. 
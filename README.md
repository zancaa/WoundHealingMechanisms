# Choosing between models and boundary descriptions in wound healing applications

This code is used to generate the simulations corresponding to Chapter 5 of my PhD thesis. This code uses Chaste, details provided below. Code in this repository is modified from https://github.com/jmosborne/TissueBoundaries (code developed by Domenic Germano, James Osborne and myself for a paper submitted to the Bulletin of Mathematical Biology).

## Chaste
The code in this repository should run using the core version of Chaste, which can be found [here](https://chaste.cs.ox.ac.uk/trac/wiki), along with documentation on how to install and run simulations. Note that Chaste is only fully supported on Linux/Unix systems, Windows and Mac OS users are encouraged to use [Docker](https://github.com/Chaste/chaste-docker). 
 
To use the code in this repository, clone this repository into the `projects` subdirectory of the Chaste code. 

## Repository structure
- `src/` additional class files required to run the simulations using Chaste
	* `BoundaryNormalForce.hpp` to apply a crawling force to nodes on the wound boundary, the direction of the force is the average of the perpendicular normals to the edges connected to the node, pointing into the wound
	* `WoundCentreForce.hpp` to apply a crawling force to the nodes on the wound boundary, additional forces all point towards the centre of the wound (the location of the centre of the wound changes over time)
	* `NagaiHondaMutationCellForce.hpp` to handle the purse string simulations that treat the wound as a cell
	* `SimpleWoundMutantTargetAreaModifier.hpp` to modify the wound cell in purse string simulations
	* `WoundCellMutationState.hpp` to define the wound cell in purse string simulations
	* `BernoulliTrialWithContactInhibitionCellCycleModel.hpp` to implement a probabilistic cell cycle model with contact inhibition  
	* `VoidAreaModifier.hpp` to track the wound area
-`test/` test files for simulation setup
	* `TestCreateWound.hpp` simulations to create a wound in a tissue and let the tissue relax to equilibrium. Must be run before all other tests (that load the output from this test)
	* `TestCompression.hpp` simulations of tissue relaxation with different levels of cellular compression
	* `TestProliferation.hpp` simulations with contact inhibited proliferation
	* `TestEdgeContractility.hpp` simulations where cell edges adjacent to the wound contract quicker than the rest of the tissue
	* `TestPurseString.hpp` simulations where the perimeter of the wound contracts at a whole-wound scale
	* `TestLocalCrawling.hpp` simulations with the local crawling force applied
	* `TestGlobalCrawling.hpp` simulations with the global crawling force applied
-`figures/` MATLAB code to generate area over time curves (all other figures were generated using Paraview)
	* `WoundAreaPlots.m` to plot wound area over time curves

## Running simulations using Chaste
To run the code (compiling using scons), the create wound test must be run first:
```
~: cd [Chaste Directory]
~/Chaste: scons b=GccOpt projects/WoundHealingMechanisms/test/TestCreateWound.hpp &> build_output.log
```
After the create wound test has been successfully run, other tests can be run, for example:
```
cd [Chaste Directory]
~/Chaste: scons b=GccOpt co=1 projects/WoundHealingMechanisms/test/TestCompression.hpp &> build_output.log
~/Chaste: ./projects/WoundHealingMechanisms/build/optimised/TestCompressionRunner &> run_output.log
```
 
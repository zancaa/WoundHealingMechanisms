#ifndef TESTCOMPRESSION_HPP_
#define TESTCOMPRESSION_HPP_


#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 
#include "CheckpointArchiveTypes.hpp"

#include "SmartPointers.hpp"

#include "ToroidalHoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "NoCellCycleModel.hpp"
#include "BernoulliTrialWithContactInhibitionCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "WoundCellMutationState.hpp"

#include "VertexBasedCellPopulation.hpp"

#include "CellVolumesWriter.hpp"
#include "BoundaryNodeWriter.hpp"
#include "WoundAreaWriter.hpp"

#include "OffLatticeSimulation.hpp"

#include "NagaiHondaForce.hpp"
#include "NagaiHondaMutationCellForce.hpp"
#include "BoundaryNormalForce.hpp"
#include "WoundCentreForce.hpp"

#include "SimpleTargetAreaModifier.hpp"
#include "SimpleWoundMutantTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "VoidAreaModifier.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp" 
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

#include "BoundaryCellWriter.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_STEADY_STATE = 20.0;
static const double M_END_TIME = 20.0;
static const double M_DT_TIME = 0.000001;
static const double M_SAMPLE_TIME = 0.1/M_DT_TIME;

// Both Width and Length must be EVEN numbers here
static const double M_DOMAIN_WIDTH = 20;
static const double M_DOMAIN_LENGTH = 20;
static const double M_DOMAIN_SCALING = 0.8;
static const double M_PERIODIC_WIDTH = 20;//M_DOMAIN_WIDTH*M_DOMAIN_SCALING;
static const double M_PERIODIC_HEIGHT = 20.0*0.5*sqrt(3);//M_DOMAIN_LENGTH*0.5*sqrt(3)*M_DOMAIN_SCALING;

static const double M_HOLEWIDTH = 2;
static const double M_HOLE_X_MIN = 7.0;//M_DOMAIN_WIDTH/2 - 3;
static const double M_HOLE_X_MAX = 13.5;//M_DOMAIN_WIDTH/2 + 3;
static const double M_HOLE_Y_MIN = 7.0*0.5*sqrt(3.0);
static const double M_HOLE_Y_MAX = 13.0*0.5*sqrt(3.0);

static const std::string M_HEAD_FOLDER = "WoundMechanisms";


class TestInternalVoid : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestCompression()
    {
        // Simulation and sweep parameters
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_param_vals"));
        unsigned num_param_vals = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-num_param_vals").c_str());

        // Set min/max target areas (proxy for compression)
        double min_cell_target_area = 0.9;
        double max_cell_target_area = 1.0;

        // Loop over parameter values
        for(unsigned sim_index=1; sim_index < num_param_vals; sim_index++)
        {    
            std::cout << " Run number " << sim_index << "... \n" << std::flush;   
            // Reseed the random number generator
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            /* 
            * == VM ==
            * 
            * Simulation internal void using the
            * Cell Vertex model.
            */

            // Compression
            {
                std::string output_directory =  M_HEAD_FOLDER + "/Pre-void/Circle";
                // Load steady state
                OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
                VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Remove the forces and boundaries - redefined here
                p_simulator->RemoveAllForces();

                std::cout << "Compression \n" << std::flush; 
                double cell_target_area = min_cell_target_area + (max_cell_target_area - min_cell_target_area) * double(sim_index) / double(num_param_vals);
                std::stringstream paramAsString;
                paramAsString << cell_target_area;
                output_directory =  M_HEAD_FOLDER + "/Circle/Compression/Cell_target_area_" + paramAsString.str();

                /* 
                * == Post-void == 
                */

                // Loop over cells to change target areas. 
                for (typename VertexMesh<2,2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
                    elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
                    ++elem_iter)
                {
                    unsigned elem_index = elem_iter->GetIndex();
                    // Get cell associated with this element
                    CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_index);

                    // Set Target Area so dont need to use a growth model in vertex simulations
                    p_cell->GetCellData()->SetItem("target area", cell_target_area);
                }

                // Add target area modifier
                MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
                p_growth_modifier->SetReferenceTargetArea(cell_target_area);
                p_growth_modifier->SetGrowthDuration(0);
                p_simulator->AddSimulationModifier(p_growth_modifier);

                // Track the area of the void
                MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
                voidarea_modifier->SetOutputDirectory(output_directory);
                p_simulator->AddSimulationModifier(voidarea_modifier);

                // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
                MAKE_PTR(NagaiHondaForce<2>, p_force);
                p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
                p_simulator->AddForce(p_force);

                // Create simulation from cell population
                p_simulator->SetDt(M_DT_TIME);
                p_simulator->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                p_simulator->SetEndTime(M_END_TIME);
                p_simulator->SetOutputDirectory(output_directory);
                
                // Run simulation
                try
                {
                    p_simulator->Solve();
                }
                catch (Exception& e)
                {
                    // If it throws then we report the error message and go to the next simulation
                    WARNING("Compression simulation didnt run" << paramAsString.str() << ".");
                    WARNING(e.GetMessage());            
                }
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
                // Tidy up
                delete p_simulator;
            }
        }
    }
};

#endif /* TESTCOMPRESSION_HPP_ */
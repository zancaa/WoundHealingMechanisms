#ifndef TESTINTERNALVOID_HPP_
#define TESTINTERNALVOID_HPP_


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
static const double M_DT_TIME = 0.001;
static const double M_SAMPLE_TIME = 100;

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

static const std::string M_HEAD_FOLDER = "WoundMethodsTest";


class TestInternalVoid : public AbstractCellBasedWithTimingsTestSuite
{
private:
    /**
    * Helper method. Smooth out edges of a vertex mesh.
    * 
    * @param rCellPopulation a cell population
    */
    void SmoothVertexMeshEdges(AbstractCellPopulation<2>& rCellPopulation)
    {
        MutableVertexMesh<2, 2>& r_mesh = static_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->rGetMesh();

        for (VertexMesh<2,2>::NodeIterator node_iter = r_mesh.GetNodeIteratorBegin();
            node_iter != r_mesh.GetNodeIteratorEnd();
            ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if (containing_element_indices.size() == 1)
            {
                // Get this element
                unsigned elem_index = (*containing_element_indices.begin());

                VertexElement<2,2>* p_element = r_mesh.GetElement(elem_index);

                // Remove node from this element and delete the node
                p_element->DeleteNode(p_element->GetNodeLocalIndex(node_index));
                r_mesh.DeleteNodePriorToReMesh(node_index);
            }
        }
        r_mesh.ReMesh();
    }

    /**
    * Helper method. Iterate over all cells and define the 'hole' by
    * killing those cells whose centres are located in a given region.
    * 
    * @param rCellPopulation a cell population
    */
    void CreateHoleInCellPopulation(AbstractCellPopulation<2>& rCellPopulation)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            // Get the coordinates of this cell centre
            c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell[0];
            double y = centre_of_cell[1];

            if ((fabs((2.0/sqrt(3.0))*y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
            {   
                cell_iter->Kill();
            }
        }
        
        /* Need to remove the deleted cells and call update note this is usually
        * performed in the Solve() method of the simulation class.
        */
        if (bool(dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation)))
        {
            dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
            dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->Update();
        }
    }

    /**
    * Helper method. Iterate over all nodes and define the 'hole' by
    * deleting nodes who are located in a given region and replacing them with 
    * a 'wound cell'.
    * 
    * @param rMesh a vertex mesh
    */
    void CreateWoundMesh(MutableVertexMesh<2,2>& rMesh)
    {
        for (VertexMesh<2, 2>::VertexElementIterator elem_iter = rMesh.GetElementIteratorBegin();
                elem_iter != rMesh.GetElementIteratorEnd();
                ++elem_iter)
        {
            // Get the coordinates of the centroid of the element
            unsigned elem_index = elem_iter->GetIndex();
            c_vector<double, 2> centroid_of_element = rMesh.GetCentroidOfElement(elem_index);
            double x = centroid_of_element[0];
            double y = centroid_of_element[1];

            if ((fabs((2.0/sqrt(3.0))*y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
            {   
                rMesh.DeleteElementPriorToReMesh(elem_index);
            }
        }
        
        /* Need to remove the deleted elements and ReMesh
        */
        rMesh.ReMesh();

        // Create the vector of nodes that defines the wound element in anti-clockwise order.
        /* For some reason, nodes aren't considered boundary nodes, so loop over nodes until find one 
         * only contained in one or two elements.
         */
        std::vector<unsigned> wound_node_indices;
        std::vector<Node<2>*> wound_nodes;
        for (VertexMesh<2,2>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
            node_iter != rMesh.GetNodeIteratorEnd();
            ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            c_vector<double,2> node_location = node_iter->rGetLocation();
            double x = node_location[0];
            double y = node_location[1];
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if ( (containing_element_indices.size() == 1 || containing_element_indices.size() == 2) &&
                    (x > 0.75) && (x < (M_PERIODIC_WIDTH - 0.75)) && (y > 0.75) && (y < (M_PERIODIC_HEIGHT - 0.75)) )
            {
                Node<2>* p_temp_node = rMesh.GetNode(node_index);
                wound_nodes.push_back(p_temp_node);
                wound_node_indices.push_back(node_index);
                break;
            }
        }

        bool repeat = false;
        unsigned node_index = wound_node_indices[0];
        while (!repeat)
        {
            std::set<unsigned> containing_elements = rMesh.GetNode(node_index)->rGetContainingElementIndices();
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                iter != containing_elements.end();
                iter++)
            {
                VertexElement<2, 2>* p_element = rMesh.GetElement(*iter);
                // Find the local index of this node in this element
                unsigned local_index = p_element->GetNodeLocalIndex(node_index);
                // Get the previous node in this element
                unsigned num_nodes_elem = p_element->GetNumNodes();
                unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                Node<2>* p_previous_node = p_element->GetNode(previous_node_local_index);
                /* If the previous (or clockwise) node is only in one or two elements, then it is the 
                 * next anti-clockwise node in the wound element.
                 */ 
                std::set<unsigned> prev_node_containing_elements = p_previous_node->rGetContainingElementIndices();
                if ( (prev_node_containing_elements.size() == 1 || prev_node_containing_elements.size() == 2) )
                {
                    node_index = p_previous_node->GetIndex();
                    // Vector should be ordered such that there are no repeats until we reach the first node again
                    if (node_index == wound_node_indices[0])
                    {
                        repeat = true;
                        break;
                    }
                    Node<2>* p_temp_node = rMesh.GetNode(node_index);
                    wound_nodes.push_back(p_temp_node);
                    wound_node_indices.push_back(node_index);
                    break;
                }
            }
        }
        // Add the wound element to the mesh
        unsigned num_elems = rMesh.GetNumElements();
        VertexElement<2,2>* p_element = new VertexElement<2,2>(num_elems, wound_nodes);
        rMesh.AddElement(p_element);
        rMesh.ReMesh();
    }

    /*
     * This is a helper method to generate cells with a cell cycle model.
     */ 

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double equilibriumVolume, double quiescentVolumeFraction, 
            double stemCellDivisionProbability, double stemCellMinimumDivisionAge)
    {

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_stem_cell_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            BernoulliTrialWithContactInhibitionCellCycleModel* p_model = new BernoulliTrialWithContactInhibitionCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetEquilibriumVolume(equilibriumVolume);
            p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            p_model->SetStemCellDivisionProbability(stemCellDivisionProbability);
            p_model->SetStemCellMinimumDivisionAge(stemCellMinimumDivisionAge);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_cell_type);
            double ave_stem_cell_cycle_duration = p_model->GetAverageStemCellCycleTime();
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * ave_stem_cell_cycle_duration;
            p_cell->SetBirthTime(birth_time);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", sqrt(3.0)/2.0);
            // p_cell->GetCellData()->SetItem("target area", 1.0);

            rCells.push_back(p_cell);
        }
    }

public:

    void TestWoundMechanisms()
    {
        // Simulation and sweep parameters
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_param_vals"));
        unsigned num_param_vals = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-num_param_vals").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_prolif_sims"));
        unsigned num_prolif_sims = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-num_prolif_sims").c_str());

        // Set min/max parameter values for each method of healing
        // Boundary adhesion
        double min_cell_boundary_adhesion = 30.0;
        double max_cell_boundary_adhesion = 60.0;
        
        // Contact inhibition levels
        double min_div_prob = 0.025;
        double max_div_prob = 0.05;

        // Wound cell deformation
        double min_wound_cell_membrane_energy = 3.0;
        double max_wound_cell_membrane_energy = 5.0;

        // Normal force
        double min_normal_force_strength = 6.0;
        double max_normal_force_strength = 8.0;

        // Wound centre force
        double min_wound_centre_force = 5.0;
        double max_wound_centre_force = 6.0;

        // Loop over parameter values
        for(unsigned sim_index=0; sim_index < num_param_vals; sim_index++)
        {    
            std::cout << " Run number " << sim_index << "... \n" << std::flush;   
            // Reseed the random number generator
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            /* 
            * == VM ==
            * 
            * Simulation internal void using the
            * Cell Vertex model.
            * Differential adhesion
            */
            // void TestWoundDifferentialAdhesion()
            // {
            //     std::string output_directory =  M_HEAD_FOLDER + "/Pre-void";
            //     // Load steady state
            //     OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            //     VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

            //     SimulationTime::Instance()->Destroy();
            //     SimulationTime::Instance()->SetStartTime(0.0);

            //     // Remove the forces and boundaries - redefined here
            //     p_simulator->RemoveAllForces();

            //     std::cout << "Differential adhesion \n" << std::flush; 
            //     p_gen->Reseed(sim_index);
            //     double cell_boundary_adhesion = min_cell_boundary_adhesion + (max_cell_boundary_adhesion - min_cell_boundary_adhesion) * double(sim_index) / double(num_param_vals);
            //     std::stringstream paramAsString;
            //     paramAsString << cell_boundary_adhesion;
            //     output_directory =  M_HEAD_FOLDER + "/DifferentialAdhesion/CellBoundaryAdhesion_" + paramAsString.str();

            //     /* 
            //     * == Post-void == 
            //     */

            //     // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
            //     MAKE_PTR(NagaiHondaForce<2>, p_force);
            //     p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
            //     p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
            //     p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
            //     p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(cell_boundary_adhesion);
            //     p_simulator->AddForce(p_force);

            //     // Track the area of the void
            //     MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
            //     voidarea_modifier->SetOutputDirectory(output_directory);
            //     p_simulator->AddSimulationModifier(voidarea_modifier);

            //     // Create simulation from cell population
            //     p_simulator->SetDt(M_DT_TIME);
            //     p_simulator->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            //     p_simulator->SetEndTime(M_END_TIME);
            //     p_simulator->SetOutputDirectory(output_directory);
                
            //     // Run simulation
            //     try
            //     {
            //         p_simulator->Solve();
            //     }
            //     catch (Exception& e)
            //     {
            //         // If it throws then we report the error message and go to the next simulation
            //         WARNING("Simulation didnt run" << paramAsString.str() << ".");
            //         WARNING(e.GetMessage());            
            //     }
            //     // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
            //     SimulationTime::Instance()->Destroy();
            //     SimulationTime::Instance()->SetStartTime(0.0);
            //     // Tidy up
            //     delete p_simulator;
            // }

            // // Cell proliferation
            // // void xTestWoundProliferation()
            // {
            //     for (unsigned index2 = 0; index2<num_prolif_sims; index2++)
            //     {
            //         std::string output_directory =  M_HEAD_FOLDER + "/Pre-void";
            //         // Load steady state
            //         OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            //         VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

            //         SimulationTime::Instance()->Destroy();
            //         SimulationTime::Instance()->SetStartTime(0.0);

            //         // Remove the forces and boundaries - redefined here
            //         p_simulator->RemoveAllForces();

            //         std::cout << "Proliferation \n" << std::flush; 
            //         p_gen->Reseed(index2);
            //         double div_prob = min_div_prob + (max_div_prob - min_div_prob) * double(sim_index) / double(num_param_vals);
            //         std::stringstream paramAsString;
            //         paramAsString << div_prob;
            //         std::stringstream runAsString;
            //         runAsString << index2;
            //         output_directory =  M_HEAD_FOLDER + "/Proliferation/Division_probability_" + paramAsString.str() + "/Run_" + runAsString.str();

            //         /* 
            //         * == Post-void == 
            //         */
            //        boost::shared_ptr<AbstractCellProperty> p_stem_cell_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

            //         // Loop over cells to set up the leading edge and proliferative hub. 
            //         for (typename VertexMesh<2,2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            //             elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            //             ++elem_iter)
            //         {
            //             unsigned elem_index = elem_iter->GetIndex();
            //             // Get cell associated with this element
            //             CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_index);

            //             // For each node we create a cell with our cell-cycle model. 
            //             BernoulliTrialWithContactInhibitionCellCycleModel* p_cc_model = static_cast<BernoulliTrialWithContactInhibitionCellCycleModel*>(p_cell->GetCellCycleModel());
            //             p_cc_model->SetDimension(2);
            //             p_cc_model->SetEquilibriumVolume(sqrt(3.0)/2.0);
            //             p_cc_model->SetQuiescentVolumeFraction(0.8);
            //             p_cc_model->SetStemCellDivisionProbability(div_prob);
            //             p_cc_model->SetStemCellMinimumDivisionAge(0.0);

            //             p_cell->SetCellProliferativeType(p_stem_cell_type);
            //             p_cell->SetBirthTime(0.0);

            //             // Set Target Area so dont need to use a growth model in vertex simulations
            //             p_cell->GetCellData()->SetItem("target area", sqrt(3.0)/2.0);
            //         }

            //         // Track the area of the void
            //         MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
            //         voidarea_modifier->SetOutputDirectory(output_directory);
            //         p_simulator->AddSimulationModifier(voidarea_modifier);

            //         // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
            //         MAKE_PTR(NagaiHondaForce<2>, p_force);
            //         p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
            //         p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
            //         p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
            //         p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
            //         p_simulator->AddForce(p_force);

            //         // Create simulation from cell population
            //         p_simulator->SetDt(M_DT_TIME);
            //         p_simulator->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            //         p_simulator->SetEndTime(M_END_TIME);
            //         p_simulator->SetOutputDirectory(output_directory);
                    
            //         // Run simulation
            //         try
            //         {
            //             p_simulator->Solve();
            //         }
            //         catch (Exception& e)
            //         {
            //             // If it throws then we report the error message and go to the next simulation
            //             WARNING("Simulation didnt run" << paramAsString.str() << ".");
            //             WARNING(e.GetMessage());            
            //         }
            //         // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
            //         SimulationTime::Instance()->Destroy();
            //         SimulationTime::Instance()->SetStartTime(0.0);
            //         // Tidy up
            //         delete p_simulator;
            //     }
            // }

            // // Wound cell
            // // void xTestWoundCell()
            // {
            //     std::string output_directory =  M_HEAD_FOLDER + "/Pre-voidCell";
            //     // Load steady state
            //     OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            //     VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

            //     SimulationTime::Instance()->Destroy();
            //     SimulationTime::Instance()->SetStartTime(0.0);

            //     // Remove the forces and boundaries - redefined here
            //     p_simulator->RemoveAllForces();

            //     std::cout << "Wound area goes to 0 \n" << std::flush; 
            //     p_gen->Reseed(sim_index);
            //     double wound_cell_membrane_energy = min_wound_cell_membrane_energy + (max_wound_cell_membrane_energy - min_wound_cell_membrane_energy) * double(sim_index)/double(num_param_vals);
            //     std::stringstream paramAsString;
            //     paramAsString << wound_cell_membrane_energy;
            //     output_directory =  M_HEAD_FOLDER + "/WoundAreaZero/WoundCellMembraneEnergy_" + paramAsString.str();

            //     /* 
            //     * == Post-void == 
            //     */

            //     // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
            //     MAKE_PTR(NagaiHondaMutationCellForce<2>, p_force);
            //     p_force->SetNagaiHondaCellDeformationEnergyParameter(50.0);
            //     p_force->SetNagaiHondaWoundDeformationEnergyParameter(0.0);
            //     p_force->SetNagaiHondaCellMembraneSurfaceEnergyParameter(1.0);
            //     p_force->SetNagaiHondaWoundMembraneSurfaceEnergyParameter(wound_cell_membrane_energy);
            //     p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
            //     p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
            //     p_force->SetNagaiHondaCellWoundAdhesionEnergyParameter(0.5);
            //     p_simulator->AddForce(p_force);

            //     // Track the area of the void
            //     MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
            //     voidarea_modifier->SetOutputDirectory(output_directory);
            //     p_simulator->AddSimulationModifier(voidarea_modifier);

            //     // Create simulation from cell population
            //     p_simulator->SetDt(M_DT_TIME);
            //     p_simulator->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            //     p_simulator->SetEndTime(M_END_TIME);
            //     p_simulator->SetOutputDirectory(output_directory);
                
            //     // Run simulation
            //     try
            //     {
            //         p_simulator->Solve();
            //     }
            //     catch (Exception& e)
            //     {
            //         // If it throws then we report the error message and go to the next simulation
            //         WARNING("Simulation didnt run" << paramAsString.str() << ".");
            //         WARNING(e.GetMessage());            
            //     }

            //     // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
            //     SimulationTime::Instance()->Destroy();
            //     SimulationTime::Instance()->SetStartTime(0.0);
            //     // Tidy up
            //     delete p_simulator;
            // }

            // Normal force
            // void xTestWoundNormalForce()
            {
                std::string output_directory =  M_HEAD_FOLDER + "/Pre-void/Circle";
                // Load steady state
                OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
                VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Remove the forces and boundaries - redefined here
                p_simulator->RemoveAllForces();

                std::cout << "Normal force \n" << std::flush; 
                p_gen->Reseed(sim_index);
                double normal_force_strength = min_normal_force_strength + (max_normal_force_strength - min_normal_force_strength) * double(sim_index) / double(num_param_vals); 
                std::stringstream paramAsString;
                paramAsString << normal_force_strength;
                output_directory =  M_HEAD_FOLDER + "Circle/NormalForce/NormalForceStrength_" + paramAsString.str();

                /* 
                * == Post-void == 
                */

                // Create Forces and pass to simulation NOTE : these are not the default ones and chosen for stability
                MAKE_PTR(NagaiHondaForce<2>, p_force);
                p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
                p_simulator->AddForce(p_force);

                MAKE_PTR(BoundaryNormalForce<2>, p_bound_force);
                p_bound_force->SetForceStrength(normal_force_strength);
                p_simulator->AddForce(p_bound_force);

                // Track the area of the void
                MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
                voidarea_modifier->SetOutputDirectory(output_directory);
                p_simulator->AddSimulationModifier(voidarea_modifier);

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
                    WARNING("Simulation didnt run" << paramAsString.str() << ".");
                    WARNING(e.GetMessage());            
                }

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
                // Tidy up
                delete p_simulator;
            }

            // // Wound centre force
            // // void xTestWoundCentreForce()
            // {
            //     std::string output_directory =  M_HEAD_FOLDER + "/Pre-void";
            //     // Load steady state
            //     OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            //     VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

            //     SimulationTime::Instance()->Destroy();
            //     SimulationTime::Instance()->SetStartTime(0.0);

            //     // Remove the forces and boundaries - redefined here
            //     p_simulator->RemoveAllForces();

            //     std::cout << "Centre force \n" << std::flush; 
            //     p_gen->Reseed(sim_index);
            //     double wound_centre_force = min_wound_centre_force + (max_wound_centre_force - min_wound_centre_force) * double(sim_index) / double(num_param_vals);
            //     std::stringstream paramAsString;
            //     paramAsString << wound_centre_force;
            //     output_directory =  M_HEAD_FOLDER + "/WoundCentreForce/WoundCentreForceStrength_" + paramAsString.str();

            //     /* 
            //     * == Post-void == 
            //     */

            //     // Create Forces and pass to simulation NOTE : these are not the default ones and chosen for stability
            //     MAKE_PTR(NagaiHondaForce<2>, p_force);
            //     p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
            //     p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
            //     p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
            //     p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
            //     p_simulator->AddForce(p_force);

            //     MAKE_PTR(WoundCentreForce<2>, p_centre_force);
            //     p_centre_force->SetForceStrength(wound_centre_force);
            //     p_simulator->AddForce(p_centre_force);

            //     // Track the area of the void
            //     MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
            //     voidarea_modifier->SetOutputDirectory(output_directory);
            //     p_simulator->AddSimulationModifier(voidarea_modifier);

            //     // Create simulation from cell population
            //     p_simulator->SetDt(M_DT_TIME);
            //     p_simulator->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            //     p_simulator->SetEndTime(M_END_TIME);
            //     p_simulator->SetOutputDirectory(output_directory);
                
            //     // // Track the area of the wound
            //     // MAKE_PTR(WoundAreaModifier<2>, wound_area_modifier);
            //     // wound_area_modifier->SetOutputDirectory(output_directory);
            //     // simulator.AddSimulationModifier(wound_area_modifier);
                
            //     // // Smooth out edges to get nice box domain
            //     // SmoothVertexMeshEdges(cell_population);
                
            //     // Run simulation
            //     try
            //     {
            //         p_simulator->Solve();
            //     }
            //     catch (Exception& e)
            //     {
            //         // If it throws then we report the error message and go to the next simulation
            //         WARNING("Simulation didnt run" << paramAsString.str() << ".");
            //         WARNING(e.GetMessage());            
            //     }

            //     // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
            //     SimulationTime::Instance()->Destroy();
            //     SimulationTime::Instance()->SetStartTime(0.0);

            //     // Tidy up
            //     delete p_simulator;

            //     // /*
            //     //  * == Smooth void == 
            //     //  */
            //     // {
            //     //     // Load steady state
            //     //     OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            //     //     VertexBasedCellPopulation<2>* p_cell_population_1 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()));

            //     //     std::string output_directory_1 =  M_HEAD_FOLDER + "/Smooth";

            //     //     // Now remove cells in a given region using a helper method
            //     //     // CreateHoleInCellPopulation(*p_cell_population_1);
            //     //     // SmoothVertexMeshEdges(*p_cell_population_1);

            //     //     // Track the area of the void
            //     //     MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_1);
            //     //     voidarea_modifier_1->SetOutputDirectory(output_directory_1);
            //     //     p_simulator_1->AddSimulationModifier(voidarea_modifier_1);

            //     //     // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            //     //     p_simulator_1->SetDt(M_DT_TIME);
            //     //     p_simulator_1->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            //     //     p_simulator_1->SetEndTime(M_END_TIME);
            //     //     p_simulator_1->SetOutputDirectory(output_directory_1);
            //     //     p_simulator_1->Solve();

            //     //     // Tidy up
            //     //     delete p_simulator_1;
            //     // }

            //     /*
            //     * == Jagged void ==
            //     */
            //     // Load steady state

            //     // {
            //     //     OffLatticeSimulation<2>* p_simulator_2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            //     //     VertexBasedCellPopulation<2>* p_cell_population_2 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator_2->rGetCellPopulation()));
                    
            //     //     std::string output_directory_2 =  M_HEAD_FOLDER + "/Vertex/Jagged";

            //     //     // Now remove cells in a given region using a helper method
            //     //     CreateHoleInCellPopulation(*p_cell_population_2);

            //     //     // Track the area of the void
            //     //     MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_2);
            //     //     voidarea_modifier_2->SetOutputDirectory(output_directory_2);
            //     //     p_simulator_2->AddSimulationModifier(voidarea_modifier_2);

            //     //     // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            //     //     p_simulator_2->SetDt(M_DT_TIME);
            //     //     p_simulator_2->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            //     //     p_simulator_2->SetEndTime(M_END_TIME);
            //     //     p_simulator_2->SetOutputDirectory(output_directory_2);
            //     //     p_simulator_2->Solve();

            //     //     // Tidy up
            //     //     delete p_simulator_2;
            //     // }

            //     /*
            //     * == Curved edges ==
            //     */
            //     // {
            //     //     OffLatticeSimulation<2>* p_simulator_2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            //     //     VertexBasedCellPopulation<2>* p_cell_population_2 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator_2->rGetCellPopulation()));
                    
            //     //     std::string output_directory_2 =  M_HEAD_FOLDER + "/Vertex/Curved";

            //     //     // Now remove cells in a given region using a helper method
            //     //     CreateHoleInCellPopulation(*p_cell_population_2);

            //     //     // Refine the edges on boundary to get smooth edges
            //     //     MAKE_PTR(VertexBoundaryRefinementModifier<2>, refinement_modifier);
            //     //     p_simulator_2->AddSimulationModifier(refinement_modifier);

            //     //     // Track the area of the void
            //     //     MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_2);
            //     //     voidarea_modifier_2->SetOutputDirectory(output_directory_2);
            //     //     p_simulator_2->AddSimulationModifier(voidarea_modifier_2);

            //     //     // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            //     //     p_simulator_2->SetDt(M_DT_TIME);
            //     //     p_simulator_2->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            //     //     p_simulator_2->SetEndTime(M_END_TIME);
            //     //     p_simulator_2->SetOutputDirectory(output_directory_2);
            //     //     p_simulator_2->Solve();

            //     //     // Tidy up
            //     //     delete p_simulator_2;
            //     // }

            // }
        }
    }
};

#endif /* TESTINTERNALVOID_HPP_ */
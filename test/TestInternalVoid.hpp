#ifndef TESTINTERNALVOID_HPP_
#define TESTINTERNALVOID_HPP_


#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 
#include "CheckpointArchiveTypes.hpp"

#include "SmartPointers.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "PeriodicNodesOnlyMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "MutableMesh.hpp"
#include "Toroidal2dMesh.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "NoCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellVolumesWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "BoundaryNodeWriter.hpp"

#include "OffLatticeSimulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"

#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "GhostNodeRemovalModifier.hpp"
#include "VoidAreaModifier.hpp"
#include "VertexBoundaryRefinementModifier.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp" 
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 *  Code based on https://github.com/jmosborne/TissueBoundaries
 */

static const double M_END_STEADY_STATE = 1.0;
static const double M_END_TIME = 0.5;
static const double M_DT_TIME = 0.001;
static const double M_SAMPLE_TIME = 0.01/M_DT_TIME;

// Both Width and Length must be EVEN numbers here
// static const double M_DOMAIN_WIDTH = 20;
// static const double M_DOMAIN_LENGTH = 20;
// static const double M_DOMAIN_SCALING = 0.8;
// static const double M_PERIODIC_WIDTH = 16;//M_DOMAIN_WIDTH*M_DOMAIN_SCALING;
// static const double M_PERIODIC_HEIGHT = 14;//M_DOMAIN_LENGTH*0.5*sqrt(3)*M_DOMAIN_SCALING;

// static const double M_HOLEWIDTH = 2.0;
// static const double M_HOLE_X_MIN = 6.0;//M_DOMAIN_WIDTH/2 - 3;
// static const double M_HOLE_X_MAX = 10.0;//M_DOMAIN_WIDTH/2 + 3;
// static const double M_HOLE_Y_MIN = 5.0;
// static const double M_HOLE_Y_MAX = 9.0;

// Both Width and Length must be EVEN numbers here
static const double M_DOMAIN_WIDTH = 14;
static const double M_DOMAIN_LENGTH = 16;
static const double M_DOMAIN_SCALING = 10.0/12.0;
static const double M_DOMAIN_SCALING_X = 12.0/14.0;
static const double M_DOMAIN_SCALING_Y = 13.0*sqrt(3.0/4.0)/(16.0*sqrt(3.0/4.0));
static const double M_PERIODIC_WIDTH = 12.0;//M_DOMAIN_WIDTH*M_DOMAIN_SCALING;
static const double M_PERIODIC_HEIGHT = 13.0*sqrt(3.0/4.0);//M_DOMAIN_LENGTH*0.5*sqrt(3)*M_DOMAIN_SCALING;

static const double M_HOLEWIDTH = 2.3;
static const double M_HOLE_X_MIN = 4.2;
static const double M_HOLE_X_MAX = 7.8;
static const double M_HOLE_Y_MIN = 4.0;
static const double M_HOLE_Y_MAX = 8.0;

// To remove 11 cells in a symmetric capsule: holewidth = 1.9; x_min = 4.2; x_max = 7.5; y_min = 4.0; y_max = 7.5

static const std::string M_HEAD_FOLDER = "InternalVoid/Capsule";


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
        double xc = M_PERIODIC_WIDTH*0.5;
        double yc;
        // Slightly differen y values in cell-centre vs vertex dynamics because of node definition
        if (bool(dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation)) || bool(dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                yc = M_PERIODIC_HEIGHT*0.5;
            }
            else if (bool(dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                // -0.25 for vertex model because shifting nodes rather than cell centres.
                yc = M_PERIODIC_HEIGHT*0.5 - 0.25;
            }
        
        if (bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation)))
        {
            std::set<unsigned> location_indices;
            std::set<unsigned> ghost_node_indices;

            for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
                cell_iter != rCellPopulation.rGetCells().end();)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double x = centre_of_cell[0];
                double y = centre_of_cell[1];
                unsigned location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
 
                // For a capsule shaped wound
                if ((fabs(y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
                // For an approximately circular wound
                // if (sqrt(pow(x-xc,2) + pow(y-yc,2)) < M_HOLEWIDTH)
                {   
                    // Delete cell and store it as a ghost node
                    rCellPopulation.RemoveCellUsingLocationIndex(location_index, (*cell_iter));
                    // Update vector of cells
                    cell_iter = rCellPopulation.rGetCells().erase(cell_iter);
 
                    // Change to chost node            
                    ghost_node_indices.insert(location_index);
                }
                else
                {
                    ++cell_iter;
                }
            }
            dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->SetGhostNodes(ghost_node_indices);
            //dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->RemoveDeadCells();
            //dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->Update();
        }
        else
        {
            for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double x = centre_of_cell[0];
                double y = centre_of_cell[1];

                // For a capsule shaped wound
                if ((fabs(y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
                // For an approximately circular wound
                // if (sqrt(pow(x-xc,2) + pow(y-yc,2)) < M_HOLEWIDTH)
                {   
                    cell_iter->Kill();
                }
            }
            
            /* Need to remove the deleted cells and call update note this is usually
            * performed in the Solve() method of the simulation class.
            */
            if (bool(dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<NodeBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<NodeBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
            else if (bool(dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<MeshBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<MeshBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
            else if (bool(dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
        }
    }

public:

    /* 
     * == OS ==
     *
     * Simulate an internal void using the
     * Overlapping Spheres model.
     *
     * Default Cut-off = 1.5
     */
    void xTestNodeBasedDefaultCutoffInternalVoid()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node/DefaultCutOff/Pre-void";
        /* 
         * == Pre-void == 
         */
         // Create simple mesh
        HoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        // p_generating_mesh->Scale(M_DOMAIN_SCALING, M_DOMAIN_SCALING);

        double cut_off_length = 1.5; //this is the default

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,2> periodic_width = zero_vector<double>(2);
        // periodic_width[0] = M_PERIODIC_WIDTH;
        // periodic_width[1] = M_PERIODIC_HEIGHT;
        periodic_width[0] = M_DOMAIN_WIDTH;
        periodic_width[1] = M_DOMAIN_LENGTH;
        PeriodicNodesOnlyMesh<2>* p_mesh = new PeriodicNodesOnlyMesh<2>(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);
        // NodesOnlyMesh<2>* p_mesh;
        // p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);
        // Convert this to a Cylindrical2dNodesOnlyMesh
        // Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(M_PERIODIC_WIDTH);
        // p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,2.0); // So factor of 16

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        // CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        // simulator.SetOutputDivisionLocations(true);
        // simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        // MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        // simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);

        // Track the area of the void
        // MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        // voidarea_modifier->SetOutputDirectory(output_directory);
        // // voidarea_modifier->SetCutoff(cut_off_length);
        // simulator.AddSimulationModifier(voidarea_modifier);

        // Run simulation
        simulator.Solve();

        // // Save simulation in steady state
		// CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // // Clear memory
        // delete p_mesh;

        // /*
        //  * == Void default cut-off == 
        //  */
        // {
        //     // Load steady state
        //     OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
        //     NodeBasedCellPopulation<2>* p_cell_population_1 = static_cast<NodeBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()));

        //     std::string output_directory_1 =  M_HEAD_FOLDER + "/Node/DefaultCutOff/Post-Void";
        //     SimulationTime::Instance()->Destroy();
        //     SimulationTime::Instance()->SetStartTime(0.0);

        //     // Now remove cells in a given region using a helper method
        //     CreateHoleInCellPopulation(*p_cell_population_1);
       
        //     // Track the area of the void
        //     MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier1);
        //     voidarea_modifier1->SetOutputDirectory(output_directory_1);
        //     // voidarea_modifier1->SetCutoff(cut_off_length);
        //     p_simulator_1->AddSimulationModifier(voidarea_modifier1);
        
        //     // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        //     p_simulator_1->SetDt(M_DT_TIME);
        //     p_simulator_1->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        //     p_simulator_1->SetEndTime(M_END_TIME);
        //     p_simulator_1->SetOutputDirectory(output_directory_1);
        //     p_simulator_1->Solve();

        //     // Tidy up
        //     delete p_simulator_1;
        // }

    }

    /*
     * == Larger cut-off ==
     * Cut-off = 2.0
     */
    void xTestNodeBasedLargeCutoffInternalVoid()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node/LargeCutoff/Pre-void";
        /* 
         * == Pre-void == 
         */
         // Create simple mesh
        HoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        p_generating_mesh->Scale(M_DOMAIN_SCALING, M_DOMAIN_SCALING);

        double cut_off_length = 2.0; 

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,2> periodic_width = zero_vector<double>(2);
        periodic_width[0] = M_PERIODIC_WIDTH;
        periodic_width[1] = M_PERIODIC_HEIGHT;
        PeriodicNodesOnlyMesh<2>* p_mesh = new PeriodicNodesOnlyMesh<2>(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        // voidarea_modifier->SetCutoff(cut_off_length);
        simulator.AddSimulationModifier(voidarea_modifier);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Clear memory
        delete p_mesh;

        /*
         * == Void default cut-off == 
         */
        {
            // Load steady state
            OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            NodeBasedCellPopulation<2>* p_cell_population_1 = static_cast<NodeBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()));

            std::string output_directory_1 =  M_HEAD_FOLDER + "/Node/LargeCutoff/Post-void";
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Now remove cells in a given region using a helper method
            CreateHoleInCellPopulation(*p_cell_population_1);
       
            // Track the area of the void
            MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier1);
            voidarea_modifier1->SetOutputDirectory(output_directory_1);
            // voidarea_modifier1->SetCutoff(cut_off_length);
            p_simulator_1->AddSimulationModifier(voidarea_modifier1);
        
            // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            p_simulator_1->SetDt(M_DT_TIME);
            p_simulator_1->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            p_simulator_1->SetEndTime(M_END_TIME);
            p_simulator_1->SetOutputDirectory(output_directory_1);
            p_simulator_1->Solve();

            // Tidy up
            delete p_simulator_1;
        }

    }

    /*
     * == Small cut-off ==
     * Cut-off = 1.0
     */
    void xTestNodeBasedSmallCutoffInternalVoid()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node/SmallCutoff/Pre-void";
        /* 
         * == Pre-void == 
         */
         // Create simple mesh
        HoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        p_generating_mesh->Scale(M_DOMAIN_SCALING, M_DOMAIN_SCALING);

        double cut_off_length = 1.0; 

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,2> periodic_width = zero_vector<double>(2);
        periodic_width[0] = M_PERIODIC_WIDTH;
        periodic_width[1] = M_PERIODIC_HEIGHT;
        PeriodicNodesOnlyMesh<2>* p_mesh = new PeriodicNodesOnlyMesh<2>(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        // voidarea_modifier->SetCutoff(cut_off_length);
        simulator.AddSimulationModifier(voidarea_modifier);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Clear memory
        delete p_mesh;

        /*
         * == Void default cut-off == 
         */
        {
            // Load steady state
            OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            NodeBasedCellPopulation<2>* p_cell_population_1 = static_cast<NodeBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()));

            std::string output_directory_1 =  M_HEAD_FOLDER + "/Node/SmallCutoff/Post-void";
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Now remove cells in a given region using a helper method
            CreateHoleInCellPopulation(*p_cell_population_1);
       
            // Track the area of the void
            MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier1);
            voidarea_modifier1->SetOutputDirectory(output_directory_1);
            // voidarea_modifier1->SetCutoff(cut_off_length);
            p_simulator_1->AddSimulationModifier(voidarea_modifier1);
        
            // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            p_simulator_1->SetDt(M_DT_TIME);
            p_simulator_1->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            p_simulator_1->SetEndTime(M_END_TIME);
            p_simulator_1->SetOutputDirectory(output_directory_1);
            p_simulator_1->Solve();

            // Tidy up
            delete p_simulator_1;
        }

    }

    /* 
     * == VT ==
     * 
     * Simulate internal voide using the
     * Voronoi Tesselation model.
     */

    /*
     * == No ghosts == 
     */
    void TestMeshBasedNoGhostsInternalVoid()
    {
        std::string output_directory =  M_HEAD_FOLDER + "/Mesh/NoGhosts/Pre-Void";

        /* // Create mesh
        ToroidalHoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, M_DOMAIN_SCALING, M_DOMAIN_SCALING);
        Toroidal2dMesh* p_mesh = generator.GetToroidalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create tissue
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Output Voroni for visualisation
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.SetWriteVtkAsPoints(true);

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking Modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        simulator.AddSimulationModifier(voidarea_modifier);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator); */

        /*
         * == No ghosts Infinite VT == 
         */
        // {
        //     // Load steady state
        //     OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
        //     MeshBasedCellPopulation<2>* p_cell_population_1 = static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()));
        //     // Remove the forces and boundaries - redefined here
        //     p_simulator_1->RemoveAllForces();

        //     std::string output_directory_1 =  M_HEAD_FOLDER + "/Mesh/NoGhosts/InfiniteVT";
        //     SimulationTime::Instance()->Destroy();
        //     SimulationTime::Instance()->SetStartTime(0.0);

        //     // Now remove cells in a given region using a helper method
        //     CreateHoleInCellPopulation(*p_cell_population_1);

        //     // Create a force law and pass it to the simulation
        //     MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        //     p_force->SetMeinekeSpringStiffness(50.0);
        //     p_force->SetCutOffLength(DBL_MAX);
        //     p_simulator_1->AddForce(p_force);

        //     // Reset end time for simulation and run for a further duration
        //     p_simulator_1->SetOutputDirectory(output_directory_1);
        //     p_simulator_1->SetEndTime(M_END_TIME);
        //     p_simulator_1->Solve();

        //     // Tidy up
        //     delete p_simulator_1;
        // }

        /*
         * == No ghosts Finite VT == 
         */
        {
            // Load steady state
            OffLatticeSimulation<2>* p_simulator_2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            MeshBasedCellPopulation<2>* p_cell_population_2 = static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator_2->rGetCellPopulation()));

            std::string output_directory_2 =  M_HEAD_FOLDER + "/Mesh/NoGhosts/FiniteVT/MidSize/Slow";
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Now remove cells in a given region using a helper method
            CreateHoleInCellPopulation(*p_cell_population_2);
            
            // Track the area of the void
            MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_2);
            voidarea_modifier_2->SetOutputDirectory(output_directory_2);
            p_simulator_2->AddSimulationModifier(voidarea_modifier_2);

            // Remove the forces - redefined here
            p_simulator_2->RemoveAllForces();
            // Create a force law and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force2);
            p_linear_force2->SetMeinekeSpringStiffness(25.0);
            p_linear_force2->SetCutOffLength(1.5);
            p_simulator_2->AddForce(p_linear_force2);
            
            // Bound the VT
            p_cell_population_2->SetBoundVoronoiTessellation(true);

            // Reset end time for simulation and run for a further duration
            p_simulator_2->SetOutputDirectory(output_directory_2);
            p_simulator_2->SetEndTime(M_END_TIME);
            p_simulator_2->Solve();

            // Tidy up
            delete p_simulator_2;
        }
    }

    /*
     * == Ghosts ==
     */
    void xTestMeshBasedGhostsInternalVoid()
    {
        std::string output_directory =  M_HEAD_FOLDER + "/Mesh/Ghosts/Pre-Void";

        // Create mesh
        ToroidalHoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, M_DOMAIN_SCALING, M_DOMAIN_SCALING);
        Toroidal2dMesh* p_mesh = generator.GetToroidalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create tissue
        double ghost_cell_spring_stiffness = 15.0; // Default
        double ghost_ghost_spring_stiffness = 1.0;
        double ghost_spring_rest_length = 1.0; // Default
         
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, std::vector<unsigned>(), false, ghost_cell_spring_stiffness, ghost_ghost_spring_stiffness, ghost_spring_rest_length);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Output Voroni for visualisation
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.SetWriteVtkAsPoints(true);

        cell_population.AddPopulationWriter<BoundaryNodeWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking Modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier_1);
        simulator.AddSimulationModifier(p_modifier_1);

        // Add volume tracking Modifier
        MAKE_PTR(GhostNodeRemovalModifier<2>, p_modifier_2);
        simulator.AddSimulationModifier(p_modifier_2);        

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(DBL_MAX);
        simulator.AddForce(p_linear_force);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        simulator.AddSimulationModifier(voidarea_modifier);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);


         /*
          * == Remove Ghosts ==
          *
          */ 
        {
            // Load steady state
            OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            MeshBasedCellPopulationWithGhostNodes<2>* p_cell_population_1 = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&(p_simulator_1->rGetCellPopulation()));

            std::string output_directory_1 =  M_HEAD_FOLDER + "/Mesh/Ghosts/Void";
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Now remove cells in a given region using a helper method
            CreateHoleInCellPopulation(*p_cell_population_1);


            // Track the area of the void
            MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_1);
            voidarea_modifier_1->SetOutputDirectory(output_directory_1);
            p_simulator_1->AddSimulationModifier(voidarea_modifier_1);

            // Reset timestep, end time for simulation and run for a further duration
            p_simulator_1->SetOutputDirectory(output_directory_1);
            p_simulator_1->SetEndTime(M_END_TIME);
            p_simulator_1->Solve();

            // Tidy up
            delete p_simulator_1;
        }
    }


    /* 
     * == VM ==
     * 
     * Simulation internal void using the
     * Cell Vertex model.
     */
    void xTestVertexBasedInternalVoid()
    {
        std::string output_directory =  M_HEAD_FOLDER + "/Vertex/Pre-void/Adhesion1.5";

        /* 
         * == Pre-void == 
         */
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->Scale(M_DOMAIN_SCALING, M_DOMAIN_SCALING);
        p_mesh->SetHeight(M_PERIODIC_HEIGHT);
        p_mesh->SetWidth(M_PERIODIC_WIDTH);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.5);
        simulator.AddForce(p_force);

        // Add target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(0);
        p_growth_modifier->SetReferenceTargetArea(0.5*sqrt(3.0));
        simulator.AddSimulationModifier(p_growth_modifier);
        
        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        simulator.AddSimulationModifier(voidarea_modifier);
        
        // Smooth out edges to get nice box domain
        SmoothVertexMeshEdges(cell_population);
        
        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        /*
         * == Smooth void == 
         */
        {
            // Load steady state
            OffLatticeSimulation<2>* p_simulator_1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            VertexBasedCellPopulation<2>* p_cell_population_1 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()));

            std::string output_directory_1 =  M_HEAD_FOLDER + "/Vertex/Smooth";
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Now remove cells in a given region using a helper method
            CreateHoleInCellPopulation(*p_cell_population_1);
            SmoothVertexMeshEdges(*p_cell_population_1);

            // Track the area of the void
            MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_1);
            voidarea_modifier_1->SetOutputDirectory(output_directory_1);
            p_simulator_1->AddSimulationModifier(voidarea_modifier_1);

            // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            p_simulator_1->SetDt(M_DT_TIME);
            p_simulator_1->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            p_simulator_1->SetEndTime(M_END_TIME);
            p_simulator_1->SetOutputDirectory(output_directory_1);
            p_simulator_1->Solve();

            // Tidy up
            delete p_simulator_1;
        }

        /*
         * == Jagged void ==
         */
         // Load steady state

        {
            OffLatticeSimulation<2>* p_simulator_2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            VertexBasedCellPopulation<2>* p_cell_population_2 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator_2->rGetCellPopulation()));
            
            std::string output_directory_2 =  M_HEAD_FOLDER + "/Vertex/Jagged";
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Now remove cells in a given region using a helper method
            CreateHoleInCellPopulation(*p_cell_population_2);

            // Track the area of the void
            MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_2);
            voidarea_modifier_2->SetOutputDirectory(output_directory_2);
            p_simulator_2->AddSimulationModifier(voidarea_modifier_2);

            // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            p_simulator_2->SetDt(M_DT_TIME);
            p_simulator_2->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            p_simulator_2->SetEndTime(M_END_TIME);
            p_simulator_2->SetOutputDirectory(output_directory_2);
            p_simulator_2->Solve();

            // Tidy up
            delete p_simulator_2;
        }

        /*
         * == Curved edges ==
         */
        {
            OffLatticeSimulation<2>* p_simulator_2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,M_END_STEADY_STATE);
            VertexBasedCellPopulation<2>* p_cell_population_2 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator_2->rGetCellPopulation()));
            
            std::string output_directory_2 =  M_HEAD_FOLDER + "/Vertex/Curved";
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Now remove cells in a given region using a helper method
            CreateHoleInCellPopulation(*p_cell_population_2);

            // Refine the edges on boundary to get smooth edges
            MAKE_PTR(VertexBoundaryRefinementModifier<2>, refinement_modifier);
            p_simulator_2->AddSimulationModifier(refinement_modifier);

            // Track the area of the void
            MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier_2);
            voidarea_modifier_2->SetOutputDirectory(output_directory_2);
            p_simulator_2->AddSimulationModifier(voidarea_modifier_2);

            // Reset timestep, sampling timestep and end time for simulation and run for a further duration
            p_simulator_2->SetDt(M_DT_TIME);
            p_simulator_2->SetSamplingTimestepMultiple(M_SAMPLE_TIME);
            p_simulator_2->SetEndTime(M_END_TIME);
            p_simulator_2->SetOutputDirectory(output_directory_2);
            p_simulator_2->Solve();

            // Tidy up
            delete p_simulator_2;
        }

    }
};

#endif /* TESTINTERNALVOID_HPP_ */
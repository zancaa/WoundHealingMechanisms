#ifndef TESTCREATEWOUND_HPP_
#define TESTCREATEWOUND_HPP_


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
static const double M_END_TIME = 5.0;
static const double M_DT_TIME = 0.001;
static const double M_SAMPLE_TIME = 100;

// Both Width and Length must be EVEN numbers here
static const double M_DOMAIN_WIDTH = 20;
static const double M_DOMAIN_LENGTH = 20;
static const double M_DOMAIN_SCALING = 0.8;
static const double M_PERIODIC_WIDTH = 20;//M_DOMAIN_WIDTH*M_DOMAIN_SCALING;
static const double M_PERIODIC_HEIGHT = 20.0*0.5*sqrt(3);//M_DOMAIN_LENGTH*0.5*sqrt(3)*M_DOMAIN_SCALING;

static const double M_HOLEWIDTH = 2.3;
static const double M_HOLE_X_MIN = 7.0;//M_DOMAIN_WIDTH/2 - 3;
static const double M_HOLE_X_MAX = 13.5;//M_DOMAIN_WIDTH/2 + 3;
static const double M_HOLE_Y_MIN = 7.0*0.5*sqrt(3.0);
static const double M_HOLE_Y_MAX = 13.0*0.5*sqrt(3.0);

static const std::string M_HEAD_FOLDER = "WoundMechanisms";


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
        double yc = M_PERIODIC_HEIGHT*0.5;
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            // Get the coordinates of this cell centre
            c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell[0];
            double y = centre_of_cell[1];

            // if ((fabs((2.0/sqrt(3.0))*y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
            if (sqrt(pow(x-xc,2) + pow(y-yc,2)) < M_HOLEWIDTH)
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
        double xc = M_PERIODIC_WIDTH*0.5;
        double yc = M_PERIODIC_HEIGHT*0.5;
        for (VertexMesh<2, 2>::VertexElementIterator elem_iter = rMesh.GetElementIteratorBegin();
                elem_iter != rMesh.GetElementIteratorEnd();
                ++elem_iter)
        {
            // Get the coordinates of the centroid of the element
            unsigned elem_index = elem_iter->GetIndex();
            c_vector<double, 2> centroid_of_element = rMesh.GetCentroidOfElement(elem_index);
            double x = centroid_of_element[0];
            double y = centroid_of_element[1];

            // if ((fabs((2.0/sqrt(3.0))*y-x)<M_HOLEWIDTH) && (x>M_HOLE_X_MIN) && (x<M_HOLE_X_MAX) && (y>M_HOLE_Y_MIN) && (y<M_HOLE_Y_MAX))
            if (sqrt(pow(x-xc,2) + pow(y-yc,2)) < M_HOLEWIDTH)
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
                p_temp_node->SetAsBoundaryNode(false);
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
                    p_temp_node->SetAsBoundaryNode(false);
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
        boost::shared_ptr<AbstractCellProperty> p_differentiated_cell_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            BernoulliTrialWithContactInhibitionCellCycleModel* p_model = new BernoulliTrialWithContactInhibitionCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetEquilibriumVolume(equilibriumVolume);
            p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            p_model->SetStemCellDivisionProbability(stemCellDivisionProbability);
            p_model->SetStemCellMinimumDivisionAge(stemCellMinimumDivisionAge);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_differentiated_cell_type);
            p_cell->SetBirthTime(0.0);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", sqrt(3.0)/2.0);
            // p_cell->GetCellData()->SetItem("target area", 1.0);

            rCells.push_back(p_cell);
        }
    }

public:
    /* 
     * == VM ==
     * 
     * Simulation internal void using the
     * Vertex dynamics model.
     */

    // Wound centre force
    void TestCreateWound()
    {
        std::string output_directory =  M_HEAD_FOLDER + "/Pre-void/Circle";

        /* 
         * == Pre-void == 
         */
         // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->SetHeight(M_PERIODIC_HEIGHT);
        p_mesh->SetWidth(M_PERIODIC_WIDTH);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        CreateHoleInCellPopulation(cell_population);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputCellVelocities(true);

        // Create Forces and pass to simulation NOTE : these are not the default ones and chosen for stability
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
        simulator.AddForce(p_force);

        MAKE_PTR(WoundCentreForce<2>, p_centre_force);
        p_centre_force->SetForceStrength(10.0);
        simulator.AddForce(p_centre_force);

        // Add target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetReferenceTargetArea(sqrt(3.0)/2.0);
        p_growth_modifier->SetGrowthDuration(0);
        simulator.AddSimulationModifier(p_growth_modifier);
        
        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        simulator.AddSimulationModifier(voidarea_modifier);
        
        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);
    }

    // Wound cell
    void TestCreateWoundCell()
    {
        // Required for purse-string mechanism
        std::string output_directory =  M_HEAD_FOLDER + "/Pre-voidCell/Circle";

        /* 
         * == Pre-void == 
         */
         // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        // p_mesh->Scale(M_DOMAIN_SCALING, M_DOMAIN_SCALING);
        p_mesh->SetHeight(M_PERIODIC_HEIGHT);
        p_mesh->SetWidth(M_PERIODIC_WIDTH);

        CreateWoundMesh(*p_mesh);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(), cells,1.0,1.0,0.0,0.0); //mature_volume = 1.0

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(WoundCellMutationState, p_wound);
        unsigned num_cells = cells.size();
        cells[num_cells-1]->SetMutationState(p_wound);
        cells[num_cells-1]->GetCellData()->SetItem("target area", 0.0);

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputCellVelocities(true);

        // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
        MAKE_PTR(NagaiHondaMutationCellForce<2>, p_force);
        p_force->SetNagaiHondaCellDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaWoundDeformationEnergyParameter(0.0);
        p_force->SetNagaiHondaCellMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaWoundMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellWoundAdhesionEnergyParameter(1.0);
        simulator.AddForce(p_force);

        // Add target area modifier
        MAKE_PTR(SimpleWoundMutantTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetReferenceTargetArea(sqrt(3.0)/2.0);
        p_growth_modifier->SetGrowthDuration(0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Track the area of the void
        MAKE_PTR(VoidAreaModifier<2>, voidarea_modifier);
        voidarea_modifier->SetOutputDirectory(output_directory);
        simulator.AddSimulationModifier(voidarea_modifier);
        
        // Run simulation
        simulator.Solve();
        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);
    }
};

#endif /* TESTCREATEWOUND_HPP_ */
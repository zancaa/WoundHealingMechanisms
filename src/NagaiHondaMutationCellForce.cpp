/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "NagaiHondaMutationCellForce.hpp"
#include "WoundCellMutationState.hpp"

template<unsigned DIM>
NagaiHondaMutationCellForce<DIM>::NagaiHondaMutationCellForce()
    : NagaiHondaForce<DIM>(),
      mNagaiHondaCellDeformationEnergyParameter(50.0),
      mNagaiHondaWoundDeformationEnergyParameter(0.0),
      mNagaiHondaCellMembraneSurfaceEnergyParameter(1.0),
      mNagaiHondaWoundMembraneSurfaceEnergyParameter(1.0),
      mNagaiHondaCellWoundAdhesionEnergyParameter(0.0)
      
{
}

template<unsigned DIM>
NagaiHondaMutationCellForce<DIM>::~NagaiHondaMutationCellForce()
{
}

template<unsigned DIM>
void NagaiHondaMutationCellForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("NagaiHondaForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of three
         * parts - a cell deformation energy, a membrane surface tension energy
         * and an adhesion energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            // Get CellPtr
            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_index);
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            c_vector<double, DIM> element_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            deformation_contribution -= 2*GetNagaiHondaDeformationEnergyParameter(p_cell)*(element_areas[elem_index] - target_areas[elem_index])*element_area_gradient;

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            double previous_edge_adhesion_parameter;
            double next_edge_adhesion_parameter;
            if (p_cell->GetMutationState()->IsType<WoundCellMutationState>())
            {
                // To avoid double counting of wound edge, no adhesion contribution
                previous_edge_adhesion_parameter = 0.0;
                next_edge_adhesion_parameter = 0.0;
            }
            else
            {
                // Compute the adhesion parameter for each of these edges
                previous_edge_adhesion_parameter = this->GetAdhesionParameter(p_previous_node, p_this_node, *p_cell_population);
                next_edge_adhesion_parameter = this->GetAdhesionParameter(p_this_node, p_next_node, *p_cell_population);
            }

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, DIM> previous_edge_gradient = -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary adhesion (note the minus sign)
            adhesion_contribution -= previous_edge_adhesion_parameter*previous_edge_gradient + next_edge_adhesion_parameter*next_edge_gradient;

            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            c_vector<double, DIM> element_perimeter_gradient;
            element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
            double cell_target_perimeter = 2*sqrt(M_PI*target_areas[elem_index]);
            membrane_surface_tension_contribution -= 2*this->GetNagaiHondaMembraneSurfaceEnergyParameter(p_cell)*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;
        }

        c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution + adhesion_contribution;
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
    }
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetNagaiHondaDeformationEnergyParameter(CellPtr pCell)
{
    if (pCell->GetMutationState()->IsType<WoundCellMutationState>())
    {
        // This cell is labelled
        return this->GetNagaiHondaWoundDeformationEnergyParameter();
    }
    else
    {
        // This cell is not labelled
        return this->GetNagaiHondaCellDeformationEnergyParameter();
    }
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter(CellPtr pCell)
{
    if (pCell->GetMutationState()->IsType<WoundCellMutationState>())
    {
        // This cell is labelled
        return this->GetNagaiHondaWoundMembraneSurfaceEnergyParameter();
    }
    else
    {
        // This cell is not labelled
        return this->GetNagaiHondaCellMembraneSurfaceEnergyParameter();
    }
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA,
                                                                      Node<DIM>* pNodeB,
                                                                      VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        unsigned element_index = *(shared_elements.begin());

        // Get cell associated with this element
        CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

        if (p_cell->GetMutationState()->IsType<WoundCellMutationState>())
        {
            // This cell is labelled
            return this->GetNagaiHondaCellWoundAdhesionEnergyParameter();
        }
        else
        {
            // This cell is not labelled
            return this->GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
        }
    }
    else
    {
        // Work out the number of labelled cells: 0,1 or 2
        unsigned num_labelled_cells = 0;
        for (std::set<unsigned>::iterator iter = shared_elements.begin();
             iter != shared_elements.end();
             ++iter)
        {
            unsigned element_index = *(iter);

            // Get cell associated with this element
            CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

            if (p_cell->GetMutationState()->IsType<WoundCellMutationState>())
            {
                num_labelled_cells++;
            }
        }

        if (num_labelled_cells == 2)
        {
            // Both cells are labelled - this should never actually happen
            return this->GetNagaiHondaCellWoundAdhesionEnergyParameter();
        }
        else if (num_labelled_cells == 1)
        {
            // One cell is labelled
            return this->GetNagaiHondaCellWoundAdhesionEnergyParameter();
        }
        else
        {
            // Neither cell is labelled
            assert(num_labelled_cells == 0);
            return this->GetNagaiHondaCellCellAdhesionEnergyParameter();
        }
    }
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetNagaiHondaCellDeformationEnergyParameter()
{
    return mNagaiHondaCellDeformationEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetNagaiHondaWoundDeformationEnergyParameter()
{
    return mNagaiHondaWoundDeformationEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetNagaiHondaCellMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaCellMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetNagaiHondaWoundMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaWoundMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaMutationCellForce<DIM>::GetNagaiHondaCellWoundAdhesionEnergyParameter()
{
    return mNagaiHondaCellWoundAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaMutationCellForce<DIM>::SetNagaiHondaCellDeformationEnergyParameter(double cellDeformationEnergyParameter)
{
    mNagaiHondaCellDeformationEnergyParameter = cellDeformationEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaMutationCellForce<DIM>::SetNagaiHondaWoundDeformationEnergyParameter(double woundDeformationEnergyParameter)
{
    mNagaiHondaWoundDeformationEnergyParameter = woundDeformationEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaMutationCellForce<DIM>::SetNagaiHondaCellMembraneSurfaceEnergyParameter(double cellMembraneSurfaceEnergyParameter)
{
    mNagaiHondaCellMembraneSurfaceEnergyParameter = cellMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaMutationCellForce<DIM>::SetNagaiHondaWoundMembraneSurfaceEnergyParameter(double woundMembraneSurfaceEnergyParameter)
{
    mNagaiHondaWoundMembraneSurfaceEnergyParameter = woundMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaMutationCellForce<DIM>::SetNagaiHondaCellWoundAdhesionEnergyParameter(double cellWoundAdhesionEnergyParameter)
{
    mNagaiHondaCellWoundAdhesionEnergyParameter = cellWoundAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaMutationCellForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variables
    *rParamsFile << "\t\t\t<NagaiHondaCellMembraneSurfaceEnergyParameter>" << mNagaiHondaCellMembraneSurfaceEnergyParameter << "</NagaiHondaCellMembraneSurfaceEnergyParameter> \n";
    *rParamsFile << "\t\t\t<NagaiHondaWoundMembraneSurfaceEnergyParameter>" << mNagaiHondaWoundMembraneSurfaceEnergyParameter << "</NagaiHondaWoundMembraneSurfaceEnergyParameter> \n";

    // Call method on direct parent class
    NagaiHondaForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class NagaiHondaMutationCellForce<1>;
template class NagaiHondaMutationCellForce<2>;
template class NagaiHondaMutationCellForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaMutationCellForce)

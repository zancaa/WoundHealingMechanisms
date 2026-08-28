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

#include "WoundAreaWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WoundCellMutationState.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
WoundAreaWriter<ELEMENT_DIM, SPACE_DIM>::WoundAreaWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("woundarea.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void WoundAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("WoundAreaWriter cannot be used with a MeshBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void WoundAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("WoundAreaWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void WoundAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("WoundAreaWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void WoundAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("WoundAreaWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void WoundAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // Helper variable that is a static cast of the cell population
    MutableVertexMesh<SPACE_DIM,SPACE_DIM>* p_mesh = static_cast<MutableVertexMesh<SPACE_DIM,SPACE_DIM>*>(&(pCellPopulation->rGetMesh()));
    double wound_area = 0.0;
    unsigned num_wounds = 0;
    std::vector<double> wound_local_areas;

    std::set<unsigned> all_wound_node_indices;
    for (typename MutableVertexMesh<SPACE_DIM,SPACE_DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
            node_iter != p_mesh->GetNodeIteratorEnd();
            ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            c_vector<double,2> node_location = node_iter->rGetLocation();
            double x = node_location[0];
            double y = node_location[1];
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if ( (containing_element_indices.size() == 1 || containing_element_indices.size() == 2) &&
                    (x > 0.75) && (x < (p_mesh->GetWidth(0) - 0.75)) && (y > 0.75) && (y < (p_mesh->GetWidth(1) - 0.75)) )
            {
                all_wound_node_indices.insert(node_index);
            }
        }

    if (all_wound_node_indices.empty())
    {
        // If there are no boundary nodes, check for wound cells.
        // Loop over cells, if cell is mutant, GetVolume() and add it to wound area.
        for (typename MutableVertexMesh<SPACE_DIM,SPACE_DIM>::VertexElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
         elem_iter != p_mesh->GetElementIteratorEnd();
         ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            CellPtr p_cell = pCellPopulation->GetCellUsingLocationIndex(elem_index);
            if (p_cell->GetMutationState()->IsType<WoundCellMutationState>())
            {
                double cell_area = pCellPopulation->rGetMesh().GetVolumeOfElement(elem_index);
                wound_area += cell_area;
                num_wounds += 1;
                wound_local_areas.push_back(cell_area);
            }
        }
    }
    else
    {
        // 'Standard' boundary node implementation
        while (!all_wound_node_indices.empty())
        {
            std::vector<unsigned> wound_node_indices;
            unsigned node_index = *(all_wound_node_indices.begin());
            wound_node_indices.push_back(node_index);
            all_wound_node_indices.erase(node_index);
            num_wounds += 1;
            double local_wound_area = 0.0;
            c_vector<double, SPACE_DIM> first_node_location = p_mesh->GetNode(node_index)->rGetLocation();
            c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

            // Find the nodes defining the wound, store in anti-clockwise orientation. Keep track of area
            bool repeat = false;
            while (!repeat)
            {
                std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
                for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    iter++)
                {
                    VertexElement<SPACE_DIM, SPACE_DIM>* p_element = p_mesh->GetElement(*iter);
                    // Find the local index of this node in this element
                    unsigned local_index = p_element->GetNodeLocalIndex(node_index);
                    // Get the previous node in this element
                    unsigned num_nodes_elem = p_element->GetNumNodes();
                    unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                    Node<SPACE_DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
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
                        wound_node_indices.push_back(node_index);
                        all_wound_node_indices.erase(node_index);
                        c_vector<double, SPACE_DIM> next_node_location = p_previous_node->rGetLocation();
                        c_vector<double, SPACE_DIM> pos_2 = p_mesh->GetVectorFromAtoB(first_node_location, next_node_location);

                        double this_x = pos_1[0];
                        double this_y = pos_1[1];
                        double next_x = pos_2[0];
                        double next_y = pos_2[1];

                        local_wound_area += 0.5 * (this_x * next_y - next_x * this_y);

                        pos_1 = pos_2;
                        break;
                    }
                }
            }
            wound_area += local_wound_area;
            wound_local_areas.push_back(local_wound_area);
        }
    }
    *this->mpOutStream << wound_area << " ";
    *this->mpOutStream << num_wounds << " ";
    for (unsigned local_index=0; local_index<wound_local_areas.size(); local_index++)
    {
        *this->mpOutStream << wound_local_areas[local_index] << " ";
    }
}

// Explicit instantiation
template class WoundAreaWriter<1,1>;
template class WoundAreaWriter<1,2>;
template class WoundAreaWriter<2,2>;
template class WoundAreaWriter<1,3>;
template class WoundAreaWriter<2,3>;
template class WoundAreaWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(WoundAreaWriter)

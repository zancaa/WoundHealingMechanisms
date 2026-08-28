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

#include "BoundaryNormalForce.hpp"
#include "MathsCustomFunctions.hpp"
#include "Exception.hpp"
#include "Debug.hpp"

template<unsigned DIM>
BoundaryNormalForce<DIM>::BoundaryNormalForce()
   : AbstractForce<DIM>(),
     mForceStrength(1.0)
{
}

template<unsigned DIM>
BoundaryNormalForce<DIM>::~BoundaryNormalForce()
{
}

template<unsigned DIM>
void BoundaryNormalForce<DIM>::SetForceStrength(double forceStrength)
{
    assert(forceStrength >= 0.0);
    mForceStrength = forceStrength;
}

template<unsigned DIM>
double BoundaryNormalForce<DIM>::GetForceStrength()
{
    return mForceStrength;
}

template<unsigned DIM>
void BoundaryNormalForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if( (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr) )
    {
        EXCEPTION("BoundaryNormalForce is to be used with VertexBasedCellPopulations only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    for (typename MutableVertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
            {
                if ( node_iter->IsBoundaryNode() )
                {
                    unsigned node_index = node_iter->GetIndex();
                    c_vector<double,DIM> node_location = node_iter->rGetLocation();
                    double x = node_location[0];
                    double y = node_location[1];
                    if ( (x > 0.75) && (x < (p_mesh->GetWidth(0) - 0.75)) && (y > 0.75) && (y < (p_mesh->GetWidth(1) - 0.75)) )
                    {
                        // Get the indices of this node's containing elements
                        std::set<unsigned> elements_containing_node = node_iter->rGetContainingElementIndices();

                        if (elements_containing_node.size() == 1)
                        {
                            // If node is only in one cell, then both of it's neighbours should also
                            // be boundary nodes.
                            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*(elements_containing_node.begin()));
                            unsigned num_nodes_elem = p_element->GetNumNodes();
                            // Find the local index of this node in this element
                            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

                            // Get the previous and next nodes in this element
                            // and check that they are in fact boundary nodes.
                            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
                            c_vector<double,DIM> prev_node_location = p_previous_node->rGetLocation();

                            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);
                            c_vector<double,DIM> next_node_location = p_next_node->rGetLocation();

                            if (!(p_previous_node->IsBoundaryNode()) || !(p_next_node->IsBoundaryNode()))
                            {
                                EXCEPTION("Have a node in one cell with non-boundary neighbours");
                            }

                            c_vector<double, DIM> clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(node_location, prev_node_location);
                            double clockwise_length = norm_2(clockwise_unit_vector);
                            clockwise_unit_vector /= norm_2(clockwise_unit_vector);

                            c_vector<double, DIM> anti_clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(node_location, next_node_location);
                            double anti_clockwise_length = norm_2(anti_clockwise_unit_vector);
                            anti_clockwise_unit_vector /= norm_2(anti_clockwise_unit_vector);

                            double sum_vector_lengths = clockwise_length + anti_clockwise_length;
                            // Calculate the outward normal at the node
                            c_vector<double, DIM> outward_normal = -clockwise_length*clockwise_unit_vector/sum_vector_lengths - anti_clockwise_length*anti_clockwise_unit_vector/sum_vector_lengths;
                            outward_normal /= norm_2(outward_normal);

                            // Apply force in outward normal direction
                            node_iter->AddAppliedForceContribution(mForceStrength*outward_normal);
                        }
                        else if (elements_containing_node.size() == 2)
                        {
                            // If node is in two elements, have to find the relevant elements and nodes
                            c_vector<double, DIM> clockwise_unit_vector;
                            double clockwise_length;
                            c_vector<double, DIM> anti_clockwise_unit_vector;
                            double anti_clockwise_length;
                            c_vector<double, DIM> outward_normal;
                            double sum_vector_lengths;
                            for (std::set<unsigned>::iterator iter = elements_containing_node.begin();
                                iter != elements_containing_node.end();
                                ++iter)
                            {
                                // For each element, check which neighbour is a boundary node
                                VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
                                unsigned num_nodes_elem = p_element->GetNumNodes();
                                // Find the local index of this node in this element
                                unsigned local_index = p_element->GetNodeLocalIndex(node_index);

                                // Get the previous and next nodes in this element and check if it is a boudary node
                                unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                                Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

                                unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                                Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

                                if ( p_previous_node->IsBoundaryNode() )
                                {
                                    c_vector<double,DIM> prev_node_location = p_previous_node->rGetLocation();
                                    clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(node_location, prev_node_location);
                                    clockwise_length = norm_2(clockwise_unit_vector);
                                    clockwise_unit_vector /= norm_2(clockwise_unit_vector);
                                }
                                else if ( p_next_node->IsBoundaryNode() )
                                {
                                    c_vector<double,DIM> next_node_location = p_next_node->rGetLocation();
                                    anti_clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(node_location, next_node_location);
                                    anti_clockwise_length = norm_2(anti_clockwise_unit_vector);
                                    anti_clockwise_unit_vector /= norm_2(anti_clockwise_unit_vector);
                                }
                                else
                                {
                                    EXCEPTION("Boundary node with no boundary node neighbours");
                                }
                            }
                            sum_vector_lengths = clockwise_length + anti_clockwise_length;
                            // Calculate the outward normal at the node
                            outward_normal = clockwise_length*clockwise_unit_vector/sum_vector_lengths + anti_clockwise_length*anti_clockwise_unit_vector/sum_vector_lengths;
                            outward_normal /= norm_2(outward_normal);

                            // Apply force in outward normal direction
                            node_iter->AddAppliedForceContribution(mForceStrength*outward_normal);
                        }
                        else
                        {
                            EXCEPTION("Have a boundary node contained in more than two cells");
                        }
                    }
                }
            }
        
}

template<unsigned DIM>
void BoundaryNormalForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ForceStrength>" << mForceStrength << "</ForceStrength>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class BoundaryNormalForce<1>;
template class BoundaryNormalForce<2>;
template class BoundaryNormalForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BoundaryNormalForce)

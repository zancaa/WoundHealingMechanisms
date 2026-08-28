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

#include "VertexBoundaryRefinementModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
VertexBoundaryRefinementModifier<DIM>::VertexBoundaryRefinementModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mMaxEdgeLength(0.25),
      mMinEdgeLength(0.0)
{
}

template<unsigned DIM>
VertexBoundaryRefinementModifier<DIM>::~VertexBoundaryRefinementModifier()
{
}

template<unsigned DIM>
void VertexBoundaryRefinementModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    RefineEdges(rCellPopulation);
}

template<unsigned DIM>
void VertexBoundaryRefinementModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We also refine at the start.
     */
    RefineEdges(rCellPopulation);
}

template<unsigned DIM>
void VertexBoundaryRefinementModifier<DIM>::RefineEdges(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("VertexBoundaryRefinementModifier is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
           
    bool recheck_edges = true;

    while (recheck_edges)
    {
        recheck_edges = false;

        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
            elem_iter != p_mesh->GetElementIteratorEnd();
            ++elem_iter)
        {
            //unsigned elem_index = elem_iter->GetIndex();

            //unsigned num_nodes = elem_iter->GetNumNodes();
            for (unsigned node_local_index = 0; node_local_index < elem_iter->GetNumNodes(); node_local_index++)
            {
                unsigned next_node_local_index = (node_local_index+1) % (elem_iter->GetNumNodes());

    //PRINT_3_VARIABLES(node_local_index,next_node_local_index,num_nodes);

                unsigned node_global_index = elem_iter->GetNodeGlobalIndex(node_local_index);
                unsigned next_node_global_index = elem_iter->GetNodeGlobalIndex(next_node_local_index);

                Node<DIM>* p_node_a = p_mesh->GetNode(node_global_index);
                Node<DIM>* p_node_b = p_mesh->GetNode(next_node_global_index);
                                
                if (p_node_a->IsBoundaryNode() && p_node_b->IsBoundaryNode())
                {
                    // Find the sets of elements containing nodes A and B
                    std::set<unsigned> node_a_elem_indices = p_node_a->rGetContainingElementIndices();
                    std::set<unsigned> node_b_elem_indices = p_node_b->rGetContainingElementIndices();

                    // Find common elements
                    std::set<unsigned> shared_elements;
                    std::set_intersection(node_a_elem_indices.begin(),
                            node_a_elem_indices.end(),
                            node_b_elem_indices.begin(),
                            node_b_elem_indices.end(),
                            std::inserter(shared_elements, shared_elements.begin()));

                    assert(shared_elements.size()>0); //otherwise not in the same element at all 

                    if(shared_elements.size() == 1)
                    {
                        // Here we have a boundary edge so add new node if needed.
                        c_vector<double,DIM> edge = p_mesh->GetVectorFromAtoB(p_node_a->rGetLocation(),p_node_b->rGetLocation());
                        if (norm_2(edge) > mMaxEdgeLength)
                        {
                            p_mesh->DivideEdge(p_node_a, p_node_b);
                            recheck_edges = true;
                        } 
                        else if (norm_2(edge) < mMinEdgeLength)
                        {
                            NEVER_REACHED; // Not implemented yet
                        }   
                    }
                }
            }
        }
    }
}


template<unsigned DIM>
double VertexBoundaryRefinementModifier<DIM>::GetMaxEdgeLength()
{
    return mMaxEdgeLength;
}

template<unsigned DIM>
void VertexBoundaryRefinementModifier<DIM>::SetMaxEdgeLength(double maxEdgeLength)
{
    assert(maxEdgeLength > 0.0);
    mMaxEdgeLength = maxEdgeLength;
}

template<unsigned DIM>
double VertexBoundaryRefinementModifier<DIM>::GetMinEdgeLength()
{
    return mMinEdgeLength;
}

template<unsigned DIM>
void VertexBoundaryRefinementModifier<DIM>::SetMinEdgeLength(double minEdgeLength)
{
    assert(minEdgeLength > 0.0);
    mMinEdgeLength = minEdgeLength;
}



template<unsigned DIM>
void VertexBoundaryRefinementModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaxEdgeLength>" << mMaxEdgeLength << "</MaxEdgeLength>\n";
    *rParamsFile << "\t\t\t<MinEdgeLength>" << mMinEdgeLength << "</MinEdgeLength>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class VertexBoundaryRefinementModifier<1>;
template class VertexBoundaryRefinementModifier<2>;
template class VertexBoundaryRefinementModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBoundaryRefinementModifier)
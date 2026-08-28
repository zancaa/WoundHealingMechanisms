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

#include "CrawlingForce.hpp"
#include "CellLabel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
CrawlingForce<DIM>::CrawlingForce()
    : AbstractForce<DIM>(),
    mActiveForceParameter(1.0)
{
}

template<unsigned DIM>
CrawlingForce<DIM>::~CrawlingForce()
{
}

template<unsigned DIM>
void CrawlingForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("CrawlingForce is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Extra force applied to cells on boundary. Cells moving towards centre of wound.
    /* Add the same y-component of force to every boundary node - this force should be the average of 
     * y-component of the forces on the boundary nodes plus the external applied force.
     */
    std::vector<double> boundary_forces;
    for (unsigned node_index=0; node_index<p_cell_population->GetNumNodes(); node_index++)
    {
        // Get pointer to this node
        Node<DIM>* p_node = p_cell_population->GetNode(node_index);
        c_vector<double, DIM> old_location = p_node->rGetLocation();
        double node_height = old_location(1);

        // Check if node is on the boundary of wound (a boundary node that is not at the external boundaries)
        if ((p_node->IsBoundaryNode()) && (node_height > 1.5))
        {
            c_vector<double, DIM> applied_force = scalar_vector<double>(DIM, DOUBLE_UNSET);
            applied_force = p_cell_population->GetNode(node_index)->rGetAppliedForce();
            boundary_forces.push_back(applied_force(1));
        }
    }

    // Calculate the average y-component of the forces on the boundary nodes
    double mean_y_boundary_force = std::accumulate(boundary_forces.begin(), boundary_forces.end(), 0.0) / boundary_forces.size();

    for (unsigned node_index=0; node_index<p_cell_population->GetNumNodes(); node_index++)
    {
        // Get pointer to this node
        Node<DIM>* p_node = p_cell_population->GetNode(node_index);
        c_vector<double, DIM> old_location = p_node->rGetLocation();
        double node_height = old_location(1);
        if ((p_node->IsBoundaryNode()) && (node_height > 1.5))
        {
            c_vector<double, DIM> force = scalar_vector<double>(DIM, DOUBLE_UNSET);
            force = p_node->rGetAppliedForce();
            double up_force = GetActiveForceParameter();
            force(1) = mean_y_boundary_force + up_force;
            p_node->ClearAppliedForce();
            p_node->AddAppliedForceContribution(force);
        }
    }
}

template<unsigned DIM>
double CrawlingForce<DIM>::GetActiveForceParameter()
{
    return mActiveForceParameter;
}

template<unsigned DIM>
void CrawlingForce<DIM>::SetActiveForceParameter(double activeForceParameter)
{
    mActiveForceParameter = activeForceParameter;
}

template<unsigned DIM>
void CrawlingForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ActiveForceParameter>" << mActiveForceParameter << "</ActiveForceParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class CrawlingForce<1>;
template class CrawlingForce<2>;
template class CrawlingForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CrawlingForce)
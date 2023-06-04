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

#include "WoundCentreForce.hpp"
#include "MathsCustomFunctions.hpp"
#include "Exception.hpp"
#include "Debug.hpp"

template<unsigned DIM>
WoundCentreForce<DIM>::WoundCentreForce()
   : AbstractForce<DIM>(),
     mForceStrength(1.0)
{
}

template<unsigned DIM>
WoundCentreForce<DIM>::~WoundCentreForce()
{
}

template<unsigned DIM>
void WoundCentreForce<DIM>::SetForceStrength(double forceStrength)
{
    assert(forceStrength >= 0.0);
    mForceStrength = forceStrength;
}

template<unsigned DIM>
double WoundCentreForce<DIM>::GetForceStrength()
{
    return mForceStrength;
}

template<unsigned DIM>
void WoundCentreForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if( (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr) )
    {
        EXCEPTION("WoundCentreForce is to be used with VertexBasedCellPopulations only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    double sum_x_locations = 0.0;
    double sum_y_locations = 0.0;

    std::set<unsigned> all_wound_node_indices;
    for (typename MutableVertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
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
                sum_x_locations += x;
                sum_y_locations += y;
            }
        }

    // Find the centre of the wound - currently not generalisable to multiple wounds
    double mean_x_boundary_location = sum_x_locations/all_wound_node_indices.size();
    double mean_y_boundary_location = sum_y_locations/all_wound_node_indices.size();

    /* Add force contribution in the direction of the wound centre to boundary nodes
     * Currently a constant force applied in the direction of the wound centre,
     * could also consider radial force whose strength depends on distance to 
     * the wound centre, for example. 
     */
    for (std::set<unsigned>::iterator iter = all_wound_node_indices.begin();
        iter != all_wound_node_indices.end();
        iter++)
    {
        c_vector<double,2> node_location = p_mesh->GetNode(*iter)->rGetLocation();
        double x = node_location[0];
        double y = node_location[1];
        c_vector<double,2> force_direction;
        force_direction(0) = -x + mean_x_boundary_location;
        force_direction(1) = -y + mean_y_boundary_location;
        p_mesh->GetNode(*iter)->AddAppliedForceContribution(mForceStrength*force_direction/norm_2(force_direction));
    }
}

template<unsigned DIM>
void WoundCentreForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ForceStrength>" << mForceStrength << "</ForceStrength>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class WoundCentreForce<1>;
template class WoundCentreForce<2>;
template class WoundCentreForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WoundCentreForce)

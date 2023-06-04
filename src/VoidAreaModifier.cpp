/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "VoidAreaModifier.hpp"

#include "NodeBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "VertexBasedCellPopulation.hpp"


#include "Debug.hpp"



template<unsigned DIM>
VoidAreaModifier<DIM>::VoidAreaModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mCutoff(1.5),
    mPixelSeparation(0.05),
    mPlotPixelContour(false)
{
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;

    OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream locationFile = output_file_handler.OpenOutputFile("voidArea.dat");
    *locationFile << "time \t";
    *locationFile << "PixelVoidArea" << "\t";
    *locationFile << "CellVoidArea" << "\t";
    *locationFile << "\n";
    locationFile->close();

    OutputFileHandler output_file_handler_2(mOutputDirectory+"/", false);
    out_stream locationFile_2 = output_file_handler_2.OpenOutputFile("voidContour.dat");
    locationFile_2->close();

    OutputFileHandler output_file_handler_3(mOutputDirectory+"/", false);
    out_stream locationFile_3 = output_file_handler_3.OpenOutputFile("voidInitialTissue.dat");
    locationFile_3->close();
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::SetCutoff(double cutoff)
{
	mCutoff = cutoff;
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::SetPixelSeparation(double pixelSeparation)
{
	mPixelSeparation = pixelSeparation;
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::SetPlotPixelContour(bool plotPixelContour)
{
	mPlotPixelContour = plotPixelContour;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
VoidAreaModifier<DIM>::~VoidAreaModifier()
{
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if(DIM == 2)
    {   
        bool plot_contour = false;
        bool plot_initial = false;

        double pixel_void_area = 0.0;
        double cell_void_area = 0.0;

        // both NodeBased and MeshBased (with no ghosts) populations are treated the same here, as both just have a cutoff radius
        // if( (bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)) 
        // || bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
        // && !(bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))) )
        OutputFileHandler output_file_handler_2(mOutputDirectory+"/", false);
        out_stream locationFile_2 = output_file_handler_2.OpenOutputFile("voidContour.dat", std::ios::app);
        SimulationTime* p_time = SimulationTime::Instance();
        *locationFile_2 << p_time->GetTime() << "\t";
        

        OutputFileHandler output_file_handler_3(mOutputDirectory+"/", false);
        out_stream locationFile_3 = output_file_handler_3.OpenOutputFile("voidInitialTissue.dat", std::ios::app);
        *locationFile_3 << p_time->GetTime() << "\t";


        
        {

            double pixel_radial_reach = 0.5*mCutoff;
            double separation_between_pixels = mPixelSeparation;
            double area_of_pixel = pow(separation_between_pixels,2);

            // minus 2 is to ensure we are in the tissue and not detecting gaps due to periodicity from each side
            unsigned pixel_tissue_width = (rCellPopulation.rGetMesh().GetWidth(0) - 2)/separation_between_pixels;
            unsigned pixel_tissue_depth = (rCellPopulation.rGetMesh().GetWidth(1) - 2)/separation_between_pixels;
            
            // unsigned pixel_grid[pixel_tissue_width][pixel_tissue_depth];
            std::vector<std::vector<unsigned>> pixel_grid(pixel_tissue_width,std::vector<unsigned>(pixel_tissue_depth));

            // int number_of_pixels_with_cell = 0;

            for(unsigned pixel_i = 0; pixel_i<pixel_tissue_width; pixel_i++)
            {
                for(unsigned pixel_j = 0; pixel_j<pixel_tissue_depth; pixel_j++)
                {
                    pixel_grid[pixel_i][pixel_j] = 0;

                    double pixel_x_coordinate = 1.0 + pixel_i*separation_between_pixels;
                    double pixel_y_coordinate = 1.0 + pixel_j*separation_between_pixels;
                    bool does_pixel_contain_cell = false;
                    // Iterate over cell population
                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter)
                    {
                        c_vector<double, DIM>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                        if( pow(cell_location[0] - pixel_x_coordinate,2) + pow(cell_location[1] - pixel_y_coordinate,2) < pow(pixel_radial_reach,2) )
                        {
                            pixel_grid[pixel_i][pixel_j] = 1;
                            does_pixel_contain_cell = true;
                            break;
                        }
                    }

                    // If pixel is within proximity of a cell, then it is in the void.
                    if(!does_pixel_contain_cell)
                    {
                        pixel_void_area = pixel_void_area + area_of_pixel;
                    }
                }

            }

            // for(unsigned pixel_i = 0; pixel_i<pixel_tissue_width; pixel_i++)
            // {
            //     for(unsigned pixel_j = 0; pixel_j<pixel_tissue_depth; pixel_j++)
            //     {
                    
            //         if(pixel_grid[pixel_i][pixel_j] == 1)
            //         {
            //             unsigned number_of_neigh_pixels = pixel_grid[pixel_i-1][pixel_j] + pixel_grid[pixel_i+1][pixel_j] + pixel_grid[pixel_i][pixel_j-1] + pixel_grid[pixel_i][pixel_j+1];
            //             if(number_of_neigh_pixels != 4)
            //             {
            //                 double pixel_x_coordinate = 1.0 + pixel_i*separation_between_pixels;
            //                 double pixel_y_coordinate = 1.0 + pixel_j*separation_between_pixels;
            //             }
            //         }
                    

            //     }
            // }

            if(plot_contour)
            {
                // Mark the boundary nodes
                if((bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))) && !(bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))) && mPlotPixelContour==false)
                {
                    MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

                    VertexMesh<DIM,DIM>* voronoi_tesselation = p_cell_population->GetVoronoiTessellation();

                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter) 
                    {

                        cell_iter->GetCellData()->SetItem("is_boundary", 0.0);
                                                
                        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                        unsigned element_index = voronoi_tesselation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
                        VertexElement<DIM,DIM>* p_element = voronoi_tesselation->GetElement(element_index);

                        std::set<unsigned> cell_neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

                        for (unsigned i=0; i<p_element->GetNumNodes(); i++)
                        {
                            unsigned node_global_index_1_a = p_element->GetNodeGlobalIndex(i);
                            unsigned node_global_index_1_b = p_element->GetNodeGlobalIndex((i+1)%(p_element->GetNumNodes()));

                            if(plot_initial)
                            {
                                // Used for initial tissue plot, but can do whole tissue if wanted...
                                Node<DIM>* p_node_a = voronoi_tesselation->GetNode(node_global_index_1_a);
                                Node<DIM>* p_node_b = voronoi_tesselation->GetNode(node_global_index_1_b);

                                double ax = p_node_a->rGetLocation()[0];
                                double ay = p_node_a->rGetLocation()[1];
                                
                                double bx = p_node_b->rGetLocation()[0];
                                double by = p_node_b->rGetLocation()[1];
                                
                                    if((ax-bx) > 0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        ax = ax - rCellPopulation.rGetMesh().GetWidth(0);
                                    }
                                    else if((ax-bx) < -0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        bx = bx - rCellPopulation.rGetMesh().GetWidth(0);
                                    }

                                    if((ay-by) > 0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        ay = ay - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                    else if((ay-by) < -0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        by = by - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                
                                *locationFile_3 << ax << "\t";
                                *locationFile_3 << ay << "\t";
                                *locationFile_3 << bx << "\t";
                                *locationFile_3 << by << "\t";
                            }

                            bool are_nodes_ab_shared_with_neighbours = false;                                

                                // Iterate over these neighbours
                            for (std::set<unsigned>::iterator neighbour_iter = cell_neighbour_indices.begin();
                                neighbour_iter != cell_neighbour_indices.end();
                                ++neighbour_iter)
                            {
                                unsigned neighbour_index = *neighbour_iter;
                                                            
                                unsigned element_index_2 = voronoi_tesselation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(neighbour_index);
                                VertexElement<DIM,DIM>* p_element_2 = voronoi_tesselation->GetElement(element_index_2);

                                for (unsigned j=0; j<p_element_2->GetNumNodes(); j++)
                                {
                                    unsigned node_global_index_2_a = p_element_2->GetNodeGlobalIndex(j);
                                    unsigned node_global_index_2_b = p_element_2->GetNodeGlobalIndex((j+1)%(p_element_2->GetNumNodes()));

                                        if(node_global_index_2_a == node_global_index_1_a)
                                        {
                                            if(node_global_index_2_b == node_global_index_1_b)
                                            {
                                                are_nodes_ab_shared_with_neighbours = true;
                                                break;
                                            }
                                        }
                                        else if(node_global_index_2_b == node_global_index_1_a)
                                        {
                                            if(node_global_index_2_a == node_global_index_1_b)
                                            {
                                                are_nodes_ab_shared_with_neighbours = true;
                                                break;
                                            }
                                        }                                

                                }
                                if(are_nodes_ab_shared_with_neighbours)
                                {
                                    break;
                                }

                            }

                            if(are_nodes_ab_shared_with_neighbours == false)
                            {
                                Node<DIM>* p_node_a = voronoi_tesselation->GetNode(node_global_index_1_a);
                                Node<DIM>* p_node_b = voronoi_tesselation->GetNode(node_global_index_1_b);

                                *locationFile_2 << p_node_a->rGetLocation()[0] << "\t";
                                *locationFile_2 << p_node_a->rGetLocation()[1] << "\t";
                                
                                *locationFile_2 << p_node_b->rGetLocation()[0] << "\t";
                                *locationFile_2 << p_node_b->rGetLocation()[1] << "\t";

                                cell_iter->GetCellData()->SetItem("is_boundary", 1.0);

                            }
                        }

                    }

                }
                else if(bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation)))
                {
                    // for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                    //     cell_iter != rCellPopulation.End();
                    //     ++cell_iter)
                    // {
                    //     cell_iter->GetCellData()->SetItem("is_boundary", 0.0);
                    // }

                    MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);
                    // MutableMesh<DIM,DIM>* p_mesh = static_cast<MutableMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

                    VertexMesh<DIM,DIM>* voronoi_tesselation = p_cell_population->GetVoronoiTessellation();

                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter) 
                    {
                        cell_iter->GetCellData()->SetItem("is_boundary", 0.0);
                        
                        // Get the location index corresponding to this cell
                        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

                        // Get the set of neighbouring location indices
                        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

                        // Iterate over these neighbours
                        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                            neighbour_iter != neighbour_indices.end();
                            ++neighbour_iter)
                        {
                            // Get the length of the edge shared with this neighbour
                            unsigned neighbour_index = *neighbour_iter;

                            if(plot_initial)
                            {
                                // Used for initial tissue plot, but can do whole tissue if wanted...

                                unsigned elementIndex1 = voronoi_tesselation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
                                unsigned elementIndex2 = voronoi_tesselation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(neighbour_index);

                                std::set<unsigned> node_indices_1;
                                for (unsigned i = 0; i < voronoi_tesselation->GetElement(elementIndex1)->GetNumNodes(); i++)
                                {
                                    node_indices_1.insert(voronoi_tesselation->GetElement(elementIndex1)->GetNodeGlobalIndex(i));
                                }
                                std::set<unsigned> node_indices_2;
                                for (unsigned i = 0; i < voronoi_tesselation->GetElement(elementIndex2)->GetNumNodes(); i++)
                                {
                                    node_indices_2.insert(voronoi_tesselation->GetElement(elementIndex2)->GetNodeGlobalIndex(i));
                                }

                                std::set<unsigned> shared_nodes;
                                std::set_intersection(node_indices_1.begin(), node_indices_1.end(),
                                                    node_indices_2.begin(), node_indices_2.end(),
                                                    std::inserter(shared_nodes, shared_nodes.begin()));
                                unsigned index1 = *(shared_nodes.begin());
                                unsigned index2 = *(++(shared_nodes.begin()));

                                c_vector<double, DIM> node_1 = voronoi_tesselation->GetNode(index1)->rGetLocation();
                                c_vector<double, DIM> node_2 = voronoi_tesselation->GetNode(index2)->rGetLocation();

                                if(plot_initial)
                                {
                                    // Used for initial tissue plot, but can do whole tissue if wanted...
                                    double ax = node_1[0];
                                    double ay = node_1[1];
                                    
                                    double bx = node_2[0];
                                    double by = node_2[1];
                                    
                                    if((ax-bx) > 0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        ax = ax - rCellPopulation.rGetMesh().GetWidth(0);
                                    }
                                    else if((ax-bx) < -0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        bx = bx - rCellPopulation.rGetMesh().GetWidth(0);
                                    }

                                    if((ay-by) > 0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        ay = ay - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                    else if((ay-by) < -0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        by = by - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                    
                                    *locationFile_3 << ax << "\t";
                                    *locationFile_3 << ay << "\t";
                                    *locationFile_3 << bx << "\t";
                                    *locationFile_3 << by << "\t";
                                }

                            }

                            if (p_cell_population->IsGhostNode(neighbour_index))
                            {
                                cell_iter->GetCellData()->SetItem("is_boundary", 1.0);

                                unsigned elementIndex1 = voronoi_tesselation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
                                unsigned elementIndex2 = voronoi_tesselation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(neighbour_index);

                                std::set<unsigned> node_indices_1;
                                for (unsigned i = 0; i < voronoi_tesselation->GetElement(elementIndex1)->GetNumNodes(); i++)
                                {
                                    node_indices_1.insert(voronoi_tesselation->GetElement(elementIndex1)->GetNodeGlobalIndex(i));
                                }
                                std::set<unsigned> node_indices_2;
                                for (unsigned i = 0; i < voronoi_tesselation->GetElement(elementIndex2)->GetNumNodes(); i++)
                                {
                                    node_indices_2.insert(voronoi_tesselation->GetElement(elementIndex2)->GetNodeGlobalIndex(i));
                                }

                                std::set<unsigned> shared_nodes;
                                std::set_intersection(node_indices_1.begin(), node_indices_1.end(),
                                                    node_indices_2.begin(), node_indices_2.end(),
                                                    std::inserter(shared_nodes, shared_nodes.begin()));
                                unsigned index1 = *(shared_nodes.begin());
                                unsigned index2 = *(++(shared_nodes.begin()));

                                c_vector<double, DIM> node_1 = voronoi_tesselation->GetNode(index1)->rGetLocation();
                                c_vector<double, DIM> node_2 = voronoi_tesselation->GetNode(index2)->rGetLocation();

                                *locationFile_2 << node_1[0] << "\t";
                                *locationFile_2 << node_1[1] << "\t";

                                *locationFile_2 << node_2[0] << "\t";
                                *locationFile_2 << node_2[1] << "\t";
                            }
                        }
                    }

                    // for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = p_mesh->EdgesBegin();
                    //     edge_iterator != p_mesh->EdgesEnd();
                    //     ++edge_iterator)
                    // {
                    //     unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
                    //     unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

                    //     if(p_cell_population->IsGhostNode(nodeA_global_index))
                    //     {
                    //         if(!(p_cell_population->IsGhostNode(nodeB_global_index)))
                    //         {
                    //             CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeB_global_index);
                    //             p_cell->GetCellData()->SetItem("is_boundary", 1.0);
                    //         }
                    //     }
                    //     else if(p_cell_population->IsGhostNode(nodeB_global_index))
                    //     {
                    //         if(!(p_cell_population->IsGhostNode(nodeA_global_index)))
                    //         {
                    //             CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeA_global_index);
                    //             p_cell->GetCellData()->SetItem("is_boundary", 1.0);
                    //         }

                    //     }
                    // }
                }
                else if(bool(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)))
                {
                    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
                    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter)
                    {
                        // bool check_node = true;
                        cell_iter->GetCellData()->SetItem("is_boundary", 0.0);

                        VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
                        for (unsigned i=0; i<p_element->GetNumNodes(); i++)
                        {
                            unsigned node_global_index_a = p_element->GetNodeGlobalIndex(i);
                            Node<DIM>* p_node_a = p_mesh->GetNode(node_global_index_a);

                            unsigned node_global_index_b = p_element->GetNodeGlobalIndex((i+1)%(p_element->GetNumNodes()) );
                            Node<DIM>* p_node_b = p_mesh->GetNode(node_global_index_b);

                            if(plot_initial)
                            {
                                // Used for initial tissue plot, but can do whole tissue if wanted...
                                double ax = p_node_a->rGetLocation()[0];
                                double ay = p_node_a->rGetLocation()[1];
                                
                                double bx = p_node_b->rGetLocation()[0];
                                double by = p_node_b->rGetLocation()[1];
                                
                                if((ax-bx) > 0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        ax = ax - rCellPopulation.rGetMesh().GetWidth(0);
                                    }
                                    else if((ax-bx) < -0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        bx = bx - rCellPopulation.rGetMesh().GetWidth(0);
                                    }

                                    if((ay-by) > 0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        ay = ay - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                    else if((ay-by) < -0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        by = by - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                
                                *locationFile_3 << ax << "\t";
                                *locationFile_3 << ay << "\t";
                                *locationFile_3 << bx << "\t";
                                *locationFile_3 << by << "\t";
                            }
                           

                            if (p_node_a->IsBoundaryNode() && p_node_b->IsBoundaryNode())
                            {
                                *locationFile_2 << p_node_a->rGetLocation()[0] << "\t";
                                *locationFile_2 << p_node_a->rGetLocation()[1] << "\t";

                                *locationFile_2 << p_node_b->rGetLocation()[0] << "\t";
                                *locationFile_2 << p_node_b->rGetLocation()[1] << "\t";

                                cell_iter->GetCellData()->SetItem("is_boundary", 1.0);
                                // check_node = false;
                            }
                        }
                    }
                }
                else if((bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))) || mPlotPixelContour==true)
                {
                // Iterate over cell population
                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter)
                    {
                        cell_iter->GetCellData()->SetItem("is_boundary", 0.0);

                        for(unsigned pixel_i = 1; pixel_i<pixel_tissue_width-1; pixel_i++)
                        {
                            bool break_out_loop = false;

                            for(unsigned pixel_j = 1; pixel_j<pixel_tissue_depth-1; pixel_j++)
                            {

                                if(pixel_grid[pixel_i][pixel_j] == 0)
                                {
                                    double pixel_x_coordinate = 1.0 + pixel_i*separation_between_pixels;
                                    double pixel_y_coordinate = 1.0 + pixel_j*separation_between_pixels;

                                    c_vector<double, DIM>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                                    
                                    if( pow(cell_location[0] - pixel_x_coordinate,2) + pow(cell_location[1] - pixel_y_coordinate,2) < pow(pixel_radial_reach+2*separation_between_pixels,2) )
                                    {
                                        cell_iter->GetCellData()->SetItem("is_boundary", 1.0);

                                        break_out_loop = true;
                                        break;
                                    }
                                }
                                if(pixel_grid[pixel_i][pixel_j] == 1)
                                {
                                    unsigned number_of_neigh_pixels = pixel_grid[pixel_i-1][pixel_j] + pixel_grid[pixel_i+1][pixel_j] + pixel_grid[pixel_i][pixel_j-1] + pixel_grid[pixel_i][pixel_j+1];
                                    if(number_of_neigh_pixels != 4)
                                    {
                                        double pixel_x_coordinate = 1.0 + pixel_i*separation_between_pixels;
                                        double pixel_y_coordinate = 1.0 + pixel_j*separation_between_pixels;

                                        *locationFile_2 << pixel_x_coordinate << "\t";
                                        *locationFile_2 << pixel_y_coordinate << "\t";
                                        
                                    }
                                }

                                if(break_out_loop)
                                {
                                    break;
                                }

                            }

                            if(break_out_loop)
                            {
                                break;
                            }

                        }
                    }
                }


                if((bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))) && mPlotPixelContour==true)
                {
                    MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

                    VertexMesh<DIM,DIM>* voronoi_tesselation = p_cell_population->GetVoronoiTessellation();

                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter) 
                    {
                                                
                        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                        unsigned element_index = voronoi_tesselation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
                        VertexElement<DIM,DIM>* p_element = voronoi_tesselation->GetElement(element_index);

                        std::set<unsigned> cell_neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

                        for (unsigned i=0; i<p_element->GetNumNodes(); i++)
                        {
                            unsigned node_global_index_1_a = p_element->GetNodeGlobalIndex(i);
                            unsigned node_global_index_1_b = p_element->GetNodeGlobalIndex((i+1)%(p_element->GetNumNodes()));

                            if(plot_initial)
                            {
                                // Used for initial tissue plot, but can do whole tissue if wanted...
                                Node<DIM>* p_node_a = voronoi_tesselation->GetNode(node_global_index_1_a);
                                Node<DIM>* p_node_b = voronoi_tesselation->GetNode(node_global_index_1_b);

                                double ax = p_node_a->rGetLocation()[0];
                                double ay = p_node_a->rGetLocation()[1];
                                
                                double bx = p_node_b->rGetLocation()[0];
                                double by = p_node_b->rGetLocation()[1];
                                
                                    if((ax-bx) > 0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        ax = ax - rCellPopulation.rGetMesh().GetWidth(0);
                                    }
                                    else if((ax-bx) < -0.5*rCellPopulation.rGetMesh().GetWidth(0))
                                    {
                                        bx = bx - rCellPopulation.rGetMesh().GetWidth(0);
                                    }

                                    if((ay-by) > 0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        ay = ay - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                    else if((ay-by) < -0.5*rCellPopulation.rGetMesh().GetWidth(1))
                                    {
                                        by = by - rCellPopulation.rGetMesh().GetWidth(1);
                                    }
                                
                                *locationFile_3 << ax << "\t";
                                *locationFile_3 << ay << "\t";
                                *locationFile_3 << bx << "\t";
                                *locationFile_3 << by << "\t";
                            }

                        }

                    }

                }
                else if((bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))) )
                {
                // Iterate over cell population
                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter)
                    {

                        c_vector<double, DIM>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                        *locationFile_3 << cell_location[0] << "\t";
                        *locationFile_3 << cell_location[1] << "\t";

                    }
                }
            }

        }
        *locationFile_3 << "\n";
        locationFile_3->close();

        *locationFile_2 << "\n";
        locationFile_2->close();

        // if(bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
        // {

        // }
        
        // We handle Vertex model a little differently as there is a very natural cell area defined
        // if(bool(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)) 
        // || bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation)))
        {
            double tissue_area = 0.0;

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
            {
                // Get the volume of this cell
                double cell_volume = rCellPopulation.GetVolumeOfCell(*cell_iter);

                tissue_area = tissue_area + cell_volume;
            }

            cell_void_area = (rCellPopulation.rGetMesh().GetWidth(0))*(rCellPopulation.rGetMesh().GetWidth(1)) - tissue_area;
            
        }

        OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
        out_stream locationFile = output_file_handler.OpenOutputFile("voidArea.dat", std::ios::app);
        // SimulationTime* p_time = SimulationTime::Instance();
        *locationFile << p_time->GetTime() << "\t";
        *locationFile << pixel_void_area << "\t";
        *locationFile << cell_void_area << "\t";
        *locationFile << "\n";
        locationFile->close();

    }

}



template<unsigned DIM>
void VoidAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void VoidAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class VoidAreaModifier<1>;
template class VoidAreaModifier<2>;
template class VoidAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VoidAreaModifier)


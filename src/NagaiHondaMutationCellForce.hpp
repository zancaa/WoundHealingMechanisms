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

#ifndef NAGAIHONDAMUTATIONCELLFORCE_HPP_
#define NAGAIHONDMUTATIONCELLFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NagaiHondaForce.hpp"

#include <iostream>

/**
 * A force class for use in vertex-based simulations, based on a model
 * model proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B
 * 81:699-719) to include different deformation energies between normal and
 * labelled cells. To include differential adhesion we override the
 * GetDeformationParameter() method.
 *
 * Each of the model parameter member variables are rescaled such that
 * mDampingConstantNormal takes the default value 1, whereas Nagai and
 * Honda (who denote the parameter by nu) take the value 0.01.
 */
template<unsigned DIM>
class NagaiHondaMutationCellForce  : public NagaiHondaForce<DIM>
{
private:

    /**
     * Deformation energy parameter for cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaCellDeformationEnergyParameter;

    /**
     * Deformation energy parameter for wound.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaWoundDeformationEnergyParameter;
    
    /**
     * Membrance surface tension energy parameter for cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaCellMembraneSurfaceEnergyParameter;

    /**
     * Membrane surface tension energy parameter for wound.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaWoundMembraneSurfaceEnergyParameter;

    /**
     * Cell wound adhesion energy parameter.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mNagaiHondaCellWoundAdhesionEnergyParameter;

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<NagaiHondaForce<DIM> >(*this);
        archive & mNagaiHondaCellDeformationEnergyParameter;
        archive & mNagaiHondaWoundDeformationEnergyParameter;
        archive & mNagaiHondaCellMembraneSurfaceEnergyParameter;
        archive & mNagaiHondaWoundMembraneSurfaceEnergyParameter;
        archive & mNagaiHondaCellWoundAdhesionEnergyParameter;
    }

public:

    /**
     * Constructor.
     */
    NagaiHondaMutationCellForce();

    /**
     * Destructor.
     */
    virtual ~NagaiHondaMutationCellForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the
     * Nagai Honda model.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden GetDeformationEnergyParameter() method.
     *
     * Get the deformation parameter for the cell. 
     *
     * @param pCell cell
     *
     * @return the deformation parameter for this cell.
     */
    virtual double GetNagaiHondaDeformationEnergyParameter(CellPtr pCell);

    /**
     * @return mNagaiHondaCellDeformationEnergyParameter
     */
    double GetNagaiHondaCellDeformationEnergyParameter();

    /**
     * @return mNagaiHondaWoundDeformationEnergyParameter
     */
    double GetNagaiHondaWoundDeformationEnergyParameter();

    /**
     * Overridden GetMembraneSurfaceEnergyParameter() method.
     *
     * Get the membrane surface tension parameter for the cell. 
     *
     * @param pCell cell
     *
     * @return the membrane surface tension parameter for this edge.
     */
    virtual double GetNagaiHondaMembraneSurfaceEnergyParameter(CellPtr pCell);

    /**
     * @return mNagaiHondaCellMembraneSurfaceEnergyParameter
     */
    double GetNagaiHondaCellMembraneSurfaceEnergyParameter();

    /**
     * @return mNagaiHondaWoundMembraneSurfaceEnergyParameter
     */
    double GetNagaiHondaWoundMembraneSurfaceEnergyParameter();

    /**
     * Overridden GetAdhesionParameter() method.
     *
     * Get the adhesion parameter for the edge between two given nodes. Depends
     * on the type of cells attached to the elements.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rVertexCellPopulation reference to the cell population
     *
     * @return the adhesion parameter for this edge.
     */
    virtual double GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

    /**
     * @return mNagaiHondaCellWoundAdhesionEnergyParameter
     */
    double GetNagaiHondaCellWoundAdhesionEnergyParameter();

    /**
     * Set mNagaiHondaCellDeformationEnergyParameter.
     *
     * @param cellDeformationEnergyParameter the new value of mNagaiHondaCellDeformationEnergyParameter
     */
    void SetNagaiHondaCellDeformationEnergyParameter(double cellDeformationEnergyParameter);

    /**
     * Set mNagaiHondaWoundDeformationEnergyParameter.
     *
     * @param woundDeformationEnergyParameter the new value of mNagaiHondaWoundDeformationEnergyParameter
     */
    void SetNagaiHondaWoundDeformationEnergyParameter(double woundDeformationEnergyParameter);

    /**
     * Set mNagaiHondaCellMembraneSurfaceEnergyParameter.
     *
     * @param cellMembraneSurfaceEnergyParameter the new value of mNagaiHondaCellMembraneSurfaceEnergyParameter
     */
    void SetNagaiHondaCellMembraneSurfaceEnergyParameter(double cellMembraneSurfaceEnergyParameter);

    /**
     * Set mNagaiHondaWoundMembraneSurfaceEnergyParameter.
     *
     * @param woundMembraneSurfaceEnergyParameter the new value of mNagaiHondaWoundMembraneSurfaceEnergyParameter
     */
    void SetNagaiHondaWoundMembraneSurfaceEnergyParameter(double woundMembraneSurfaceEnergyParameter);

    /**
     * Set mNagaiHondaCellWoundAdhesionEnergyParameter.
     *
     * @param cellWoundAdhesionEnergyParameter the new value of mNagaiHondaCellWoundAdhesionEnergyParameter
     */
    void SetNagaiHondaCellWoundAdhesionEnergyParameter(double cellWoundAdhesionEnergyParameter);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaMutationCellForce)

#endif /*NAGAIHONDAMUTATIONCELLFORCE_HPP_*/

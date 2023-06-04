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

#include "BernoulliTrialWithContactInhibitionCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

BernoulliTrialWithContactInhibitionCellCycleModel::BernoulliTrialWithContactInhibitionCellCycleModel()
    : AbstractCellCycleModel(),
      mQuiescentVolumeFraction(DOUBLE_UNSET),
      mEquilibriumVolume(DOUBLE_UNSET),
      mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
      mCurrentQuiescentDuration(0.0),
      mStemCellDivisionProbability(0.1),
      mStemCellMinimumDivisionAge(1.0),
      mTransitCellDivisionProbability(0.1),
      mTransitCellMinimumDivisionAge(2.0)
{
}

BernoulliTrialWithContactInhibitionCellCycleModel::BernoulliTrialWithContactInhibitionCellCycleModel(const BernoulliTrialWithContactInhibitionCellCycleModel& rModel)
   : AbstractCellCycleModel(rModel),
      mQuiescentVolumeFraction(rModel.mQuiescentVolumeFraction),
      mEquilibriumVolume(rModel.mEquilibriumVolume),
      mCurrentQuiescentOnsetTime(rModel.mCurrentQuiescentOnsetTime),
      mCurrentQuiescentDuration(rModel.mCurrentQuiescentDuration),
      mStemCellDivisionProbability(rModel.mStemCellDivisionProbability),
      mStemCellMinimumDivisionAge(rModel.mStemCellMinimumDivisionAge),
      mTransitCellDivisionProbability(rModel.mTransitCellDivisionProbability),
      mTransitCellMinimumDivisionAge(rModel.mTransitCellMinimumDivisionAge)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

bool BernoulliTrialWithContactInhibitionCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
    {
        // Check that cells are large enough to divide. 
        double cell_volume = mpCell->GetCellData()->GetItem("volume");
        double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

        if (cell_volume < quiescent_volume)
        {
            // Update the duration of the current period of contact inhibition.
            mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;

            // Update labels
            mpCell->RemoveCellProperty<CellLabel>();
            boost::shared_ptr<AbstractCellProperty> p_label =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
            mpCell->AddCellProperty(p_label);

            // Increase cell cycle duration
            double dt = SimulationTime::Instance()->GetTimeStep();
            mStemCellMinimumDivisionAge += dt;
            mTransitCellMinimumDivisionAge += dt;
        }
        else
        {
            // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
            mCurrentQuiescentDuration = 0.0;
            mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
            mpCell->RemoveCellProperty<CellLabel>();
        }
    }

    if (!mReadyToDivide)
    {
        double dt = SimulationTime::Instance()->GetTimeStep();
        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>() && !(mpCell->HasCellProperty<CellLabel>()))
        {
            if (GetAge() > mStemCellMinimumDivisionAge)
            {
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                if (p_gen->ranf() < mStemCellDivisionProbability*dt)
                {
                    mReadyToDivide = true;
                }
            }
        }
        if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && !(mpCell->HasCellProperty<CellLabel>()))
        {
            if (GetAge() > mTransitCellMinimumDivisionAge)
            {
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                if (p_gen->ranf() < mTransitCellDivisionProbability*dt)
                {
                    mReadyToDivide = true;
                }
            }
        }
    }

    return mReadyToDivide;
}

AbstractCellCycleModel* BernoulliTrialWithContactInhibitionCellCycleModel::CreateCellCycleModel()
{
    return new BernoulliTrialWithContactInhibitionCellCycleModel(*this);
}

void BernoulliTrialWithContactInhibitionCellCycleModel::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetQuiescentVolumeFraction() const
{
    return mQuiescentVolumeFraction;
}

void BernoulliTrialWithContactInhibitionCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetEquilibriumVolume() const
{
    return mEquilibriumVolume;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetCurrentQuiescentDuration() const
{
    return mCurrentQuiescentDuration;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetCurrentQuiescentOnsetTime() const
{
    return mCurrentQuiescentOnsetTime;
}

void BernoulliTrialWithContactInhibitionCellCycleModel::SetStemCellDivisionProbability(double stemCellDivisionProbability)
{
    mStemCellDivisionProbability = stemCellDivisionProbability;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetStemCellDivisionProbability()
{
    return mStemCellDivisionProbability;
}

void BernoulliTrialWithContactInhibitionCellCycleModel::SetStemCellMinimumDivisionAge(double stemCellMinimumDivisionAge)
{
    mStemCellMinimumDivisionAge = stemCellMinimumDivisionAge;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetStemCellMinimumDivisionAge()
{
    return mStemCellMinimumDivisionAge;
}

void BernoulliTrialWithContactInhibitionCellCycleModel::SetTransitCellDivisionProbability(double transitCellDivisionProbability)
{
    mTransitCellDivisionProbability = transitCellDivisionProbability;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetTransitCellDivisionProbability()
{
    return mTransitCellDivisionProbability;
}

void BernoulliTrialWithContactInhibitionCellCycleModel::SetTransitCellMinimumDivisionAge(double transitCellMinimumDivisionAge)
{
    mTransitCellMinimumDivisionAge = transitCellMinimumDivisionAge;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetTransitCellMinimumDivisionAge()
{
    return mTransitCellMinimumDivisionAge;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetAverageTransitCellCycleTime()
{
    return mTransitCellMinimumDivisionAge + 1.0/mTransitCellDivisionProbability;
}

double BernoulliTrialWithContactInhibitionCellCycleModel::GetAverageStemCellCycleTime()
{
    return mStemCellMinimumDivisionAge + 1.0/mStemCellDivisionProbability;
}

void BernoulliTrialWithContactInhibitionCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<QuiescentVolumeFraction>" << mQuiescentVolumeFraction << "</QuiescentVolumeFraction>\n";
    *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BernoulliTrialWithContactInhibitionCellCycleModel)

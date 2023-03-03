// EnergyPlus, Copyright (c) 1996-2023, The Board of Trustees of the University of Illinois,
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
// National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
// contributors. All rights reserved.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without the U.S. Department of Energy's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// C++ Headers
#include <cmath>

// ObjexxFCL Headers
#include <ObjexxFCL/Array.functions.hh>
#include <ObjexxFCL/Fmath.hh>

// EnergyPlus Headers
#include <EnergyPlus/Autosizing/HeatingCapacitySizing.hh>
#include <EnergyPlus/BranchNodeConnections.hh>
#include <EnergyPlus/Data/EnergyPlusData.hh>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/DataHeatBalFanSys.hh>
#include <EnergyPlus/DataHeatBalSurface.hh>
#include <EnergyPlus/DataHeatBalance.hh>
#include <EnergyPlus/DataIPShortCuts.hh>
#include <EnergyPlus/DataLoopNode.hh>
#include <EnergyPlus/DataSizing.hh>
#include <EnergyPlus/DataSurfaces.hh>
#include <EnergyPlus/DataZoneEnergyDemands.hh>
#include <EnergyPlus/DataZoneEquipment.hh>
#include <EnergyPlus/FluidProperties.hh>
#include <EnergyPlus/General.hh>
#include <EnergyPlus/GeneralRoutines.hh>
#include <EnergyPlus/GlobalNames.hh>
#include <EnergyPlus/HWBaseboardRadiator.hh>
#include <EnergyPlus/HeatBalanceIntRadExchange.hh>
#include <EnergyPlus/HeatBalanceSurfaceManager.hh>
#include <EnergyPlus/InputProcessing/InputProcessor.hh>
#include <EnergyPlus/NodeInputManager.hh>
#include <EnergyPlus/OutputProcessor.hh>
#include <EnergyPlus/Plant/DataPlant.hh>
#include <EnergyPlus/PlantUtilities.hh>
#include <EnergyPlus/Psychrometrics.hh>
#include <EnergyPlus/ScheduleManager.hh>
#include <EnergyPlus/UtilityRoutines.hh>

namespace EnergyPlus {

namespace HWBaseboardRadiator {

    // Module -- (ref: Object: ZoneHVAC:Baseboard:RadiantConvective:Water)

    // Module containing the routines dealing with the hot water baseboard heaters

    // MODULE INFORMATION:
    //       AUTHOR         Daeho Kang
    //       DATE WRITTEN   Aug 2007
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS MODULE:
    // The purpose of this module is to simulate hot water baseboard heaters.

    // METHODOLOGY EMPLOYED:

    // REFERENCES:
    // 1. I=B=R Ratings for Baseboards, Baseboard Radiation,
    //   Finned Tube (Commercial) Radiation, and Indirect Fired Water Heaters, January 2007 Edition
    // 2. Incropera and DeWitt, Fundamentals of Heat and Mass Transfer, Chapter 11.3 and 11.4,
    //   eq. 11.15, 11.17, and 11.33

    // OTHER NOTES:
    // na

    // USE STATEMENTS:
    // Use statements for data only modules
    // Using/Aliasing
    using DataHVACGlobals::SmallLoad;
    using FluidProperties::GetDensityGlycol;
    using FluidProperties::GetSpecificHeatGlycol;
    using Psychrometrics::PsyCpAirFnW;
    using Psychrometrics::PsyRhoAirFnPbTdbW;

    std::string const cCMO_BBRadiator_Water("ZoneHVAC:Baseboard:RadiantConvective:Water");
    std::string const cCMO_BBRadiator_Water_Design("ZoneHVAC:Baseboard:RadiantConvective:Water:Design");

    void SimHWBaseboard(EnergyPlusData &state,
                        std::string const &EquipName,
                        int const ControlledZoneNum,
                        bool const FirstHVACIteration,
                        Real64 &PowerMet,
                        int &CompIndex)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Russ Taylor
        //       DATE WRITTEN   Nov 1997

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine simulates the Baseboard Radiators.

        // Using/Aliasing
        using ScheduleManager::GetCurrentScheduleValue;

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int BaseboardNum; // Index of unit in baseboard array
        auto &GetInputFlag = state.dataHWBaseboardRad->GetInputFlag;
        Real64 QZnReq; // Zone load not yet satisfied
        Real64 MaxWaterFlow;
        Real64 MinWaterFlow;

        if (GetInputFlag) {
            GetHWBaseboardInput(state);
            GetInputFlag = false;
        }

        auto &HWBaseboard = state.dataHWBaseboardRad->HWBaseboard;
        auto &CheckEquipName = state.dataHWBaseboardRad->CheckEquipName;
        auto &HWBaseboardDesignObject = state.dataHWBaseboardRad->HWBaseboardDesignObject;

        // Find the correct Baseboard Equipment
        if (CompIndex == 0) {
            BaseboardNum = UtilityRoutines::FindItemInList(EquipName, HWBaseboard, &HWBaseboardParams::EquipID);
            if (BaseboardNum == 0) {
                ShowFatalError(state, format("SimHWBaseboard: Unit not found={}", EquipName));
            }
            CompIndex = BaseboardNum;
        } else {
            BaseboardNum = CompIndex;
            if (BaseboardNum > state.dataHWBaseboardRad->NumHWBaseboards || BaseboardNum < 1) {
                ShowFatalError(state,
                               format("SimHWBaseboard:  Invalid CompIndex passed={}, Number of Units={}, Entered Unit name={}",
                                      BaseboardNum,
                                      state.dataHWBaseboardRad->NumHWBaseboards,
                                      EquipName));
            }
            if (CheckEquipName(BaseboardNum)) {
                if (EquipName != HWBaseboard(BaseboardNum).EquipID) {
                    ShowFatalError(state,
                                   format("SimHWBaseboard: Invalid CompIndex passed={}, Unit name={}, stored Unit Name for that index={}",
                                          BaseboardNum,
                                          EquipName,
                                          HWBaseboard(BaseboardNum).EquipID));
                }
                CheckEquipName(BaseboardNum) = false;
            }
        }

        if (CompIndex > 0) {
            HWBaseboardDesignData HWBaseboardDesignDataObject{
                HWBaseboardDesignObject(HWBaseboard(BaseboardNum).DesignObjectPtr)}; // Contains the data for variable flow hydronic systems

            InitHWBaseboard(state, BaseboardNum, ControlledZoneNum, FirstHVACIteration);

            QZnReq = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(ControlledZoneNum).RemainingOutputReqToHeatSP;

            // On the first HVAC iteration the system values are given to the controller, but after that
            // the demand limits are in place and there needs to be feedback to the Zone Equipment
            if (FirstHVACIteration) {
                MaxWaterFlow = HWBaseboard(BaseboardNum).WaterMassFlowRateMax;
                MinWaterFlow = 0.0;
            } else {
                MaxWaterFlow = state.dataLoopNodes->Node(HWBaseboard(BaseboardNum).WaterInletNode).MassFlowRateMaxAvail;
                MinWaterFlow = state.dataLoopNodes->Node(HWBaseboard(BaseboardNum).WaterInletNode).MassFlowRateMinAvail;
            }

            switch (HWBaseboard(BaseboardNum).EquipType) {
            case DataPlant::PlantEquipmentType::Baseboard_Rad_Conv_Water: { // 'ZoneHVAC:Baseboard:RadiantConvective:Water'
                ControlCompOutput(state,
                                  HWBaseboard(BaseboardNum).EquipID,
                                  cCMO_BBRadiator_Water,
                                  BaseboardNum,
                                  FirstHVACIteration,
                                  QZnReq,
                                  HWBaseboard(BaseboardNum).WaterInletNode,
                                  MaxWaterFlow,
                                  MinWaterFlow,
                                  HWBaseboardDesignDataObject.Offset,
                                  HWBaseboard(BaseboardNum).ControlCompTypeNum,
                                  HWBaseboard(BaseboardNum).CompErrIndex,
                                  _,
                                  _,
                                  _,
                                  _,
                                  _,
                                  HWBaseboard(BaseboardNum).plantLoc);
            } break;
            default: {
                ShowSevereError(state, format("SimBaseboard: Errors in Baseboard={}", HWBaseboard(BaseboardNum).EquipID));
                ShowContinueError(state, format("Invalid or unimplemented equipment type={}", HWBaseboard(BaseboardNum).EquipType));
                ShowFatalError(state, "Preceding condition causes termination.");
            } break;
            }

            PowerMet = HWBaseboard(BaseboardNum).TotPower;

            UpdateHWBaseboard(state, BaseboardNum);

            ReportHWBaseboard(state, BaseboardNum);

        } else {
            ShowFatalError(state, format("SimHWBaseboard: Unit not found={}", EquipName));
        }
    }

    void GetHWBaseboardInput(EnergyPlusData &state)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Daeho Kang
        //       DATE WRITTEN   Aug 2007
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine gets the input for the baseboard units.

        // METHODOLOGY EMPLOYED:
        // Standard input processor calls.

        // Using/Aliasing
        using BranchNodeConnections::TestCompSet;
        using DataLoopNode::ObjectIsNotParent;
        using DataSizing::AutoSize;
        using DataSizing::CapacityPerFloorArea;
        using DataSizing::FractionOfAutosizedHeatingCapacity;
        using DataSizing::HeatingDesignCapacity;
        using GlobalNames::VerifyUniqueBaseboardName;
        using NodeInputManager::GetOnlySingleNode;
        using ScheduleManager::GetCurrentScheduleValue;
        using ScheduleManager::GetScheduleIndex;

        // SUBROUTINE PARAMETER DEFINITIONS:
        static constexpr std::string_view RoutineName("GetHWBaseboardInput:");
        Real64 constexpr MaxFraction(1.0);
        Real64 constexpr MinFraction(0.0);
        Real64 constexpr MaxWaterTempAvg(150.0);        // Maximum limit of average water temperature in degree C
        Real64 constexpr MinWaterTempAvg(20.0);         // Minimum limit of average water temperature in degree C
        Real64 constexpr HighWaterMassFlowRate(10.0);   // Maximum limit of water mass flow rate in kg/s
        Real64 constexpr LowWaterMassFlowRate(0.00001); // Minimum limit of water mass flow rate in kg/s
        Real64 constexpr MaxWaterFlowRate(10.0);        // Maximum limit of water volume flow rate in m3/s
        Real64 constexpr MinWaterFlowRate(0.00001);     // Minimum limit of water volume flow rate in m3/s
        Real64 constexpr WaterMassFlowDefault(0.063);   // Default water mass flow rate in kg/s
        //    INTEGER, PARAMETER   :: MaxDistribSurfaces    = 20         ! Maximum number of surfaces that a baseboard heater can radiate to
        int constexpr MinDistribSurfaces(1);                  // Minimum number of surfaces that a baseboard heater can radiate to
        int constexpr iHeatCAPMAlphaNum(2);                   // get input index to HW baseboard heating capacity sizing method
        int constexpr iHeatDesignCapacityNumericNum(3);       // get input index to HW baseboard heating capacity
        int constexpr iHeatCapacityPerFloorAreaNumericNum(1); // get input index to HW baseboard heating capacity per floor area sizing
        int constexpr iHeatFracOfAutosizedCapacityNumericNum(
            2); //  get input index to HW baseboard heating capacity sizing as fraction of autosized heating capacity

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 AllFracsSummed; // Sum of the fractions radiant
        int NumAlphas;  // Number of Alphas for each GetobjectItem call
        int NumNumbers; // Number of Numbers for each GetobjectItem call
        int SurfNum;    // Surface number Do loop counter
        int IOStat;
        bool ErrorsFound(false); // If errors detected in input

        auto &HWBaseboard = state.dataHWBaseboardRad->HWBaseboard;
        auto &HWBaseboardDesignObject = state.dataHWBaseboardRad->HWBaseboardDesignObject;
        auto &CheckEquipName = state.dataHWBaseboardRad->CheckEquipName;
        auto &HWBaseboardNumericFields = state.dataHWBaseboardRad->HWBaseboardNumericFields;
        auto &HWBaseboardDesignNumericFields = state.dataHWBaseboardRad->HWBaseboardDesignNumericFields;
        auto &HWBaseboardDesignNames = state.dataHWBaseboardRad->HWBaseboardDesignNames;

        state.dataHWBaseboardRad->NumHWBaseboards = state.dataInputProcessing->inputProcessor->getNumObjectsFound(state, cCMO_BBRadiator_Water);
        state.dataHWBaseboardRad->NumHWBaseboardDesignObjs = state.dataInputProcessing->inputProcessor->getNumObjectsFound(state, cCMO_BBRadiator_Water_Design);

        // Count total number of baseboard units

        HWBaseboard.allocate(state.dataHWBaseboardRad->NumHWBaseboards);
        HWBaseboardDesignObject.allocate(state.dataHWBaseboardRad->NumHWBaseboardDesignObjs);
        CheckEquipName.allocate(state.dataHWBaseboardRad->NumHWBaseboards);
        HWBaseboardNumericFields.allocate(state.dataHWBaseboardRad->NumHWBaseboards);
        HWBaseboardDesignNumericFields.allocate(state.dataHWBaseboardRad->NumHWBaseboardDesignObjs);
        HWBaseboardDesignNames.allocate(state.dataHWBaseboardRad->NumHWBaseboardDesignObjs);
        CheckEquipName = true;

        // Get the data from the user input related to design data for baseboard heaters
        for (int BaseboardDesignNum = 1; BaseboardDesignNum <= state.dataHWBaseboardRad->NumHWBaseboardDesignObjs; ++BaseboardDesignNum) {
            state.dataInputProcessing->inputProcessor->getObjectItem(state,
                                                                     cCMO_BBRadiator_Water_Design,
                                                                     BaseboardDesignNum,
                                                                     state.dataIPShortCut->cAlphaArgs,
                                                                     NumAlphas,
                                                                     state.dataIPShortCut->rNumericArgs,
                                                                     NumNumbers,
                                                                     IOStat,
                                                                     state.dataIPShortCut->lNumericFieldBlanks,
                                                                     state.dataIPShortCut->lAlphaFieldBlanks,
                                                                     state.dataIPShortCut->cAlphaFieldNames,
                                                                     state.dataIPShortCut->cNumericFieldNames);

            HWBaseboardDesignNumericFields(BaseboardDesignNum).FieldNames.allocate(NumNumbers);
            HWBaseboardDesignNumericFields(BaseboardDesignNum).FieldNames = "";
            HWBaseboardDesignNumericFields(BaseboardDesignNum).FieldNames = state.dataIPShortCut->cNumericFieldNames;

	    auto &thisHWBaseboardDesignObject = state.dataHWBaseboardRad->HWBaseboardDesignObject(BaseboardDesignNum);

	    // ErrorsFound will be set to True if problem was found, left untouched otherwise
            VerifyUniqueBaseboardName(
                state, cCMO_BBRadiator_Water_Design, state.dataIPShortCut->cAlphaArgs(1), ErrorsFound, cCMO_BBRadiator_Water_Design + " Name");

            HWBaseboardDesignObject(BaseboardDesignNum).designName = state.dataIPShortCut->cAlphaArgs(1); // Name of this baseboard design object
            HWBaseboardDesignNames(BaseboardDesignNum) = state.dataIPShortCut->cAlphaArgs(1);

            // Determine HW radiant baseboard heating design capacity sizing method
            if (UtilityRoutines::SameString(state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum), "HeatingDesignCapacity")) {
                thisHWBaseboardDesignObject.HeatingCapMethod = HeatingDesignCapacity;
            } else if (UtilityRoutines::SameString(state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum), "CapacityPerFloorArea")) {
                thisHWBaseboardDesignObject.HeatingCapMethod = CapacityPerFloorArea;
                if (!state.dataIPShortCut->lNumericFieldBlanks(iHeatCapacityPerFloorAreaNumericNum)) {
                    thisHWBaseboardDesignObject.ScaledHeatingCapacity =
                        state.dataIPShortCut->rNumericArgs(iHeatCapacityPerFloorAreaNumericNum);
                    if (thisHWBaseboardDesignObject.ScaledHeatingCapacity <= 0.0) {
                        ShowSevereError(
                            state,
                            format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboardDesignObject.designName));
                        ShowContinueError(state,
                                          format("Input for {} = {}",
                                                 state.dataIPShortCut->cAlphaFieldNames(iHeatCAPMAlphaNum),
                                                 state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum)));
                        ShowContinueError(state,
                                          format("Illegal {} = {:.7T}",
                                                 state.dataIPShortCut->cNumericFieldNames(iHeatCapacityPerFloorAreaNumericNum),
                                                 state.dataIPShortCut->rNumericArgs(iHeatCapacityPerFloorAreaNumericNum)));
                        ErrorsFound = true;
                    } else if (thisHWBaseboardDesignObject.ScaledHeatingCapacity == AutoSize) {
                        ShowSevereError(
                            state,
                            format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboardDesignObject.designName));
                        ShowContinueError(state,
                                          format("Input for {} = {}",
                                                 state.dataIPShortCut->cAlphaFieldNames(iHeatCAPMAlphaNum),
                                                 state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum)));
                        ShowContinueError(
                            state, format("Illegal {} = Autosize", state.dataIPShortCut->cNumericFieldNames(iHeatCapacityPerFloorAreaNumericNum)));
                        ErrorsFound = true;
                    }
                } else {
                    ShowSevereError(
                        state, format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboardDesignObject.designName));
                    ShowContinueError(state,
                                      format("Input for {} = {}",
                                             state.dataIPShortCut->cAlphaFieldNames(iHeatCAPMAlphaNum),
                                             state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum)));
                    ShowContinueError(
                        state,
                        format("Blank field not allowed for {}", state.dataIPShortCut->cNumericFieldNames(iHeatCapacityPerFloorAreaNumericNum)));
                    ErrorsFound = true;
                }
            } else if (UtilityRoutines::SameString(state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum), "FractionOfAutosizedHeatingCapacity")) {
                thisHWBaseboardDesignObject.HeatingCapMethod = FractionOfAutosizedHeatingCapacity;
                if (!state.dataIPShortCut->lNumericFieldBlanks(iHeatFracOfAutosizedCapacityNumericNum)) {
                    thisHWBaseboardDesignObject.ScaledHeatingCapacity =
                        state.dataIPShortCut->rNumericArgs(iHeatFracOfAutosizedCapacityNumericNum);
                    if (thisHWBaseboardDesignObject.ScaledHeatingCapacity < 0.0) {
                        ShowSevereError(
                            state,
                            format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboardDesignObject.designName));
                        ShowContinueError(state,
                                          format("Illegal {} = {:.7T}",
                                                 state.dataIPShortCut->cNumericFieldNames(iHeatFracOfAutosizedCapacityNumericNum),
                                                 state.dataIPShortCut->rNumericArgs(iHeatFracOfAutosizedCapacityNumericNum)));
                        ErrorsFound = true;
                    }
                } else {
                    ShowSevereError(
                        state, format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboardDesignObject.designName));
                    ShowContinueError(state,
                                      format("Input for {} = {}",
                                             state.dataIPShortCut->cAlphaFieldNames(iHeatCAPMAlphaNum),
                                             state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum)));
                    ShowContinueError(
                        state,
                        format("Blank field not allowed for {}", state.dataIPShortCut->cNumericFieldNames(iHeatFracOfAutosizedCapacityNumericNum)));
                    ErrorsFound = true;
                }
            }

            thisHWBaseboardDesignObject.Offset = state.dataIPShortCut->rNumericArgs(3);
            // Set default convergence tolerance
            if (thisHWBaseboardDesignObject.Offset <= 0.0) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was less than the allowable minimum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water_Design,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(3)));
                ShowContinueError(state, format("...reset to a default value=[{:.2R}].", MaxFraction));
                thisHWBaseboardDesignObject.Offset = 0.001;
            }

            thisHWBaseboardDesignObject.FracRadiant = state.dataIPShortCut->rNumericArgs(4);
            if (thisHWBaseboardDesignObject.FracRadiant < MinFraction) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was lower than the allowable minimum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(4)));
                ShowContinueError(state, format("...reset to minimum value=[{:.2R}].", MinFraction));
                thisHWBaseboardDesignObject.FracRadiant = MinFraction;
            }
            if (thisHWBaseboardDesignObject.FracRadiant > MaxFraction) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was higher than the allowable maximum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(4)));
                ShowContinueError(state, format("...reset to maximum value=[{:.2R}].", MaxFraction));
                thisHWBaseboardDesignObject.FracRadiant = MaxFraction;
            }

            thisHWBaseboardDesignObject.FracDistribPerson = state.dataIPShortCut->rNumericArgs(5);
            if (thisHWBaseboardDesignObject.FracDistribPerson < MinFraction) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was lower than the allowable minimum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(5)));
                ShowContinueError(state, format("...reset to minimum value=[{:.3R}].", MinFraction));
                thisHWBaseboardDesignObject.FracDistribPerson = MinFraction;
            }
            if (thisHWBaseboardDesignObject.FracDistribPerson > MaxFraction) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was higher than the allowable maximum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(5)));
                ShowContinueError(state, format("...reset to maximum value=[{:.3R}].", MaxFraction));
                thisHWBaseboardDesignObject.FracDistribPerson = MaxFraction;
            }
        }

        // Get the data from the user input related to baseboard heaters
        for (int BaseboardNum = 1; BaseboardNum <= state.dataHWBaseboardRad->NumHWBaseboards; ++BaseboardNum) {

            state.dataInputProcessing->inputProcessor->getObjectItem(state,
                                                                     cCMO_BBRadiator_Water,
                                                                     BaseboardNum,
                                                                     state.dataIPShortCut->cAlphaArgs,
                                                                     NumAlphas,
                                                                     state.dataIPShortCut->rNumericArgs,
                                                                     NumNumbers,
                                                                     IOStat,
                                                                     state.dataIPShortCut->lNumericFieldBlanks,
                                                                     state.dataIPShortCut->lAlphaFieldBlanks,
                                                                     state.dataIPShortCut->cAlphaFieldNames,
                                                                     state.dataIPShortCut->cNumericFieldNames);

            HWBaseboardNumericFields(BaseboardNum).FieldNames.allocate(NumNumbers);
            HWBaseboardNumericFields(BaseboardNum).FieldNames = "";
            HWBaseboardNumericFields(BaseboardNum).FieldNames = state.dataIPShortCut->cNumericFieldNames;

            // ErrorsFound will be set to True if problem was found, left untouched otherwise
            VerifyUniqueBaseboardName(
                state, cCMO_BBRadiator_Water, state.dataIPShortCut->cAlphaArgs(1), ErrorsFound, cCMO_BBRadiator_Water + " Name");

	    auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
	    
            thisHWBaseboard.EquipID = state.dataIPShortCut->cAlphaArgs(1); // Name of this baseboard
            thisHWBaseboard.EquipType =
                DataPlant::PlantEquipmentType::Baseboard_Rad_Conv_Water; //'ZoneHVAC:Baseboard:RadiantConvective:Water'

            thisHWBaseboard.designObjectName = state.dataIPShortCut->cAlphaArgs(2); // Name of the design object for this baseboard
            thisHWBaseboard.DesignObjectPtr =
                UtilityRoutines::FindItemInList(thisHWBaseboard.designObjectName, HWBaseboardDesignNames);
            HWBaseboardDesignData HWBaseboardDesignDataObject{
                HWBaseboardDesignObject(thisHWBaseboard.DesignObjectPtr)}; // Contains the data for the design object

            // Get schedule
            thisHWBaseboard.Schedule = state.dataIPShortCut->cAlphaArgs(3);
            if (state.dataIPShortCut->lAlphaFieldBlanks(3)) {
                thisHWBaseboard.SchedPtr = DataGlobalConstants::ScheduleAlwaysOn;
            } else {
                thisHWBaseboard.SchedPtr = GetScheduleIndex(state, state.dataIPShortCut->cAlphaArgs(3));
                if (thisHWBaseboard.SchedPtr == 0) {
                    ShowSevereError(state,
                                    format("{}{}=\"{}\", {}=\"{}\" not found.",
                                           RoutineName,
                                           cCMO_BBRadiator_Water,
                                           state.dataIPShortCut->cAlphaArgs(1),
                                           state.dataIPShortCut->cAlphaFieldNames(3),
                                           state.dataIPShortCut->cAlphaArgs(3)));
                    ErrorsFound = true;
                }
            }

            // Get inlet node number
            thisHWBaseboard.WaterInletNode = GetOnlySingleNode(state,
                                                                         state.dataIPShortCut->cAlphaArgs(4),
                                                                         ErrorsFound,
                                                                         DataLoopNode::ConnectionObjectType::ZoneHVACBaseboardRadiantConvectiveWater,
                                                                         state.dataIPShortCut->cAlphaArgs(1),
                                                                         DataLoopNode::NodeFluidType::Water,
                                                                         DataLoopNode::ConnectionType::Inlet,
                                                                         NodeInputManager::CompFluidStream::Primary,
                                                                         ObjectIsNotParent);

            // Get outlet node number
            thisHWBaseboard.WaterOutletNode = GetOnlySingleNode(state,
                                                                          state.dataIPShortCut->cAlphaArgs(5),
                                                                          ErrorsFound,
                                                                          DataLoopNode::ConnectionObjectType::ZoneHVACBaseboardRadiantConvectiveWater,
                                                                          state.dataIPShortCut->cAlphaArgs(1),
                                                                          DataLoopNode::NodeFluidType::Water,
                                                                          DataLoopNode::ConnectionType::Outlet,
                                                                          NodeInputManager::CompFluidStream::Primary,
                                                                          ObjectIsNotParent);
            TestCompSet(state,
                        cCMO_BBRadiator_Water,
                        state.dataIPShortCut->cAlphaArgs(1),
                        state.dataIPShortCut->cAlphaArgs(4),
                        state.dataIPShortCut->cAlphaArgs(5),
                        "Hot Water Nodes");

            thisHWBaseboard.WaterTempAvg = state.dataIPShortCut->rNumericArgs(1);
            if (thisHWBaseboard.WaterTempAvg > MaxWaterTempAvg + 0.001) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was higher than the allowable maximum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(1)));
                ShowContinueError(state, format("...reset to maximum value=[{:.2R}].", MaxWaterTempAvg));
                thisHWBaseboard.WaterTempAvg = MaxWaterTempAvg;
            } else if (thisHWBaseboard.WaterTempAvg < MinWaterTempAvg - 0.001) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was lower than the allowable minimum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(1)));
                ShowContinueError(state, format("...reset to minimum value=[{:.2R}].", MinWaterTempAvg));
                thisHWBaseboard.WaterTempAvg = MinWaterTempAvg;
            }

            thisHWBaseboard.WaterMassFlowRateStd = state.dataIPShortCut->rNumericArgs(2);
            if (thisHWBaseboard.WaterMassFlowRateStd < LowWaterMassFlowRate - 0.0001 ||
                thisHWBaseboard.WaterMassFlowRateStd > HighWaterMassFlowRate + 0.0001) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} is an invalid Standard Water mass flow rate.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(2)));
                ShowContinueError(state, format("...reset to a default value=[{:.1R}].", WaterMassFlowDefault));
                thisHWBaseboard.WaterMassFlowRateStd = WaterMassFlowDefault;
            }

            // Determine HW radiant baseboard heating design capacity sizing method
            thisHWBaseboard.HeatingCapMethod = HWBaseboardDesignDataObject.HeatingCapMethod;
            if (thisHWBaseboard.HeatingCapMethod == HeatingDesignCapacity) {
                if (!state.dataIPShortCut->lNumericFieldBlanks(iHeatDesignCapacityNumericNum)) {
                    thisHWBaseboard.ScaledHeatingCapacity = state.dataIPShortCut->rNumericArgs(iHeatDesignCapacityNumericNum);
                    if (thisHWBaseboard.ScaledHeatingCapacity < 0.0 && thisHWBaseboard.ScaledHeatingCapacity != AutoSize) {
                        ShowSevereError(state, format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboard.EquipID));
                        ShowContinueError(state,
                                          format("Illegal {} = {:.7T}",
                                                 state.dataIPShortCut->cNumericFieldNames(iHeatDesignCapacityNumericNum),
                                                 state.dataIPShortCut->rNumericArgs(iHeatDesignCapacityNumericNum)));
                        ErrorsFound = true;
                    }
                } else {
                    ShowSevereError(state, format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboard.EquipID));
                    ShowContinueError(state,
                                      format("Input for {} = {}",
                                             state.dataIPShortCut->cAlphaFieldNames(iHeatCAPMAlphaNum),
                                             state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum)));
                    ShowContinueError(
                        state, format("Blank field not allowed for {}", state.dataIPShortCut->cNumericFieldNames(iHeatDesignCapacityNumericNum)));
                    ErrorsFound = true;
                }
            } else if (thisHWBaseboard.HeatingCapMethod == CapacityPerFloorArea) {
                thisHWBaseboard.ScaledHeatingCapacity = HWBaseboardDesignDataObject.ScaledHeatingCapacity;

            } else if (thisHWBaseboard.HeatingCapMethod == FractionOfAutosizedHeatingCapacity) {
                thisHWBaseboard.ScaledHeatingCapacity = HWBaseboardDesignDataObject.ScaledHeatingCapacity;

            } else {
                ShowSevereError(state, format("{} = {}", state.dataIPShortCut->cCurrentModuleObject, thisHWBaseboard.EquipID));
                ShowContinueError(state,
                                  format("Illegal {} = {}",
                                         state.dataIPShortCut->cAlphaFieldNames(iHeatCAPMAlphaNum),
                                         state.dataIPShortCut->cAlphaArgs(iHeatCAPMAlphaNum)));
                ErrorsFound = true;
            }

            thisHWBaseboard.WaterVolFlowRateMax = state.dataIPShortCut->rNumericArgs(4);
            if (std::abs(thisHWBaseboard.WaterVolFlowRateMax) <= MinWaterFlowRate) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was less than the allowable minimum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(4)));
                ShowContinueError(state, format("...reset to minimum value=[{:.2R}].", MinWaterFlowRate));
                thisHWBaseboard.WaterVolFlowRateMax = MinWaterFlowRate;
            } else if (thisHWBaseboard.WaterVolFlowRateMax > MaxWaterFlowRate) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", {} was higher than the allowable maximum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1),
                                        state.dataIPShortCut->cNumericFieldNames(4)));
                ShowContinueError(state, format("...reset to maximum value=[{:.2R}].", MaxWaterFlowRate));
                thisHWBaseboard.WaterVolFlowRateMax = MaxWaterFlowRate;
            }

            // Remaining fraction is added to the zone as convective heat transfer
            AllFracsSummed = HWBaseboardDesignDataObject.FracDistribPerson;
            if (AllFracsSummed > MaxFraction) {
                ShowWarningError(state,
                                 format("{}{}=\"{}\", Fraction Radiant was higher than the allowable maximum.",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1)));
                HWBaseboardDesignDataObject.FracRadiant = MaxFraction;
                thisHWBaseboard.FracConvect = 0.0;
            } else {
                thisHWBaseboard.FracConvect = 1.0 - AllFracsSummed;
            }

            thisHWBaseboard.TotSurfToDistrib = NumNumbers - 4;
            //      IF (thisHWBaseboard%TotSurfToDistrib > MaxDistribSurfaces) THEN
            //        CALL ShowWarningError(state, RoutineName//cCMO_BBRadiator_Water//'="'//TRIM(state.dataIPShortCut->cAlphaArgs(1))// &
            //          '", the number of surface/radiant fraction groups entered was higher than the allowable maximum.')
            //        CALL ShowContinueError(state, '...only the maximum value=['//TRIM(RoundSigDigits(MaxDistribSurfaces))// &
            //           '] will be processed.')
            //        thisHWBaseboard%TotSurfToDistrib = MaxDistribSurfaces
            //      END IF
            if ((thisHWBaseboard.TotSurfToDistrib < MinDistribSurfaces) && (HWBaseboardDesignDataObject.FracRadiant > MinFraction)) {
                ShowSevereError(state,
                                std::string{RoutineName} + cCMO_BBRadiator_Water + "=\"" + state.dataIPShortCut->cAlphaArgs(1) +
                                    "\", the number of surface/radiant fraction groups entered was less than the allowable minimum.");
                ShowContinueError(state, format("...the minimum that must be entered=[{}].", MinDistribSurfaces));
                ErrorsFound = true;
                thisHWBaseboard.TotSurfToDistrib = 0; // error
            }

            thisHWBaseboard.SurfaceName.allocate(thisHWBaseboard.TotSurfToDistrib);
            thisHWBaseboard.SurfaceName = "";
            thisHWBaseboard.SurfacePtr.allocate(thisHWBaseboard.TotSurfToDistrib);
            thisHWBaseboard.SurfacePtr = 0;
            thisHWBaseboard.FracDistribToSurf.allocate(thisHWBaseboard.TotSurfToDistrib);
            thisHWBaseboard.FracDistribToSurf = 0.0;

            thisHWBaseboard.ZonePtr =
                DataZoneEquipment::GetZoneEquipControlledZoneNum(state, DataZoneEquipment::ZoneEquip::BBWater, thisHWBaseboard.EquipID);

            AllFracsSummed = HWBaseboardDesignDataObject.FracDistribPerson;
            for (SurfNum = 1; SurfNum <= thisHWBaseboard.TotSurfToDistrib; ++SurfNum) {
                thisHWBaseboard.SurfaceName(SurfNum) = state.dataIPShortCut->cAlphaArgs(SurfNum + 5);
                thisHWBaseboard.SurfacePtr(SurfNum) =
                    HeatBalanceIntRadExchange::GetRadiantSystemSurface(state,
                                                                       cCMO_BBRadiator_Water,
                                                                       thisHWBaseboard.EquipID,
                                                                       thisHWBaseboard.ZonePtr,
                                                                       thisHWBaseboard.SurfaceName(SurfNum),
                                                                       ErrorsFound);
                thisHWBaseboard.FracDistribToSurf(SurfNum) = state.dataIPShortCut->rNumericArgs(SurfNum + 4);
                if (thisHWBaseboard.FracDistribToSurf(SurfNum) > MaxFraction) {
                    ShowWarningError(state,
                                     format("{}{}=\"{}\", {}was greater than the allowable maximum.",
                                            RoutineName,
                                            cCMO_BBRadiator_Water,
                                            state.dataIPShortCut->cAlphaArgs(1),
                                            state.dataIPShortCut->cNumericFieldNames(SurfNum + 4)));
                    ShowContinueError(state, format("...reset to maximum value=[{:.2R}].", MaxFraction));
                    thisHWBaseboard.TotSurfToDistrib = MaxFraction;
                }
                if (thisHWBaseboard.FracDistribToSurf(SurfNum) < MinFraction) {
                    ShowWarningError(state,
                                     format("{}{}=\"{}\", {}was less than the allowable minimum.",
                                            RoutineName,
                                            cCMO_BBRadiator_Water,
                                            state.dataIPShortCut->cAlphaArgs(1),
                                            state.dataIPShortCut->cNumericFieldNames(SurfNum + 4)));
                    ShowContinueError(state, format("...reset to maximum value=[{:.2R}].", MinFraction));
                    thisHWBaseboard.TotSurfToDistrib = MinFraction;
                }
                if (thisHWBaseboard.SurfacePtr(SurfNum) != 0) {
                    state.dataSurface->SurfIntConvSurfGetsRadiantHeat(thisHWBaseboard.SurfacePtr(SurfNum)) = true;
                }

                AllFracsSummed += thisHWBaseboard.FracDistribToSurf(SurfNum);
            } // Surfaces

            if (AllFracsSummed > (MaxFraction + 0.01)) {
                ShowSevereError(state,
                                format("{}{}=\"{}\", Summed radiant fractions for people + surface groups > 1.0",
                                       RoutineName,
                                       cCMO_BBRadiator_Water,
                                       state.dataIPShortCut->cAlphaArgs(1)));
                ErrorsFound = true;
            }
            if ((AllFracsSummed < (MaxFraction - 0.01)) &&
                (HWBaseboardDesignDataObject.FracRadiant >
                 MinFraction)) { // User didn't distribute all of the | radiation warn that some will be lost
                ShowWarningError(state,
                                 format("{}{}=\"{}\", Summed radiant fractions for people + surface groups < 1.0",
                                        RoutineName,
                                        cCMO_BBRadiator_Water,
                                        state.dataIPShortCut->cAlphaArgs(1)));
                ShowContinueError(state, "The rest of the radiant energy delivered by the baseboard heater will be lost");
            }
        }

        if (ErrorsFound) {
            ShowFatalError(state, format("{}{}Errors found getting input. Program terminates.", RoutineName, cCMO_BBRadiator_Water));
        }

        // Setup Report variables for the Coils
        for (int BaseboardNum = 1; BaseboardNum <= state.dataHWBaseboardRad->NumHWBaseboards; ++BaseboardNum) {
            // CurrentModuleObject='ZoneHVAC:Baseboard:RadiantConvective:Water'
            auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
            SetupOutputVariable(state,
                                "Baseboard Total Heating Rate",
                                OutputProcessor::Unit::W,
                                thisHWBaseboard.TotPower,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);

            SetupOutputVariable(state,
                                "Baseboard Convective Heating Rate",
                                OutputProcessor::Unit::W,
                                thisHWBaseboard.ConvPower,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Radiant Heating Rate",
                                OutputProcessor::Unit::W,
                                thisHWBaseboard.RadPower,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Total Heating Energy",
                                OutputProcessor::Unit::J,
                                thisHWBaseboard.TotEnergy,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Summed,
                                thisHWBaseboard.EquipID,
                                _,
                                "ENERGYTRANSFER",
                                "BASEBOARD",
                                _,
                                "System");

            SetupOutputVariable(state,
                                "Baseboard Convective Heating Energy",
                                OutputProcessor::Unit::J,
                                thisHWBaseboard.ConvEnergy,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Summed,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Radiant Heating Energy",
                                OutputProcessor::Unit::J,
                                thisHWBaseboard.RadEnergy,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Summed,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Hot Water Energy",
                                OutputProcessor::Unit::J,
                                thisHWBaseboard.Energy,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Summed,
                                thisHWBaseboard.EquipID,
                                _,
                                "PLANTLOOPHEATINGDEMAND",
                                "BASEBOARD",
                                _,
                                "System");
            SetupOutputVariable(state,
                                "Baseboard Hot Water Mass Flow Rate",
                                OutputProcessor::Unit::kg_s,
                                thisHWBaseboard.WaterMassFlowRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Air Mass Flow Rate",
                                OutputProcessor::Unit::kg_s,
                                thisHWBaseboard.AirMassFlowRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Air Inlet Temperature",
                                OutputProcessor::Unit::C,
                                thisHWBaseboard.AirInletTemp,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Air Outlet Temperature",
                                OutputProcessor::Unit::C,
                                thisHWBaseboard.AirOutletTemp,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Water Inlet Temperature",
                                OutputProcessor::Unit::C,
                                thisHWBaseboard.WaterInletTemp,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
            SetupOutputVariable(state,
                                "Baseboard Water Outlet Temperature",
                                OutputProcessor::Unit::C,
                                thisHWBaseboard.WaterOutletTemp,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisHWBaseboard.EquipID);
        }
    }

    void InitHWBaseboard(EnergyPlusData &state, int const BaseboardNum, int const ControlledZoneNum, bool const FirstHVACIteration)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Russ Taylor
        //                      Rick Strand
        //       DATE WRITTEN   Nov 1997
        //                      Feb 2001
        //       MODIFIED       Aug 2007 Daeho Kang (Add radiant component)
        //                      Sept 2010 Brent Griffith (plant interactions)

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine initializes the baseboard units, and determines the UA values during simulation.

        // METHODOLOGY EMPLOYED:
        // The initialization subrotines both in high temperature radiant radiator
        // and convective only baseboard radiator are combined and modified. In addition,
        // an UA value calculation by LMTD method is added.
        // The heater is assumed to be crossflow with both fluids unmixed.

        // REFERENCES:
        // 1. Incropera and DeWitt, Fundamentals of Heat and Mass Transfer
        // Chapter 11.3, p. 510, eq. 11.15 and 11.17
        // 2. I=B=R Ratings for Baseboards, Baseboard Radiation, Finned Tube (Commercial) Radiation,
        // and Indirect Fired Water Heaters, January 2007 Edition

        // Using/Aliasing
        using PlantUtilities::InitComponentNodes;
        using PlantUtilities::ScanPlantLoopsForObject;

        Real64 constexpr Constant(0.0062); // Constant of linear equation for air mass flow rate
        Real64 constexpr Coeff(0.0000275); // Correlation coefficient to capacity
        static constexpr std::string_view RoutineName("BaseboardRadiatorWater:InitHWBaseboard");

        auto &MyOneTimeFlag = state.dataHWBaseboardRad->MyOneTimeFlag;
        auto &MyEnvrnFlag = state.dataHWBaseboardRad->MyEnvrnFlag;

        int WaterInletNode;
        Real64 RhoAirStdInit;
        Real64 rho; // local fluid density
        Real64 Cp;  // local fluid specific heat
        bool errFlag;

        auto &MySizeFlag = state.dataHWBaseboardRad->MySizeFlag;
        auto &ZeroSourceSumHATsurf = state.dataHWBaseboardRad->ZeroSourceSumHATsurf;
        auto &QBBRadSource = state.dataHWBaseboardRad->QBBRadSource;
        auto &QBBRadSrcAvg = state.dataHWBaseboardRad->QBBRadSrcAvg;
        auto &LastQBBRadSrc = state.dataHWBaseboardRad->LastQBBRadSrc;
        auto &LastSysTimeElapsed = state.dataHWBaseboardRad->LastSysTimeElapsed;
        auto &LastTimeStepSys = state.dataHWBaseboardRad->LastTimeStepSys;
        auto &SetLoopIndexFlag = state.dataHWBaseboardRad->SetLoopIndexFlag;
        auto &HWBaseboard = state.dataHWBaseboardRad->HWBaseboard;

        // Do the one time initializations
        if (MyOneTimeFlag) {

            // Initialize the environment and sizing flags
            MyEnvrnFlag.allocate(state.dataHWBaseboardRad->NumHWBaseboards);
            MySizeFlag.allocate(state.dataHWBaseboardRad->NumHWBaseboards);
            ZeroSourceSumHATsurf.dimension(state.dataGlobal->NumOfZones, 0.0);
            QBBRadSource.dimension(state.dataHWBaseboardRad->NumHWBaseboards, 0.0);
            QBBRadSrcAvg.dimension(state.dataHWBaseboardRad->NumHWBaseboards, 0.0);
            LastQBBRadSrc.dimension(state.dataHWBaseboardRad->NumHWBaseboards, 0.0);
            LastSysTimeElapsed.dimension(state.dataHWBaseboardRad->NumHWBaseboards, 0.0);
            LastTimeStepSys.dimension(state.dataHWBaseboardRad->NumHWBaseboards, 0.0);
            SetLoopIndexFlag.allocate(state.dataHWBaseboardRad->NumHWBaseboards);
            MyEnvrnFlag = true;
            MySizeFlag = true;
            MyOneTimeFlag = false;
            SetLoopIndexFlag = true;
            for (int Loop = 1; Loop <= state.dataHWBaseboardRad->NumHWBaseboards; ++Loop) {
                // Air mass flow rate is obtained from the following linear equation (reset if autosize is used)
                // m_dot = 0.0062 + 2.75e-05*q
                HWBaseboard(Loop).AirMassFlowRateStd = Constant + Coeff * HWBaseboard(Loop).RatedCapacity;
            }
        }

        auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);

        if (SetLoopIndexFlag(BaseboardNum)) {
            if (allocated(state.dataPlnt->PlantLoop)) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisHWBaseboard.EquipID,
                                        thisHWBaseboard.EquipType,
                                        thisHWBaseboard.plantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
                if (errFlag) {
                    ShowFatalError(state, "InitHWBaseboard: Program terminated for previous conditions.");
                }
                SetLoopIndexFlag(BaseboardNum) = false;
            }
        }

        if (!state.dataGlobal->SysSizingCalc && MySizeFlag(BaseboardNum) && !SetLoopIndexFlag(BaseboardNum)) {
            // For each coil, do the sizing once
            SizeHWBaseboard(state, BaseboardNum);
            MySizeFlag(BaseboardNum) = false;
        }

        // Do the Begin Environment initializations
        if (state.dataGlobal->BeginEnvrnFlag && MyEnvrnFlag(BaseboardNum)) {
            // Initialize
            RhoAirStdInit = state.dataEnvrn->StdRhoAir;
            WaterInletNode = thisHWBaseboard.WaterInletNode;

            rho = GetDensityGlycol(state,
                                   state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                   DataGlobalConstants::HWInitConvTemp,
                                   state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                   RoutineName);

            thisHWBaseboard.WaterMassFlowRateMax = rho * thisHWBaseboard.WaterVolFlowRateMax;

            InitComponentNodes(state,
                               0.0,
                               thisHWBaseboard.WaterMassFlowRateMax,
                               thisHWBaseboard.WaterInletNode,
                               thisHWBaseboard.WaterOutletNode);

            state.dataLoopNodes->Node(WaterInletNode).Temp = 60.0;

            Cp = GetSpecificHeatGlycol(state,
                                       state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                       state.dataLoopNodes->Node(WaterInletNode).Temp,
                                       state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                       RoutineName);

            state.dataLoopNodes->Node(WaterInletNode).Enthalpy = Cp * state.dataLoopNodes->Node(WaterInletNode).Temp;
            state.dataLoopNodes->Node(WaterInletNode).Quality = 0.0;
            state.dataLoopNodes->Node(WaterInletNode).Press = 0.0;
            state.dataLoopNodes->Node(WaterInletNode).HumRat = 0.0;

            ZeroSourceSumHATsurf = 0.0;
            QBBRadSource = 0.0;
            QBBRadSrcAvg = 0.0;
            LastQBBRadSrc = 0.0;
            LastSysTimeElapsed = 0.0;
            LastTimeStepSys = 0.0;

            MyEnvrnFlag(BaseboardNum) = false;
        }

        if (!state.dataGlobal->BeginEnvrnFlag) {
            MyEnvrnFlag(BaseboardNum) = true;
        }

        if (state.dataGlobal->BeginTimeStepFlag && FirstHVACIteration) {
            int ZoneNum = thisHWBaseboard.ZonePtr;
            ZeroSourceSumHATsurf(ZoneNum) = state.dataHeatBal->Zone(ZoneNum).sumHATsurf(state);
            QBBRadSrcAvg(BaseboardNum) = 0.0;
            LastQBBRadSrc(BaseboardNum) = 0.0;
            LastSysTimeElapsed(BaseboardNum) = 0.0;
            LastTimeStepSys(BaseboardNum) = 0.0;
        }

        // Do the every time step initializations
        WaterInletNode = thisHWBaseboard.WaterInletNode;
        int ZoneNode = state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).ZoneNode;
        thisHWBaseboard.WaterMassFlowRate = state.dataLoopNodes->Node(WaterInletNode).MassFlowRate;
        thisHWBaseboard.WaterInletTemp = state.dataLoopNodes->Node(WaterInletNode).Temp;
        thisHWBaseboard.WaterInletEnthalpy = state.dataLoopNodes->Node(WaterInletNode).Enthalpy;
        thisHWBaseboard.AirInletTemp = state.dataLoopNodes->Node(ZoneNode).Temp;
        thisHWBaseboard.AirInletHumRat = state.dataLoopNodes->Node(ZoneNode).HumRat;

        thisHWBaseboard.TotPower = 0.0;
        thisHWBaseboard.Power = 0.0;
        thisHWBaseboard.ConvPower = 0.0;
        thisHWBaseboard.RadPower = 0.0;
        thisHWBaseboard.TotEnergy = 0.0;
        thisHWBaseboard.Energy = 0.0;
        thisHWBaseboard.ConvEnergy = 0.0;
        thisHWBaseboard.RadEnergy = 0.0;
    }

    void SizeHWBaseboard(EnergyPlusData &state, int const BaseboardNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Fred Buhl
        //       DATE WRITTEN   February 2002
        //       MODIFIED       August 2009 Daeho Kang (Add UA autosizing by LMTD)
        //                      Aug 2013 Daeho Kang, add component sizing table entries
        //                      July 2014, B.Nigusse, added scalable sizing
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for sizing hot water baseboard components

        // METHODOLOGY EMPLOYED:
        // Obtains flow rates from the zone sizing arrays and plant sizing data.

        // REFERENCES:
        // na

        // Using/Aliasing
        using namespace DataSizing;
        using DataHVACGlobals::HeatingCapacitySizing;

        using PlantUtilities::RegisterPlantCompDesignFlow;

        // Locals
        // SUBROUTINE ARGUMENT DEFINITIONS:

        // SUBROUTINE PARAMETER DEFINITIONS:
        Real64 constexpr AirInletTempStd(18.0); // I=B=R rating document
        Real64 constexpr CPAirStd(1005.0);      // Average specific heat of air at between 25C and 40C in J/kg-k
        Real64 constexpr Constant(0.0062);      // Constant of linear equation for air mass flow rate
        Real64 constexpr Coeff(0.0000275);      // Correlation coefficient to capacity
        static constexpr std::string_view RoutineName("SizeHWBaseboard");
        static constexpr std::string_view RoutineNameFull("BaseboardRadiatorWater:SizeHWBaseboard");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int PltSizNum;     // do loop index for plant sizing
        int PltSizHeatNum; // index of plant sizing object for 1st heating loop
        Real64 DesCoilLoad;
        Real64 WaterInletTempStd;
        Real64 WaterOutletTempStd;
        Real64 AirOutletTempStd;
        Real64 DeltaT1;
        Real64 DeltaT2;
        Real64 LMTD;
        Real64 AirMassFlowRate;
        Real64 WaterMassFlowRateStd;
        Real64 rho;                     // local fluid density
        Real64 Cp;                      // local fluid specific heat
        bool ErrorsFound;               // If errors detected in input
        bool FlowAutoSize;              // Indicator to autosize for maximum water vloume flow
        bool CapAutoSize;               // Indicator to autosize for capacity
        Real64 WaterVolFlowRateMaxDes;  // Design maximum water volume flow for reproting
        Real64 WaterVolFlowRateMaxUser; // User hard-sized maximum water volume flow for reproting
        Real64 RatedCapacityDes;        // Design rated capacity for reproting

        std::string CompName;     // component name
        std::string CompType;     // component type
        std::string SizingString; // input field sizing description (e.g., Nominal Capacity)
        Real64 TempSize;          // autosized value of coil input field
        int FieldNum = 1;         // IDD numeric field number where input field description is found
        int SizingMethod; // Integer representation of sizing method name (e.g., CoolingAirflowSizing, HeatingAirflowSizing, CoolingCapacitySizing,
                          // HeatingCapacitySizing, etc.)
        bool PrintFlag;   // TRUE when sizing information is reported in the eio file
        int CapSizingMethod(0); // capacity sizing methods (HeatingDesignCapacity, CapacityPerFloorArea, FractionOfAutosizedCoolingCapacity, and
                                // FractionOfAutosizedHeatingCapacity )

        PltSizHeatNum = 0;
        PltSizNum = 0;
        DesCoilLoad = 0.0;
        ErrorsFound = false;
        FlowAutoSize = false;
        CapAutoSize = false;
        WaterVolFlowRateMaxDes = 0.0;
        WaterVolFlowRateMaxUser = 0.0;
        RatedCapacityDes = 0.0;
        state.dataSize->DataScalableCapSizingON = false;

        auto &ZoneEqSizing(state.dataSize->ZoneEqSizing);
        auto &HWBaseboardNumericFields = state.dataHWBaseboardRad->HWBaseboardNumericFields;

        auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
	
        if (state.dataSize->CurZoneEqNum > 0) {

            CompType = cCMO_BBRadiator_Water;
            CompName = thisHWBaseboard.EquipID;
            state.dataSize->DataHeatSizeRatio = 1.0;
            state.dataSize->DataFracOfAutosizedHeatingCapacity = 1.0;
            state.dataSize->DataZoneNumber = thisHWBaseboard.ZonePtr;
            SizingMethod = HeatingCapacitySizing;
            FieldNum = 3;
            PrintFlag = false;
            SizingString = HWBaseboardNumericFields(BaseboardNum).FieldNames(FieldNum) + " [W]";
            CapSizingMethod = thisHWBaseboard.HeatingCapMethod;
            ZoneEqSizing(state.dataSize->CurZoneEqNum).SizingMethod(SizingMethod) = CapSizingMethod;
            if (CapSizingMethod == HeatingDesignCapacity || CapSizingMethod == CapacityPerFloorArea ||
                CapSizingMethod == FractionOfAutosizedHeatingCapacity) {
                if (CapSizingMethod == HeatingDesignCapacity) {
                    if (thisHWBaseboard.ScaledHeatingCapacity == AutoSize) {
                        CheckZoneSizing(state, CompType, CompName);
                        ZoneEqSizing(state.dataSize->CurZoneEqNum).HeatingCapacity = true;
                        ZoneEqSizing(state.dataSize->CurZoneEqNum).DesHeatingLoad =
                            state.dataSize->FinalZoneSizing(state.dataSize->CurZoneEqNum).NonAirSysDesHeatLoad;
                    }
                    TempSize = thisHWBaseboard.ScaledHeatingCapacity;

                } else if (CapSizingMethod == CapacityPerFloorArea) {
                    ZoneEqSizing(state.dataSize->CurZoneEqNum).HeatingCapacity = true;
                    ZoneEqSizing(state.dataSize->CurZoneEqNum).DesHeatingLoad =
                        thisHWBaseboard.ScaledHeatingCapacity * state.dataHeatBal->Zone(state.dataSize->DataZoneNumber).FloorArea;
                    TempSize = ZoneEqSizing(state.dataSize->CurZoneEqNum).DesHeatingLoad;
                    state.dataSize->DataScalableCapSizingON = true;
                } else if (CapSizingMethod == FractionOfAutosizedHeatingCapacity) {
                    CheckZoneSizing(state, CompType, CompName);
                    ZoneEqSizing(state.dataSize->CurZoneEqNum).HeatingCapacity = true;
                    state.dataSize->DataFracOfAutosizedHeatingCapacity = thisHWBaseboard.ScaledHeatingCapacity;
                    ZoneEqSizing(state.dataSize->CurZoneEqNum).DesHeatingLoad =
                        state.dataSize->FinalZoneSizing(state.dataSize->CurZoneEqNum).NonAirSysDesHeatLoad;
                    TempSize = AutoSize;
                    state.dataSize->DataScalableCapSizingON = true;
                } else {
                    TempSize = thisHWBaseboard.ScaledHeatingCapacity;
                }
                bool errorsFound = false;
                HeatingCapacitySizer sizerHeatingCapacity;
                sizerHeatingCapacity.overrideSizingString(SizingString);
                sizerHeatingCapacity.initializeWithinEP(state, CompType, CompName, PrintFlag, RoutineName);
                TempSize = sizerHeatingCapacity.size(state, TempSize, errorsFound);
                if (thisHWBaseboard.ScaledHeatingCapacity == AutoSize) {
                    thisHWBaseboard.RatedCapacity = AutoSize;
                } else {
                    thisHWBaseboard.RatedCapacity = TempSize;
                }
                RatedCapacityDes = TempSize;
                state.dataSize->DataScalableCapSizingON = false;
            }
        }

        // find the appropriate heating Plant Sizing object
        PltSizHeatNum = state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).PlantSizNum;

        if (PltSizHeatNum > 0) {
            if (state.dataSize->CurZoneEqNum > 0) {

                if (thisHWBaseboard.WaterVolFlowRateMax == AutoSize) {
                    FlowAutoSize = true;
                }
                if (!FlowAutoSize && !state.dataSize->ZoneSizingRunDone) { // Simulation continue
                    if (thisHWBaseboard.WaterVolFlowRateMax > 0.0) {
                        BaseSizer::reportSizerOutput(state,
                                                     cCMO_BBRadiator_Water,
                                                     thisHWBaseboard.EquipID,
                                                     "User-Specified Maximum Water Flow Rate [m3/s]",
                                                     thisHWBaseboard.WaterVolFlowRateMax);
                    }
                } else {
                    CheckZoneSizing(state, cCMO_BBRadiator_Water, thisHWBaseboard.EquipID);
                    DesCoilLoad = RatedCapacityDes;
                    if (DesCoilLoad >= SmallLoad) {
                        Cp = GetSpecificHeatGlycol(state,
                                                   state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                                   DataGlobalConstants::HWInitConvTemp,
                                                   state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                                   RoutineName);
                        rho = GetDensityGlycol(state,
                                               state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                               DataGlobalConstants::HWInitConvTemp,
                                               state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                               RoutineName);
                        WaterVolFlowRateMaxDes = DesCoilLoad / (state.dataSize->PlantSizData(PltSizHeatNum).DeltaT * Cp * rho);
                    } else {
                        WaterVolFlowRateMaxDes = 0.0;
                    }

                    if (FlowAutoSize) {
                        thisHWBaseboard.WaterVolFlowRateMax = WaterVolFlowRateMaxDes;
                        BaseSizer::reportSizerOutput(state,
                                                     cCMO_BBRadiator_Water,
                                                     thisHWBaseboard.EquipID,
                                                     "Design Size Maximum Water Flow Rate [m3/s]",
                                                     WaterVolFlowRateMaxDes);
                    } else { // Hard-sized with sizing data
                        if (thisHWBaseboard.WaterVolFlowRateMax > 0.0 && WaterVolFlowRateMaxDes > 0.0) {
                            WaterVolFlowRateMaxUser = thisHWBaseboard.WaterVolFlowRateMax;
                            BaseSizer::reportSizerOutput(state,
                                                         cCMO_BBRadiator_Water,
                                                         thisHWBaseboard.EquipID,
                                                         "Design Size Maximum Water Flow Rate [m3/s]",
                                                         WaterVolFlowRateMaxDes,
                                                         "User-Specified Maximum Water Flow Rate [m3/s]",
                                                         WaterVolFlowRateMaxUser);
                            if (state.dataGlobal->DisplayExtraWarnings) {
                                if ((std::abs(WaterVolFlowRateMaxDes - WaterVolFlowRateMaxUser) / WaterVolFlowRateMaxUser) >
                                    state.dataSize->AutoVsHardSizingThreshold) {
                                    ShowMessage(state,
                                                format("SizeHWBaseboard: Potential issue with equipment sizing for "
                                                       "ZoneHVAC:Baseboard:RadiantConvective:Water=\"{}\".",
                                                       thisHWBaseboard.EquipID));
                                    ShowContinueError(state,
                                                      format("User-Specified Maximum Water Flow Rate of {:.5R} [m3/s]", WaterVolFlowRateMaxUser));
                                    ShowContinueError(
                                        state, format("differs from Design Size Maximum Water Flow Rate of {:.5R} [m3/s]", WaterVolFlowRateMaxDes));
                                    ShowContinueError(state, "This may, or may not, indicate mismatched component sizes.");
                                    ShowContinueError(state, "Verify that the value entered is intended and is consistent with other components.");
                                }
                            }
                        }
                    }
                }
                if (thisHWBaseboard.WaterTempAvg > 0.0 && thisHWBaseboard.WaterMassFlowRateStd > 0.0 &&
                    thisHWBaseboard.RatedCapacity > 0.0) {
                    DesCoilLoad = thisHWBaseboard.RatedCapacity;
                    WaterMassFlowRateStd = thisHWBaseboard.WaterMassFlowRateStd;
                } else if (thisHWBaseboard.RatedCapacity == AutoSize || thisHWBaseboard.RatedCapacity == 0.0) {
                    DesCoilLoad = RatedCapacityDes;
                    rho = GetDensityGlycol(state,
                                           state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                           DataGlobalConstants::HWInitConvTemp,
                                           state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                           RoutineNameFull);
                    WaterMassFlowRateStd = thisHWBaseboard.WaterVolFlowRateMax * rho;
                }
                if (DesCoilLoad >= SmallLoad) {
                    // Calculate UA value
                    // Air mass flow rate is obtained from the following linear equation
                    // m_dot = 0.0062 + 2.75e-05*q
                    AirMassFlowRate = Constant + Coeff * DesCoilLoad;
                    Cp = GetSpecificHeatGlycol(state,
                                               state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                               thisHWBaseboard.WaterTempAvg,
                                               state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                               RoutineName);
                    WaterInletTempStd = (DesCoilLoad / (2.0 * WaterMassFlowRateStd * Cp)) + thisHWBaseboard.WaterTempAvg;
                    WaterOutletTempStd = std::abs((2.0 * thisHWBaseboard.WaterTempAvg) - WaterInletTempStd);
                    AirOutletTempStd = (DesCoilLoad / (AirMassFlowRate * CPAirStd)) + AirInletTempStd;
                    thisHWBaseboard.AirMassFlowRateStd = AirMassFlowRate;
                    // Check Ta,out < Tw,in
                    if (AirOutletTempStd >= WaterInletTempStd) {
                        ShowSevereError(
                            state, format("SizeHWBaseboard: ZoneHVAC:Baseboard:RadiantConvective:Water=\"{}\".", thisHWBaseboard.EquipID));
                        ShowContinueError(state, "...Air Outlet temperature must be below the Water Inlet temperature");
                        ShowContinueError(
                            state,
                            format("...Air Outlet Temperature=[{:.2R}], Water Inlet Temperature=[{:.2R}].", AirOutletTempStd, WaterInletTempStd));
                        AirOutletTempStd = WaterInletTempStd - 0.01;
                        ShowContinueError(state, format("...Air Outlet Temperature set to [{:.2R}].", AirOutletTempStd));
                    }
                    // Check Tw,out < Ta,in
                    if (AirInletTempStd >= WaterOutletTempStd) {
                        ShowSevereError(
                            state, format("SizeHWBaseboard: ZoneHVAC:Baseboard:RadiantConvective:Water=\"{}\".", thisHWBaseboard.EquipID));
                        ShowContinueError(state, "...Water Outlet temperature must be below the Air Inlet temperature");
                        ShowContinueError(
                            state,
                            format("...Air Inlet Temperature=[{:.2R}], Water Outlet Temperature=[{:.2R}].", AirInletTempStd, WaterOutletTempStd));
                        WaterOutletTempStd = AirInletTempStd + 0.01;
                        ShowContinueError(state, format("...Water Outlet Temperature set to [{:.2R}].", WaterOutletTempStd));
                    }
                    // LMTD calculation
                    DeltaT1 = WaterInletTempStd - AirOutletTempStd;
                    DeltaT2 = WaterOutletTempStd - AirInletTempStd;
                    LMTD = (DeltaT1 - DeltaT2) / (std::log(DeltaT1 / DeltaT2));
                    thisHWBaseboard.UA = DesCoilLoad / LMTD;
                } else {
                    thisHWBaseboard.UA = 0.0;
                }
                // Report an UA value
                BaseSizer::reportSizerOutput(
                    state, cCMO_BBRadiator_Water, thisHWBaseboard.EquipID, "U-Factor times Area [W/C]", thisHWBaseboard.UA);
            }
        } else {
            // if there is no heating Sizing:Plant object and autosizng was requested, issue an error message
            if (thisHWBaseboard.WaterVolFlowRateMax == AutoSize || thisHWBaseboard.RatedCapacity == AutoSize ||
                thisHWBaseboard.RatedCapacity == 0.0) {
                ShowSevereError(state, "Autosizing of hot water baseboard requires a heating loop Sizing:Plant object");
                ShowContinueError(state, format("Occurs in Hot Water Baseboard Heater={}", thisHWBaseboard.EquipID));
                ErrorsFound = true;
            }
            // calculate UA from rated capacities
            thisHWBaseboard.RatedCapacity = RatedCapacityDes;
            DesCoilLoad = RatedCapacityDes;

            if (DesCoilLoad >= SmallLoad) {
                WaterMassFlowRateStd = thisHWBaseboard.WaterMassFlowRateStd;
                // m_dot = 0.0062 + 2.75e-05*q
                AirMassFlowRate = Constant + Coeff * DesCoilLoad;
                Cp = GetSpecificHeatGlycol(state,
                                           state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                           thisHWBaseboard.WaterTempAvg,
                                           state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                           RoutineName);
                WaterInletTempStd = (DesCoilLoad / (2.0 * WaterMassFlowRateStd * Cp)) + thisHWBaseboard.WaterTempAvg;
                WaterOutletTempStd = std::abs((2.0 * thisHWBaseboard.WaterTempAvg) - WaterInletTempStd);
                AirOutletTempStd = (DesCoilLoad / (AirMassFlowRate * CPAirStd)) + AirInletTempStd;
                thisHWBaseboard.AirMassFlowRateStd = AirMassFlowRate;

                // Check Ta,out < Tw,in
                if (AirOutletTempStd >= WaterInletTempStd) {
                    ShowSevereError(state,
                                    format("SizeHWBaseboard: ZoneHVAC:Baseboard:RadiantConvective:Water=\"{}\".", thisHWBaseboard.EquipID));
                    ShowContinueError(state, "...Air Outlet temperature must be below the Water Inlet temperature");
                    ShowContinueError(
                        state, format("...Air Outlet Temperature=[{:.2R}], Water Inlet Temperature=[{:.2R}].", AirOutletTempStd, WaterInletTempStd));
                    AirOutletTempStd = WaterInletTempStd - 0.01;
                    ShowContinueError(state, format("...Air Outlet Temperature set to [{:.2R}].", AirOutletTempStd));
                }
                // Check Tw,out < Ta,in
                if (AirInletTempStd >= WaterOutletTempStd) {
                    ShowSevereError(state,
                                    format("SizeHWBaseboard: ZoneHVAC:Baseboard:RadiantConvective:Water=\"{}\".", thisHWBaseboard.EquipID));
                    ShowContinueError(state, "...Water Outlet temperature must be below the Air Inlet temperature");
                    ShowContinueError(
                        state, format("...Air Inlet Temperature=[{:.2R}], Water Outlet Temperature=[{:.2R}].", AirInletTempStd, WaterOutletTempStd));
                    WaterOutletTempStd = AirInletTempStd + 0.01;
                    ShowContinueError(state, format("...Water Outlet Temperature set to [{:.2R}].", WaterOutletTempStd));
                }
                // LMTD calculation
                DeltaT1 = WaterInletTempStd - AirOutletTempStd;
                DeltaT2 = WaterOutletTempStd - AirInletTempStd;
                LMTD = (DeltaT1 - DeltaT2) / (std::log(DeltaT1 / DeltaT2));
                thisHWBaseboard.UA = DesCoilLoad / LMTD;
            } else {
                thisHWBaseboard.UA = 0.0;
            }
            // Report an UA value
            BaseSizer::reportSizerOutput(
                state, cCMO_BBRadiator_Water, thisHWBaseboard.EquipID, "U-Factor times Area [W/C]", thisHWBaseboard.UA);
        }
        // save the design water flow rate for use by the water loop sizing algorithms
        RegisterPlantCompDesignFlow(state, thisHWBaseboard.WaterInletNode, thisHWBaseboard.WaterVolFlowRateMax);

        if (ErrorsFound) {
            ShowFatalError(state, "Preceding sizing errors cause program termination");
        }
    }

    void CalcHWBaseboard(EnergyPlusData &state, int &BaseboardNum, Real64 &LoadMet)
    {
        // SUBROUTINE INFORMATION:
        //       AUTHOR         Russ Taylor
        //       DATE WRITTEN   Nov 1997
        //       MODIFIED       May 2000 Fred Buhl
        //                      Aug 2007 Daeho Kang (Add the calculation of radiant heat source)
        //                      Sep 2011 LKL/BG - resimulate only zones needing it for Radiant systems
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine calculates both the convective and radiant heat transfer rate
        // in a hot water baseboard heater.  The heater is assumed to be crossflowwith
        // both fluids unmixed.  The air flow is buoyancy driven and a constant airflow
        // and a constant airflow velocity of 0.5m/s is assumed.

        // METHODOLOGY EMPLOYED:
        // na

        // REFERENCES:
        // Incropera and DeWitt, Fundamentals of Heat and Mass Transfer
        // Chapter 11.4, p. 523, eq. 11.33

        // Using/Aliasing
        using namespace DataSizing;
        using PlantUtilities::SetActuatedBranchFlowRate;
        using ScheduleManager::GetCurrentScheduleValue;

        // Locals
        // SUBROUTINE ARGUMENT DEFINITIONS:

        // SUBROUTINE PARAMETER DEFINITIONS:
        Real64 constexpr MinFrac(0.0005); // Minimum fraction that delivers radiant heats to surfaces
        static constexpr std::string_view RoutineName("CalcHWBaseboard");

        // INTERFACE BLOCK SPECIFICATIONS
        // na

        // DERIVED TYPE DEFINITIONS
        // na

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int ZoneNum;
        Real64 RadHeat;
        Real64 BBHeat;
        Real64 AirInletTemp;
        Real64 AirOutletTemp;
        Real64 WaterInletTemp;
        Real64 WaterOutletTemp;
        Real64 WaterMassFlowRate;
        Real64 AirMassFlowRate;
        Real64 CapacitanceAir;
        Real64 CapacitanceWater;
        Real64 CapacitanceMax;
        Real64 CapacitanceMin;
        Real64 CapacityRatio;
        Real64 NTU;
        Real64 Effectiveness;
        Real64 AA;
        Real64 BB;
        Real64 CC;
        Real64 QZnReq;
        Real64 Cp;

        auto &HWBaseboardDesignObject = state.dataHWBaseboardRad->HWBaseboardDesignObject;
        auto &QBBRadSource = state.dataHWBaseboardRad->QBBRadSource;
        auto &ZeroSourceSumHATsurf = state.dataHWBaseboardRad->ZeroSourceSumHATsurf;

        auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
        ZoneNum = thisHWBaseboard.ZonePtr;
        QZnReq = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(ZoneNum).RemainingOutputReqToHeatSP;
        AirInletTemp = thisHWBaseboard.AirInletTemp;
        AirOutletTemp = AirInletTemp;
        WaterInletTemp = thisHWBaseboard.WaterInletTemp;
        WaterOutletTemp = WaterInletTemp;
        WaterMassFlowRate = state.dataLoopNodes->Node(thisHWBaseboard.WaterInletNode).MassFlowRate;
        HWBaseboardDesignData HWBaseboardDesignDataObject{
            HWBaseboardDesignObject(thisHWBaseboard.DesignObjectPtr)}; // Contains the data for the design object

        if (QZnReq > SmallLoad && !state.dataZoneEnergyDemand->CurDeadBandOrSetback(ZoneNum) &&
            (GetCurrentScheduleValue(state, thisHWBaseboard.SchedPtr) > 0) && (WaterMassFlowRate > 0.0)) {
            // Calculate air mass flow rate
            AirMassFlowRate = thisHWBaseboard.AirMassFlowRateStd * (WaterMassFlowRate / thisHWBaseboard.WaterMassFlowRateMax);
            CapacitanceAir = PsyCpAirFnW(thisHWBaseboard.AirInletHumRat) * AirMassFlowRate;
            Cp = GetSpecificHeatGlycol(state,
                                       state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidName,
                                       WaterInletTemp,
                                       state.dataPlnt->PlantLoop(thisHWBaseboard.plantLoc.loopNum).FluidIndex,
                                       RoutineName);

            CapacitanceWater = Cp * WaterMassFlowRate;
            CapacitanceMax = max(CapacitanceAir, CapacitanceWater);
            CapacitanceMin = min(CapacitanceAir, CapacitanceWater);
            CapacityRatio = CapacitanceMin / CapacitanceMax;
            NTU = thisHWBaseboard.UA / CapacitanceMin;

            // The effectiveness is given by the following formula:
            // Effectiveness = 1. - EXP((1./CapacityRatio)*(NTU)**0.22*(EXP(-CapacityRatio*(NTU)**0.78)-1.))
            // To prevent possible underflows (numbers smaller than the computer can handle) we must break
            // the calculation up into steps and check the size of the exponential arguments.
            AA = -CapacityRatio * std::pow(NTU, 0.78);
            if (AA < -20.0) {
                BB = 0.0;
            } else {
                BB = std::exp(AA);
            }
            CC = (1.0 / CapacityRatio) * std::pow(NTU, 0.22) * (BB - 1.0);
            if (CC < -20.0) {
                Effectiveness = 1.0;
            } else {
                Effectiveness = 1.0 - std::exp(CC);
            }

            AirOutletTemp = AirInletTemp + Effectiveness * CapacitanceMin * (WaterInletTemp - AirInletTemp) / CapacitanceAir;
            WaterOutletTemp = WaterInletTemp - CapacitanceAir * (AirOutletTemp - AirInletTemp) / CapacitanceWater;
            BBHeat = CapacitanceWater * (WaterInletTemp - WaterOutletTemp);
            RadHeat = BBHeat * HWBaseboardDesignDataObject.FracRadiant;
            QBBRadSource(BaseboardNum) = RadHeat;

            if (HWBaseboardDesignDataObject.FracRadiant <= MinFrac) {
                LoadMet = BBHeat;
            } else {

                // Now, distribute the radiant energy of all systems to the appropriate surfaces, to people, and the air
                DistributeBBRadGains(state);
                // Now "simulate" the system by recalculating the heat balances
                HeatBalanceSurfaceManager::CalcHeatBalanceOutsideSurf(state, ZoneNum);

                HeatBalanceSurfaceManager::CalcHeatBalanceInsideSurf(state, ZoneNum);

                // Here an assumption is made regarding radiant heat transfer to people.
                // While the radiant heat transfer to people array will be used by the thermal comfort
                // routines, the energy transfer to people would get lost from the perspective
                // of the heat balance.  So, to avoid this net loss of energy which clearly
                // gets added to the zones, we must account for it somehow.  This assumption
                // that all energy radiated to people is converted to convective energy is
                // not very precise, but at least it conserves energy. The system impact to heat balance
                // should include this.
                LoadMet = (state.dataHeatBal->Zone(ZoneNum).sumHATsurf(state) - ZeroSourceSumHATsurf(ZoneNum)) +
                          (BBHeat * thisHWBaseboard.FracConvect) + (RadHeat * HWBaseboardDesignDataObject.FracDistribPerson);
            }
            thisHWBaseboard.WaterOutletEnthalpy = thisHWBaseboard.WaterInletEnthalpy - BBHeat / WaterMassFlowRate;
        } else {
            CapacitanceWater = 0.0;
            CapacitanceMax = 0.0;
            CapacitanceMin = 0.0;
            NTU = 0.0;
            Effectiveness = 0.0;
            AirOutletTemp = AirInletTemp;
            WaterOutletTemp = WaterInletTemp;
            BBHeat = 0.0;
            LoadMet = 0.0;
            RadHeat = 0.0;
            WaterMassFlowRate = 0.0;
            AirMassFlowRate = 0.0;
            QBBRadSource(BaseboardNum) = 0.0;
            thisHWBaseboard.WaterOutletEnthalpy = thisHWBaseboard.WaterInletEnthalpy;
            SetActuatedBranchFlowRate(state, WaterMassFlowRate, thisHWBaseboard.WaterInletNode, thisHWBaseboard.plantLoc, false);
        }

        thisHWBaseboard.WaterOutletTemp = WaterOutletTemp;
        thisHWBaseboard.AirOutletTemp = AirOutletTemp;
        thisHWBaseboard.WaterMassFlowRate = WaterMassFlowRate;
        thisHWBaseboard.AirMassFlowRate = AirMassFlowRate;
        thisHWBaseboard.TotPower = LoadMet;
        thisHWBaseboard.Power = BBHeat;
        thisHWBaseboard.ConvPower = BBHeat - RadHeat;
        thisHWBaseboard.RadPower = RadHeat;
    }

    void UpdateHWBaseboard(EnergyPlusData &state, int const BaseboardNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Russ Taylor
        //                      Rick Strand
        //       DATE WRITTEN   Nov 1997
        //                      February 2001
        //       MODIFIED       Aug 2007 Daeho Kang (Add the update of radiant source)
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:

        // METHODOLOGY EMPLOYED:
        // The update subrotines both in high temperature radiant radiator
        // and convective only baseboard radiator are combined and modified.

        // REFERENCES:
        // na

        // Using/Aliasing
        using PlantUtilities::SafeCopyPlantNode;

        // Locals
        // SUBROUTINE ARGUMENT DEFINITIONS:

        // SUBROUTINE PARAMETER DEFINITIONS:
        // na

        // INTERFACE BLOCK SPECIFICATIONS
        // na

        // DERIVED TYPE DEFINITIONS
        // na

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int WaterInletNode;
        int WaterOutletNode;
        auto &Iter = state.dataHWBaseboardRad->Iter;
        auto &MyEnvrnFlag2 = state.dataHWBaseboardRad->MyEnvrnFlag2;
        auto &LastSysTimeElapsed = state.dataHWBaseboardRad->LastSysTimeElapsed;
        auto &QBBRadSrcAvg = state.dataHWBaseboardRad->QBBRadSrcAvg;
        auto &LastQBBRadSrc = state.dataHWBaseboardRad->LastQBBRadSrc;
        auto &QBBRadSource = state.dataHWBaseboardRad->QBBRadSource;
        auto &LastTimeStepSys = state.dataHWBaseboardRad->LastTimeStepSys;

        if (state.dataGlobal->BeginEnvrnFlag && MyEnvrnFlag2) {
            Iter = 0;
            MyEnvrnFlag2 = false;
        }
        if (!state.dataGlobal->BeginEnvrnFlag) {
            MyEnvrnFlag2 = true;
        }

        // First, update the running average if necessary...
        if (LastSysTimeElapsed(BaseboardNum) == state.dataHVACGlobal->SysTimeElapsed) {
            QBBRadSrcAvg(BaseboardNum) -= LastQBBRadSrc(BaseboardNum) * LastTimeStepSys(BaseboardNum) / state.dataGlobal->TimeStepZone;
        }
        // Update the running average and the "last" values with the current values of the appropriate variables
        QBBRadSrcAvg(BaseboardNum) += QBBRadSource(BaseboardNum) * state.dataHVACGlobal->TimeStepSys / state.dataGlobal->TimeStepZone;

        LastQBBRadSrc(BaseboardNum) = QBBRadSource(BaseboardNum);
        LastSysTimeElapsed(BaseboardNum) = state.dataHVACGlobal->SysTimeElapsed;
        LastTimeStepSys(BaseboardNum) = state.dataHVACGlobal->TimeStepSys;

        auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
        WaterInletNode = thisHWBaseboard.WaterInletNode;
        WaterOutletNode = thisHWBaseboard.WaterOutletNode;

        // Set the outlet air nodes of the Baseboard
        // Set the outlet water nodes for the Coil
        SafeCopyPlantNode(state, WaterInletNode, WaterOutletNode);
        state.dataLoopNodes->Node(WaterOutletNode).Temp = thisHWBaseboard.WaterOutletTemp;
        state.dataLoopNodes->Node(WaterOutletNode).Enthalpy = thisHWBaseboard.WaterOutletEnthalpy;
    }

    void UpdateBBRadSourceValAvg(EnergyPlusData &state, bool &HWBaseboardSysOn) // .TRUE. if the radiant system has run this zone time step
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Rick Strand
        //       DATE WRITTEN   February 2001
        //       MODIFIED       Aug 2007 Daeho Kang (Modification only for baseboard)
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // To transfer the average value of the heat source over the entire
        // zone time step back to the heat balance routines so that the heat
        // balance algorithms can simulate one last time with the average source
        // to maintain some reasonable amount of continuity and energy balance
        // in the temperature and flux histories.

        // METHODOLOGY EMPLOYED:
        // All of the record keeping for the average term is done in the Update
        // routine so the only other thing that this subroutine does is check to
        // see if the system was even on.  If any average term is non-zero, then
        // one or more of the radiant systems was running.

        HWBaseboardSysOn = false;

        // If this was never allocated, then there are no radiant systems in this input file (just RETURN)
        if (!allocated(state.dataHWBaseboardRad->QBBRadSrcAvg)) return;

        // If it was allocated, then we have to check to see if this was running at all...
        for (int BaseboardNum = 1; BaseboardNum <= state.dataHWBaseboardRad->NumHWBaseboards; ++BaseboardNum) {
            if (state.dataHWBaseboardRad->QBBRadSrcAvg(BaseboardNum) != 0.0) {
                HWBaseboardSysOn = true;
                break; // DO loop
            }
        }

        state.dataHWBaseboardRad->QBBRadSource = state.dataHWBaseboardRad->QBBRadSrcAvg;

        DistributeBBRadGains(state); // QBBRadSource has been modified so we need to redistribute gains
    }

    void DistributeBBRadGains(EnergyPlusData &state)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Rick Strand
        //       DATE WRITTEN   February 2001
        //       MODIFIED       Aug. 2007 Daeho Kang (Modification only for baseboard)
        //                      April 2010 Brent Griffith, max limit to protect surface temperature calcs
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // To distribute the gains from the hot water basebaord heater
        // as specified in the user input file.  This includes distribution
        // of long wavelength radiant gains to surfaces and "people."

        // METHODOLOGY EMPLOYED:
        // We must cycle through all of the radiant systems because each
        // surface could feel the effect of more than one radiant system.
        // Note that the energy radiated to people is assumed to affect them
        // but them it is assumed to be convected to the air.

        // Using/Aliasing
        using DataHeatBalFanSys::MaxRadHeatFlux;

        // SUBROUTINE PARAMETER DEFINITIONS:
        Real64 constexpr SmallestArea(0.001); // Smallest area in meters squared (to avoid a divide by zero)

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int RadSurfNum;           // Counter for surfaces receiving radiation from radiant heater
        int SurfNum;              // Pointer to the Surface derived type
        int ZoneNum;              // Pointer to the Zone derived type
        Real64 ThisSurfIntensity; // temporary for W/m2 term for rad on a surface

        auto &HWBaseboardDesignObject = state.dataHWBaseboardRad->HWBaseboardDesignObject;
        auto &QBBRadSource = state.dataHWBaseboardRad->QBBRadSource;

        // Initialize arrays
        state.dataHeatBalFanSys->SurfQHWBaseboard = 0.0;
        state.dataHeatBalFanSys->ZoneQHWBaseboardToPerson = 0.0;

        for (int BaseboardNum = 1; BaseboardNum <= state.dataHWBaseboardRad->NumHWBaseboards; ++BaseboardNum) {

            auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
            HWBaseboardDesignData HWBaseboardDesignDataObject{
                HWBaseboardDesignObject(thisHWBaseboard.DesignObjectPtr)}; // Contains the data for the design object
            ZoneNum = thisHWBaseboard.ZonePtr;
            if (ZoneNum <= 0) continue;
            state.dataHeatBalFanSys->ZoneQHWBaseboardToPerson(ZoneNum) += QBBRadSource(BaseboardNum) * HWBaseboardDesignDataObject.FracDistribPerson;

            for (RadSurfNum = 1; RadSurfNum <= thisHWBaseboard.TotSurfToDistrib; ++RadSurfNum) {
                SurfNum = thisHWBaseboard.SurfacePtr(RadSurfNum);
                if (state.dataSurface->Surface(SurfNum).Area > SmallestArea) {
                    ThisSurfIntensity = (QBBRadSource(BaseboardNum) * thisHWBaseboard.FracDistribToSurf(RadSurfNum) /
                                         state.dataSurface->Surface(SurfNum).Area);
                    state.dataHeatBalFanSys->SurfQHWBaseboard(SurfNum) += ThisSurfIntensity;
                    state.dataHeatBalSurf->AnyRadiantSystems = true;
                    // CR 8074, trap for excessive intensity (throws off surface balance )
                    if (ThisSurfIntensity > MaxRadHeatFlux) {
                        ShowSevereError(state, "DistributeBBRadGains:  excessive thermal radiation heat flux intensity detected");
                        ShowContinueError(state, format("Surface = {}", state.dataSurface->Surface(SurfNum).Name));
                        ShowContinueError(state, format("Surface area = {:.3R} [m2]", state.dataSurface->Surface(SurfNum).Area));
                        ShowContinueError(state, format("Occurs in {} = {}", cCMO_BBRadiator_Water, thisHWBaseboard.EquipID));
                        ShowContinueError(state, format("Radiation intensity = {:.2R} [W/m2]", ThisSurfIntensity));
                        ShowContinueError(state, format("Assign a larger surface area or more surfaces in {}", cCMO_BBRadiator_Water));
                        ShowFatalError(state, "DistributeBBRadGains:  excessive thermal radiation heat flux intensity detected");
                    }
                } else {
                    ShowSevereError(state, "DistributeBBRadGains:  surface not large enough to receive thermal radiation heat flux");
                    ShowContinueError(state, format("Surface = {}", state.dataSurface->Surface(SurfNum).Name));
                    ShowContinueError(state, format("Surface area = {:.3R} [m2]", state.dataSurface->Surface(SurfNum).Area));
                    ShowContinueError(state, format("Occurs in {} = {}", cCMO_BBRadiator_Water, thisHWBaseboard.EquipID));
                    ShowContinueError(state, format("Assign a larger surface area or more surfaces in {}", cCMO_BBRadiator_Water));
                    ShowFatalError(state, "DistributeBBRadGains:  surface not large enough to receive thermal radiation heat flux");
                }
            }
        }
    }

    void ReportHWBaseboard(EnergyPlusData &state, int const BaseboardNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Daeho Kang
        //       DATE WRITTEN   Aug 2007
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
        thisHWBaseboard.TotEnergy = thisHWBaseboard.TotPower * state.dataHVACGlobal->TimeStepSys * DataGlobalConstants::SecInHour;
        thisHWBaseboard.Energy = thisHWBaseboard.Power * state.dataHVACGlobal->TimeStepSys * DataGlobalConstants::SecInHour;
        thisHWBaseboard.ConvEnergy =
            thisHWBaseboard.ConvPower * state.dataHVACGlobal->TimeStepSys * DataGlobalConstants::SecInHour;
        thisHWBaseboard.RadEnergy = thisHWBaseboard.RadPower * state.dataHVACGlobal->TimeStepSys * DataGlobalConstants::SecInHour;
    }

    void UpdateHWBaseboardPlantConnection(EnergyPlusData &state,
                                          int const BaseboardTypeNum,                                  // type index
                                          std::string const &BaseboardName,                            // component name
                                          [[maybe_unused]] int const EquipFlowCtrl,                    // Flow control mode for the equipment
                                          [[maybe_unused]] int const LoopNum,                          // Plant loop index for where called from
                                          [[maybe_unused]] const DataPlant::LoopSideLocation LoopSide, // Plant loop side index for where called from
                                          int &CompIndex,                                              // Chiller number pointer
                                          [[maybe_unused]] bool const FirstHVACIteration,
                                          bool &InitLoopEquip // If not zero, calculate the max load for operating conditions
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Brent Griffith
        //       DATE WRITTEN   Sept. 2010
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // update sim routine called from plant

        // METHODOLOGY EMPLOYED:
        // check input, provide comp index, call utility routines

        // Using/Aliasing
        using DataPlant::PlantEquipTypeNames;

        using PlantUtilities::PullCompInterconnectTrigger;

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        auto &HWBaseboard = state.dataHWBaseboardRad->HWBaseboard;

        int BaseboardNum;

        // Find the correct baseboard
        if (CompIndex == 0) {
            BaseboardNum = UtilityRoutines::FindItemInList(BaseboardName, HWBaseboard, &HWBaseboardParams::EquipID);
            if (BaseboardNum == 0) {
                ShowFatalError(state, format("UpdateHWBaseboardPlantConnection: Specified baseboard not valid ={}", BaseboardName));
            }
            CompIndex = BaseboardNum;
        } else {
            BaseboardNum = CompIndex;
            if (BaseboardNum > state.dataHWBaseboardRad->NumHWBaseboards || BaseboardNum < 1) {
                ShowFatalError(
                    state,
                    format("UpdateHWBaseboardPlantConnection:  Invalid CompIndex passed={}, Number of baseboards={}, Entered baseboard name={}",
                           BaseboardNum,
                           state.dataHWBaseboardRad->NumHWBaseboards,
                           BaseboardName));
            }
            auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);

            if (state.dataGlobal->KickOffSimulation) {
                if (BaseboardName != thisHWBaseboard.EquipID) {
                    ShowFatalError(state,
                                   format("UpdateHWBaseboardPlantConnection: Invalid CompIndex passed={}, baseboard name={}, stored baseboard Name "
                                          "for that index={}",
                                          BaseboardNum,
                                          BaseboardName,
                                          thisHWBaseboard.EquipID));
                }
                if (BaseboardTypeNum != static_cast<int>(DataPlant::PlantEquipmentType::Baseboard_Rad_Conv_Water)) {
                    ShowFatalError(state,
                                   format("UpdateHWBaseboardPlantConnection: Invalid CompIndex passed={}, baseboard name={}, stored baseboard Name "
                                          "for that index={}",
                                          BaseboardNum,
                                          BaseboardName,
                                          PlantEquipTypeNames[BaseboardTypeNum]));
                }
            }
        }

        if (InitLoopEquip) {
            return;
        }

        auto &thisHWBaseboard = state.dataHWBaseboardRad->HWBaseboard(BaseboardNum);
        PullCompInterconnectTrigger(state,
                                    thisHWBaseboard.plantLoc,
                                    thisHWBaseboard.BBLoadReSimIndex,
                                    thisHWBaseboard.plantLoc,
                                    DataPlant::CriteriaType::HeatTransferRate,
                                    thisHWBaseboard.Power);

        PullCompInterconnectTrigger(state,
                                    thisHWBaseboard.plantLoc,
                                    thisHWBaseboard.BBMassFlowReSimIndex,
                                    thisHWBaseboard.plantLoc,
                                    DataPlant::CriteriaType::MassFlowRate,
                                    thisHWBaseboard.WaterMassFlowRate);

        PullCompInterconnectTrigger(state,
                                    thisHWBaseboard.plantLoc,
                                    thisHWBaseboard.BBInletTempFlowReSimIndex,
                                    thisHWBaseboard.plantLoc,
                                    DataPlant::CriteriaType::Temperature,
                                    thisHWBaseboard.WaterOutletTemp);
    }

    //*****************************************************************************************

} // namespace HWBaseboardRadiator

} // namespace EnergyPlus

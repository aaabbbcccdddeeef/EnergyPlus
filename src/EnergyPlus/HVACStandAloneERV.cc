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
#include <EnergyPlus/Autosizing/SystemAirFlowSizing.hh>
#include <EnergyPlus/BranchNodeConnections.hh>
#include <EnergyPlus/CurveManager.hh>
#include <EnergyPlus/Data/EnergyPlusData.hh>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/DataHeatBalance.hh>
#include <EnergyPlus/DataLoopNode.hh>
#include <EnergyPlus/DataSizing.hh>
#include <EnergyPlus/DataZoneControls.hh>
#include <EnergyPlus/DataZoneEquipment.hh>
#include <EnergyPlus/Fans.hh>
#include <EnergyPlus/General.hh>
#include <EnergyPlus/GeneralRoutines.hh>
#include <EnergyPlus/GlobalNames.hh>
#include <EnergyPlus/HVACFan.hh>
#include <EnergyPlus/HVACStandAloneERV.hh>
#include <EnergyPlus/HeatRecovery.hh>
#include <EnergyPlus/InputProcessing/InputProcessor.hh>
#include <EnergyPlus/MixedAir.hh>
#include <EnergyPlus/NodeInputManager.hh>
#include <EnergyPlus/OutAirNodeManager.hh>
#include <EnergyPlus/OutputProcessor.hh>
#include <EnergyPlus/ScheduleManager.hh>
#include <EnergyPlus/UtilityRoutines.hh>

namespace EnergyPlus::HVACStandAloneERV {

// Module containing the routines dealing with stand alone energy recovery ventilators (ERVs)

// MODULE INFORMATION:
//       AUTHOR         Richard Raustad, FSEC
//       DATE WRITTEN   June 2003
//       MODIFIED       na
//       RE-ENGINEERED  na

// PURPOSE OF THIS MODULE:
// To encapsulate the data and algorithms needed to simulate stand alone
// energy recovery ventilators that condition outdoor ventilation air and
// supply that air directly to a zone.

// METHODOLOGY EMPLOYED:
// These units are modeled as a collection of components: air-to-air generic heat exchanger,
// supply air fan, exhaust air fan and an optional controller to avoid overheating
// of the supply air (economizer or free cooling operation).

// REFERENCES: none

// OTHER NOTES: none

// USE STATEMENTS:
// Use statements for data only modules
// Using/Aliasing
using namespace DataLoopNode;
using namespace DataHVACGlobals;
using Fans::GetFanVolFlow;
using ScheduleManager::GetCurrentScheduleValue;
using ScheduleManager::GetScheduleIndex;

void SimStandAloneERV(EnergyPlusData &state,
                      std::string_view CompName,     // name of the Stand Alone ERV unit
                      int const ZoneNum,             // number of zone being served unused1208
                      bool const FirstHVACIteration, // TRUE if 1st HVAC simulation of system timestep
                      Real64 &SensLoadMet,           // net sensible load supplied by the ERV unit to the zone (W)
                      Real64 &LatLoadMet,            // net latent load supplied by ERV unit to the zone (kg/s),
                      int &CompIndex                 // pointer to correct component
)
{

    // SUBROUTINE INFORMATION:
    //       AUTHOR         Richard Raustad, FSEC
    //       DATE WRITTEN   June 2003
    //       MODIFIED       Don Shirey, Aug 2009 (LatLoadMet)
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // Manages the simulation of a Stand Alone ERV unit. Called from SimZoneEquipment

    // Using/Aliasing

    // Locals
    // SUBROUTINE ARGUMENT DEFINITIONS:
    // ZoneNum not used at this time, future modifications may require zone information
    // dehumid = negative

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    int StandAloneERVNum; // index of Stand Alone ERV unit being simulated

    // First time SimStandAloneERV is called, get the input for all Stand Alone ERV units
    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    // Find the correct Stand Alone ERV unit index
    if (CompIndex == 0) {
        StandAloneERVNum = UtilityRoutines::FindItem(CompName, state.dataHVACStandAloneERV->StandAloneERV);
        if (StandAloneERVNum == 0) {
            ShowFatalError(state, format("SimStandAloneERV: Unit not found={}", CompName));
        }
        CompIndex = StandAloneERVNum;
    } else {
        StandAloneERVNum = CompIndex;
        if (StandAloneERVNum > state.dataHVACStandAloneERV->NumStandAloneERVs || StandAloneERVNum < 1) {
            ShowFatalError(state,
                           format("SimStandAloneERV:  Invalid CompIndex passed={}, Number of Units={}, Entered Unit name={}",
                                  StandAloneERVNum,
                                  state.dataHVACStandAloneERV->NumStandAloneERVs,
                                  CompName));
        }
        if (state.dataHVACStandAloneERV->CheckEquipName(StandAloneERVNum)) {
            if (CompName != state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum).Name) {
                ShowFatalError(state,
                               format("SimStandAloneERV: Invalid CompIndex passed={}, Unit name={}, stored Unit Name for that index={}",
                                      StandAloneERVNum,
                                      CompName,
                                      state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum).Name));
            }
            state.dataHVACStandAloneERV->CheckEquipName(StandAloneERVNum) = false;
        }
    }

    // Initialize the Stand Alone ERV unit
    InitStandAloneERV(state, StandAloneERVNum, ZoneNum, FirstHVACIteration);

    CalcStandAloneERV(state, StandAloneERVNum, FirstHVACIteration, SensLoadMet, LatLoadMet);

    ReportStandAloneERV(state, StandAloneERVNum);
}

void GetStandAloneERV(EnergyPlusData &state)
{

    // SUBROUTINE INFORMATION:
    //       AUTHOR         Richard Raustad
    //       DATE WRITTEN   June 2003
    //       MODIFIED       July 2012, Chandan Sharma - FSEC: Added zone sys avail managers
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // Obtains input data for Stand Alone ERV units and stores it in the Stand Alone ERV data structure

    // METHODOLOGY EMPLOYED:
    // Uses "Get" routines to read in data.

    // Using/Aliasing
    using BranchNodeConnections::SetUpCompSets;
    using DataSizing::AutoSize;
    using Fans::GetFanAvailSchPtr;
    using Fans::GetFanDesignVolumeFlowRate;
    using Fans::GetFanIndex;
    using Fans::GetFanOutletNode;
    using Fans::GetFanType;

    using NodeInputManager::GetOnlySingleNode;
    using HeatRecovery::GetHeatExchangerObjectTypeNum;
    using Curve::GetCurveIndex;
    using OutAirNodeManager::CheckOutAirNodeNumber;

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    int StandAloneERVIndex;  // loop index
    int StandAloneERVNum;    // current Stand Alone ERV number
    Array1D_string Alphas;   // Alpha items for object
    Array1D<Real64> Numbers; // Numeric items for object
    Array1D_string cAlphaFields;
    Array1D_string cNumericFields;
    Array1D_bool lAlphaBlanks;
    Array1D_bool lNumericBlanks;
    std::string CompSetSupplyFanInlet;
    std::string CompSetSupplyFanOutlet;
    std::string CompSetExhaustFanInlet;
    std::string CompSetExhaustFanOutlet;
    std::string CurrentModuleObject; // Object type for getting and error messages
    int SAFanTypeNum;                // Integer equivalent to fan type
    int EAFanTypeNum;                // Integer equivalent to fan type
    int NumArg;
    int NumAlphas;             // Number of Alphas for each GetObjectItem call
    int NumNumbers;            // Number of Numbers for each GetObjectItem call
    int MaxAlphas;             // Max between the two objects gotten here
    int MaxNumbers;            // Max between the two objects gotten here
    int IOStatus;              // Used in GetObjectItem
    bool ErrorsFound(false);   // Set to true if errors in input, fatal at end of routine
    Real64 AirFlowRate;        // used to find zone with humidistat
    int NodeNumber;            // used to find zone with humidistat
    int HStatZoneNum;          // used to find zone with humidistat
    int NumHstatZone;          // index to humidity controlled zones
    bool ZoneNodeFound(false); // used to find zone with humidistat
    bool HStatFound(false);    // used to find zone with humidistat
    bool errFlag;              // Error flag used in mining calls
    Real64 SAFanVolFlowRate;   // supply air fan volumetric flow rate [m3/s]
    Real64 EAFanVolFlowRate;   // exhaust air fan volumetric flow rate [m3/s]
    Real64 HXSupAirFlowRate;   // HX supply air flow rate [m3/s]
    Real64 HighRHOARatio;      // local variable for HighRHOAFlowRatio
    bool ZoneInletNodeFound;   // used for warning when zone node not listed in equipment connections
    bool ZoneExhaustNodeFound; // used for warning when zone node not listed in equipment connections
    int ZoneInletCZN;          // used for warning when zone node not listed in equipment connections
    int ZoneExhaustCZN;        // used for warning when zone node not listed in equipment connections

    state.dataInputProcessing->inputProcessor->getObjectDefMaxArgs(state, "ZoneHVAC:EnergyRecoveryVentilator", NumArg, NumAlphas, NumNumbers);
    MaxAlphas = NumAlphas;
    MaxNumbers = NumNumbers;
    state.dataInputProcessing->inputProcessor->getObjectDefMaxArgs(
        state, "ZoneHVAC:EnergyRecoveryVentilator:Controller", NumArg, NumAlphas, NumNumbers);
    MaxAlphas = max(MaxAlphas, NumAlphas);
    MaxNumbers = max(MaxNumbers, NumNumbers);

    Alphas.allocate(MaxAlphas);
    Numbers.dimension(MaxNumbers, 0.0);
    cAlphaFields.allocate(MaxAlphas);
    cNumericFields.allocate(MaxNumbers);
    lNumericBlanks.dimension(MaxNumbers, false);
    lAlphaBlanks.dimension(MaxAlphas, false);

    state.dataHVACStandAloneERV->GetERVInputFlag = false;

    // find the number of each type of Stand Alone ERV unit
    CurrentModuleObject = "ZoneHVAC:EnergyRecoveryVentilator";

    state.dataHVACStandAloneERV->NumStandAloneERVs = state.dataInputProcessing->inputProcessor->getNumObjectsFound(state, CurrentModuleObject);

    // allocate the data structures
    state.dataHVACStandAloneERV->StandAloneERV.allocate(state.dataHVACStandAloneERV->NumStandAloneERVs);
    state.dataHVACStandAloneERV->HeatExchangerUniqueNames.reserve(static_cast<unsigned>(state.dataHVACStandAloneERV->NumStandAloneERVs));
    state.dataHVACStandAloneERV->SupplyAirFanUniqueNames.reserve(static_cast<unsigned>(state.dataHVACStandAloneERV->NumStandAloneERVs));
    state.dataHVACStandAloneERV->ExhaustAirFanUniqueNames.reserve(static_cast<unsigned>(state.dataHVACStandAloneERV->NumStandAloneERVs));
    state.dataHVACStandAloneERV->ControllerUniqueNames.reserve(static_cast<unsigned>(state.dataHVACStandAloneERV->NumStandAloneERVs));
    state.dataHVACStandAloneERV->CheckEquipName.dimension(state.dataHVACStandAloneERV->NumStandAloneERVs, true);

    // loop over Stand Alone ERV units; get and load the input data
    for (StandAloneERVIndex = 1; StandAloneERVIndex <= state.dataHVACStandAloneERV->NumStandAloneERVs; ++StandAloneERVIndex) {

        state.dataInputProcessing->inputProcessor->getObjectItem(state,
                                                                 CurrentModuleObject,
                                                                 StandAloneERVIndex,
                                                                 Alphas,
                                                                 NumAlphas,
                                                                 Numbers,
                                                                 NumNumbers,
                                                                 IOStatus,
                                                                 lNumericBlanks,
                                                                 lAlphaBlanks,
                                                                 cAlphaFields,
                                                                 cNumericFields);
        StandAloneERVNum = StandAloneERVIndex; // separate variables in case other objects read by this module at some point later

	auto &thisStandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum);
        thisStandAloneERV.Name = Alphas(1);
        thisStandAloneERV.UnitType = CurrentModuleObject;

        if (lAlphaBlanks(2)) {
            thisStandAloneERV.SchedPtr = DataGlobalConstants::ScheduleAlwaysOn;
        } else {
            thisStandAloneERV.SchedPtr =
                GetScheduleIndex(state, Alphas(2)); // convert schedule name to pointer
            if (thisStandAloneERV.SchedPtr == 0) {
                ShowSevereError(state,
                                format("{}, \"{}\" {} not found = {}",
                                       CurrentModuleObject,
                                       thisStandAloneERV.Name,
                                       cAlphaFields(2),
                                       Alphas(2)));
                ErrorsFound = true;
            }
        }

        GlobalNames::IntraObjUniquenessCheck(
            state, Alphas(3), CurrentModuleObject, cAlphaFields(3), state.dataHVACStandAloneERV->HeatExchangerUniqueNames, ErrorsFound);
        thisStandAloneERV.HeatExchangerName = Alphas(3);
        errFlag = false;
        thisStandAloneERV.HeatExchangerTypeNum =
            GetHeatExchangerObjectTypeNum(state, thisStandAloneERV.HeatExchangerName, errFlag);
        if (errFlag) {
            ShowContinueError(
                state, format("... occurs in {} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ErrorsFound = true;
        }

        errFlag = false;
        HXSupAirFlowRate =
            HeatRecovery::GetSupplyAirFlowRate(state, thisStandAloneERV.HeatExchangerName, errFlag);
        if (errFlag) {
            ShowContinueError(
                state, format("... occurs in {} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ErrorsFound = true;
        }
        thisStandAloneERV.DesignHXVolFlowRate = HXSupAirFlowRate;

        thisStandAloneERV.SupplyAirFanName = Alphas(4);
        GlobalNames::IntraObjUniquenessCheck(
            state, Alphas(4), CurrentModuleObject, cAlphaFields(4), state.dataHVACStandAloneERV->SupplyAirFanUniqueNames, ErrorsFound);

        errFlag = false;
        if (HVACFan::checkIfFanNameIsAFanSystem(state,
                                                thisStandAloneERV
                                                    .SupplyAirFanName)) { // no object type in input, so check if Fan:SystemModel
            thisStandAloneERV.SupplyAirFanType_Num = DataHVACGlobals::FanType_SystemModelObject;
            state.dataHVACFan->fanObjs.emplace_back(new HVACFan::FanSystem(state, thisStandAloneERV.SupplyAirFanName)); // call constructor
            thisStandAloneERV.SupplyAirFanIndex = HVACFan::getFanObjectVectorIndex(state, thisStandAloneERV.SupplyAirFanName);
            thisStandAloneERV.SupplyAirFanSchPtr = state.dataHVACFan->fanObjs[thisStandAloneERV.SupplyAirFanIndex]->availSchedIndex;
            thisStandAloneERV.DesignSAFanVolFlowRate = state.dataHVACFan->fanObjs[thisStandAloneERV.SupplyAirFanIndex]->designAirVolFlowRate;
            thisStandAloneERV.SupplyAirOutletNode = state.dataHVACFan->fanObjs[thisStandAloneERV.SupplyAirFanIndex]->outletNodeNum;
        } else {
            GetFanType(state,
                       thisStandAloneERV.SupplyAirFanName,
                       SAFanTypeNum,
                       errFlag,
                       CurrentModuleObject,
                       thisStandAloneERV.Name);
            if (errFlag) {
                ErrorsFound = true;
            }
            thisStandAloneERV.SupplyAirFanType_Num = SAFanTypeNum;

            errFlag = false;
            thisStandAloneERV.SupplyAirFanSchPtr = GetFanAvailSchPtr(
                state, cFanTypes(SAFanTypeNum), thisStandAloneERV.SupplyAirFanName, errFlag);
            if (errFlag) {
                ShowContinueError(
                    state, format("... occurs in {} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
                ErrorsFound = true;
            }

            GetFanIndex(state,
                        thisStandAloneERV.SupplyAirFanName,
                        thisStandAloneERV.SupplyAirFanIndex,
                        errFlag,
                        CurrentModuleObject + " \"" + thisStandAloneERV.Name + "\"");

            // Set the SA Design Fan Volume Flow Rate
            // get from fan module
            errFlag = false;
            SAFanVolFlowRate = GetFanDesignVolumeFlowRate(
                state, cFanTypes(SAFanTypeNum), thisStandAloneERV.SupplyAirFanName, errFlag);
            if (errFlag) {
                ShowContinueError(
                    state, format("... occurs in {} ={}", CurrentModuleObject, thisStandAloneERV.Name));
                ErrorsFound = true;
            }
            thisStandAloneERV.DesignSAFanVolFlowRate = SAFanVolFlowRate;
            errFlag = false;
            thisStandAloneERV.SupplyAirOutletNode = GetFanOutletNode(
                state, cFanTypes(SAFanTypeNum), thisStandAloneERV.SupplyAirFanName, errFlag);
        }

        thisStandAloneERV.ExhaustAirFanName = Alphas(5);
        GlobalNames::IntraObjUniquenessCheck(
            state, Alphas(5), CurrentModuleObject, cAlphaFields(5), state.dataHVACStandAloneERV->ExhaustAirFanUniqueNames, ErrorsFound);
        errFlag = false;
        if (HVACFan::checkIfFanNameIsAFanSystem(state,
                                                thisStandAloneERV
                                                    .ExhaustAirFanName)) { // no object type in input, so check if Fan:SystemModel
            thisStandAloneERV.ExhaustAirFanType_Num = DataHVACGlobals::FanType_SystemModelObject;
            state.dataHVACFan->fanObjs.emplace_back(new HVACFan::FanSystem(state, thisStandAloneERV.ExhaustAirFanName)); // call constructor

            thisStandAloneERV.ExhaustAirFanIndex = HVACFan::getFanObjectVectorIndex(state, thisStandAloneERV.ExhaustAirFanName);
            thisStandAloneERV.ExhaustAirFanSchPtr = state.dataHVACFan->fanObjs[thisStandAloneERV.ExhaustAirFanIndex]->availSchedIndex;
            thisStandAloneERV.DesignEAFanVolFlowRate = state.dataHVACFan->fanObjs[thisStandAloneERV.ExhaustAirFanIndex]->designAirVolFlowRate;
            thisStandAloneERV.ExhaustAirOutletNode = state.dataHVACFan->fanObjs[thisStandAloneERV.ExhaustAirFanIndex]->outletNodeNum;

        } else {
            GetFanType(state,
                       thisStandAloneERV.ExhaustAirFanName,
                       EAFanTypeNum,
                       errFlag,
                       CurrentModuleObject,
                       thisStandAloneERV.Name);
            if (!errFlag) {
                thisStandAloneERV.ExhaustAirFanType_Num = EAFanTypeNum;
                // error for fan availability schedule?
                thisStandAloneERV.ExhaustAirFanSchPtr = GetFanAvailSchPtr(
                    state, cFanTypes(EAFanTypeNum), thisStandAloneERV.ExhaustAirFanName, errFlag);
                GetFanIndex(state,
                            thisStandAloneERV.ExhaustAirFanName,
                            thisStandAloneERV.ExhaustAirFanIndex,
                            errFlag,
                            CurrentModuleObject + " \"" + thisStandAloneERV.Name + "\"");
            } else {
                ErrorsFound = true;
            }

            // Set the EA Design Fan Volume Flow Rate
            // get from fan module
            errFlag = false;
            EAFanVolFlowRate = GetFanDesignVolumeFlowRate(
                state, cFanTypes(EAFanTypeNum), thisStandAloneERV.ExhaustAirFanName, errFlag);
            if (errFlag) {
                ShowContinueError(
                    state, format("... occurs in {} ={}", CurrentModuleObject, thisStandAloneERV.Name));
                ErrorsFound = true;
            }
            thisStandAloneERV.DesignEAFanVolFlowRate = EAFanVolFlowRate;

            thisStandAloneERV.ExhaustAirOutletNode = GetFanOutletNode(
                state, cFanTypes(EAFanTypeNum), thisStandAloneERV.ExhaustAirFanName, errFlag);
            if (errFlag) {
                ShowContinueError(
                    state, format("... occurs in {} ={}", CurrentModuleObject, thisStandAloneERV.Name));
                ErrorsFound = true;
            }
        }

        errFlag = false;
        thisStandAloneERV.SupplyAirInletNode = HeatRecovery::GetSupplyInletNode(state, thisStandAloneERV.HeatExchangerName, errFlag);
        thisStandAloneERV.ExhaustAirInletNode = HeatRecovery::GetSecondaryInletNode(state, thisStandAloneERV.HeatExchangerName, errFlag);
        if (errFlag) {
            ShowContinueError(state,
                              format("... occurs in {} ={}", CurrentModuleObject, thisStandAloneERV.Name));
            ErrorsFound = true;
        }
        thisStandAloneERV.SupplyAirInletNode =
            GetOnlySingleNode(state,
                              state.dataLoopNodes->NodeID(thisStandAloneERV.SupplyAirInletNode),
                              ErrorsFound,
                              DataLoopNode::ConnectionObjectType::ZoneHVACEnergyRecoveryVentilator,
                              Alphas(1),
                              DataLoopNode::NodeFluidType::Air,
                              DataLoopNode::ConnectionType::Inlet,
                              NodeInputManager::CompFluidStream::Primary,
                              ObjectIsParent);
        thisStandAloneERV.SupplyAirOutletNode =
            GetOnlySingleNode(state,
                              state.dataLoopNodes->NodeID(thisStandAloneERV.SupplyAirOutletNode),
                              ErrorsFound,
                              DataLoopNode::ConnectionObjectType::ZoneHVACEnergyRecoveryVentilator,
                              Alphas(1),
                              DataLoopNode::NodeFluidType::Air,
                              DataLoopNode::ConnectionType::Outlet,
                              NodeInputManager::CompFluidStream::Primary,
                              ObjectIsParent);
        thisStandAloneERV.ExhaustAirInletNode =
            GetOnlySingleNode(state,
                              state.dataLoopNodes->NodeID(thisStandAloneERV.ExhaustAirInletNode),
                              ErrorsFound,
                              DataLoopNode::ConnectionObjectType::ZoneHVACEnergyRecoveryVentilator,
                              Alphas(1),
                              DataLoopNode::NodeFluidType::Air,
                              DataLoopNode::ConnectionType::Inlet,
                              NodeInputManager::CompFluidStream::Secondary,
                              ObjectIsParent);
        thisStandAloneERV.ExhaustAirOutletNode =
            GetOnlySingleNode(state,
                              state.dataLoopNodes->NodeID(thisStandAloneERV.ExhaustAirOutletNode),
                              ErrorsFound,
                              DataLoopNode::ConnectionObjectType::ZoneHVACEnergyRecoveryVentilator,
                              Alphas(1),
                              DataLoopNode::NodeFluidType::Air,
                              DataLoopNode::ConnectionType::ReliefAir,
                              NodeInputManager::CompFluidStream::Secondary,
                              ObjectIsParent);

        //   Check that supply air inlet node is an OA node
        if (!CheckOutAirNodeNumber(state, thisStandAloneERV.SupplyAirInletNode)) {
            ShowSevereError(state, format("For {} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ShowContinueError(state,
                              format(" Node name of supply air inlet node not valid Outdoor Air Node = {}",
                                     state.dataLoopNodes->NodeID(thisStandAloneERV.SupplyAirInletNode)));
            ShowContinueError(state, "...does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node.");
            ErrorsFound = true;
        }

        //   Check to make sure inlet and exhaust nodes are listed in a ZoneHVAC:EquipmentConnections object
        ZoneInletNodeFound = false;
        ZoneExhaustNodeFound = false;
        for (int ControlledZoneNum = 1; ControlledZoneNum <= state.dataGlobal->NumOfZones; ++ControlledZoneNum) {
            if (!ZoneInletNodeFound) {
                for (int NodeNumber = 1; NodeNumber <= state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).NumInletNodes; ++NodeNumber) {
                    if (state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).InletNode(NodeNumber) ==
                        thisStandAloneERV.SupplyAirOutletNode) {
                        ZoneInletNodeFound = true;
                        ZoneInletCZN = ControlledZoneNum;
                        break; // found zone inlet node
                    }
                }
            }
            if (!ZoneExhaustNodeFound) {
                for (int NodeNumber = 1; NodeNumber <= state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).NumExhaustNodes; ++NodeNumber) {
                    if (state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).ExhaustNode(NodeNumber) ==
                        thisStandAloneERV.ExhaustAirInletNode) {
                        ZoneExhaustNodeFound = true;
                        ZoneExhaustCZN = ControlledZoneNum;
                        break; // found zone exhaust node
                    }
                }
            }
        }
        if (!ZoneInletNodeFound) {
            ShowSevereError(state, format("For {} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ShowContinueError(state, "... Node name of supply air outlet node does not appear in a ZoneHVAC:EquipmentConnections object.");
            ShowContinueError(state,
                              format("... Supply air outlet node = {}",
                                     state.dataLoopNodes->NodeID(thisStandAloneERV.SupplyAirOutletNode)));
            ErrorsFound = true;
        }
        if (!ZoneExhaustNodeFound) {
            ShowSevereError(state, format("For {} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ShowContinueError(state, "... Node name of exhaust air inlet node does not appear in a ZoneHVAC:EquipmentConnections object.");
            ShowContinueError(state,
                              format("... Exhaust air inlet node = {}",
                                     state.dataLoopNodes->NodeID(thisStandAloneERV.ExhaustAirInletNode)));
            ErrorsFound = true;
        }
        //   If nodes are found, make sure they are in the same zone
        if (ZoneInletNodeFound && ZoneExhaustNodeFound) {
            if (ZoneInletCZN != ZoneExhaustCZN) {
                ShowSevereError(state,
                                format("For {} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
                ShowContinueError(state,
                                  "... Node name of supply air outlet node and exhasut air inlet node must appear in the same "
                                  "ZoneHVAC:EquipmentConnections object.");
                ShowContinueError(
                    state,
                    format("... Supply air outlet node = {}",
                           state.dataLoopNodes->NodeID(thisStandAloneERV.SupplyAirOutletNode)));
                ShowContinueError(
                    state, format("... ZoneHVAC:EquipmentConnections Zone Name = {}", state.dataZoneEquip->ZoneEquipConfig(ZoneInletCZN).ZoneName));
                ShowContinueError(
                    state,
                    format("... Exhaust air inlet node = {}",
                           state.dataLoopNodes->NodeID(thisStandAloneERV.ExhaustAirInletNode)));
                ShowContinueError(
                    state, format("... ZoneHVAC:EquipmentConnections Zone Name = {}", state.dataZoneEquip->ZoneEquipConfig(ZoneExhaustCZN).ZoneName));
                ErrorsFound = true;
            }
        }

        thisStandAloneERV.ControllerName = Alphas(6);
        // If controller name is blank the ERV unit will operate with no controller
        if (lAlphaBlanks(6)) {
            thisStandAloneERV.ControllerName = "xxxxx";
            thisStandAloneERV.ControllerNameDefined = false;
        } else {
            // Verify controller name in Stand Alone ERV object matches name of valid controller object
            GlobalNames::IntraObjUniquenessCheck(
                state, Alphas(6), CurrentModuleObject, cAlphaFields(6), state.dataHVACStandAloneERV->ControllerUniqueNames, ErrorsFound);
            thisStandAloneERV.ControllerNameDefined = true;
            if (ErrorsFound) {
                thisStandAloneERV.ControllerNameDefined = false;
            }

            if (state.dataInputProcessing->inputProcessor->getObjectItemNum(
                    state,
                    "ZoneHVAC:EnergyRecoveryVentilator:Controller",
                    thisStandAloneERV.ControllerName) <= 0) {
                ShowSevereError(
                    state, format("{} controller type ZoneHVAC:EnergyRecoveryVentilator:Controller not found = {}", CurrentModuleObject, Alphas(6)));
                ErrorsFound = true;
                thisStandAloneERV.ControllerNameDefined = false;
            }
        }

        if (!lAlphaBlanks(7)) {
            thisStandAloneERV.AvailManagerListName = Alphas(7);
        }

        // Read supply and exhaust air flow rates
        thisStandAloneERV.SupplyAirVolFlow = Numbers(1);
        thisStandAloneERV.ExhaustAirVolFlow = Numbers(2);

        // Read ventilation rate per floor area for autosizing HX and fans
        thisStandAloneERV.AirVolFlowPerFloorArea = Numbers(3);
        thisStandAloneERV.AirVolFlowPerOccupant = Numbers(4);

        if (thisStandAloneERV.SupplyAirVolFlow == AutoSize &&
            thisStandAloneERV.DesignSAFanVolFlowRate != AutoSize) {
            ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ShowContinueError(state,
                              format("... When autosizing ERV, supply air fan = {} \"{}\" must also be autosized.",
                                     cFanTypes(SAFanTypeNum),
                                     thisStandAloneERV.SupplyAirFanName));
        }

        if (thisStandAloneERV.ExhaustAirVolFlow == AutoSize &&
            thisStandAloneERV.DesignEAFanVolFlowRate != AutoSize) {
            ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ShowContinueError(state,
                              format("... When autosizing ERV, exhaust air fan = {} \"{}\" must also be autosized.",
                                     cFanTypes(EAFanTypeNum),
                                     thisStandAloneERV.ExhaustAirFanName));
        }

        if (thisStandAloneERV.SupplyAirVolFlow == AutoSize && HXSupAirFlowRate != AutoSize) {
            ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ShowContinueError(
                state,
                format("... When autosizing ERV {}, nominal supply air flow rate for heat exchanger with name = {} must also be autosized.",
                       cNumericFields(1),
                       thisStandAloneERV.HeatExchangerName));
        }

        if (thisStandAloneERV.ExhaustAirVolFlow == AutoSize && HXSupAirFlowRate != AutoSize) {
            ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
            ShowContinueError(
                state,
                format("... When autosizing ERV {}, nominal supply air flow rate for heat exchanger with name = {} must also be autosized.",
                       cNumericFields(2),
                       thisStandAloneERV.HeatExchangerName));
        }

        // Compare the ERV SA flow rates to SA fan object.
        if (thisStandAloneERV.DesignSAFanVolFlowRate != AutoSize &&
            thisStandAloneERV.SupplyAirVolFlow != AutoSize) {
            if (thisStandAloneERV.SupplyAirVolFlow >
                thisStandAloneERV.DesignSAFanVolFlowRate) {
                ShowWarningError(state,
                                 format("{} = {} has a {} > Max Volume Flow Rate defined in the associated fan object, should be <=",
                                        CurrentModuleObject,
                                        thisStandAloneERV.Name,
                                        cNumericFields(1)));
                ShowContinueError(state,
                                  format("... Entered value={:.2R}... Fan [{} \"{}\"] Max Value = {:.2R}",
                                         thisStandAloneERV.SupplyAirVolFlow,
                                         cFanTypes(SAFanTypeNum),
                                         thisStandAloneERV.SupplyAirFanName,
                                         thisStandAloneERV.DesignSAFanVolFlowRate));
                ShowContinueError(state,
                                  format(" The ERV {} is reset to the supply air fan flow rate and the simulation continues.", cNumericFields(1)));
                thisStandAloneERV.SupplyAirVolFlow =
                    thisStandAloneERV.DesignSAFanVolFlowRate;
            }
        }
        if (thisStandAloneERV.SupplyAirVolFlow != AutoSize) {
            if (thisStandAloneERV.SupplyAirVolFlow <= 0.0) {
                ShowSevereError(state,
                                format("{} = {} has a {} <= 0.0, it must be >0.0",
                                       CurrentModuleObject,
                                       thisStandAloneERV.Name,
                                       cNumericFields(1)));
                ShowContinueError(state,
                                  format("... Entered value={:.2R}", thisStandAloneERV.SupplyAirVolFlow));
                ErrorsFound = true;
            }
        } else {
            if (thisStandAloneERV.AirVolFlowPerFloorArea == 0.0 &&
                thisStandAloneERV.AirVolFlowPerOccupant == 0.0) {
                ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
                ShowContinueError(
                    state,
                    format("... Autosizing {} requires at least one input for {} or {}.", cNumericFields(1), cNumericFields(3), cNumericFields(4)));
                ErrorsFound = true;
            }
            // both inputs must be autosized
            if (thisStandAloneERV.ExhaustAirVolFlow != AutoSize) {
                ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
                ShowContinueError(state, format("... When autosizing, {} and {} must both be autosized.", cNumericFields(1), cNumericFields(2)));
                ErrorsFound = true;
            }
        }

        // Compare the ERV EA flow rates to EA fan object.
        if (thisStandAloneERV.DesignEAFanVolFlowRate != AutoSize &&
            thisStandAloneERV.ExhaustAirVolFlow != AutoSize) {
            if (thisStandAloneERV.ExhaustAirVolFlow >
                thisStandAloneERV.DesignEAFanVolFlowRate) {
                ShowWarningError(state,
                                 format("{} = {} has an {} > Max Volume Flow Rate defined in the associated fan object, should be <=",
                                        CurrentModuleObject,
                                        thisStandAloneERV.Name,
                                        cNumericFields(2)));
                ShowContinueError(state,
                                  format("... Entered value={:.2R}... Fan [{}:{}] Max Value = {:.2R}",
                                         thisStandAloneERV.ExhaustAirVolFlow,
                                         cFanTypes(EAFanTypeNum),
                                         thisStandAloneERV.ExhaustAirFanName,
                                         thisStandAloneERV.DesignEAFanVolFlowRate));
                ShowContinueError(state,
                                  format(" The ERV {} is reset to the exhaust air fan flow rate and the simulation continues.", cNumericFields(2)));
                thisStandAloneERV.ExhaustAirVolFlow =
                    thisStandAloneERV.DesignEAFanVolFlowRate;
            }
        }
        if (thisStandAloneERV.ExhaustAirVolFlow != AutoSize) {
            if (thisStandAloneERV.ExhaustAirVolFlow <= 0.0) {
                ShowSevereError(state,
                                format("{} = {} has an {} <= 0.0, it must be >0.0",
                                       CurrentModuleObject,
                                       thisStandAloneERV.Name,
                                       cNumericFields(2)));
                ShowContinueError(state,
                                  format("... Entered value={:.2R}", thisStandAloneERV.ExhaustAirVolFlow));
                ErrorsFound = true;
            }
        } else {
            if (thisStandAloneERV.AirVolFlowPerFloorArea == 0.0 &&
                thisStandAloneERV.AirVolFlowPerOccupant == 0.0) {
                ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
                ShowContinueError(
                    state,
                    format("... Autosizing {} requires at least one input for {} or {}.", cNumericFields(2), cNumericFields(3), cNumericFields(4)));
                ErrorsFound = true;
            }
            if (thisStandAloneERV.SupplyAirVolFlow != AutoSize) {
                ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, thisStandAloneERV.Name));
                ShowContinueError(state, format("... When autosizing, {} and {} must both be autosized.", cNumericFields(1), cNumericFields(2)));
                ErrorsFound = true;
            }
        }

        // Add supply fan to component sets array
        CompSetSupplyFanInlet = "UNDEFINED";
        CompSetSupplyFanOutlet = state.dataLoopNodes->NodeID(thisStandAloneERV.SupplyAirOutletNode);

        // Add exhaust fan to component sets array
        CompSetExhaustFanInlet = "UNDEFINED";
        CompSetExhaustFanOutlet = state.dataLoopNodes->NodeID(thisStandAloneERV.ExhaustAirOutletNode);

        // Add HX to component sets array
        SetUpCompSets(state,
                      thisStandAloneERV.UnitType,
                      thisStandAloneERV.Name,
                      "UNDEFINED",
                      thisStandAloneERV.HeatExchangerName,
                      "UNDEFINED",
                      "UNDEFINED");

        // Add supply fan to component sets array
        SetUpCompSets(state,
                      thisStandAloneERV.UnitType,
                      thisStandAloneERV.Name,
                      "UNDEFINED",
                      thisStandAloneERV.SupplyAirFanName,
                      CompSetSupplyFanInlet,
                      CompSetSupplyFanOutlet);

        // Add exhaust fan to component sets array
        SetUpCompSets(state,
                      thisStandAloneERV.UnitType,
                      thisStandAloneERV.Name,
                      "UNDEFINED",
                      thisStandAloneERV.ExhaustAirFanName,
                      CompSetExhaustFanInlet,
                      CompSetExhaustFanOutlet);

        // Verify HX name in Stand Alone ERV object matches name of valid HX object
        if (state.dataInputProcessing->inputProcessor->getObjectItemNum(
                state, "HeatExchanger:AirToAir:SensibleAndLatent", thisStandAloneERV.HeatExchangerName) <=
            0) {
            ShowSevereError(state,
                            format("{} heat exchanger type HeatExchanger:AirToAir:SensibleAndLatent not found = {}",
                                   CurrentModuleObject,
                                   thisStandAloneERV.HeatExchangerName));
            ErrorsFound = true;
        }
        // Verify supply air fan name in Stand Alone ERV object matches name of valid fan object
        if (thisStandAloneERV.SupplyAirFanType_Num != DataHVACGlobals::FanType_SystemModelObject) {
            if (state.dataInputProcessing->inputProcessor->getObjectItemNum(
                    state, "Fan:OnOff", thisStandAloneERV.SupplyAirFanName) <= 0) {
                ShowSevereError(state,
                                format("{} supply fan type Fan:OnOff not found = {}",
                                       CurrentModuleObject,
                                       thisStandAloneERV.SupplyAirFanName));
                ErrorsFound = true;
            }
        } else {
            if (state.dataInputProcessing->inputProcessor->getObjectItemNum(
                    state, "Fan:SystemModel", thisStandAloneERV.SupplyAirFanName) <= 0) {
                ShowSevereError(state,
                                format("{} supply fan type Fan:SystemModel not found = {}",
                                       CurrentModuleObject,
                                       thisStandAloneERV.SupplyAirFanName));
                ErrorsFound = true;
            }
        }

        // Verify exhaust air fan name in Stand Alone ERV object matches name of valid fan object
        if (thisStandAloneERV.ExhaustAirFanType_Num != DataHVACGlobals::FanType_SystemModelObject) {
            if (state.dataInputProcessing->inputProcessor->getObjectItemNum(
                    state, "Fan:OnOff", thisStandAloneERV.ExhaustAirFanName) <= 0) {
                ShowSevereError(state,
                                format("{} exhaust fan type Fan:OnOff not found = {}",
                                       CurrentModuleObject,
                                       thisStandAloneERV.ExhaustAirFanName));
                ErrorsFound = true;
            }
        } else {
            if (state.dataInputProcessing->inputProcessor->getObjectItemNum(
                    state, "Fan:SystemModel", thisStandAloneERV.ExhaustAirFanName) <= 0) {
                ShowSevereError(state,
                                format("{} exhaust fan type Fan:SystemModel not found = {}",
                                       CurrentModuleObject,
                                       thisStandAloneERV.ExhaustAirFanName));
                ErrorsFound = true;
            }
        }
    }

    int OutAirNum = 0;
    CurrentModuleObject = "ZoneHVAC:EnergyRecoveryVentilator:Controller";
    int NumERVCtrlrs = state.dataInputProcessing->inputProcessor->getNumObjectsFound(state, CurrentModuleObject);

    for (int ERVControllerNum = 1; ERVControllerNum <= NumERVCtrlrs; ++ERVControllerNum) {
        state.dataInputProcessing->inputProcessor->getObjectItem(state,
                                                                 CurrentModuleObject,
                                                                 ERVControllerNum,
                                                                 Alphas,
                                                                 NumAlphas,
                                                                 Numbers,
                                                                 NumNumbers,
                                                                 IOStatus,
                                                                 lNumericBlanks,
                                                                 lAlphaBlanks,
                                                                 cAlphaFields,
                                                                 cNumericFields);
        MixedAir::CheckOAControllerName(state, Alphas(1), CurrentModuleObject, cAlphaFields(1), ErrorsFound);
        ++OutAirNum;
        auto &thisOAController(state.dataMixedAir->OAController(OutAirNum));

        thisOAController.Name = Alphas(1);
        thisOAController.ControllerType = CurrentModuleObject;
        thisOAController.ControllerType_Num = MixedAir::MixedAirControllerType::ControllerStandAloneERV;
        int WhichERV = UtilityRoutines::FindItemInList(thisOAController.Name, state.dataHVACStandAloneERV->StandAloneERV, &StandAloneERVData::ControllerName);
        if (WhichERV != 0) {
            AirFlowRate = state.dataHVACStandAloneERV->StandAloneERV(WhichERV).SupplyAirVolFlow;
            state.dataHVACStandAloneERV->StandAloneERV(WhichERV).ControllerIndex = OutAirNum;
        } else {
            ShowSevereError(
                state, format("GetERVController: Could not find ZoneHVAC:EnergyRecoveryVentilator with {} = \"{}\"", cAlphaFields(1), Alphas(1)));
            ErrorsFound = true;
            AirFlowRate = -1000.0;
        }
        thisOAController.MaxOA = AirFlowRate;
        thisOAController.MinOA = AirFlowRate;
        //    OAController(OutAirNum)%TempLim = Numbers(1)
        if (lNumericBlanks(1)) {
            thisOAController.TempLim = BlankNumeric;
        } else {
            thisOAController.TempLim = Numbers(1);
        }
        //    OAController(OutAirNum)%TempLowLim = Numbers(2)
        if (lNumericBlanks(2)) {
            thisOAController.TempLowLim = BlankNumeric;
        } else {
            thisOAController.TempLowLim = Numbers(2);
        }
        //    OAController(OutAirNum)%EnthLim = Numbers(3)
        if (lNumericBlanks(3)) {
            thisOAController.EnthLim = BlankNumeric;
        } else {
            thisOAController.EnthLim = Numbers(3);
        }
        //    OAController(OutAirNum)%DPTempLim = Numbers(4)
        if (lNumericBlanks(4)) {
            thisOAController.DPTempLim = BlankNumeric;
        } else {
            thisOAController.DPTempLim = Numbers(4);
        }

        if (WhichERV != 0) {
            NodeNumber = state.dataHVACStandAloneERV->StandAloneERV(WhichERV).SupplyAirInletNode;
        } else {
            NodeNumber = 0;
        }
        thisOAController.OANode = NodeNumber;
        // set the inlet node to also equal the OA node because this is a special controller for economizing stand alone ERV
        // with the assumption that equipment is bypassed....(moved from module MixedAir)
        thisOAController.InletNode = NodeNumber;

        if (WhichERV != 0) {
            NodeNumber = state.dataHVACStandAloneERV->StandAloneERV(WhichERV).ExhaustAirInletNode;
        } else {
            NodeNumber = 0;
        }
        thisOAController.RetNode = NodeNumber;

        if (!lAlphaBlanks(2)) {
            thisOAController.EnthalpyCurvePtr = GetCurveIndex(state, Alphas(2));
            if (GetCurveIndex(state, Alphas(2)) == 0) {
                ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                ShowContinueError(state, format("...{} not found:{}", cAlphaFields(2), Alphas(2)));
                ErrorsFound = true;
            } else {
                // Verify Curve Object, only legal types are Quadratic and Cubic
                ErrorsFound |= Curve::CheckCurveDims(state,
                                                     thisOAController.EnthalpyCurvePtr, // Curve index
                                                     {1},                               // Valid dimensions
                                                     "GetStandAloneERV: ",              // Routine name
                                                     CurrentModuleObject,               // Object Type
                                                     thisOAController.Name,             // Object Name
                                                     cAlphaFields(2));                  // Field Name
            }
        }

        // Changed by AMIT for new implementation of the controller:outside air
        if (Alphas(3) == "EXHAUSTAIRTEMPERATURELIMIT" && Alphas(4) == "EXHAUSTAIRENTHALPYLIMIT") {
            thisOAController.Econo = MixedAir::EconoOp::DifferentialDryBulbAndEnthalpy;
        } else if (Alphas(3) == "EXHAUSTAIRTEMPERATURELIMIT" && Alphas(4) == "NOEXHAUSTAIRENTHALPYLIMIT") {
            thisOAController.Econo = MixedAir::EconoOp::DifferentialDryBulb;
        } else if (Alphas(3) == "NOEXHAUSTAIRTEMPERATURELIMIT" && Alphas(4) == "EXHAUSTAIRENTHALPYLIMIT") {
            thisOAController.Econo = MixedAir::EconoOp::DifferentialEnthalpy;
        } else if (Alphas(3) == "NOEXHAUSTAIRTEMPERATURELIMIT" && Alphas(4) == "NOEXHAUSTAIRENTHALPYLIMIT") {
            if ((!lNumericBlanks(1)) || (!lNumericBlanks(3)) || (!lNumericBlanks(4)) || (!lAlphaBlanks(2))) {
                // This means that any of the FIXED DRY BULB, FIXED ENTHALPY, FIXED DEW POINT AND DRY BULB OR
                // ELECTRONIC ENTHALPY ECONOMIZER STRATEGY is present
                thisOAController.Econo = MixedAir::EconoOp::FixedDryBulb;
            }
        } else if ((!lAlphaBlanks(3)) && (!lAlphaBlanks(4))) {
            if ((lNumericBlanks(1)) && (lNumericBlanks(3)) && (lNumericBlanks(4)) && lAlphaBlanks(2)) {
                ShowWarningError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                ShowContinueError(state, format("... Invalid {}{} = {}{}", cAlphaFields(3), cAlphaFields(4), Alphas(3), Alphas(4)));
                ShowContinueError(state, "... Assumed NO EXHAUST AIR TEMP LIMIT and NO EXHAUST AIR ENTHALPY LIMIT.");
                thisOAController.Econo = MixedAir::EconoOp::NoEconomizer;
            } else {
                // This means that any of the FIXED DRY BULB, FIXED ENTHALPY, FIXED DEW POINT AND DRY BULB OR
                // ELECTRONIC ENTHALPY ECONOMIZER STRATEGY is present
                thisOAController.Econo = MixedAir::EconoOp::FixedDryBulb;
            }
        } else if ((lAlphaBlanks(3)) && (!lAlphaBlanks(4))) {
            if ((lNumericBlanks(1)) && (lNumericBlanks(3)) && (lNumericBlanks(4)) && lAlphaBlanks(2)) {
                ShowWarningError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                ShowContinueError(state, format("... Invalid {} = {}", cAlphaFields(4), Alphas(4)));
                ShowContinueError(state, "... Assumed  NO EXHAUST AIR ENTHALPY LIMIT.");
                thisOAController.Econo = MixedAir::EconoOp::NoEconomizer;
            } else {
                // This means that any of the FIXED DRY BULB, FIXED ENTHALPY, FIXED DEW POINT AND DRY BULB OR
                // ELECTRONIC ENTHALPY ECONOMIZER STRATEGY is present
                thisOAController.Econo = MixedAir::EconoOp::FixedDryBulb;
            }
        } else if ((!lAlphaBlanks(3)) && (lAlphaBlanks(4))) {
            if ((lNumericBlanks(1)) && (lNumericBlanks(3)) && (lNumericBlanks(4)) && lAlphaBlanks(2)) {
                ShowWarningError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                ShowContinueError(state, format("... Invalid {} = {}", cAlphaFields(3), Alphas(3)));
                ShowContinueError(state, "... Assumed NO EXHAUST AIR TEMP LIMIT ");
                thisOAController.Econo = MixedAir::EconoOp::NoEconomizer;
            } else {
                // This means that any of the FIXED DRY BULB, FIXED ENTHALPY, FIXED DEW POINT AND DRY BULB OR
                // ELECTRONIC ENTHALPY ECONOMIZER STRATEGY is present
                thisOAController.Econo = MixedAir::EconoOp::FixedDryBulb;
            }
        } else { // NO Economizer
            thisOAController.Econo = MixedAir::EconoOp::NoEconomizer;
        }

        thisOAController.FixedMin = false;
        thisOAController.EconBypass = true;

        //   Initialize to one in case high humidity control is NOT used
        HighRHOARatio = 1.0;
        //   READ Modify Air Flow Data
        //   High humidity control option is YES, read in additional data
        if (UtilityRoutines::SameString(Alphas(6), "Yes")) {

            HStatZoneNum = UtilityRoutines::FindItemInList(Alphas(7), state.dataHeatBal->Zone);
            thisOAController.HumidistatZoneNum = HStatZoneNum;

            // Get the node number for the zone with the humidistat
            if (HStatZoneNum > 0) {
                ZoneNodeFound = false;
                HStatFound = false;
                if (state.dataZoneEquip->ZoneEquipConfig(HStatZoneNum).IsControlled) {
                    //         Find the controlled zone number for the specified humidistat location
                    thisOAController.NodeNumofHumidistatZone = state.dataZoneEquip->ZoneEquipConfig(HStatZoneNum).ZoneNode;
                    ZoneNodeFound = true;
                }
                if (!ZoneNodeFound) {
                    ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                    ShowContinueError(state, "... Did not find Air Node (Zone with Humidistat)");
                    ShowContinueError(state, format("... Specified {} = {}", cAlphaFields(7), Alphas(7)));
                    ShowContinueError(state, "... A ZoneHVAC:EquipmentConnections object must be specified for this zone.");
                    ErrorsFound = true;
                } else {
                    for (NumHstatZone = 1; NumHstatZone <= state.dataZoneCtrls->NumHumidityControlZones; ++NumHstatZone) {
                        if (state.dataZoneCtrls->HumidityControlZone(NumHstatZone).ActualZoneNum != HStatZoneNum) continue;
                        HStatFound = true;
                        break;
                    }
                    if (!HStatFound) {
                        ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                        ShowContinueError(state, "... Did not find zone humidistat");
                        ShowContinueError(state, "... A ZoneControl:Humidistat object must be specified for this zone.");
                        ErrorsFound = true;
                    }
                }
            } else {
                ShowSevereError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                ShowContinueError(state, "... Did not find Air Node (Zone with Humidistat)");
                ShowContinueError(state, "... A ZoneHVAC:EquipmentConnections object must be specified for this zone.");
                ErrorsFound = true;
            }

            if (Numbers(5) <= 0.0 && NumNumbers > 4) {

                ShowWarningError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                ShowContinueError(state, format("... {} must be greater than 0.", cNumericFields(5)));
                ShowContinueError(state, format("... {} is reset to 1 and the simulation continues.", cNumericFields(5)));

                HighRHOARatio = 1.0;

            } else if (NumNumbers > 4) {

                HighRHOARatio = Numbers(5);

            } else {

                HighRHOARatio = 1.0;
            }

            if (UtilityRoutines::SameString(Alphas(8), "Yes")) {
                thisOAController.ModifyDuringHighOAMoisture = false;
            } else {
                thisOAController.ModifyDuringHighOAMoisture = true;
            }

        } else if (!UtilityRoutines::SameString(Alphas(6), "No") && NumAlphas > 4 && (!lAlphaBlanks(5))) {
            ShowWarningError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
            ShowContinueError(state, format("... Invalid {} = {}", cAlphaFields(6), Alphas(6)));
            ShowContinueError(state, format("... {} is assumed to be \"No\" and the simulation continues.", cAlphaFields(6)));
        } // IF(UtilityRoutines::SameString(Alphas(6),'Yes'))THEN

        thisOAController.HighRHOAFlowRatio = HighRHOARatio;
        if (WhichERV != 0) {
            state.dataHVACStandAloneERV->StandAloneERV(WhichERV).HighRHOAFlowRatio = HighRHOARatio;
        }

        //   Check for a time of day outside air schedule
        thisOAController.EconomizerOASchedPtr = GetScheduleIndex(state, Alphas(5));

        if (WhichERV != 0) {
            auto &thisStandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(WhichERV);
	    thisStandAloneERV.EconomizerOASchedPtr = GetScheduleIndex(state, Alphas(5));

            // Compare the ERV SA fan flow rates to modified air flow rate.
            if (HighRHOARatio > 1.0 && thisStandAloneERV.SupplyAirVolFlow != AutoSize &&
                thisStandAloneERV.DesignSAFanVolFlowRate != AutoSize) {
                if (thisStandAloneERV.SupplyAirVolFlow * HighRHOARatio >
                    thisStandAloneERV.DesignSAFanVolFlowRate) {
                    ShowWarningError(state, format("{} \"{}\"", CurrentModuleObject, Alphas(1)));
                    ShowContinueError(state, format("... A {} was entered as {:.4R}", cNumericFields(5), HighRHOARatio));
                    ShowContinueError(state,
                                      "... This flow ratio results in a Supply Air Volume Flow Rate through the ERV which is greater than the "
                                      "Max Volume specified in the supply air fan object.");
                    ShowContinueError(state,
                                      format("... Associated fan object = {} \"{}\"",
                                             cFanTypes(SAFanTypeNum),
                                             thisStandAloneERV.SupplyAirFanName));
                    ShowContinueError(state,
                                      format("... Modified value                   = {:.2R}",
                                             thisStandAloneERV.SupplyAirVolFlow * HighRHOARatio));
                    ShowContinueError(state,
                                      format(" ... Supply Fan Max Volume Flow Rate = {:.2R}",
                                             thisStandAloneERV.DesignSAFanVolFlowRate));
                    ShowContinueError(state, "... The ERV supply air fan will limit the air flow through the ERV and the simulation continues.");
                }
            }

            // Compare the ERV EA fan flow rates to modified air flow rate.
            if (HighRHOARatio > 1.0 && thisStandAloneERV.ExhaustAirVolFlow != AutoSize &&
                thisStandAloneERV.DesignEAFanVolFlowRate != AutoSize) {
                if (thisStandAloneERV.ExhaustAirVolFlow * HighRHOARatio >
                    thisStandAloneERV.DesignEAFanVolFlowRate) {
                    ShowWarningError(state, format("ZoneHVAC:EnergyRecoveryVentilator:Controller \"{}\"", Alphas(1)));
                    ShowContinueError(state, format("... A {} was entered as {:.4R}", cNumericFields(5), HighRHOARatio));
                    ShowContinueError(state,
                                      "... This flow ratio results in an Exhaust Air Volume Flow Rate through the ERV which is greater than the "
                                      "Max Volume specified in the exhaust air fan object.");
                    ShowContinueError(state,
                                      format("... Associated fan object = {} \"{}\"",
                                             cFanTypes(EAFanTypeNum),
                                             thisStandAloneERV.ExhaustAirFanName));
                    ShowContinueError(state,
                                      format("... Modified value                    = {:.2R}",
                                             thisStandAloneERV.ExhaustAirVolFlow * HighRHOARatio));
                    ShowContinueError(state,
                                      format(" ... Exhaust Fan Max Volume Flow Rate = {:.2R}",
                                             thisStandAloneERV.DesignEAFanVolFlowRate));
                    ShowContinueError(state, "... The ERV exhaust air fan will limit the air flow through the ERV and the simulation continues.");
                }
            }
        } // IF(WhichERV /= 0)THEN
    }

    if (ErrorsFound) {
        ShowFatalError(state, "Errors found in getting ZoneHVAC:EnergyRecoveryVentilator input.");
    }

    // Setup report variables for the stand alone ERVs
    for (int StandAloneERVIndex = 1; StandAloneERVIndex <= state.dataHVACStandAloneERV->NumStandAloneERVs; ++StandAloneERVIndex) {
	auto &thisStandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVIndex);
        SetupOutputVariable(state,
                            "Zone Ventilator Sensible Cooling Rate",
                            OutputProcessor::Unit::W,
                            thisStandAloneERV.SensCoolingRate,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Sensible Cooling Energy",
                            OutputProcessor::Unit::J,
                            thisStandAloneERV.SensCoolingEnergy,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Summed,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Latent Cooling Rate",
                            OutputProcessor::Unit::W,
                            thisStandAloneERV.LatCoolingRate,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Latent Cooling Energy",
                            OutputProcessor::Unit::J,
                            thisStandAloneERV.LatCoolingEnergy,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Summed,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Total Cooling Rate",
                            OutputProcessor::Unit::W,
                            thisStandAloneERV.TotCoolingRate,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Total Cooling Energy",
                            OutputProcessor::Unit::J,
                            thisStandAloneERV.TotCoolingEnergy,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Summed,
                            thisStandAloneERV.Name);

        SetupOutputVariable(state,
                            "Zone Ventilator Sensible Heating Rate",
                            OutputProcessor::Unit::W,
                            thisStandAloneERV.SensHeatingRate,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Sensible Heating Energy",
                            OutputProcessor::Unit::J,
                            thisStandAloneERV.SensHeatingEnergy,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Summed,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Latent Heating Rate",
                            OutputProcessor::Unit::W,
                            thisStandAloneERV.LatHeatingRate,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Latent Heating Energy",
                            OutputProcessor::Unit::J,
                            thisStandAloneERV.LatHeatingEnergy,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Summed,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Total Heating Rate",
                            OutputProcessor::Unit::W,
                            thisStandAloneERV.TotHeatingRate,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Total Heating Energy",
                            OutputProcessor::Unit::J,
                            thisStandAloneERV.TotHeatingEnergy,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Summed,
                            thisStandAloneERV.Name);

        SetupOutputVariable(state,
                            "Zone Ventilator Electricity Rate",
                            OutputProcessor::Unit::W,
                            thisStandAloneERV.ElecUseRate,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Electricity Energy",
                            OutputProcessor::Unit::J,
                            thisStandAloneERV.ElecUseEnergy,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Summed,
                            thisStandAloneERV.Name);
        SetupOutputVariable(state,
                            "Zone Ventilator Supply Fan Availability Status",
                            OutputProcessor::Unit::None,
                            thisStandAloneERV.AvailStatus,
                            OutputProcessor::SOVTimeStepType::System,
                            OutputProcessor::SOVStoreType::Average,
                            thisStandAloneERV.Name);
    }

    Alphas.deallocate();
    Numbers.deallocate();
    cAlphaFields.deallocate();
    cNumericFields.deallocate();
    lNumericBlanks.deallocate();
    lAlphaBlanks.deallocate();
}

void InitStandAloneERV(EnergyPlusData &state,
                       int const StandAloneERVNum,   // number of the current Stand Alone ERV unit being simulated
                       int const ZoneNum,            // number of zone being served unused1208
                       bool const FirstHVACIteration // TRUE if first HVAC iteration
)
{

    // SUBROUTINE INFORMATION:
    //       AUTHOR         Richard Raustad, FSEC
    //       DATE WRITTEN   June 2003
    //       MODIFIED       July 2012, Chandan Sharma - FSEC: Added zone sys avail managers
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine is for initializations of the Stand Alone ERV unit information.

    // METHODOLOGY EMPLOYED:
    // Uses the status flags to trigger initializations.

    using DataZoneEquipment::CheckZoneEquipmentList;
    using MixedAir::SimOAController;

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    int SupInNode;    // supply air inlet node number
    int ExhInNode;    // exhaust air inlet node number
    int SupInletNode; // supply air inlet node number for Stand Alone ERV 'StandAloneERVNum'

    auto &Node(state.dataLoopNodes->Node);

    // Do the one time initializations
    if (state.dataHVACStandAloneERV->MyOneTimeFlag) {

        state.dataHVACStandAloneERV->MyEnvrnFlag.allocate(state.dataHVACStandAloneERV->NumStandAloneERVs);
        state.dataHVACStandAloneERV->MySizeFlag_InitStandAloneERV.allocate(state.dataHVACStandAloneERV->NumStandAloneERVs);
        state.dataHVACStandAloneERV->MyZoneEqFlag.allocate(state.dataHVACStandAloneERV->NumStandAloneERVs);
        state.dataHVACStandAloneERV->MyEnvrnFlag = true;
        state.dataHVACStandAloneERV->MySizeFlag_InitStandAloneERV = true;
        state.dataHVACStandAloneERV->MyZoneEqFlag = true;
        state.dataHVACStandAloneERV->MyOneTimeFlag = false;
    }

    auto &thisStandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum);
    
    if (allocated(state.dataHVACGlobal->ZoneComp)) {
        if (state.dataHVACStandAloneERV->MyZoneEqFlag(StandAloneERVNum)) { // initialize the name of each availability manager list and zone number
            state.dataHVACGlobal->ZoneComp(DataZoneEquipment::ZoneEquip::ERVStandAlone).ZoneCompAvailMgrs(StandAloneERVNum).AvailManagerListName =
                thisStandAloneERV.AvailManagerListName;
            state.dataHVACGlobal->ZoneComp(DataZoneEquipment::ZoneEquip::ERVStandAlone).ZoneCompAvailMgrs(StandAloneERVNum).ZoneNum = ZoneNum;
            state.dataHVACStandAloneERV->MyZoneEqFlag(StandAloneERVNum) = false;
        }
        thisStandAloneERV.AvailStatus =
            state.dataHVACGlobal->ZoneComp(DataZoneEquipment::ZoneEquip::ERVStandAlone).ZoneCompAvailMgrs(StandAloneERVNum).AvailStatus;
    }

    // need to check all units to see if they are on Zone Equipment List or issue warning
    if (!state.dataHVACStandAloneERV->ZoneEquipmentListChecked && state.dataZoneEquip->ZoneEquipInputsFilled) {
        state.dataHVACStandAloneERV->ZoneEquipmentListChecked = true;
        for (int Loop = 1; Loop <= state.dataHVACStandAloneERV->NumStandAloneERVs; ++Loop) {
            if (CheckZoneEquipmentList(
                    state, state.dataHVACStandAloneERV->StandAloneERV(Loop).UnitType, state.dataHVACStandAloneERV->StandAloneERV(Loop).Name))
                continue;
            ShowSevereError(state,
                            format("InitStandAloneERV: Unit=[{},{}] is not on any ZoneHVAC:EquipmentList.  It will not be simulated.",
                                   state.dataHVACStandAloneERV->StandAloneERV(Loop).UnitType,
                                   state.dataHVACStandAloneERV->StandAloneERV(Loop).Name));
        }
    }

    if (!state.dataGlobal->SysSizingCalc && state.dataHVACStandAloneERV->MySizeFlag_InitStandAloneERV(StandAloneERVNum)) {
        SizeStandAloneERV(state, StandAloneERVNum);
        state.dataHVACStandAloneERV->MySizeFlag_InitStandAloneERV(StandAloneERVNum) = false;
    }

    // Do the Begin Environment initializations
    if (state.dataGlobal->BeginEnvrnFlag && state.dataHVACStandAloneERV->MyEnvrnFlag(StandAloneERVNum)) {
        SupInNode = thisStandAloneERV.SupplyAirInletNode;
        ExhInNode = thisStandAloneERV.ExhaustAirInletNode;
        // set the mass flow rates from the input volume flow rates
        thisStandAloneERV.MaxSupAirMassFlow =
            state.dataEnvrn->StdRhoAir * thisStandAloneERV.SupplyAirVolFlow;
        thisStandAloneERV.MaxExhAirMassFlow =
            state.dataEnvrn->StdRhoAir * thisStandAloneERV.ExhaustAirVolFlow;
        thisStandAloneERV.DesignSAFanMassFlowRate =
            state.dataEnvrn->StdRhoAir * thisStandAloneERV.DesignSAFanVolFlowRate;
        thisStandAloneERV.DesignEAFanMassFlowRate =
            state.dataEnvrn->StdRhoAir * thisStandAloneERV.DesignEAFanVolFlowRate;
        // set the node max and min mass flow rates
        Node(SupInNode).MassFlowRateMax = thisStandAloneERV.MaxSupAirMassFlow;
        Node(SupInNode).MassFlowRateMin = 0.0;
        Node(ExhInNode).MassFlowRateMax = thisStandAloneERV.MaxExhAirMassFlow;
        Node(ExhInNode).MassFlowRateMin = 0.0;
        state.dataHVACStandAloneERV->MyEnvrnFlag(StandAloneERVNum) = false;
        //   Initialize OA Controller on BeginEnvrnFlag
        if (thisStandAloneERV.ControllerNameDefined) {
            SimOAController(state,
                            thisStandAloneERV.ControllerName,
                            thisStandAloneERV.ControllerIndex,
                            FirstHVACIteration,
                            0);
        }
    } // end one time inits

    if (!state.dataGlobal->BeginEnvrnFlag) {
        state.dataHVACStandAloneERV->MyEnvrnFlag(StandAloneERVNum) = true;
    }

    // These initializations are done every iteration
    thisStandAloneERV.ElecUseRate = 0.0;
    thisStandAloneERV.SensCoolingRate = 0.0;
    thisStandAloneERV.LatCoolingRate = 0.0;
    thisStandAloneERV.TotCoolingRate = 0.0;
    thisStandAloneERV.SensHeatingRate = 0.0;
    thisStandAloneERV.LatHeatingRate = 0.0;
    thisStandAloneERV.TotHeatingRate = 0.0;
    SupInletNode = thisStandAloneERV.SupplyAirInletNode;
    ExhInNode = thisStandAloneERV.ExhaustAirInletNode;

    // Set the inlet node mass flow rate
    if (GetCurrentScheduleValue(state, thisStandAloneERV.SchedPtr) > 0.0) {

        //   IF optional ControllerName is defined SimOAController ONLY to set economizer and Modifyairflow flags
        if (thisStandAloneERV.ControllerNameDefined) {
            //     Initialize a flow rate for controller
            Node(SupInletNode).MassFlowRate = thisStandAloneERV.MaxSupAirMassFlow;
            SimOAController(state,
                            thisStandAloneERV.ControllerName,
                            thisStandAloneERV.ControllerIndex,
                            FirstHVACIteration,
                            0);
        }

        if (GetCurrentScheduleValue(state, thisStandAloneERV.SupplyAirFanSchPtr) > 0 ||
            (state.dataHVACGlobal->ZoneCompTurnFansOn && !state.dataHVACGlobal->ZoneCompTurnFansOff)) {
            if (thisStandAloneERV.ControllerNameDefined) {
                if (state.dataMixedAir->OAController(thisStandAloneERV.ControllerIndex)
                        .HighHumCtrlActive) {
                    Node(SupInletNode).MassFlowRate = min(thisStandAloneERV.DesignSAFanMassFlowRate,
                                                          thisStandAloneERV.MaxSupAirMassFlow *
                                                              thisStandAloneERV.HighRHOAFlowRatio);
                } else {
                    Node(SupInletNode).MassFlowRate = min(thisStandAloneERV.DesignSAFanMassFlowRate,
                                                          thisStandAloneERV.MaxSupAirMassFlow);
                }
            } else {
                Node(SupInletNode).MassFlowRate = min(thisStandAloneERV.DesignSAFanMassFlowRate,
                                                      thisStandAloneERV.MaxSupAirMassFlow);
            }
        } else {
            Node(SupInletNode).MassFlowRate = 0.0;
        }
        Node(SupInletNode).MassFlowRateMaxAvail = Node(SupInletNode).MassFlowRate;
        Node(SupInletNode).MassFlowRateMinAvail = Node(SupInletNode).MassFlowRate;

        if (GetCurrentScheduleValue(state, thisStandAloneERV.ExhaustAirFanSchPtr) > 0) {
            if (thisStandAloneERV.ControllerNameDefined) {
                if (state.dataMixedAir->OAController(thisStandAloneERV.ControllerIndex)
                        .HighHumCtrlActive) {
                    Node(ExhInNode).MassFlowRate = min(thisStandAloneERV.DesignEAFanMassFlowRate,
                                                       thisStandAloneERV.MaxExhAirMassFlow *
                                                           thisStandAloneERV.HighRHOAFlowRatio);
                } else {
                    Node(ExhInNode).MassFlowRate = min(thisStandAloneERV.DesignEAFanMassFlowRate,
                                                       thisStandAloneERV.MaxExhAirMassFlow);
                }
            } else {
                Node(ExhInNode).MassFlowRate = min(thisStandAloneERV.DesignEAFanMassFlowRate,
                                                   thisStandAloneERV.MaxExhAirMassFlow);
            }
        } else {
            Node(ExhInNode).MassFlowRate = 0.0;
        }
        Node(ExhInNode).MassFlowRateMaxAvail = Node(ExhInNode).MassFlowRate;
        Node(ExhInNode).MassFlowRateMinAvail = Node(ExhInNode).MassFlowRate;
    } else {
        Node(SupInletNode).MassFlowRate = 0.0;
        Node(SupInletNode).MassFlowRateMaxAvail = 0.0;
        Node(SupInletNode).MassFlowRateMinAvail = 0.0;
        Node(ExhInNode).MassFlowRate = 0.0;
        Node(ExhInNode).MassFlowRateMaxAvail = 0.0;
        Node(ExhInNode).MassFlowRateMinAvail = 0.0;
    }
}

void SizeStandAloneERV(EnergyPlusData &state, int const StandAloneERVNum)
{

    // SUBROUTINE INFORMATION:
    //       AUTHOR         Richard Raustad
    //       DATE WRITTEN   October 2007
    //       MODIFIED       August 2013 Daeho Kang, add component sizing table entries
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine is for sizing Stand Alone ERV Components for which flow rates have not been
    // specified in the input.

    // METHODOLOGY EMPLOYED:
    // Obtains flow rates from the zone or system sizing arrays.

    // Using/Aliasing
    using Fans::SetFanData;
    using Fans::SimulateFanComponents;
    using HeatRecovery::SetHeatExchangerData;
    using ScheduleManager::GetScheduleMaxValue;

    static constexpr std::string_view RoutineName("SizeStandAloneERV: ");

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    std::string ZoneName;              // Name of zone
    Real64 ZoneMult;                   // Zone multiplier
    int PeopleNum;                     // Index to people object
    Real64 NumberOfPeople;             // Maximum number of people in zone
    int PeopleSchPtr;                  // Pointer to people schedule
    Real64 MaxPeopleSch;               // maximum people schedule value
    Real64 FloorArea;                  // Floor area of zone (m2)
    bool IsAutoSize;                   // Indicator to autosize
    Real64 SupplyAirVolFlowDes;        // Autosized supply air flow for reporting
    Real64 SupplyAirVolFlowUser;       // Hardsized supply air flow for reporting
    Real64 DesignSAFanVolFlowRateDes;  // Autosized supply air fan flow for reporting
    Real64 DesignSAFanVolFlowRateUser; // Hardsized supply air fan flow for reporting
    Real64 ExhaustAirVolFlowDes;       // Autosized exhaust air flow for reporting
    Real64 ExhaustAirVolFlowUser;      // Hardsized exhaust air flow for reporting

    IsAutoSize = false;
    SupplyAirVolFlowDes = 0.0;
    SupplyAirVolFlowUser = 0.0;
    DesignSAFanVolFlowRateDes = 0.0;
    DesignSAFanVolFlowRateUser = 0.0;
    ExhaustAirVolFlowDes = 0.0;
    ExhaustAirVolFlowUser = 0.0;
    constexpr std::string_view CompType("ZoneHVAC:EnergyRecoveryVentilator");

    auto &thisStandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum);
    std::string CompName(thisStandAloneERV.Name);
    bool PrintFlag = true;
    bool ErrorsFound = false;

    auto &ZoneEqSizing(state.dataSize->ZoneEqSizing);

    if (thisStandAloneERV.SupplyAirVolFlow == DataSizing::AutoSize) {
        IsAutoSize = true;
    }

    if (state.dataSize->CurZoneEqNum > 0) {

        //      Sizing objects are not required for stand alone ERV
        //      CALL CheckZoneSizing('ZoneHVAC:EnergyRecoveryVentilator',StandAloneERV(StandAloneERVNum)%Name)
        ZoneName = state.dataZoneEquip->ZoneEquipConfig(state.dataSize->CurZoneEqNum).ZoneName;
        int ZoneNum = state.dataSize->CurZoneEqNum;
        ZoneMult = state.dataHeatBal->Zone(ZoneNum).Multiplier * state.dataHeatBal->Zone(ZoneNum).ListMultiplier;
        FloorArea = 0.0;
        FloorArea = state.dataHeatBal->Zone(ZoneNum).FloorArea;
        NumberOfPeople = 0.0;
        MaxPeopleSch = 0.0;
        for (PeopleNum = 1; PeopleNum <= state.dataHeatBal->TotPeople; ++PeopleNum) {
            if (ZoneNum != state.dataHeatBal->People(PeopleNum).ZonePtr) continue;
            PeopleSchPtr = state.dataHeatBal->People(PeopleNum).NumberOfPeoplePtr;
            MaxPeopleSch = GetScheduleMaxValue(state, PeopleSchPtr);
            NumberOfPeople = NumberOfPeople + (state.dataHeatBal->People(PeopleNum).NumberOfPeople * MaxPeopleSch);
        }
        SupplyAirVolFlowDes = FloorArea * thisStandAloneERV.AirVolFlowPerFloorArea +
                              NumberOfPeople * thisStandAloneERV.AirVolFlowPerOccupant;
        SupplyAirVolFlowDes = ZoneMult * SupplyAirVolFlowDes;

        if (SupplyAirVolFlowDes < SmallAirVolFlow) {
            SupplyAirVolFlowDes = 0.0;
        }

        // Size ERV supply flow rate
        Real64 TempSize = thisStandAloneERV.SupplyAirVolFlow;
        std::string SizingString = "Supply Air Flow Rate [m3/s]";
        if (IsAutoSize) {
            state.dataSize->DataConstantUsedForSizing = SupplyAirVolFlowDes;
            state.dataSize->DataFractionUsedForSizing = 1.0;
            TempSize = SupplyAirVolFlowDes;
            if (thisStandAloneERV.ControllerNameDefined) {
                state.dataMixedAir->OAController(thisStandAloneERV.ControllerIndex).MaxOA =
                    SupplyAirVolFlowDes * thisStandAloneERV.HighRHOAFlowRatio;
                state.dataMixedAir->OAController(thisStandAloneERV.ControllerIndex).MinOA =
                    SupplyAirVolFlowDes;
            }
        } else {
            state.dataSize->DataConstantUsedForSizing = thisStandAloneERV.SupplyAirVolFlow;
            state.dataSize->DataFractionUsedForSizing = 1.0;
        }
        if (TempSize > 0.0) {
            SystemAirFlowSizer sizerSystemAirFlow;
            sizerSystemAirFlow.overrideSizingString(SizingString);
            sizerSystemAirFlow.initializeWithinEP(state, CompType, CompName, PrintFlag, RoutineName);
            TempSize = sizerSystemAirFlow.size(state, TempSize, ErrorsFound);
        }
        thisStandAloneERV.SupplyAirVolFlow = TempSize;
    }

    // Size ERV exhaust flow rate
    state.dataSize->DataFractionUsedForSizing = 1.0;
    IsAutoSize = false;
    if (thisStandAloneERV.ExhaustAirVolFlow == DataSizing::AutoSize) {
        IsAutoSize = true;
    }

    if (state.dataSize->CurZoneEqNum > 0) {

        ExhaustAirVolFlowDes = SupplyAirVolFlowDes;

        if (ExhaustAirVolFlowDes < SmallAirVolFlow) {
            ExhaustAirVolFlowDes = 0.0;
        }

        if (ExhaustAirVolFlowDes > thisStandAloneERV.SupplyAirVolFlow) {
            ExhaustAirVolFlowDes = thisStandAloneERV.SupplyAirVolFlow;
        }

        std::string SizingString = "Exhaust Air Flow Rate [m3/s]";
        Real64 TempSize = thisStandAloneERV.ExhaustAirVolFlow;
        if (IsAutoSize) {
            TempSize = ExhaustAirVolFlowDes;
            state.dataSize->DataConstantUsedForSizing = ExhaustAirVolFlowDes;
        } else {
            state.dataSize->DataConstantUsedForSizing = thisStandAloneERV.ExhaustAirVolFlow;
        }
        state.dataSize->DataFractionUsedForSizing = 1.0;
        if (TempSize > 0.0) {
            SystemAirFlowSizer sizerSystemAirFlow;
            sizerSystemAirFlow.overrideSizingString(SizingString);
            sizerSystemAirFlow.initializeWithinEP(state, CompType, CompName, PrintFlag, RoutineName);
            TempSize = sizerSystemAirFlow.size(state, TempSize, ErrorsFound);
        }
        thisStandAloneERV.ExhaustAirVolFlow = TempSize;
        thisStandAloneERV.DesignEAFanVolFlowRate =
            TempSize * thisStandAloneERV.HighRHOAFlowRatio;
    }

    // Set Zone equipment sizing data for autosizing the fans and heat exchanger
    ZoneEqSizing(state.dataSize->CurZoneEqNum).AirVolFlow = thisStandAloneERV.SupplyAirVolFlow *
                                                            thisStandAloneERV.HighRHOAFlowRatio;
    ZoneEqSizing(state.dataSize->CurZoneEqNum).OAVolFlow = thisStandAloneERV.SupplyAirVolFlow;
    ZoneEqSizing(state.dataSize->CurZoneEqNum).SystemAirFlow = true;
    ZoneEqSizing(state.dataSize->CurZoneEqNum).DesignSizeFromParent = true;

    // Check supply fan flow rate or set flow rate if autosized in fan object
    IsAutoSize = false;
    if (thisStandAloneERV.DesignSAFanVolFlowRate == DataSizing::AutoSize) {
        IsAutoSize = true;
    }
    DesignSAFanVolFlowRateDes = thisStandAloneERV.SupplyAirVolFlow *
                                thisStandAloneERV.HighRHOAFlowRatio;
    if (IsAutoSize) {
        thisStandAloneERV.DesignSAFanVolFlowRate = DesignSAFanVolFlowRateDes;
    } else {
        if (thisStandAloneERV.DesignSAFanVolFlowRate > 0.0 && DesignSAFanVolFlowRateDes > 0.0) {
            DesignSAFanVolFlowRateUser = thisStandAloneERV.DesignSAFanVolFlowRate;
            if (state.dataGlobal->DisplayExtraWarnings) {
                if ((std::abs(DesignSAFanVolFlowRateDes - DesignSAFanVolFlowRateUser) / DesignSAFanVolFlowRateUser) >
                    state.dataSize->AutoVsHardSizingThreshold) {
                    ShowMessage(state,
                                format("SizeStandAloneERV: Potential issue with equipment sizing for ZoneHVAC:EnergyRecoveryVentilator {} {}",
                                       cFanTypes(thisStandAloneERV.SupplyAirFanType_Num),
                                       thisStandAloneERV.SupplyAirFanName));
                    ShowContinueError(state, format("User-Specified Supply Fan Maximum Flow Rate of {:.5R} [m3/s]", DesignSAFanVolFlowRateUser));
                    ShowContinueError(state, format("differs from the ERV Supply Air Flow Rate of {:.5R} [m3/s]", DesignSAFanVolFlowRateDes));
                    ShowContinueError(state, "This may, or may not, indicate mismatched component sizes.");
                    ShowContinueError(state, "Verify that the value entered is intended and is consistent with other components.");
                }
            }
        }
    }

    // simulate the fan to size using the flow rate specified above
    // (i.e., ZoneEqSizing( CurZoneEqNum ).AirVolFlow = StandAloneERV( StandAloneERVNum ).SupplyAirVolFlow * StandAloneERV( StandAloneERVNum
    // ).HighRHOAFlowRatio;)
    if (!(thisStandAloneERV.SupplyAirFanType_Num == DataHVACGlobals::FanType_SystemModelObject)) {
        SimulateFanComponents(state,
                              thisStandAloneERV.SupplyAirFanName,
                              true,
                              thisStandAloneERV.SupplyAirFanIndex);
    } else {
        state.dataHVACFan->fanObjs[thisStandAloneERV.SupplyAirFanIndex]->simulate(
            state, _, state.dataHVACGlobal->ZoneCompTurnFansOn, state.dataHVACGlobal->ZoneCompTurnFansOff, _);
    }
    if (!(thisStandAloneERV.ExhaustAirFanType_Num == DataHVACGlobals::FanType_SystemModelObject)) {
        SimulateFanComponents(state,
                              thisStandAloneERV.ExhaustAirFanName,
                              true,
                              thisStandAloneERV.ExhaustAirFanIndex);
    } else {
        state.dataHVACFan->fanObjs[thisStandAloneERV.ExhaustAirFanIndex]->simulate(
            state, _, state.dataHVACGlobal->ZoneCompTurnFansOn, state.dataHVACGlobal->ZoneCompTurnFansOff, _);
    }

    // now reset the ZoneEqSizing variable to NOT use the multiplier for HighRHOAFlowRatio for sizing HXs
    ZoneEqSizing(state.dataSize->CurZoneEqNum).AirVolFlow = thisStandAloneERV.SupplyAirVolFlow;
}

void CalcStandAloneERV(EnergyPlusData &state,
                       int const StandAloneERVNum,    // Unit index in ERV data structure
                       bool const FirstHVACIteration, // flag for 1st HVAC iteration in the time step
                       Real64 &SensLoadMet,           // sensible zone load met by unit (W)
                       Real64 &LatentMassLoadMet      // latent zone load met by unit (kg/s), dehumid = negative
)
{

    // SUBROUTINE INFORMATION:
    //       AUTHOR         Richard Raustad, FSEC
    //       DATE WRITTEN   June 2003
    //       MODIFIED       Don Shirey, Aug 2009 (LatentMassLoadMet)
    //                      July 2012, Chandan Sharma - FSEC: Added zone sys avail managers
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // Simulate the components making up the Stand Alone ERV unit.

    // METHODOLOGY EMPLOYED:
    // Simulates the unit components sequentially in the air flow direction.

    // Using/Aliasing
    using Fans::SimulateFanComponents;
    using HeatRecovery::SimHeatRecovery;

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    Real64 AirMassFlow;   // total mass flow through supply side of the ERV (supply air outlet node)
    // (so enthalpy routines work without error)
    Real64 TotLoadMet;           // total zone load met by unit (W)
    Real64 LatLoadMet;           // latent zone load met by unit (W)
    bool HXUnitOn;               // flag to operate heat exchanger heat recovery
    bool EconomizerFlag;         // economizer signal from OA controller
    bool HighHumCtrlFlag;        // high humditiy control signal from OA controller
    Real64 TotalExhaustMassFlow; // total exhaust air mass flow rate in controlled zone
    Real64 TotalSupplyMassFlow;  // total supply air mass flow rate in controlled zone

    auto &thisStandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum);

    int SupInletNode = thisStandAloneERV.SupplyAirInletNode;
    int SupOutletNode = thisStandAloneERV.SupplyAirOutletNode;
    int ExhaustInletNode = thisStandAloneERV.ExhaustAirInletNode;

    // Stand alone ERV's HX is ON by default
    HXUnitOn = true;

    // Get stand alone ERV's controller economizer and high humidity control status
    if (thisStandAloneERV.ControllerNameDefined) {
        EconomizerFlag = state.dataMixedAir->OAController(thisStandAloneERV.ControllerIndex).EconoActive;
        HighHumCtrlFlag =
            state.dataMixedAir->OAController(thisStandAloneERV.ControllerIndex).HighHumCtrlActive;
    } else {
        EconomizerFlag = false;
        HighHumCtrlFlag = false;
    }

    SimHeatRecovery(state,
                    thisStandAloneERV.HeatExchangerName,
                    FirstHVACIteration,
                    thisStandAloneERV.HeatExchangerIndex,
                    ContFanCycCoil,
                    _,
                    HXUnitOn,
                    _,
                    _,
                    EconomizerFlag,
                    HighHumCtrlFlag);
    thisStandAloneERV.ElecUseRate = state.dataHVACGlobal->AirToAirHXElecPower;

    if (thisStandAloneERV.SupplyAirFanType_Num != DataHVACGlobals::FanType_SystemModelObject) {
        SimulateFanComponents(state,
                              thisStandAloneERV.SupplyAirFanName,
                              FirstHVACIteration,
                              thisStandAloneERV.SupplyAirFanIndex,
                              _,
                              state.dataHVACGlobal->ZoneCompTurnFansOn,
                              state.dataHVACGlobal->ZoneCompTurnFansOff);
        thisStandAloneERV.ElecUseRate +=
            Fans::GetFanPower(state, thisStandAloneERV.SupplyAirFanIndex);
    } else {
        state.dataHVACFan->fanObjs[thisStandAloneERV.SupplyAirFanIndex]->simulate(
            state, _, state.dataHVACGlobal->ZoneCompTurnFansOn, state.dataHVACGlobal->ZoneCompTurnFansOff, _);
        thisStandAloneERV.ElecUseRate +=
            state.dataHVACFan->fanObjs[thisStandAloneERV.SupplyAirFanIndex]->fanPower();
    }

    if (thisStandAloneERV.ExhaustAirFanType_Num != DataHVACGlobals::FanType_SystemModelObject) {
        SimulateFanComponents(state,
                              thisStandAloneERV.ExhaustAirFanName,
                              FirstHVACIteration,
                              thisStandAloneERV.ExhaustAirFanIndex); // why no Turn on off flags here?
        thisStandAloneERV.ElecUseRate +=
            Fans::GetFanPower(state, thisStandAloneERV.ExhaustAirFanIndex);
    } else {
        state.dataHVACFan->fanObjs[thisStandAloneERV.ExhaustAirFanIndex]->simulate(
            state, _, state.dataHVACGlobal->ZoneCompTurnFansOn, state.dataHVACGlobal->ZoneCompTurnFansOff, _);
        thisStandAloneERV.ElecUseRate +=
            state.dataHVACFan->fanObjs[thisStandAloneERV.ExhaustAirFanIndex]->fanPower();
    }

    AirMassFlow = state.dataLoopNodes->Node(SupOutletNode).MassFlowRate;
    CalcZoneSensibleLatentOutput(AirMassFlow,
                                 state.dataLoopNodes->Node(SupOutletNode).Temp,
                                 state.dataLoopNodes->Node(SupOutletNode).HumRat,
                                 state.dataLoopNodes->Node(ExhaustInletNode).Temp,
                                 state.dataLoopNodes->Node(ExhaustInletNode).HumRat,
                                 SensLoadMet,
                                 LatLoadMet,
                                 TotLoadMet);
    LatentMassLoadMet = AirMassFlow * (state.dataLoopNodes->Node(SupOutletNode).HumRat -
                                       state.dataLoopNodes->Node(ExhaustInletNode).HumRat); // kg/s, dehumidification = negative

    if (SensLoadMet < 0.0) {
        thisStandAloneERV.SensCoolingRate = std::abs(SensLoadMet);
        thisStandAloneERV.SensHeatingRate = 0.0;
    } else {
        thisStandAloneERV.SensCoolingRate = 0.0;
        thisStandAloneERV.SensHeatingRate = SensLoadMet;
    }
    if (TotLoadMet < 0.0) {
        thisStandAloneERV.TotCoolingRate = std::abs(TotLoadMet);
        thisStandAloneERV.TotHeatingRate = 0.0;
    } else {
        thisStandAloneERV.TotCoolingRate = 0.0;
        thisStandAloneERV.TotHeatingRate = TotLoadMet;
    }
    if (LatLoadMet < 0.0) {
        thisStandAloneERV.LatCoolingRate = std::abs(LatLoadMet);
        thisStandAloneERV.LatHeatingRate = 0.0;
    } else {
        thisStandAloneERV.LatCoolingRate = 0.0;
        thisStandAloneERV.LatHeatingRate = LatLoadMet;
    }

    // Provide a one time message when exhaust flow rate is greater than supply flow rate
    if (thisStandAloneERV.FlowError && !state.dataGlobal->WarmupFlag) {
        TotalExhaustMassFlow = state.dataLoopNodes->Node(ExhaustInletNode).MassFlowRate;
        TotalSupplyMassFlow = state.dataLoopNodes->Node(SupInletNode).MassFlowRate;
        if (TotalExhaustMassFlow > TotalSupplyMassFlow && !state.dataHeatBal->ZoneAirMassFlow.EnforceZoneMassBalance) {
            ShowWarningError(state,
                             format("For {} \"{}\" there is unbalanced exhaust air flow.",
                                    thisStandAloneERV.UnitType,
                                    thisStandAloneERV.Name));
            ShowContinueError(state, format("... The exhaust air mass flow rate = {:.6R}", state.dataLoopNodes->Node(ExhaustInletNode).MassFlowRate));
            ShowContinueError(state, format("... The  supply air mass flow rate = {:.6R}", state.dataLoopNodes->Node(SupInletNode).MassFlowRate));
            ShowContinueErrorTimeStamp(state, "");
            ShowContinueError(state, "... Unless there is balancing infiltration / ventilation air flow, this will result in");
            ShowContinueError(state, "... load due to induced outside air being neglected in the simulation.");
            thisStandAloneERV.FlowError = false;
        }
    }
}

void ReportStandAloneERV(EnergyPlusData &state, int const StandAloneERVNum) // number of the current Stand Alone ERV being simulated
{

    // SUBROUTINE INFORMATION:
    //       AUTHOR         Richard Raustad, FSEC
    //       DATE WRITTEN   June 2003
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // Fill remaining report variables

    Real64 ReportingConstant = state.dataHVACGlobal->TimeStepSys * DataGlobalConstants::SecInHour;

    auto &thisStandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum);
    
    thisStandAloneERV.ElecUseEnergy = thisStandAloneERV.ElecUseRate * ReportingConstant;
    thisStandAloneERV.SensCoolingEnergy = thisStandAloneERV.SensCoolingRate * ReportingConstant;
    thisStandAloneERV.LatCoolingEnergy = thisStandAloneERV.LatCoolingRate * ReportingConstant;
    thisStandAloneERV.TotCoolingEnergy = thisStandAloneERV.TotCoolingRate * ReportingConstant;
    thisStandAloneERV.SensHeatingEnergy = thisStandAloneERV.SensHeatingRate * ReportingConstant;
    thisStandAloneERV.LatHeatingEnergy = thisStandAloneERV.LatHeatingRate * ReportingConstant;
    thisStandAloneERV.TotHeatingEnergy = thisStandAloneERV.TotHeatingRate * ReportingConstant;

    if (thisStandAloneERV.FirstPass) { // reset sizing flags so other zone equipment can size normally
        if (!state.dataGlobal->SysSizingCalc) {
            DataSizing::resetHVACSizingGlobals(state, state.dataSize->CurZoneEqNum, 0, thisStandAloneERV.FirstPass);
        }
    }
}

//        End of Reporting subroutines for the Module

//        Utility subroutines/functions for the HeatingCoil Module

Real64 GetSupplyAirFlowRate(EnergyPlusData &state,
                            std::string const &ERVType,     // must be "ZoneHVAC:EnergyRecoveryVentilator"
                            std::string const &ERVCtrlName, // must match a controller name in the ERV data structure
                            bool &ErrorsFound               // set to true if problem
)
{

    // FUNCTION INFORMATION:
    //       AUTHOR         Linda Lawrie
    //       DATE WRITTEN   October 2006
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS FUNCTION:
    // This function looks up the ERVCtrlName in the ERV Stand Alone list and returns the
    // Supply Air Flow rate, if found.  If incorrect name is given, ErrorsFound is returned as true
    // and supply air flow rate as negative.

    // Return value
    Real64 AirFlowRate; // returned supply air flow rate of the ERV unit

    // FUNCTION LOCAL VARIABLE DECLARATIONS:
    int WhichERV;

    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    if (UtilityRoutines::SameString(ERVType, "ZoneHVAC:EnergyRecoveryVentilator")) {
        WhichERV = UtilityRoutines::FindItem(ERVCtrlName, state.dataHVACStandAloneERV->StandAloneERV, &StandAloneERVData::ControllerName);
        if (WhichERV != 0) {
            AirFlowRate = state.dataHVACStandAloneERV->StandAloneERV(WhichERV).SupplyAirVolFlow;
        }
    } else {
        WhichERV = 0;
    }

    if (WhichERV == 0) {
        ShowSevereError(state, format("Could not find ZoneHVAC:EnergyRecoveryVentilator with Controller Name=\"{}\"", ERVCtrlName));
        ErrorsFound = true;
        AirFlowRate = -1000.0;
    }

    return AirFlowRate;
}

int GetSupplyAirInletNode(EnergyPlusData &state,
                          std::string const &ERVType,     // must be "ZoneHVAC:EnergyRecoveryVentilator"
                          std::string const &ERVCtrlName, // must match a controller name in the ERV data structure
                          bool &ErrorsFound               // set to true if problem
)
{

    // FUNCTION INFORMATION:
    //       AUTHOR         Linda Lawrie
    //       DATE WRITTEN   October 2006
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS FUNCTION:
    // This function looks up the ERVCtrlName in the ERV Stand Alone list and returns the
    // Supply Air Inlet Node Number, if found.  If incorrect name is given, ErrorsFound is returned as true
    // and Supply Air Inlet Node Number as zero.

    // Return value
    int AirInletNode(0); // returned air inlet node number of the ERV unit

    // FUNCTION LOCAL VARIABLE DECLARATIONS:
    int WhichERV;

    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    if (UtilityRoutines::SameString(ERVType, "ZoneHVAC:EnergyRecoveryVentilator")) {
        WhichERV = UtilityRoutines::FindItem(ERVCtrlName, state.dataHVACStandAloneERV->StandAloneERV, &StandAloneERVData::ControllerName);
        if (WhichERV != 0) {
            AirInletNode = state.dataHVACStandAloneERV->StandAloneERV(WhichERV).SupplyAirInletNode;
        }
    } else {
        WhichERV = 0;
    }

    if (WhichERV == 0) {
        ShowSevereError(state, format("Could not find ZoneHVAC:EnergyRecoveryVentilator with Controller Name=\"{}\"", ERVCtrlName));
        ErrorsFound = true;
        AirInletNode = 0;
    }

    return AirInletNode;
}

int GetExhaustAirInletNode(EnergyPlusData &state,
                           std::string const &ERVType,     // must be "ZoneHVAC:EnergyRecoveryVentilator"
                           std::string const &ERVCtrlName, // must match a controller name in the ERV data structure
                           bool &ErrorsFound               // set to true if problem
)
{

    // FUNCTION INFORMATION:
    //       AUTHOR         Linda Lawrie
    //       DATE WRITTEN   October 2006
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS FUNCTION:
    // This function looks up the ERVCtrlName in the ERV Stand Alone list and returns the
    // Exhaust Air Inlet Node Number, if found.  If incorrect name is given, ErrorsFound is returned as true
    // and Exhaust Air Inlet Node Number as zero.

    // Return value
    int AirInletNode(0); // returned air inlet node number of the ERV unit

    // FUNCTION LOCAL VARIABLE DECLARATIONS:
    int WhichERV;

    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    if (UtilityRoutines::SameString(ERVType, "ZoneHVAC:EnergyRecoveryVentilator")) {
        WhichERV = UtilityRoutines::FindItem(ERVCtrlName, state.dataHVACStandAloneERV->StandAloneERV, &StandAloneERVData::ControllerName);
        if (WhichERV != 0) {
            AirInletNode = state.dataHVACStandAloneERV->StandAloneERV(WhichERV).ExhaustAirInletNode;
        }
    } else {
        WhichERV = 0;
    }

    if (WhichERV == 0) {
        ShowSevereError(state, format("Could not find ZoneHVAC:EnergyRecoveryVentilator with Controller Name=\"{}\"", ERVCtrlName));
        ErrorsFound = true;
        AirInletNode = 0;
    }

    return AirInletNode;
}

int GetStandAloneERVOutAirNode(EnergyPlusData &state, int const StandAloneERVNum)
{
    // FUNCTION INFORMATION:
    //       AUTHOR         B Griffith
    //       DATE WRITTEN   Dec  2006
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS FUNCTION:
    // lookup function for OA inlet node for ventilation rate reporting

    // METHODOLOGY EMPLOYED:
    // <description>

    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    if (StandAloneERVNum > 0 && StandAloneERVNum <= state.dataHVACStandAloneERV->NumStandAloneERVs) {
	    return state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum).SupplyAirInletNode;
    } else {
        return 0;
    }
}

int GetStandAloneERVZoneInletAirNode(EnergyPlusData &state, int const StandAloneERVNum)
{
    // FUNCTION INFORMATION:
    //       AUTHOR         B Griffith
    //       DATE WRITTEN   Dec  2006
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS FUNCTION:
    // lookup function for OA inlet node for ventilation rate reporting

    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    if (StandAloneERVNum > 0 && StandAloneERVNum <= state.dataHVACStandAloneERV->NumStandAloneERVs) {
        return state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum).SupplyAirOutletNode;
    } else {
        return 0;
    }
}

int GetStandAloneERVReturnAirNode(EnergyPlusData &state, int const StandAloneERVNum)
{
    // FUNCTION INFORMATION:
    //       AUTHOR         B Griffith
    //       DATE WRITTEN   Dec  2006
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS FUNCTION:
    // lookup function for OA inlet node for ventilation rate reporting

    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    if (StandAloneERVNum > 0 && StandAloneERVNum <= state.dataHVACStandAloneERV->NumStandAloneERVs) {
        return state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVNum).ExhaustAirInletNode;
    } else {
        return 0;
    }
}

bool GetStandAloneERVNodeNumber(EnergyPlusData &state, int const NodeNumber)
{
    // PURPOSE OF THIS FUNCTION:
    // Check if a node is used by a stand alone ERV
    // and can be excluded from an airflow network.

    // Return value
    bool StandSloneERVAFNException;

    if (state.dataHVACStandAloneERV->GetERVInputFlag) {
        GetStandAloneERV(state);
        state.dataHVACStandAloneERV->GetERVInputFlag = false;
    }

    StandSloneERVAFNException = false;

    for (int StandAloneERVIndex = 1; StandAloneERVIndex <= state.dataHVACStandAloneERV->NumStandAloneERVs; ++StandAloneERVIndex) {

        auto &StandAloneERV = state.dataHVACStandAloneERV->StandAloneERV(StandAloneERVIndex);
        bool ErrorsFound{false};
        int SupplyFanInletNodeIndex(0);
        int SupplyFanOutletNodeIndex(0);
        int ExhaustFanInletNodeIndex(0);
        int ExhaustFanOutletNodeIndex(0);
        Real64 SupplyFanAirFlow;
        Real64 ExhaustFanAirFlow;

        // Get supply air fan inlet and outlet node index and air flow
        // ZoneHVAC:EnergyRecoveryVentilator only accepts Fan:SystemModel or Fan:OnOff
        if (StandAloneERV.SupplyAirFanType_Num == DataHVACGlobals::FanType_SystemModelObject) {
            // Fan:SystemModel
            SupplyFanInletNodeIndex = state.dataHVACFan->fanObjs[StandAloneERV.SupplyAirFanIndex]->inletNodeNum;
            SupplyFanOutletNodeIndex = state.dataHVACFan->fanObjs[StandAloneERV.SupplyAirFanIndex]->outletNodeNum;
            SupplyFanAirFlow = state.dataHVACFan->fanObjs[StandAloneERV.SupplyAirFanIndex]->designAirVolFlowRate;
        } else {
            // Fan:OnOff
            SupplyFanInletNodeIndex = Fans::GetFanInletNode(state, "Fan:OnOff", StandAloneERV.SupplyAirFanName, ErrorsFound);
            SupplyFanOutletNodeIndex = Fans::GetFanOutletNode(state, "Fan:OnOff", StandAloneERV.SupplyAirFanName, ErrorsFound);
            GetFanVolFlow(state, StandAloneERV.SupplyAirFanIndex, SupplyFanAirFlow);
            if (ErrorsFound) {
                ShowWarningError(state, format("Could not retrieve fan outlet node for this unit=\"{}\".", StandAloneERV.Name));
                ErrorsFound = true;
            }
        }
        // Get exhaust air fan inlet and outlet node index and air flow
        if (StandAloneERV.ExhaustAirFanType_Num == DataHVACGlobals::FanType_SystemModelObject) {
            // Fan:SystemModel
            ExhaustFanInletNodeIndex = state.dataHVACFan->fanObjs[StandAloneERV.ExhaustAirFanIndex]->inletNodeNum;
            ExhaustFanOutletNodeIndex = state.dataHVACFan->fanObjs[StandAloneERV.ExhaustAirFanIndex]->outletNodeNum;
            ExhaustFanAirFlow = state.dataHVACFan->fanObjs[StandAloneERV.ExhaustAirFanIndex]->designAirVolFlowRate;
        } else {
            // Fan:OnOff
            ExhaustFanInletNodeIndex = Fans::GetFanInletNode(state, "Fan:OnOff", StandAloneERV.ExhaustAirFanName, ErrorsFound);
            ExhaustFanOutletNodeIndex = Fans::GetFanOutletNode(state, "Fan:OnOff", StandAloneERV.ExhaustAirFanName, ErrorsFound);
            GetFanVolFlow(state, StandAloneERV.ExhaustAirFanIndex, ExhaustFanAirFlow);
            if (ErrorsFound) {
                ShowWarningError(state, format("Could not retrieve fan outlet node for this unit=\"{}\".", StandAloneERV.Name));
                ErrorsFound = true;
            }
        }

        // If a standalone ERV's airflow is unbalanced it shouldn't be model along with an AFN
        if (std::abs(SupplyFanAirFlow - ExhaustFanAirFlow) >= 1E-20 ||
            std::abs(StandAloneERV.DesignSAFanVolFlowRate - StandAloneERV.DesignEAFanVolFlowRate) >= 1E-20) {
            break;
        }

        // Supply air fan nodes
        if (NodeNumber == SupplyFanInletNodeIndex || NodeNumber == SupplyFanOutletNodeIndex || NodeNumber == ExhaustFanInletNodeIndex ||
            NodeNumber == ExhaustFanOutletNodeIndex) {
            StandSloneERVAFNException = true;
            break;
        }

        // Supply air inlet node
        if (NodeNumber == StandAloneERV.SupplyAirInletNode) {
            StandSloneERVAFNException = true;
            break;
        }
    }

    return StandSloneERVAFNException;
}

} // namespace EnergyPlus::HVACStandAloneERV

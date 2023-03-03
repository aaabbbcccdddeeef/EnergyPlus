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
#include <EnergyPlus/Autosizing/Base.hh>
#include <EnergyPlus/BranchNodeConnections.hh>
#include <EnergyPlus/Data/EnergyPlusData.hh>
#include <EnergyPlus/DataDefineEquip.hh>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/DataLoopNode.hh>
#include <EnergyPlus/DataSizing.hh>
#include <EnergyPlus/DataZoneEnergyDemands.hh>
#include <EnergyPlus/DataZoneEquipment.hh>
#include <EnergyPlus/FluidProperties.hh>
#include <EnergyPlus/General.hh>
#include <EnergyPlus/GeneralRoutines.hh>
#include <EnergyPlus/HVACSingleDuctInduc.hh>
#include <EnergyPlus/InputProcessing/InputProcessor.hh>
#include <EnergyPlus/MixerComponent.hh>
#include <EnergyPlus/NodeInputManager.hh>
#include <EnergyPlus/OutputProcessor.hh>
#include <EnergyPlus/Plant/DataPlant.hh>
#include <EnergyPlus/PlantUtilities.hh>
#include <EnergyPlus/Psychrometrics.hh>
#include <EnergyPlus/ScheduleManager.hh>
#include <EnergyPlus/UtilityRoutines.hh>
#include <EnergyPlus/WaterCoils.hh>

namespace EnergyPlus {

namespace HVACSingleDuctInduc {

    // Module containing routines dealing terminal 4 pipe induction terminal units

    // MODULE INFORMATION:
    //       AUTHOR         Fred Buhl
    //       DATE WRITTEN   June 15 2004
    //       MODIFIED       Brent Griffith, Sept 2010, plant upgrades, fluid props
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS MODULE:
    // To encapsulate the data and algorithms needed to simulate 4 pipe induction terminal units

    // METHODOLOGY EMPLOYED:
    // The terminal boxes are modeled as compound components: heating coil, cooling coil and
    // mixer. The combined components are controlled to meet the zone load.

    // Using/Aliasing
    using namespace DataLoopNode;
    using namespace ScheduleManager;
    using DataHVACGlobals::SmallAirVolFlow;
    using DataHVACGlobals::SmallLoad;
    using DataHVACGlobals::SmallMassFlow;
    using Psychrometrics::PsyCpAirFnW;
    using Psychrometrics::PsyHFnTdbW;
    using Psychrometrics::PsyRhoAirFnPbTdbW;

    void SimIndUnit(EnergyPlusData &state,
                    std::string_view CompName,     // name of the terminal unit
                    bool const FirstHVACIteration, // TRUE if first HVAC iteration in time step
                    int const ZoneNum,             // index of zone served by the terminal unit
                    int const ZoneNodeNum,         // zone node number of zone served by the terminal unit
                    int &CompIndex                 // which terminal unit in data structure
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Fred Buhl
        //       DATE WRITTEN   June 18 2004
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // Manages the simulation of a passive (no fan) induction terminal unit.
        // Called from SimZoneAirLoopEquipment in module ZoneAirLoopEquipmentManager.

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int IUNum; // index of terminal unit being simulated

        // First time SimIndUnit is called, get the input for all the passive terminal induction units
        if (state.dataHVACSingleDuctInduc->GetIUInputFlag) {
            GetIndUnits(state);
            state.dataHVACSingleDuctInduc->GetIUInputFlag = false;
        }

        auto &IndUnit = state.dataHVACSingleDuctInduc->IndUnit;

        // Get the induction unit index
        if (CompIndex == 0) {
            IUNum = UtilityRoutines::FindItemInList(CompName, state.dataHVACSingleDuctInduc->IndUnit);
            if (IUNum == 0) {
                ShowFatalError(state, format("SimIndUnit: Induction Unit not found={}", CompName));
            }
            CompIndex = IUNum;
        } else {
            IUNum = CompIndex;
            if (IUNum > state.dataHVACSingleDuctInduc->NumIndUnits || IUNum < 1) {
                ShowFatalError(state,
                               format("SimIndUnit: Invalid CompIndex passed={}, Number of Induction Units={}, System name={}",
                                      CompIndex,
                                      state.dataHVACSingleDuctInduc->NumIndUnits,
                                      CompName));
            }
            if (state.dataHVACSingleDuctInduc->CheckEquipName(IUNum)) {
                if (CompName != IndUnit(IUNum).Name) {
                    ShowFatalError(state,
                                   format("SimIndUnit: Invalid CompIndex passed={}, Induction Unit name={}, stored Induction Unit for that index={}",
                                          CompIndex,
                                          CompName,
                                          state.dataHVACSingleDuctInduc->IndUnit(IUNum).Name));
                }
                state.dataHVACSingleDuctInduc->CheckEquipName(IUNum) = false;
            }
        }

        state.dataSize->CurTermUnitSizingNum =
            state.dataDefineEquipment->AirDistUnit(state.dataHVACSingleDuctInduc->IndUnit(IUNum).ADUNum).TermUnitSizingNum;
        // initialize the unit
        InitIndUnit(state, IUNum, FirstHVACIteration);

        state.dataSize->TermUnitIU = true;

        // Select the correct unit type
        switch (state.dataHVACSingleDuctInduc->IndUnit(IUNum).UnitType_Num) {
        case SingleDuct_CV::FourPipeInduc: {
            SimFourPipeIndUnit(state, IUNum, ZoneNum, ZoneNodeNum, FirstHVACIteration);
        } break;
        default: {
            ShowSevereError(state, format("Illegal Induction Unit Type used={}", state.dataHVACSingleDuctInduc->IndUnit(IUNum).UnitType));
            ShowContinueError(state, format("Occurs in Induction Unit={}", state.dataHVACSingleDuctInduc->IndUnit(IUNum).Name));
            ShowFatalError(state, "Preceding condition causes termination.");
        } break;
        }

        state.dataSize->TermUnitIU = false;

        // the tasks usually done by the Update and Report routines are not required in a compound terminal unit.

        // Update the current unit's outlet nodes. No update needed

        // Fill the report variables. There are no report variables
        state.dataHVACSingleDuctInduc->IndUnit(IUNum).ReportIndUnit(state);
    }

    void GetIndUnits(EnergyPlusData &state)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Fred Buhl
        //       DATE WRITTEN   June 15 2004
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // Obtains input data for passive induction air terminal units and stores it in the
        // induction terminal unit data structures

        // METHODOLOGY EMPLOYED:
        // Uses "Get" routines to read in data.

        // Using/Aliasing
        using BranchNodeConnections::SetUpCompSets;
        using BranchNodeConnections::TestCompSet;
        using NodeInputManager::GetOnlySingleNode;
        using namespace DataSizing;

        using MixerComponent::GetZoneMixerIndex;
        using WaterCoils::GetCoilWaterInletNode;

        // SUBROUTINE PARAMETER DEFINITIONS:
        static constexpr std::string_view RoutineName("GetIndUnits "); // include trailing blank space

        int IUNum;                       // current fan coil number
        std::string CurrentModuleObject; // for ease in getting objects
        Array1D_string Alphas;           // Alpha input items for object
        Array1D_string cAlphaFields;     // Alpha field names
        Array1D_string cNumericFields;   // Numeric field names
        Array1D<Real64> Numbers;         // Numeric input items for object
        Array1D_bool lAlphaBlanks;       // Logical array, alpha field input BLANK = .TRUE.
        Array1D_bool lNumericBlanks;     // Logical array, numeric field input BLANK = .TRUE.
        int NumAlphas(0);                // Number of Alphas for each GetObjectItem call
        int NumNumbers(0);               // Number of Numbers for each GetObjectItem call
        int TotalArgs(0);                // Total number of alpha and numeric arguments (max) for a
        //  certain object in the input file
        int IOStatus;            // Used in GetObjectItem
        bool ErrorsFound(false); // Set to true if errors in input, fatal at end of routine
        bool IsNotOK;            // Flag to verify name
        bool AirNodeFound;
        int ADUNum;
        bool errFlag;

        // find the number of each type of induction unit
        CurrentModuleObject = "AirTerminal:SingleDuct:ConstantVolume:FourPipeInduction";
        state.dataHVACSingleDuctInduc->NumFourPipes = state.dataInputProcessing->inputProcessor->getNumObjectsFound(state, CurrentModuleObject);
        state.dataHVACSingleDuctInduc->NumIndUnits = state.dataHVACSingleDuctInduc->NumFourPipes;
        // allocate the data structures
        state.dataHVACSingleDuctInduc->IndUnit.allocate(state.dataHVACSingleDuctInduc->NumIndUnits);
        state.dataHVACSingleDuctInduc->CheckEquipName.dimension(state.dataHVACSingleDuctInduc->NumIndUnits, true);

        state.dataInputProcessing->inputProcessor->getObjectDefMaxArgs(state, CurrentModuleObject, TotalArgs, NumAlphas, NumNumbers);

        Alphas.allocate(NumAlphas);
        cAlphaFields.allocate(NumAlphas);
        cNumericFields.allocate(NumNumbers);
        Numbers.dimension(NumNumbers, 0.0);
        lAlphaBlanks.dimension(NumAlphas, true);
        lNumericBlanks.dimension(NumNumbers, true);

        // loop over Series PIUs; get and load the input data
        for (int IUIndex = 1; IUIndex <= state.dataHVACSingleDuctInduc->NumFourPipes; ++IUIndex) {

            state.dataInputProcessing->inputProcessor->getObjectItem(state,
                                                                     CurrentModuleObject,
                                                                     IUIndex,
                                                                     Alphas,
                                                                     NumAlphas,
                                                                     Numbers,
                                                                     NumNumbers,
                                                                     IOStatus,
                                                                     lNumericBlanks,
                                                                     lAlphaBlanks,
                                                                     cAlphaFields,
                                                                     cNumericFields);

            IUNum = IUIndex;

	    auto &thisIndUnit = state.dataHVACSingleDuctInduc->IndUnit(IUNum);
	    
            thisIndUnit.Name = Alphas(1);
            thisIndUnit.UnitType = CurrentModuleObject;
            thisIndUnit.UnitType_Num = SingleDuct_CV::FourPipeInduc;
            thisIndUnit.Sched = Alphas(2);
            if (lAlphaBlanks(2)) {
                thisIndUnit.SchedPtr = DataGlobalConstants::ScheduleAlwaysOn;
            } else {
                thisIndUnit.SchedPtr = GetScheduleIndex(state, Alphas(2)); // convert schedule name to pointer
                if (thisIndUnit.SchedPtr == 0) {
                    ShowSevereError(state,
                                    format("{}{}: invalid {} entered ={} for {}={}",
                                           RoutineName,
                                           CurrentModuleObject,
                                           cAlphaFields(2),
                                           thisIndUnit.Sched,
                                           cAlphaFields(1),
                                           thisIndUnit.Name));
                    ErrorsFound = true;
                }
            }
            thisIndUnit.MaxTotAirVolFlow = Numbers(1);

            thisIndUnit.InducRatio = lNumericBlanks(2) ? 2.5 : Numbers(2);

            thisIndUnit.PriAirInNode =
                GetOnlySingleNode(state,
                                  Alphas(3),
                                  ErrorsFound,
                                  DataLoopNode::ConnectionObjectType::AirTerminalSingleDuctConstantVolumeFourPipeInduction,
                                  thisIndUnit.Name,
                                  DataLoopNode::NodeFluidType::Air,
                                  DataLoopNode::ConnectionType::Inlet,
                                  NodeInputManager::CompFluidStream::Primary,
                                  ObjectIsParent,
                                  cAlphaFields(3));
            thisIndUnit.SecAirInNode =
                GetOnlySingleNode(state,
                                  Alphas(4),
                                  ErrorsFound,
                                  DataLoopNode::ConnectionObjectType::AirTerminalSingleDuctConstantVolumeFourPipeInduction,
                                  thisIndUnit.Name,
                                  DataLoopNode::NodeFluidType::Air,
                                  DataLoopNode::ConnectionType::Inlet,
                                  NodeInputManager::CompFluidStream::Primary,
                                  ObjectIsParent,
                                  cAlphaFields(4));
            thisIndUnit.OutAirNode =
                GetOnlySingleNode(state,
                                  Alphas(5),
                                  ErrorsFound,
                                  DataLoopNode::ConnectionObjectType::AirTerminalSingleDuctConstantVolumeFourPipeInduction,
                                  thisIndUnit.Name,
                                  DataLoopNode::NodeFluidType::Air,
                                  DataLoopNode::ConnectionType::Outlet,
                                  NodeInputManager::CompFluidStream::Primary,
                                  ObjectIsParent,
                                  cAlphaFields(5));

            thisIndUnit.HCoilType = Alphas(6); // type (key) of heating coil
            if (UtilityRoutines::SameString(thisIndUnit.HCoilType, "Coil:Heating:Water")) {
                thisIndUnit.HeatingCoilType = DataPlant::PlantEquipmentType::CoilWaterSimpleHeating;
            }

            thisIndUnit.HCoil = Alphas(7); // name of heating coil object
            IsNotOK = false;
            thisIndUnit.HWControlNode = GetCoilWaterInletNode(
                state, thisIndUnit.HCoilType, thisIndUnit.HCoil, IsNotOK);
            if (IsNotOK) {
                ShowContinueError(state, format("In {} = {}", CurrentModuleObject, thisIndUnit.Name));
                ShowContinueError(state, "..Only Coil:Heating:Water is allowed.");
                ErrorsFound = true;
            }
            thisIndUnit.MaxVolHotWaterFlow = Numbers(3);
            thisIndUnit.MinVolHotWaterFlow = Numbers(4);
            thisIndUnit.HotControlOffset = Numbers(5);

            thisIndUnit.CCoilType = Alphas(8); // type (key) of cooling coil

            if (UtilityRoutines::SameString(thisIndUnit.CCoilType, "Coil:Cooling:Water")) {
                thisIndUnit.CoolingCoilType = DataPlant::PlantEquipmentType::CoilWaterCooling;
            } else if (UtilityRoutines::SameString(thisIndUnit.CCoilType, "Coil:Cooling:Water:DetailedGeometry")) {
                thisIndUnit.CoolingCoilType = DataPlant::PlantEquipmentType::CoilWaterDetailedFlatCooling;
            }

            thisIndUnit.CCoil = Alphas(9); // name of cooling coil object
            IsNotOK = false;
            thisIndUnit.CWControlNode = GetCoilWaterInletNode(
                state, thisIndUnit.CCoilType, thisIndUnit.CCoil, IsNotOK);
            if (IsNotOK) {
                ShowContinueError(state, format("In {} = {}", CurrentModuleObject, thisIndUnit.Name));
                ShowContinueError(state, "..Only Coil:Cooling:Water or Coil:Cooling:Water:DetailedGeometry is allowed.");
                ErrorsFound = true;
            }
            thisIndUnit.MaxVolColdWaterFlow = Numbers(6);
            thisIndUnit.MinVolColdWaterFlow = Numbers(7);
            thisIndUnit.ColdControlOffset = Numbers(8);

            // Get the Zone Mixer name and check that it is OK
            errFlag = false;
            thisIndUnit.MixerName = Alphas(10);
            GetZoneMixerIndex(state,
                              thisIndUnit.MixerName,
                              thisIndUnit.Mixer_Num,
                              errFlag,
                              CurrentModuleObject);
            if (errFlag) {
                ShowContinueError(state, format("...specified in {} = {}", CurrentModuleObject, thisIndUnit.Name));
                ErrorsFound = true;
            }

            // Add heating coil to component sets array
            SetUpCompSets(state,
                          thisIndUnit.UnitType,
                          thisIndUnit.Name,
                          thisIndUnit.HCoilType,
                          thisIndUnit.HCoil,
                          Alphas(4),
                          "UNDEFINED");
            // Add cooling coil to component sets array
            SetUpCompSets(state,
                          thisIndUnit.UnitType,
                          thisIndUnit.Name,
                          thisIndUnit.CCoilType,
                          thisIndUnit.CCoil,
                          "UNDEFINED",
                          "UNDEFINED");

            // Register component set data
            TestCompSet(state,
                        thisIndUnit.UnitType,
                        thisIndUnit.Name,
                        state.dataLoopNodes->NodeID(thisIndUnit.PriAirInNode),
                        state.dataLoopNodes->NodeID(thisIndUnit.OutAirNode),
                        "Air Nodes");

            AirNodeFound = false;
            for (ADUNum = 1; ADUNum <= (int)state.dataDefineEquipment->AirDistUnit.size(); ++ADUNum) {
                if (thisIndUnit.OutAirNode == state.dataDefineEquipment->AirDistUnit(ADUNum).OutletNodeNum) {
                    thisIndUnit.ADUNum = ADUNum;
                }
            }
            // one assumes if there isn't one assigned, it's an error?
            if (thisIndUnit.ADUNum == 0) {
                ShowSevereError(state,
                                format("{}No matching Air Distribution Unit, for Unit = [{},{}].",
                                       RoutineName,
                                       thisIndUnit.UnitType,
                                       thisIndUnit.Name));
                ShowContinueError(
                    state,
                    format("...should have outlet node={}", state.dataLoopNodes->NodeID(thisIndUnit.OutAirNode)));
                ErrorsFound = true;
            } else {
                // Fill the Zone Equipment data with the supply air inlet node number of this unit.
                for (int CtrlZone = 1; CtrlZone <= state.dataGlobal->NumOfZones; ++CtrlZone) {
	            auto &thisZoneEquipConfig = state.dataZoneEquip->ZoneEquipConfig(CtrlZone);
                    if (!thisZoneEquipConfig.IsControlled) continue;
                    for (int SupAirIn = 1; SupAirIn <= thisZoneEquipConfig.NumInletNodes; ++SupAirIn) {
                        if (thisIndUnit.OutAirNode == thisZoneEquipConfig.InletNode(SupAirIn)) {
                            if (thisZoneEquipConfig.AirDistUnitCool(SupAirIn).OutNode > 0) {
                                ShowSevereError(state, "Error in connecting a terminal unit to a zone");
                                ShowContinueError(state,
                                                  format("{} already connects to another zone",
                                                         state.dataLoopNodes->NodeID(thisIndUnit.OutAirNode)));
                                ShowContinueError(state,
                                                  format("Occurs for terminal unit {} = {}",
                                                         thisIndUnit.UnitType,
                                                         thisIndUnit.Name));
                                ShowContinueError(state, "Check terminal unit node names for errors");
                                ErrorsFound = true;
                            } else {
                                thisZoneEquipConfig.AirDistUnitCool(SupAirIn).InNode = thisIndUnit.PriAirInNode;
                                thisZoneEquipConfig.AirDistUnitCool(SupAirIn).OutNode = thisIndUnit.OutAirNode;
                                state.dataDefineEquipment->AirDistUnit(thisIndUnit.ADUNum).TermUnitSizingNum =
                                    thisZoneEquipConfig.AirDistUnitCool(SupAirIn).TermUnitSizingIndex;
                                state.dataDefineEquipment->AirDistUnit(thisIndUnit.ADUNum).ZoneEqNum = CtrlZone;
                                thisIndUnit.CtrlZoneNum = CtrlZone;
                            }
                            thisIndUnit.CtrlZoneInNodeIndex = SupAirIn;
                            AirNodeFound = true;
                            break;
                        }
                    }
                }
                if (!AirNodeFound) {
                    ShowSevereError(
                        state,
                        format("The outlet air node from the {} = {}", CurrentModuleObject, thisIndUnit.Name));
                    ShowContinueError(state, format("did not have a matching Zone Equipment Inlet Node, Node ={}", Alphas(3)));
                    ErrorsFound = true;
                }
            }
            // report variable for all single duct air terminals
            SetupOutputVariable(state,
                                "Zone Air Terminal Outdoor Air Volume Flow Rate",
                                OutputProcessor::Unit::m3_s,
                                thisIndUnit.OutdoorAirFlowRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisIndUnit.Name);
        }

        Alphas.deallocate();
        cAlphaFields.deallocate();
        cNumericFields.deallocate();
        Numbers.deallocate();
        lAlphaBlanks.deallocate();
        lNumericBlanks.deallocate();
        if (ErrorsFound) {
            ShowFatalError(state, format("{}Errors found in getting input. Preceding conditions cause termination.", RoutineName));
        }
    }

    void InitIndUnit(EnergyPlusData &state,
                     int const IUNum,              // number of the current induction unit being simulated
                     bool const FirstHVACIteration // TRUE if first air loop solution this HVAC step
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Fred Buhl
        //       DATE WRITTEN   June 21 2004
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for initialization of the passive induction
        // terminal boxes

        // METHODOLOGY EMPLOYED:
        // Uses the status flags to trigger initializations.

        // Using/Aliasing

        using DataZoneEquipment::CheckZoneEquipmentList;
        using FluidProperties::GetDensityGlycol;
        using PlantUtilities::InitComponentNodes;
        using PlantUtilities::ScanPlantLoopsForObject;

        // SUBROUTINE PARAMETER DEFINITIONS:
        static constexpr std::string_view RoutineName("InitIndUnit");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int PriNode;     // primary air inlet node number
        int SecNode;     // secondary air inlet node number
        int OutletNode;  // unit air outlet node
        int HotConNode;  // hot water control node number
        int ColdConNode; // cold water control node  number
        Real64 IndRat;   // unit induction ratio
        Real64 RhoAir;   // air density at outside pressure and standard temperature and humidity

        int Loop;         // Loop checking control variable
        Real64 rho;       // local fluid density
        int HWOutletNode; // local node index for hot water coil's outlet node
        int CWOutletNode; // local node index for cold water coil's outlet node
        bool errFlag(false);

        auto &ZoneEquipmentListChecked = state.dataHVACSingleDuctInduc->ZoneEquipmentListChecked;

        // Do the one time initializations
        if (state.dataHVACSingleDuctInduc->MyOneTimeFlag) {

            state.dataHVACSingleDuctInduc->MyEnvrnFlag.allocate(state.dataHVACSingleDuctInduc->NumIndUnits);
            state.dataHVACSingleDuctInduc->MySizeFlag.allocate(state.dataHVACSingleDuctInduc->NumIndUnits);
            state.dataHVACSingleDuctInduc->MyPlantScanFlag.allocate(state.dataHVACSingleDuctInduc->NumIndUnits);
            state.dataHVACSingleDuctInduc->MyAirDistInitFlag.allocate(state.dataHVACSingleDuctInduc->NumIndUnits);
            state.dataHVACSingleDuctInduc->MyEnvrnFlag = true;
            state.dataHVACSingleDuctInduc->MySizeFlag = true;
            state.dataHVACSingleDuctInduc->MyPlantScanFlag = true;
            state.dataHVACSingleDuctInduc->MyAirDistInitFlag = true;
            state.dataHVACSingleDuctInduc->MyOneTimeFlag = false;
        }


	auto &thisIndUnit = state.dataHVACSingleDuctInduc->IndUnit(IUNum);
	
        if (state.dataHVACSingleDuctInduc->MyPlantScanFlag(IUNum) && allocated(state.dataPlnt->PlantLoop)) {
            if (thisIndUnit.HeatingCoilType == DataPlant::PlantEquipmentType::CoilWaterSimpleHeating) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisIndUnit.HCoil,
                                        thisIndUnit.HeatingCoilType,
                                        thisIndUnit.HWPlantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
            }
            if (errFlag) {
                ShowContinueError(state,
                                  format("Reference Unit=\"{}\", type={}",
                                         thisIndUnit.Name,
                                         thisIndUnit.UnitType));
            }
            if (thisIndUnit.CoolingCoilType == DataPlant::PlantEquipmentType::CoilWaterCooling ||
                thisIndUnit.CoolingCoilType == DataPlant::PlantEquipmentType::CoilWaterDetailedFlatCooling) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisIndUnit.CCoil,
                                        thisIndUnit.CoolingCoilType,
                                        thisIndUnit.CWPlantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
            }
            if (errFlag) {
                ShowContinueError(state,
                                  format("Reference Unit=\"{}\", type={}",
                                         thisIndUnit.Name,
                                         thisIndUnit.UnitType));
                ShowFatalError(state, "InitIndUnit: Program terminated for previous conditions.");
            }
            state.dataHVACSingleDuctInduc->MyPlantScanFlag(IUNum) = false;
        } else if (state.dataHVACSingleDuctInduc->MyPlantScanFlag(IUNum) && !state.dataGlobal->AnyPlantInModel) {
            state.dataHVACSingleDuctInduc->MyPlantScanFlag(IUNum) = false;
        }

        if (state.dataHVACSingleDuctInduc->MyAirDistInitFlag(IUNum)) {
            // save the induction ratio in the term unit sizing array for use in the system sizing calculation
            if (state.dataSize->CurTermUnitSizingNum > 0) {
                state.dataSize->TermUnitSizing(state.dataSize->CurTermUnitSizingNum).InducRat =
                    thisIndUnit.InducRatio;
            }
            if (thisIndUnit.AirLoopNum == 0) {
                if ((thisIndUnit.CtrlZoneNum > 0) &&
                    (thisIndUnit.CtrlZoneInNodeIndex > 0)) {
                    thisIndUnit.AirLoopNum =
                        state.dataZoneEquip->ZoneEquipConfig(thisIndUnit.CtrlZoneNum)
                            .InletNodeAirLoopNum(thisIndUnit.CtrlZoneInNodeIndex);
                    state.dataDefineEquipment->AirDistUnit(thisIndUnit.ADUNum).AirLoopNum =
                        thisIndUnit.AirLoopNum;
                }
            } else {
                state.dataHVACSingleDuctInduc->MyAirDistInitFlag(IUNum) = false;
            }
        }
        if (!ZoneEquipmentListChecked && state.dataZoneEquip->ZoneEquipInputsFilled) {
            ZoneEquipmentListChecked = true;
            // Check to see if there is a Air Distribution Unit on the Zone Equipment List
            for (Loop = 1; Loop <= state.dataHVACSingleDuctInduc->NumIndUnits; ++Loop) {
                if (state.dataHVACSingleDuctInduc->IndUnit(Loop).ADUNum == 0) continue;
                if (CheckZoneEquipmentList(state,
                                           "ZONEHVAC:AIRDISTRIBUTIONUNIT",
                                           state.dataDefineEquipment->AirDistUnit(state.dataHVACSingleDuctInduc->IndUnit(Loop).ADUNum).Name))
                    continue;
                ShowSevereError(state,
                                format("InitIndUnit: ADU=[Air Distribution Unit,{}] is not on any ZoneHVAC:EquipmentList.",
                                       state.dataDefineEquipment->AirDistUnit(state.dataHVACSingleDuctInduc->IndUnit(Loop).ADUNum).Name));
                ShowContinueError(state,
                                  format("...Unit=[{},{}] will not be simulated.",
                                         state.dataHVACSingleDuctInduc->IndUnit(Loop).UnitType,
                                         state.dataHVACSingleDuctInduc->IndUnit(Loop).Name));
            }
        }

        if (!state.dataGlobal->SysSizingCalc && state.dataHVACSingleDuctInduc->MySizeFlag(IUNum)) {

            SizeIndUnit(state, IUNum);
            state.dataHVACSingleDuctInduc->MySizeFlag(IUNum) = false;
        }

        // Do the Begin Environment initializations
        if (state.dataGlobal->BeginEnvrnFlag && state.dataHVACSingleDuctInduc->MyEnvrnFlag(IUNum)) {
            RhoAir = state.dataEnvrn->StdRhoAir;
            PriNode = thisIndUnit.PriAirInNode;
            SecNode = thisIndUnit.SecAirInNode;
            OutletNode = thisIndUnit.OutAirNode;
            IndRat = thisIndUnit.InducRatio;
            // set the mass flow rates from the input volume flow rates
            if (UtilityRoutines::SameString(thisIndUnit.UnitType,
                                            "AirTerminal:SingleDuct:ConstantVolume:FourPipeInduction")) {
                thisIndUnit.MaxTotAirMassFlow = RhoAir * thisIndUnit.MaxTotAirVolFlow;
                thisIndUnit.MaxPriAirMassFlow = thisIndUnit.MaxTotAirMassFlow / (1.0 + IndRat);
                thisIndUnit.MaxSecAirMassFlow = IndRat * thisIndUnit.MaxTotAirMassFlow / (1.0 + IndRat);
                state.dataLoopNodes->Node(PriNode).MassFlowRateMax = thisIndUnit.MaxPriAirMassFlow;
                state.dataLoopNodes->Node(PriNode).MassFlowRateMin = thisIndUnit.MaxPriAirMassFlow;
                state.dataLoopNodes->Node(SecNode).MassFlowRateMax = thisIndUnit.MaxSecAirMassFlow;
                state.dataLoopNodes->Node(SecNode).MassFlowRateMin = thisIndUnit.MaxSecAirMassFlow;
                state.dataLoopNodes->Node(OutletNode).MassFlowRateMax = thisIndUnit.MaxTotAirMassFlow;
            }

            HotConNode = thisIndUnit.HWControlNode;
            if (HotConNode > 0 && !state.dataHVACSingleDuctInduc->MyPlantScanFlag(IUNum)) {

                rho = GetDensityGlycol(state,
                                       state.dataPlnt->PlantLoop(thisIndUnit.HWPlantLoc.loopNum).FluidName,
                                       DataGlobalConstants::HWInitConvTemp,
                                       state.dataPlnt->PlantLoop(thisIndUnit.HWPlantLoc.loopNum).FluidIndex,
                                       RoutineName);
                thisIndUnit.MaxHotWaterFlow = rho * thisIndUnit.MaxVolHotWaterFlow;
                thisIndUnit.MinHotWaterFlow = rho * thisIndUnit.MinVolHotWaterFlow;
                // get component outlet node from plant structure
                HWOutletNode = DataPlant::CompData::getPlantComponent(state, thisIndUnit.HWPlantLoc).NodeNumOut;
                InitComponentNodes(state,
                                   thisIndUnit.MinHotWaterFlow,
                                   thisIndUnit.MaxHotWaterFlow,
                                   HotConNode,
                                   HWOutletNode);
            }

            ColdConNode = thisIndUnit.CWControlNode;
            if (ColdConNode > 0) {
                rho = GetDensityGlycol(state,
                                       state.dataPlnt->PlantLoop(thisIndUnit.CWPlantLoc.loopNum).FluidName,
                                       DataGlobalConstants::CWInitConvTemp,
                                       state.dataPlnt->PlantLoop(thisIndUnit.CWPlantLoc.loopNum).FluidIndex,
                                       RoutineName);
                thisIndUnit.MaxColdWaterFlow = rho * thisIndUnit.MaxVolColdWaterFlow;
                thisIndUnit.MinColdWaterFlow = rho * thisIndUnit.MinVolColdWaterFlow;

                CWOutletNode = DataPlant::CompData::getPlantComponent(state, thisIndUnit.CWPlantLoc).NodeNumOut;
                InitComponentNodes(state,
                                   thisIndUnit.MinColdWaterFlow,
                                   thisIndUnit.MaxColdWaterFlow,
                                   ColdConNode,
                                   CWOutletNode);
            }

            state.dataHVACSingleDuctInduc->MyEnvrnFlag(IUNum) = false;
        } // end one time inits

        if (!state.dataGlobal->BeginEnvrnFlag) {
            state.dataHVACSingleDuctInduc->MyEnvrnFlag(IUNum) = true;
        }

        auto &thisPriNode = state.dataLoopNodes->Node(thisIndUnit.PriAirInNode);
        auto &thisSecNode = state.dataLoopNodes->Node(thisIndUnit.SecAirInNode);

        // Do the start of HVAC time step initializations
        if (FirstHVACIteration) {
            // check for upstream zero flow. If nonzero and schedule ON, set primary flow to max
            if (GetCurrentScheduleValue(state, thisIndUnit.SchedPtr) > 0.0 &&
                thisPriNode.MassFlowRate > 0.0) {
                if (UtilityRoutines::SameString(thisIndUnit.UnitType,
                                                "AirTerminal:SingleDuct:ConstantVolume:FourPipeInduction")) {
                    thisPriNode.MassFlowRate = thisIndUnit.MaxPriAirMassFlow;
                    thisSecNode.MassFlowRate = thisIndUnit.MaxSecAirMassFlow;
                }
            } else {
                thisPriNode.MassFlowRate = 0.0;
                thisSecNode.MassFlowRate = 0.0;
            }
            // reset the max and min avail flows
            if (GetCurrentScheduleValue(state, thisIndUnit.SchedPtr) > 0.0 &&
                thisPriNode.MassFlowRateMaxAvail > 0.0) {
                if (UtilityRoutines::SameString(thisIndUnit.UnitType,
                                                "AirTerminal:SingleDuct:ConstantVolume:FourPipeInduction")) {
                    thisPriNode.MassFlowRateMaxAvail = thisIndUnit.MaxPriAirMassFlow;
                    thisPriNode.MassFlowRateMinAvail = thisIndUnit.MaxPriAirMassFlow;
                    thisSecNode.MassFlowRateMaxAvail = thisIndUnit.MaxSecAirMassFlow;
                    thisSecNode.MassFlowRateMinAvail = thisIndUnit.MaxSecAirMassFlow;
                }
            } else {
                thisPriNode.MassFlowRateMaxAvail = 0.0;
                thisPriNode.MassFlowRateMinAvail = 0.0;
                thisSecNode.MassFlowRateMaxAvail = 0.0;
                thisSecNode.MassFlowRateMinAvail = 0.0;
            }
        }
    }

    void SizeIndUnit(EnergyPlusData &state, int const IUNum)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Fred Buhl
        //       DATE WRITTEN   June 22 2004
        //       MODIFIED       August 2013 Daeho Kang, add component sizing table entries
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for sizing induction terminal units for which flow rates have not been
        // specified in the input

        // METHODOLOGY EMPLOYED:
        // Accesses zone sizing array for air flow rates and zone and plant sizing arrays to
        // calculate coil water flow rates.

        // Using/Aliasing
        using namespace DataSizing;
        using FluidProperties::GetDensityGlycol;
        using FluidProperties::GetSpecificHeatGlycol;

        using PlantUtilities::MyPlantSizingIndex;
        using WaterCoils::GetCoilWaterInletNode;
        using WaterCoils::GetCoilWaterOutletNode;
        using WaterCoils::SetCoilDesFlow;

        // SUBROUTINE PARAMETER DEFINITIONS:
        static constexpr std::string_view RoutineName("SizeIndUnit");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int PltSizHeatNum; // index of plant sizing object for 1st heating loop
        int PltSizCoolNum; // index of plant sizing object for 1st cooling loop
        Real64 DesCoilLoad;
        Real64 DesPriVolFlow;
        Real64 RhoAir;
        Real64 CpAir;
        int CoilWaterInletNode(0);
        int CoilWaterOutletNode(0);
        bool ErrorsFound;
        Real64 Cp;  // local fluid specific heat
        Real64 rho; // local fluid density
        bool IsAutoSize;
        Real64 MaxTotAirVolFlowDes;     // Desing size maximum air volume flow for reproting
        Real64 MaxTotAirVolFlowUser;    // User hard-sized maximum air volume flow for reporting
        Real64 MaxVolHotWaterFlowDes;   // Desing size maximum hot water flow for reproting
        Real64 MaxVolHotWaterFlowUser;  // User hard-sized maximum hot water flow for reporting
        Real64 MaxVolColdWaterFlowDes;  // Desing size maximum cold water flow for reproting
        Real64 MaxVolColdWaterFlowUser; // User hard-sized maximum cold water flow for reporting

        PltSizHeatNum = 0;
        PltSizCoolNum = 0;
        DesPriVolFlow = 0.0;
        CpAir = 0.0;
        RhoAir = state.dataEnvrn->StdRhoAir;
        ErrorsFound = false;
        IsAutoSize = false;
        MaxTotAirVolFlowDes = 0.0;
        MaxTotAirVolFlowUser = 0.0;
        MaxVolHotWaterFlowDes = 0.0;
        MaxVolHotWaterFlowUser = 0.0;
        MaxVolColdWaterFlowDes = 0.0;
        MaxVolColdWaterFlowUser = 0.0;

        auto &TermUnitSizing(state.dataSize->TermUnitSizing);

        auto &thisIndUnit = state.dataHVACSingleDuctInduc->IndUnit(IUNum);
	
        if (thisIndUnit.MaxTotAirVolFlow == AutoSize) {
            IsAutoSize = true;
        }

        if (state.dataSize->CurZoneEqNum > 0) {
            if (!IsAutoSize && !state.dataSize->ZoneSizingRunDone) { // simulation continue
                if (thisIndUnit.MaxTotAirVolFlow > 0.0) {
                    BaseSizer::reportSizerOutput(state,
                                                 thisIndUnit.UnitType,
                                                 thisIndUnit.Name,
                                                 "User-Specified Maximum Total Air Flow Rate [m3/s]",
                                                 thisIndUnit.MaxTotAirVolFlow);
                }
            } else {
                CheckZoneSizing(state, thisIndUnit.UnitType, thisIndUnit.Name);
                if (state.dataSize->CurTermUnitSizingNum > 0) {
                    MaxTotAirVolFlowDes = max(state.dataSize->TermUnitFinalZoneSizing(state.dataSize->CurTermUnitSizingNum).DesCoolVolFlow,
                                              state.dataSize->TermUnitFinalZoneSizing(state.dataSize->CurTermUnitSizingNum).DesHeatVolFlow);
                } else {
                    MaxTotAirVolFlowDes = 0.0;
                }
                if (MaxTotAirVolFlowDes < SmallAirVolFlow) {
                    MaxTotAirVolFlowDes = 0.0;
                }
                if (IsAutoSize) {
                    thisIndUnit.MaxTotAirVolFlow = MaxTotAirVolFlowDes;
                    BaseSizer::reportSizerOutput(state,
                                                 thisIndUnit.UnitType,
                                                 thisIndUnit.Name,
                                                 "Design Size Maximum Total Air Flow Rate [m3/s]",
                                                 MaxTotAirVolFlowDes);
                } else {
                    if (thisIndUnit.MaxTotAirVolFlow > 0.0 && MaxTotAirVolFlowDes > 0.0) {
                        MaxTotAirVolFlowUser = thisIndUnit.MaxTotAirVolFlow;
                        BaseSizer::reportSizerOutput(state,
                                                     thisIndUnit.UnitType,
                                                     thisIndUnit.Name,
                                                     "Design Size Maximum Total Air Flow Rate [m3/s]",
                                                     MaxTotAirVolFlowDes,
                                                     "User-Specified Maximum Total Air Flow Rate [m3/s]",
                                                     MaxTotAirVolFlowUser);
                        if (state.dataGlobal->DisplayExtraWarnings) {
                            if ((std::abs(MaxTotAirVolFlowDes - MaxTotAirVolFlowUser) / MaxTotAirVolFlowUser) >
                                state.dataSize->AutoVsHardSizingThreshold) {
                                ShowMessage(state,
                                            format("SizeHVACSingleDuctInduction: Potential issue with equipment sizing for {} = \"{}\".",
                                                   thisIndUnit.UnitType,
                                                   thisIndUnit.Name));
                                ShowContinueError(state, format("User-Specified Maximum Total Air Flow Rate of {:.5R} [m3/s]", MaxTotAirVolFlowUser));
                                ShowContinueError(
                                    state, format("differs from Design Size Maximum Total Air Flow Rate of {:.5R} [m3/s]", MaxTotAirVolFlowDes));
                                ShowContinueError(state, "This may, or may not, indicate mismatched component sizes.");
                                ShowContinueError(state, "Verify that the value entered is intended and is consistent with other components.");
                            }
                        }
                    }
                }
            }
        }

        IsAutoSize = false;
        if (thisIndUnit.MaxVolHotWaterFlow == AutoSize) {
            IsAutoSize = true;
        }
        if ((state.dataSize->CurZoneEqNum > 0) && (state.dataSize->CurTermUnitSizingNum > 0)) {
            if (!IsAutoSize && !state.dataSize->ZoneSizingRunDone) { // simulation continue
                if (thisIndUnit.MaxVolHotWaterFlow > 0.0) {
                    BaseSizer::reportSizerOutput(state,
                                                 thisIndUnit.UnitType,
                                                 thisIndUnit.Name,
                                                 "User-Specified Maximum Hot Water Flow Rate [m3/s]",
                                                 thisIndUnit.MaxVolHotWaterFlow);
                }
            } else {
                CheckZoneSizing(state, thisIndUnit.UnitType, thisIndUnit.Name);

                if (UtilityRoutines::SameString(thisIndUnit.HCoilType, "Coil:Heating:Water")) {

                    CoilWaterInletNode =
                        GetCoilWaterInletNode(state, "Coil:Heating:Water", thisIndUnit.HCoil, ErrorsFound);
                    CoilWaterOutletNode =
                        GetCoilWaterOutletNode(state, "Coil:Heating:Water", thisIndUnit.HCoil, ErrorsFound);
                    if (IsAutoSize) {
                        PltSizHeatNum = MyPlantSizingIndex(state,
                                                           "Coil:Heating:Water",
                                                           thisIndUnit.HCoil,
                                                           CoilWaterInletNode,
                                                           CoilWaterOutletNode,
                                                           ErrorsFound);
                        if (PltSizHeatNum > 0) {

                            auto &thisTermUnitFinalZoneSizing = state.dataSize->TermUnitFinalZoneSizing(state.dataSize->CurTermUnitSizingNum); 				
                            if (thisTermUnitFinalZoneSizing.DesHeatMassFlow >= SmallAirVolFlow) {
                                DesPriVolFlow = thisIndUnit.MaxTotAirVolFlow /
                                                (1.0 + thisIndUnit.InducRatio);
                                CpAir = PsyCpAirFnW(thisTermUnitFinalZoneSizing.HeatDesHumRat);
                                // the design heating coil load is the zone load minus whatever the central system does. Note that
                                // DesHeatCoilInTempTU is really the primary air inlet temperature for the unit.
                                if (thisTermUnitFinalZoneSizing.ZoneTempAtHeatPeak > 0.0) {
                                    DesCoilLoad =
                                        thisTermUnitFinalZoneSizing.NonAirSysDesHeatLoad -
                                        CpAir * RhoAir * DesPriVolFlow *
                                            (thisTermUnitFinalZoneSizing.DesHeatCoilInTempTU - thisTermUnitFinalZoneSizing.ZoneTempAtHeatPeak);
                                } else {
                                    DesCoilLoad = CpAir * RhoAir * DesPriVolFlow *
                                                  (state.dataSize->ZoneSizThermSetPtLo(state.dataSize->CurZoneEqNum) -
                                                   thisTermUnitFinalZoneSizing.DesHeatCoilInTempTU);
                                }
                                thisIndUnit.DesHeatingLoad = DesCoilLoad;
                                Cp = GetSpecificHeatGlycol(
                                    state,
                                    state.dataPlnt->PlantLoop(thisIndUnit.HWPlantLoc.loopNum).FluidName,
                                    DataGlobalConstants::HWInitConvTemp,
                                    state.dataPlnt->PlantLoop(thisIndUnit.HWPlantLoc.loopNum).FluidIndex,
                                    RoutineName);

                                rho = GetDensityGlycol(
                                    state,
                                    state.dataPlnt->PlantLoop(thisIndUnit.HWPlantLoc.loopNum).FluidName,
                                    DataGlobalConstants::HWInitConvTemp,
                                    state.dataPlnt->PlantLoop(thisIndUnit.HWPlantLoc.loopNum).FluidIndex,
                                    RoutineName);

                                MaxVolHotWaterFlowDes = DesCoilLoad / (state.dataSize->PlantSizData(PltSizHeatNum).DeltaT * Cp * rho);
                                MaxVolHotWaterFlowDes = max(MaxVolHotWaterFlowDes, 0.0);
                            } else {
                                MaxVolHotWaterFlowDes = 0.0;
                            }
                        } else {
                            ShowSevereError(state, "Autosizing of water flow requires a heating loop Sizing:Plant object");
                            ShowContinueError(state,
                                              format("Occurs in{} Object={}",
                                                     thisIndUnit.UnitType,
                                                     thisIndUnit.Name));
                            ErrorsFound = true;
                        }
                    }
                    if (IsAutoSize) {
                        thisIndUnit.MaxVolHotWaterFlow = MaxVolHotWaterFlowDes;
                        auto &thisTermUnitFinalZoneSizing = state.dataSize->TermUnitFinalZoneSizing(state.dataSize->CurTermUnitSizingNum);
			BaseSizer::reportSizerOutput(state,
                                                     thisIndUnit.UnitType,
                                                     thisIndUnit.Name,
                                                     "Design Size Maximum Hot Water Flow Rate [m3/s]",
                                                     MaxVolHotWaterFlowDes);

			BaseSizer::reportSizerOutput(
                            state,
                            thisIndUnit.UnitType,
                            thisIndUnit.Name,
                            "Design Size Inlet Air Temperature [C]",
                            thisTermUnitFinalZoneSizing.DesHeatCoilInTempTU);
                        BaseSizer::reportSizerOutput(
                            state,
                            thisIndUnit.UnitType,
                            thisIndUnit.Name,
                            "Design Size Inlet Air Humidity Ratio [kgWater/kgDryAir]",
                            thisTermUnitFinalZoneSizing.DesHeatCoilInHumRatTU);
                    } else {
                        if (thisIndUnit.MaxVolHotWaterFlow > 0.0 && MaxVolHotWaterFlowDes > 0.0) {
                            MaxVolHotWaterFlowUser = thisIndUnit.MaxVolHotWaterFlow;
                            BaseSizer::reportSizerOutput(state,
                                                         thisIndUnit.UnitType,
                                                         thisIndUnit.Name,
                                                         "Design Size Maximum Hot Water Flow Rate [m3/s]",
                                                         MaxVolHotWaterFlowDes,
                                                         "User-Specified Maximum Hot Water Flow Rate [m3/s]",
                                                         MaxVolHotWaterFlowUser);
                            if (state.dataGlobal->DisplayExtraWarnings) {
                                if ((std::abs(MaxVolHotWaterFlowDes - MaxVolHotWaterFlowUser) / MaxVolHotWaterFlowUser) >
                                    state.dataSize->AutoVsHardSizingThreshold) {
                                    ShowMessage(state,
                                                format("SizeHVACSingleDuctInduction: Potential issue with equipment sizing for {} = \"{}\".",
                                                       thisIndUnit.UnitType,
                                                       thisIndUnit.Name));
                                    ShowContinueError(state,
                                                      format("User-Specified Maximum Hot Water Flow Rate of {:.5R} [m3/s]", MaxVolHotWaterFlowUser));
                                    ShowContinueError(
                                        state,
                                        format("differs from Design Size Maximum Hot Water Flow Rate of {:.5R} [m3/s]", MaxVolHotWaterFlowDes));
                                    ShowContinueError(state, "This may, or may not, indicate mismatched component sizes.");
                                    ShowContinueError(state, "Verify that the value entered is intended and is consistent with other components.");
                                }
                            }
                        }
                    }
                } else {
                    thisIndUnit.MaxVolHotWaterFlow = 0.0;
                }
            }
        }

        IsAutoSize = false;
        if (thisIndUnit.MaxVolColdWaterFlow == AutoSize) {
            IsAutoSize = true;
        }
        if ((state.dataSize->CurZoneEqNum > 0) && (state.dataSize->CurTermUnitSizingNum > 0)) {
            if (!IsAutoSize && !state.dataSize->ZoneSizingRunDone) { // simulation continue
                if (thisIndUnit.MaxVolColdWaterFlow > 0.0) {
                    BaseSizer::reportSizerOutput(state,
                                                 thisIndUnit.UnitType,
                                                 thisIndUnit.Name,
                                                 "User-Specified Maximum Cold Water Flow Rate [m3/s]",
                                                 thisIndUnit.MaxVolColdWaterFlow);
                }
            } else {
                CheckZoneSizing(state, thisIndUnit.UnitType, thisIndUnit.Name);

                if (UtilityRoutines::SameString(thisIndUnit.CCoilType, "Coil:Cooling:Water") ||
                    UtilityRoutines::SameString(thisIndUnit.CCoilType, "Coil:Cooling:Water:DetailedGeometry")) {

                    CoilWaterInletNode = GetCoilWaterInletNode(state,
                                                               thisIndUnit.CCoilType,
                                                               thisIndUnit.CCoil,
                                                               ErrorsFound);
                    CoilWaterOutletNode = GetCoilWaterOutletNode(state,
                                                                 thisIndUnit.CCoilType,
                                                                 thisIndUnit.CCoil,
                                                                 ErrorsFound);
                    if (IsAutoSize) {
                        PltSizCoolNum = MyPlantSizingIndex(state,
                                                           thisIndUnit.CCoilType,
                                                           thisIndUnit.CCoil,
                                                           CoilWaterInletNode,
                                                           CoilWaterOutletNode,
                                                           ErrorsFound);
                        if (PltSizCoolNum > 0) {
                            auto &thisTermUnitFinalZoneSizing = state.dataSize->TermUnitFinalZoneSizing(state.dataSize->CurTermUnitSizingNum);
                            if (thisTermUnitFinalZoneSizing.DesCoolMassFlow >= SmallAirVolFlow) {
                                DesPriVolFlow = thisIndUnit.MaxTotAirVolFlow / (1.0 + thisIndUnit.InducRatio);
                                CpAir = PsyCpAirFnW(thisTermUnitFinalZoneSizing.CoolDesHumRat);
                                // the design cooling coil load is the zone load minus whatever the central system does. Note that
                                // DesCoolCoilInTempTU is really the primary air inlet temperature for the unit.
                                if (thisTermUnitFinalZoneSizing.ZoneTempAtCoolPeak > 0.0) {
                                    DesCoilLoad =
                                        thisTermUnitFinalZoneSizing.NonAirSysDesCoolLoad -
                                        CpAir * RhoAir * DesPriVolFlow *
                                            (thisTermUnitFinalZoneSizing.ZoneTempAtCoolPeak - thisTermUnitFinalZoneSizing.DesCoolCoilInTempTU);
                                } else {
                                    DesCoilLoad = CpAir * RhoAir * DesPriVolFlow *
                                                  (thisTermUnitFinalZoneSizing.DesCoolCoilInTempTU -
                                                   state.dataSize->ZoneSizThermSetPtHi(state.dataSize->CurZoneEqNum));
                                }
                                thisIndUnit.DesCoolingLoad = DesCoilLoad;
                                Cp = GetSpecificHeatGlycol(
                                    state,
                                    state.dataPlnt->PlantLoop(thisIndUnit.CWPlantLoc.loopNum).FluidName,
                                    5.0,
                                    state.dataPlnt->PlantLoop(thisIndUnit.CWPlantLoc.loopNum).FluidIndex,
                                    RoutineName);

                                rho = GetDensityGlycol(
                                    state,
                                    state.dataPlnt->PlantLoop(thisIndUnit.CWPlantLoc.loopNum).FluidName,
                                    5.0,
                                    state.dataPlnt->PlantLoop(thisIndUnit.CWPlantLoc.loopNum).FluidIndex,
                                    RoutineName);

                                MaxVolColdWaterFlowDes = DesCoilLoad / (state.dataSize->PlantSizData(PltSizCoolNum).DeltaT * Cp * rho);
                                MaxVolColdWaterFlowDes = max(MaxVolColdWaterFlowDes, 0.0);
                            } else {
                                MaxVolColdWaterFlowDes = 0.0;
                            }
                        } else {
                            ShowSevereError(state, "Autosizing of water flow requires a cooling loop Sizing:Plant object");
                            ShowContinueError(state,
                                              format("Occurs in{} Object={}",
                                                     thisIndUnit.UnitType,
                                                     thisIndUnit.Name));
                            ErrorsFound = true;
                        }
                    }
                    if (IsAutoSize) {
                        thisIndUnit.MaxVolColdWaterFlow = MaxVolColdWaterFlowDes;
                        BaseSizer::reportSizerOutput(state,
                                                     thisIndUnit.UnitType,
                                                     thisIndUnit.Name,
                                                     "Design Size Maximum Cold Water Flow Rate [m3/s]",
                                                     MaxVolColdWaterFlowDes);
                    } else {
                        if (thisIndUnit.MaxVolColdWaterFlow > 0.0 && MaxVolColdWaterFlowDes > 0.0) {
                            MaxVolColdWaterFlowUser = thisIndUnit.MaxVolColdWaterFlow;
                            BaseSizer::reportSizerOutput(state,
                                                         thisIndUnit.UnitType,
                                                         thisIndUnit.Name,
                                                         "Design Size Maximum Cold Water Flow Rate [m3/s]",
                                                         MaxVolColdWaterFlowDes,
                                                         "User-Specified Maximum Cold Water Flow Rate [m3/s]",
                                                         MaxVolColdWaterFlowUser);
                            if (state.dataGlobal->DisplayExtraWarnings) {
                                if ((std::abs(MaxVolColdWaterFlowDes - MaxVolColdWaterFlowUser) / MaxVolColdWaterFlowUser) >
                                    state.dataSize->AutoVsHardSizingThreshold) {
                                    ShowMessage(state,
                                                format("SizeHVACSingleDuctInduction: Potential issue with equipment sizing for {} = \"{}\".",
                                                       thisIndUnit.UnitType,
                                                       thisIndUnit.Name));
                                    ShowContinueError(
                                        state, format("User-Specified Maximum Cold Water Flow Rate of {:.5R} [m3/s]", MaxVolColdWaterFlowUser));
                                    ShowContinueError(
                                        state,
                                        format("differs from Design Size Maximum Cold Water Flow Rate of {:.5R} [m3/s]", MaxVolColdWaterFlowDes));
                                    ShowContinueError(state, "This may, or may not, indicate mismatched component sizes.");
                                    ShowContinueError(state, "Verify that the value entered is intended and is consistent with other components.");
                                }
                            }
                        }
                    }
                } else {
                    thisIndUnit.MaxVolColdWaterFlow = 0.0;
                }
            }
        }

        if (state.dataSize->CurTermUnitSizingNum > 0) {
	    auto &thisTermUnitSizing = TermUnitSizing(state.dataSize->CurTermUnitSizingNum);
		// note we save the induced air flow for use by the hw and cw coil sizing routines
            thisTermUnitSizing.AirVolFlow = thisIndUnit.MaxTotAirVolFlow * thisIndUnit.InducRatio / (1.0 + thisIndUnit.InducRatio);
            // save the max hot and cold water flows for use in coil sizing
            thisTermUnitSizing.MaxHWVolFlow = thisIndUnit.MaxVolHotWaterFlow;
            thisTermUnitSizing.MaxCWVolFlow = thisIndUnit.MaxVolColdWaterFlow;
            // save the design load used for reporting
            thisTermUnitSizing.DesCoolingLoad = thisIndUnit.DesCoolingLoad;
            thisTermUnitSizing.DesHeatingLoad = thisIndUnit.DesHeatingLoad;
            // save the induction ratio for use in subsequent sizing calcs
            thisTermUnitSizing.InducRat = thisIndUnit.InducRatio;
            if (UtilityRoutines::SameString(thisIndUnit.HCoilType, "Coil:Heating:Water")) {
                SetCoilDesFlow(state,
                               thisIndUnit.HCoilType,
                               thisIndUnit.HCoil,
                               thisTermUnitSizing.AirVolFlow,
                               ErrorsFound);
            }
            if (UtilityRoutines::SameString(thisIndUnit.CCoilType, "Coil:Cooling:Water:DetailedGeometry")) {
                SetCoilDesFlow(state,
                               thisIndUnit.CCoilType,
                               thisIndUnit.CCoil,
                               thisTermUnitSizing.AirVolFlow,
                               ErrorsFound);
            }
        }
    }

    void SimFourPipeIndUnit(EnergyPlusData &state,
                            int const IUNum,              // number of the current unit being simulated
                            int const ZoneNum,            // number of zone being served
                            int const ZoneNodeNum,        // zone node number
                            bool const FirstHVACIteration // TRUE if 1st HVAC simulation of system timestep
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Fred Buhl
        //       DATE WRITTEN   June 23 2004
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // Simulate a 4 pipe induction unit; adjust its heating or cooling
        // coil outputs to match the zone load.

        // METHODOLOGY EMPLOYED:
        // (1) From the zone load and the primary air inlet conditions calculate the coil load
        //     in the secondary air stream
        // (2) If there is a cooling coil load, set the heating coil off and control the cooling
        //     coil to meet the coil load
        // (3) If there is a heating coil load, control the heating coil to meet the load and keep
        //     the cooling coil off.

        // Using/Aliasing
        using namespace DataZoneEnergyDemands;

        using General::SolveRoot;
        using PlantUtilities::SetComponentFlowRate;

        // SUBROUTINE PARAMETER DEFINITIONS:
        int constexpr SolveMaxIter(50);

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 QZnReq;           // heating or cooling needed by zone [Watts]
        Real64 QToHeatSetPt;     // [W]  remaining load to heating setpoint
        Real64 QToCoolSetPt;     // [W]  remaining load to cooling setpoint
        Real64 PowerMet;         // power supplied
        bool UnitOn;             // TRUE if unit is on
        Real64 MaxHotWaterFlow;  // maximum water flow for heating [kg/s]
        Real64 MinHotWaterFlow;  // minimum water flow for heating [kg/s]
        Real64 MaxColdWaterFlow; // maximum water flow for cooling [kg/s]
        Real64 MinColdWaterFlow; // minimum water flow for cooling [kg/s]
        Real64 HWFlow;           // hot water flow [kg/s]
        Real64 CWFlow;           // cold water flow [kg/s]
        int PriNode;             // unit primary air inlet node
        int SecNode;             // unit secondary air inlet node
        int OutletNode;          // unit air outlet node
        int HotControlNode;      // hot water coil inlet node
        int ColdControlNode;     // cold water coil inlet node
        Real64 QPriOnly;         // unit output with no zone coils active
        Real64 PriAirMassFlow;   // primary air mass flow rate [kg/s]
        Real64 SecAirMassFlow;   // secondary air mass flow rate [kg/s]
        Real64 InducRat;         // Induction Ratio
        int SolFlag;
        Real64 ErrTolerance;
        int HWOutletNode;
        int CWOutletNode;
        UnitOn = true;
        PowerMet = 0.0;

        auto &thisIndUnit = state.dataHVACSingleDuctInduc->IndUnit(IUNum);

	InducRat = thisIndUnit.InducRatio;
        PriNode = thisIndUnit.PriAirInNode;
        SecNode = thisIndUnit.SecAirInNode;
        OutletNode = thisIndUnit.OutAirNode;
        HotControlNode = thisIndUnit.HWControlNode;
        HWOutletNode = DataPlant::CompData::getPlantComponent(state, thisIndUnit.HWPlantLoc).NodeNumOut;
        ColdControlNode = thisIndUnit.CWControlNode;
        CWOutletNode = DataPlant::CompData::getPlantComponent(state, thisIndUnit.CWPlantLoc).NodeNumOut;
        PriAirMassFlow = state.dataLoopNodes->Node(PriNode).MassFlowRateMaxAvail;
        SecAirMassFlow = InducRat * PriAirMassFlow;
        QZnReq = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(ZoneNum).RemainingOutputRequired;
        QToHeatSetPt = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(ZoneNum).RemainingOutputReqToHeatSP;
        QToCoolSetPt = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(ZoneNum).RemainingOutputReqToCoolSP;
        // On the first HVAC iteration the system values are given to the controller, but after that
        // the demand limits are in place and there needs to be feedback to the Zone Equipment

        MaxHotWaterFlow = thisIndUnit.MaxHotWaterFlow;
        SetComponentFlowRate(state, MaxHotWaterFlow, HotControlNode, HWOutletNode, thisIndUnit.HWPlantLoc);

        MinHotWaterFlow = thisIndUnit.MinHotWaterFlow;
        SetComponentFlowRate(state, MinHotWaterFlow, HotControlNode, HWOutletNode, thisIndUnit.HWPlantLoc);

        MaxColdWaterFlow = thisIndUnit.MaxColdWaterFlow;
        SetComponentFlowRate(state, MaxColdWaterFlow, ColdControlNode, CWOutletNode, thisIndUnit.CWPlantLoc);

        MinColdWaterFlow = thisIndUnit.MinColdWaterFlow;
        SetComponentFlowRate(state, MinColdWaterFlow, ColdControlNode, CWOutletNode, thisIndUnit.CWPlantLoc);

        if (GetCurrentScheduleValue(state, thisIndUnit.SchedPtr) <= 0.0) UnitOn = false;
        if (PriAirMassFlow <= SmallMassFlow) UnitOn = false;

        // Set the unit's air inlet nodes mass flow rates
        state.dataLoopNodes->Node(PriNode).MassFlowRate = PriAirMassFlow;
        state.dataLoopNodes->Node(SecNode).MassFlowRate = SecAirMassFlow;
        // initialize the water inlet nodes to minimum
        // fire the unit at min water flow
        CalcFourPipeIndUnit(state, IUNum, FirstHVACIteration, ZoneNodeNum, MinHotWaterFlow, MinColdWaterFlow, QPriOnly);
        // the load to be met by the secondary air stream coils is QZnReq-PowerMet

        if (UnitOn) {

            if (QToHeatSetPt - QPriOnly > SmallLoad) {
                // heating coil
                // check that it can meet the load
                CalcFourPipeIndUnit(state, IUNum, FirstHVACIteration, ZoneNodeNum, MaxHotWaterFlow, MinColdWaterFlow, PowerMet);
                if (PowerMet > QToHeatSetPt + SmallLoad) {
                    ErrTolerance = thisIndUnit.HotControlOffset;
                    auto f = // (THIS_AUTO_OK)
                        [&state, IUNum, FirstHVACIteration, ZoneNodeNum, MinColdWaterFlow, QToHeatSetPt, QPriOnly, PowerMet](Real64 const HWFlow) {
                            Real64 UnitOutput;
                            CalcFourPipeIndUnit(state, IUNum, FirstHVACIteration, ZoneNodeNum, HWFlow, MinColdWaterFlow, UnitOutput);
                            return (QToHeatSetPt - UnitOutput) / (PowerMet - QPriOnly);
                        };
                    SolveRoot(state, ErrTolerance, SolveMaxIter, SolFlag, HWFlow, f, MinHotWaterFlow, MaxHotWaterFlow);
                    if (SolFlag == -1) {
                        if (thisIndUnit.HWCoilFailNum1 == 0) {
                            ShowWarningMessage(state,
                                               format("SimFourPipeIndUnit: Hot water coil control failed for {}=\"{}\"",
                                                      thisIndUnit.UnitType,
                                                      thisIndUnit.Name));
                            ShowContinueErrorTimeStamp(state, "");
                            ShowContinueError(state, format("  Iteration limit [{}] exceeded in calculating hot water mass flow rate", SolveMaxIter));
                        }
                        ShowRecurringWarningErrorAtEnd(
                            state,
                            format("SimFourPipeIndUnit: Hot water coil control failed (iteration limit [{}]) for {}=\"{}\"",
                                   SolveMaxIter,
                                   thisIndUnit.UnitType,
                                   thisIndUnit.Name),
                            thisIndUnit.HWCoilFailNum1);
                    } else if (SolFlag == -2) {
                        if (thisIndUnit.HWCoilFailNum2 == 0) {
                            ShowWarningMessage(state,
                                               format("SimFourPipeIndUnit: Hot water coil control failed (maximum flow limits) for {}=\"{}\"",
                                                      thisIndUnit.UnitType,
                                                      thisIndUnit.Name));
                            ShowContinueErrorTimeStamp(state, "");
                            ShowContinueError(state, "...Bad hot water maximum flow rate limits");
                            ShowContinueError(state, format("...Given minimum water flow rate={:.3R} kg/s", MinHotWaterFlow));
                            ShowContinueError(state, format("...Given maximum water flow rate={:.3R} kg/s", MaxHotWaterFlow));
                        }
                        ShowRecurringWarningErrorAtEnd(state,
                                                       "SimFourPipeIndUnit: Hot water coil control failed (flow limits) for " +
                                                           thisIndUnit.UnitType + "=\"" +
                                                           thisIndUnit.Name + "\"",
                                                       thisIndUnit.HWCoilFailNum2,
                                                       MaxHotWaterFlow,
                                                       MinHotWaterFlow,
                                                       _,
                                                       "[kg/s]",
                                                       "[kg/s]");
                    }
                }
            } else if (QToCoolSetPt - QPriOnly < -SmallLoad) {
                // cooling coil
                // check that it can meet the load
                CalcFourPipeIndUnit(state, IUNum, FirstHVACIteration, ZoneNodeNum, MinHotWaterFlow, MaxColdWaterFlow, PowerMet);
                if (PowerMet < QToCoolSetPt - SmallLoad) {
                    ErrTolerance = thisIndUnit.ColdControlOffset;
                    auto f = // (THIS_AUTO_OK)
                        [&state, IUNum, FirstHVACIteration, ZoneNodeNum, MinHotWaterFlow, QToCoolSetPt, QPriOnly, PowerMet](Real64 const CWFlow) {
                            Real64 UnitOutput;
                            CalcFourPipeIndUnit(state, IUNum, FirstHVACIteration, ZoneNodeNum, MinHotWaterFlow, CWFlow, UnitOutput);
                            return (QToCoolSetPt - UnitOutput) / (PowerMet - QPriOnly);
                        };
                    SolveRoot(state, ErrTolerance, SolveMaxIter, SolFlag, CWFlow, f, MinColdWaterFlow, MaxColdWaterFlow);
                    if (SolFlag == -1) {
                        if (thisIndUnit.CWCoilFailNum1 == 0) {
                            ShowWarningMessage(state,
                                               format("SimFourPipeIndUnit: Cold water coil control failed for {}=\"{}\"",
                                                      thisIndUnit.UnitType,
                                                      thisIndUnit.Name));
                            ShowContinueErrorTimeStamp(state, "");
                            ShowContinueError(state,
                                              format("  Iteration limit [{}] exceeded in calculating cold water mass flow rate", SolveMaxIter));
                        }
                        ShowRecurringWarningErrorAtEnd(state,
                                                       format("SimFourPipeIndUnit: Cold water coil control failed (iteration limit [{}]) for {}=\"{}",
                                                              SolveMaxIter,
                                                              thisIndUnit.UnitType,
                                                              thisIndUnit.Name),
                                                       thisIndUnit.CWCoilFailNum1);
                    } else if (SolFlag == -2) {
                        if (thisIndUnit.CWCoilFailNum2 == 0) {
                            ShowWarningMessage(state,
                                               format("SimFourPipeIndUnit: Cold water coil control failed (maximum flow limits) for {}=\"{}\"",
                                                      thisIndUnit.UnitType,
                                                      thisIndUnit.Name));
                            ShowContinueErrorTimeStamp(state, "");
                            ShowContinueError(state, "...Bad cold water maximum flow rate limits");
                            ShowContinueError(state, format("...Given minimum water flow rate={:.3R} kg/s", MinColdWaterFlow));
                            ShowContinueError(state, format("...Given maximum water flow rate={:.3R} kg/s", MaxColdWaterFlow));
                        }
                        ShowRecurringWarningErrorAtEnd(state,
                                                       "SimFourPipeIndUnit: Cold water coil control failed (flow limits) for " +
                                                           thisIndUnit.UnitType + "=\"" +
                                                           thisIndUnit.Name + "\"",
                                                       thisIndUnit.CWCoilFailNum2,
                                                       MaxColdWaterFlow,
                                                       MinColdWaterFlow,
                                                       _,
                                                       "[kg/s]",
                                                       "[kg/s]");
                    }
                }
            } else {
                CalcFourPipeIndUnit(state, IUNum, FirstHVACIteration, ZoneNodeNum, MinHotWaterFlow, MinColdWaterFlow, PowerMet);
            }

        } else {
            // unit off
            CalcFourPipeIndUnit(state, IUNum, FirstHVACIteration, ZoneNodeNum, MinHotWaterFlow, MinColdWaterFlow, PowerMet);
        }
        state.dataLoopNodes->Node(OutletNode).MassFlowRateMax = thisIndUnit.MaxTotAirMassFlow;

        // At this point we are done. There is no output to report or pass back up: the output provided is calculated
        // one level up in the calling routine SimZoneAirLoopEquipment. All the inlet and outlet flow rates and
        // conditions have been set by CalcFourPipeIndUnit either explicitly or as a result of the simple component calls.
    }

    void CalcFourPipeIndUnit(EnergyPlusData &state,
                             int const IUNum,               // Unit index
                             bool const FirstHVACIteration, // flag for 1st HVAV iteration in the time step
                             int const ZoneNode,            // zone node number
                             Real64 const HWFlow,           // hot water flow (kg/s)
                             Real64 const CWFlow,           // cold water flow (kg/s)
                             Real64 &LoadMet                // load met by unit (watts)
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Fred Buhl
        //       DATE WRITTEN   June 2004
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // Simulate the components making up the 4 pipe induction unit.

        // METHODOLOGY EMPLOYED:
        // Simulates the unit components sequentially in the air flow direction.

        // REFERENCES:
        // na

        // Using/Aliasing
        using MixerComponent::SimAirMixer;
        using PlantUtilities::SetComponentFlowRate;
        using WaterCoils::SimulateWaterCoilComponents;

        // Locals
        // SUBROUTINE ARGUMENT DEFINITIONS:

        // SUBROUTINE PARAMETER DEFINITIONS:
        // na

        // INTERFACE BLOCK SPECIFICATIONS
        // na

        // DERIVED TYPE DEFINITIONS
        // na

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int OutletNode;        // unit air outlet node
        int PriNode;           // unit primary air inlet node
        int HotControlNode;    // the hot water inlet node
        int ColdControlNode;   // the cold water inlet node
        Real64 PriAirMassFlow; // primary air mass flow rate [kg/s]
        Real64 SecAirMassFlow; // secondary air mass flow rate [kg/s]
        Real64 TotAirMassFlow; // total air mass flow rate [kg/s]
        Real64 InducRat;       // induction ratio
        Real64 mdotHW;         // local temporary hot water flow rate [kg/s]
        Real64 mdotCW;         // local temporary cold water flow rate [kg/s]
        int HWOutletNode;
        int CWOutletNode;

        auto &thisIndUnit = state.dataHVACSingleDuctInduc->IndUnit(IUNum);

        PriNode = thisIndUnit.PriAirInNode;
        OutletNode = thisIndUnit.OutAirNode;
        PriAirMassFlow = state.dataLoopNodes->Node(PriNode).MassFlowRateMaxAvail;
        InducRat = thisIndUnit.InducRatio;
        SecAirMassFlow = InducRat * PriAirMassFlow;
        TotAirMassFlow = PriAirMassFlow + SecAirMassFlow;
        HotControlNode = thisIndUnit.HWControlNode;
        HWOutletNode = DataPlant::CompData::getPlantComponent(state, thisIndUnit.HWPlantLoc).NodeNumOut;

        ColdControlNode = thisIndUnit.CWControlNode;
        CWOutletNode = DataPlant::CompData::getPlantComponent(state, thisIndUnit.CWPlantLoc).NodeNumOut;

        mdotHW = HWFlow;
        SetComponentFlowRate(state, mdotHW, HotControlNode, HWOutletNode, thisIndUnit.HWPlantLoc);

        //  Node(HotControlNode)%MassFlowRate = HWFlow

        mdotCW = CWFlow;
        SetComponentFlowRate(state, mdotCW, ColdControlNode, CWOutletNode, thisIndUnit.CWPlantLoc);
        //  Node(ColdControlNode)%MassFlowRate = CWFlow

        SimulateWaterCoilComponents(state, thisIndUnit.HCoil, FirstHVACIteration, thisIndUnit.HCoil_Num);
        SimulateWaterCoilComponents(state, thisIndUnit.CCoil, FirstHVACIteration, thisIndUnit.CCoil_Num);
        SimAirMixer(state, thisIndUnit.MixerName, thisIndUnit.Mixer_Num);
        LoadMet = TotAirMassFlow * Psychrometrics::PsyDeltaHSenFnTdb2W2Tdb1W1(state.dataLoopNodes->Node(OutletNode).Temp,
                                                                              state.dataLoopNodes->Node(OutletNode).HumRat,
                                                                              state.dataLoopNodes->Node(ZoneNode).Temp,
                                                                              state.dataLoopNodes->Node(ZoneNode).HumRat);
    }

    bool FourPipeInductionUnitHasMixer(EnergyPlusData &state, std::string_view CompName) // component (mixer) name
    {

        // FUNCTION INFORMATION:
        //       AUTHOR         Linda Lawrie
        //       DATE WRITTEN   September 2011
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS FUNCTION:
        // Given a mixer name, this routine determines if that mixer is found on
        // PIUnits.

        if (state.dataHVACSingleDuctInduc->GetIUInputFlag) {
            GetIndUnits(state);
            state.dataHVACSingleDuctInduc->GetIUInputFlag = false;
        }

        if (state.dataHVACSingleDuctInduc->NumIndUnits == 0)
            return false;
	
        return (UtilityRoutines::FindItemInList(CompName, state.dataHVACSingleDuctInduc->IndUnit, &IndUnitData::MixerName) > 0);
    }

    void IndUnitData::ReportIndUnit(EnergyPlusData &state)
    {
        // Purpose: this subroutine for reporting

        // set zone OA volume flow rate
        this->CalcOutdoorAirVolumeFlowRate(state);
    }

    void IndUnitData::CalcOutdoorAirVolumeFlowRate(EnergyPlusData &state)
    {
        // calculates zone outdoor air volume flow rate using the supply air flow rate and OA fraction
        if (this->AirLoopNum > 0) {
            this->OutdoorAirFlowRate = (state.dataLoopNodes->Node(this->PriAirInNode).MassFlowRate / state.dataEnvrn->StdRhoAir) *
                                       state.dataAirLoop->AirLoopFlow(this->AirLoopNum).OAFrac;
        } else {
            this->OutdoorAirFlowRate = 0.0;
        }
    }

} // namespace HVACSingleDuctInduc

} // namespace EnergyPlus

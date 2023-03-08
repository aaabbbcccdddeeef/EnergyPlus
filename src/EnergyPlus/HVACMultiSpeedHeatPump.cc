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
#include <AirflowNetwork/Solver.hpp>
#include <EnergyPlus/Autosizing/Base.hh>
#include <EnergyPlus/BranchNodeConnections.hh>
#include <EnergyPlus/DXCoils.hh>
#include <EnergyPlus/Data/EnergyPlusData.hh>
#include <EnergyPlus/DataAirSystems.hh>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/DataLoopNode.hh>
#include <EnergyPlus/DataSizing.hh>
#include <EnergyPlus/DataZoneControls.hh>
#include <EnergyPlus/DataZoneEnergyDemands.hh>
#include <EnergyPlus/DataZoneEquipment.hh>
#include <EnergyPlus/EMSManager.hh>
#include <EnergyPlus/Fans.hh>
#include <EnergyPlus/FluidProperties.hh>
#include <EnergyPlus/General.hh>
#include <EnergyPlus/GeneralRoutines.hh>
#include <EnergyPlus/HVACMultiSpeedHeatPump.hh>
#include <EnergyPlus/HeatingCoils.hh>
#include <EnergyPlus/InputProcessing/InputProcessor.hh>
#include <EnergyPlus/NodeInputManager.hh>
#include <EnergyPlus/OutputProcessor.hh>
#include <EnergyPlus/Plant/DataPlant.hh>
#include <EnergyPlus/PlantUtilities.hh>
#include <EnergyPlus/Psychrometrics.hh>
#include <EnergyPlus/ScheduleManager.hh>
#include <EnergyPlus/SteamCoils.hh>
#include <EnergyPlus/UtilityRoutines.hh>
#include <EnergyPlus/WaterCoils.hh>
#include <EnergyPlus/ZoneTempPredictorCorrector.hh>

namespace EnergyPlus {

namespace HVACMultiSpeedHeatPump {

    // Module containing the Multi Speed Heat Pump simulation routines

    // MODULE INFORMATION:
    //       AUTHOR         Lixing Gu, Florida Solar Energy Center
    //       DATE WRITTEN   June 2007
    //       MODIFIED       Bereket Nigusse, FSEC, June 2010 - deprecated supply air flow fraction through controlled
    //                      zone from the furnace object input field. Now, the flow fraction is calculated internally
    //                      Brent Griffith, NREL, Dec 2010 -- upgrade to new plant for heat recovery, general fluid props.
    //                      Bereket Nigusse, FSEC, Jan. 2012 -- added hot water and steam heating coil

    //       RE-ENGINEERED  na

    // PURPOSE OF THIS MODULE:
    // To encapsulate the data and algorithms required to simulate Multi Speed Heat Pump in
    // EnergyPlus.

    // Module currently models air-cooled or evap-cooled direct expansion systems
    // (split or packaged) with mulptiple speeds. Air-side performance is modeled to determine
    // coil discharge air conditions. The module also determines the DX unit's energy
    // usage. Neither the air-side performance nor the energy usage includes the effect
    // of supply air fan heat/energy usage. The supply air fan is modeled by other modules.

    // USE STATEMENTS:
    // Use statements for data only modules
    // Using/Aliasing
    using namespace DataLoopNode;
    using DataHVACGlobals::BlowThru;
    using DataHVACGlobals::Coil_HeatingElectric;
    using DataHVACGlobals::Coil_HeatingElectric_MultiStage;
    using DataHVACGlobals::Coil_HeatingGas_MultiStage;
    using DataHVACGlobals::Coil_HeatingGasOrOtherFuel;
    using DataHVACGlobals::Coil_HeatingSteam;
    using DataHVACGlobals::Coil_HeatingWater;
    using DataHVACGlobals::ContFanCycCoil;
    using DataHVACGlobals::CycFanCycCoil;
    using DataHVACGlobals::DrawThru;
    using DataHVACGlobals::SmallAirVolFlow;
    using DataHVACGlobals::SmallLoad;
    using DataHVACGlobals::SmallMassFlow;
    using namespace Psychrometrics;
    using namespace ScheduleManager;

    // Curve Types
    enum class CurveType
    {
        Invalid = -1,
        Linear,      // Linear curve type
        BiLinear,    // Bi-linear curve type
        Quadratic,   // Quadratic curve type
        BiQuadratic, // Bi-quadratic curve type
        Cubic,       // Cubic curve type
        Num
    };

    static constexpr std::string_view fluidNameSteam("STEAM");

    void SimMSHeatPump(EnergyPlusData &state,
                       std::string_view CompName,     // Name of the unitary engine driven heat pump system
                       bool const FirstHVACIteration, // TRUE if 1st HVAC simulation of system time step
                       int const AirLoopNum,          // air loop index
                       int &CompIndex                 // Index to changeover-bypass VAV system
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Lixing Gu, Florida Solar Energy Center
        //       DATE WRITTEN   June. 2007
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // Manages the simulation of multispeed heat pump.

        // Using/Aliasing

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int MSHeatPumpNum;        // index of fan coil unit being simulated
        Real64 OnOffAirFlowRatio; // Ratio of compressor ON airflow to average airflow over timestep
        Real64 QZnLoad;           // Zone load required by all zones served by this air loop system
        Real64 QSensUnitOut;      // MSHP sensible capacity output [W]

        // First time SimMSHeatPump is called, get the input
        if (state.dataHVACMultiSpdHP->GetInputFlag) {
            GetMSHeatPumpInput(state);
            state.dataHVACMultiSpdHP->GetInputFlag = false; // Set GetInputFlag false so you don't get coil inputs again
        }

        if (CompIndex == 0) {
            MSHeatPumpNum = UtilityRoutines::FindItemInList(CompName, state.dataHVACMultiSpdHP->MSHeatPump);
            if (MSHeatPumpNum == 0) {
                ShowFatalError(state, format("MultiSpeed Heat Pump is not found={}", CompName));
            }
            CompIndex = MSHeatPumpNum;
        } else {
            MSHeatPumpNum = CompIndex;
            if (MSHeatPumpNum > state.dataHVACMultiSpdHP->NumMSHeatPumps || MSHeatPumpNum < 1) {
                ShowFatalError(state,
                               format("SimMSHeatPump: Invalid CompIndex passed={}, Number of MultiSpeed Heat Pumps={}, Heat Pump name={}",
                                      MSHeatPumpNum,
                                      state.dataHVACMultiSpdHP->NumMSHeatPumps,
                                      CompName));
            }
            if (state.dataHVACMultiSpdHP->CheckEquipName(MSHeatPumpNum)) {
                if (CompName != state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum).Name) {
                    ShowFatalError(state,
                                   format("SimMSHeatPump: Invalid CompIndex passed={}, Heat Pump name={}{}",
                                          MSHeatPumpNum,
                                          CompName,
                                          state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum).Name));
                }
                state.dataHVACMultiSpdHP->CheckEquipName(MSHeatPumpNum) = false;
            }
        }

        OnOffAirFlowRatio = 0.0;

        // Initialize the engine driven heat pump
        InitMSHeatPump(state, MSHeatPumpNum, FirstHVACIteration, AirLoopNum, QZnLoad, OnOffAirFlowRatio);

        SimMSHP(state, MSHeatPumpNum, FirstHVACIteration, AirLoopNum, QSensUnitOut, QZnLoad, OnOffAirFlowRatio);

        // Update the unit outlet nodes
        UpdateMSHeatPump(state, MSHeatPumpNum);

        // Report the result of the simulation
        ReportMSHeatPump(state, MSHeatPumpNum);
    }

    //******************************************************************************

    void SimMSHP(EnergyPlusData &state,
                 int const MSHeatPumpNum,       // number of the current engine driven Heat Pump being simulated
                 bool const FirstHVACIteration, // TRUE if 1st HVAC simulation of system timestep
                 int const AirLoopNum,          // air loop index
                 Real64 &QSensUnitOut,          // cooling/heating deliveded to zones [W]
                 Real64 const QZnReq,           // required zone load
                 Real64 &OnOffAirFlowRatio      // ratio of compressor ON airflow to AVERAGE airflow over timestep
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Lixing Gu
        //       DATE WRITTEN   June 2007
        //       MODIFIED       na
        //       RE-ENGINEERED  Revised based on SimPTHP

        // PURPOSE OF THIS SUBROUTINE:
        // Simulate a multispeed heat pump; adjust its output to match the
        // required system load.

        // METHODOLOGY EMPLOYED:
        // Calls ControlMSHPOutput to obtain the desired unit output

        // Using/Aliasing
        using namespace DataZoneEnergyDemands;
        using DataHVACGlobals::SmallLoad;
        using DataHVACGlobals::SmallMassFlow;

        // Locals
        Real64 SupHeaterLoad;

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 PartLoadFrac; // compressor part load fraction
        Real64 SpeedRatio;   // compressor speed ratio
        Real64 QTotUnitOut;
        int SpeedNum;                                      // Speed number
        Real64 SaveMassFlowRate;                           // saved inlet air mass flow rate [kg/s]

        // zero the fan, DX coils, and supplemental electric heater electricity consumption
        state.dataHVACGlobal->DXElecHeatingPower = 0.0;
        state.dataHVACGlobal->DXElecCoolingPower = 0.0;
        state.dataHVACMultiSpdHP->SaveCompressorPLR = 0.0;
        state.dataHVACGlobal->ElecHeatingCoilPower = 0.0;
        state.dataHVACGlobal->SuppHeatingCoilPower = 0.0;

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        // initialize local variables
        bool UnitOn = true;
        int OutletNode = thisMSHeatPump.AirOutletNodeNum;
        int InletNode = thisMSHeatPump.AirInletNodeNum;
        Real64 AirMassFlow = state.dataLoopNodes->Node(InletNode).MassFlowRate;
        int OpMode = thisMSHeatPump.OpMode;
        int ZoneNum = thisMSHeatPump.ControlZoneNum;
        DataHVACGlobals::CompressorOperation CompressorOp = DataHVACGlobals::CompressorOperation::On;

        // set the on/off flags
        if (thisMSHeatPump.OpMode == CycFanCycCoil) {
            // cycling unit only runs if there is a cooling or heating load.
            if (std::abs(QZnReq) < SmallLoad || AirMassFlow < SmallMassFlow || state.dataZoneEnergyDemand->CurDeadBandOrSetback(ZoneNum)) {
                UnitOn = false;
            }
        } else if (thisMSHeatPump.OpMode == ContFanCycCoil) {
            // continuous unit: fan runs if scheduled on; coil runs only if there is a cooling or heating load
            if (AirMassFlow < SmallMassFlow) {
                UnitOn = false;
            }
        }

        state.dataHVACGlobal->OnOffFanPartLoadFraction = 1.0;

        SaveMassFlowRate = state.dataLoopNodes->Node(InletNode).MassFlowRate;
        if (thisMSHeatPump.EMSOverrideCoilSpeedNumOn) {
            Real64 SpeedVal = thisMSHeatPump.EMSOverrideCoilSpeedNumValue;

            if (!FirstHVACIteration && thisMSHeatPump.OpMode == CycFanCycCoil && QZnReq < 0.0 &&
                state.dataAirLoop->AirLoopControlInfo(AirLoopNum).EconoActive) {
                CompressorOp = DataHVACGlobals::CompressorOperation::Off;
                ControlMSHPOutputEMS(state,
                                     MSHeatPumpNum,
                                     FirstHVACIteration,
                                     CompressorOp,
                                     OpMode,
                                     QZnReq,
                                     SpeedVal,
                                     SpeedNum,
                                     SpeedRatio,
                                     PartLoadFrac,
                                     OnOffAirFlowRatio,
                                     SupHeaterLoad);
                if (ceil(SpeedVal) == thisMSHeatPump.NumOfSpeedCooling && SpeedRatio == 1.0) {
                    state.dataLoopNodes->Node(InletNode).MassFlowRate = SaveMassFlowRate;
                    CompressorOp = DataHVACGlobals::CompressorOperation::On;
                    ControlMSHPOutputEMS(state,
                                         MSHeatPumpNum,
                                         FirstHVACIteration,
                                         CompressorOp,
                                         OpMode,
                                         QZnReq,
                                         SpeedVal,
                                         SpeedNum,
                                         SpeedRatio,
                                         PartLoadFrac,
                                         OnOffAirFlowRatio,
                                         SupHeaterLoad);
                }
            } else {
                ControlMSHPOutputEMS(state,
                                     MSHeatPumpNum,
                                     FirstHVACIteration,
                                     CompressorOp,
                                     OpMode,
                                     QZnReq,
                                     SpeedVal,
                                     SpeedNum,
                                     SpeedRatio,
                                     PartLoadFrac,
                                     OnOffAirFlowRatio,
                                     SupHeaterLoad);
            }
        } else {
            if (!FirstHVACIteration && thisMSHeatPump.OpMode == CycFanCycCoil && QZnReq < 0.0 &&
                state.dataAirLoop->AirLoopControlInfo(AirLoopNum).EconoActive) {
                // for cycling fan, cooling load, check whether furnace can meet load with compressor off
                CompressorOp = DataHVACGlobals::CompressorOperation::Off;
                ControlMSHPOutput(state,
                                  MSHeatPumpNum,
                                  FirstHVACIteration,
                                  CompressorOp,
                                  OpMode,
                                  QZnReq,
                                  ZoneNum,
                                  SpeedNum,
                                  SpeedRatio,
                                  PartLoadFrac,
                                  OnOffAirFlowRatio,
                                  SupHeaterLoad);
                if (SpeedNum == thisMSHeatPump.NumOfSpeedCooling && SpeedRatio == 1.0) {
                    // compressor on (reset inlet air mass flow rate to starting value)
                    state.dataLoopNodes->Node(InletNode).MassFlowRate = SaveMassFlowRate;
                    CompressorOp = DataHVACGlobals::CompressorOperation::On;
                    ControlMSHPOutput(state,
                                      MSHeatPumpNum,
                                      FirstHVACIteration,
                                      CompressorOp,
                                      OpMode,
                                      QZnReq,
                                      ZoneNum,
                                      SpeedNum,
                                      SpeedRatio,
                                      PartLoadFrac,
                                      OnOffAirFlowRatio,
                                      SupHeaterLoad);
                }
            } else {
                // compressor on
                ControlMSHPOutput(state,
                                  MSHeatPumpNum,
                                  FirstHVACIteration,
                                  CompressorOp,
                                  OpMode,
                                  QZnReq,
                                  ZoneNum,
                                  SpeedNum,
                                  SpeedRatio,
                                  PartLoadFrac,
                                  OnOffAirFlowRatio,
                                  SupHeaterLoad);
            }
        }

        if (thisMSHeatPump.HeatCoilType != MultiSpeedHeatingCoil) {
            state.dataHVACMultiSpdHP->SaveCompressorPLR = PartLoadFrac;
        } else {
            if (SpeedNum > 1) {
                state.dataHVACMultiSpdHP->SaveCompressorPLR = 1.0;
            }

            if (PartLoadFrac == 1.0 && state.dataHVACMultiSpdHP->SaveCompressorPLR < 1.0 && (!thisMSHeatPump.Staged)) {
                PartLoadFrac = state.dataHVACMultiSpdHP->SaveCompressorPLR;
            }
        }

        CalcMSHeatPump(state,
                       MSHeatPumpNum,
                       FirstHVACIteration,
                       CompressorOp,
                       SpeedNum,
                       SpeedRatio,
                       PartLoadFrac,
                       QSensUnitOut,
                       QZnReq,
                       OnOffAirFlowRatio,
                       SupHeaterLoad);

        // calculate delivered capacity
        AirMassFlow = state.dataLoopNodes->Node(InletNode).MassFlowRate;

        QTotUnitOut = AirMassFlow * (state.dataLoopNodes->Node(OutletNode).Enthalpy -
                                     state.dataLoopNodes->Node(thisMSHeatPump.NodeNumOfControlledZone).Enthalpy);

        // report variables
        thisMSHeatPump.CompPartLoadRatio = state.dataHVACMultiSpdHP->SaveCompressorPLR;
        if (thisMSHeatPump.OpMode == CycFanCycCoil) {
            if (SupHeaterLoad > 0.0) {
                thisMSHeatPump.FanPartLoadRatio = 1.0;
            } else {
                if (SpeedNum < 2) {
                    thisMSHeatPump.FanPartLoadRatio = PartLoadFrac;
                } else {
                    thisMSHeatPump.FanPartLoadRatio = 1.0;
                }
            }
        } else {
            if (UnitOn) {
                thisMSHeatPump.FanPartLoadRatio = 1.0;
            } else {
                if (SpeedNum < 2) {
                    thisMSHeatPump.FanPartLoadRatio = PartLoadFrac;
                } else {
                    thisMSHeatPump.FanPartLoadRatio = 1.0;
                }
            }
        }

        if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::HeatingMode) {
            thisMSHeatPump.TotHeatEnergyRate = std::abs(max(0.0, QTotUnitOut));
            thisMSHeatPump.SensHeatEnergyRate = std::abs(max(0.0, QSensUnitOut));
            thisMSHeatPump.LatHeatEnergyRate = std::abs(max(0.0, (QTotUnitOut - QSensUnitOut)));
            thisMSHeatPump.TotCoolEnergyRate = 0.0;
            thisMSHeatPump.SensCoolEnergyRate = 0.0;
            thisMSHeatPump.LatCoolEnergyRate = 0.0;
        }
        if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::CoolingMode) {
            thisMSHeatPump.TotCoolEnergyRate = std::abs(min(0.0, QTotUnitOut));
            thisMSHeatPump.SensCoolEnergyRate = std::abs(min(0.0, QSensUnitOut));
            thisMSHeatPump.LatCoolEnergyRate = std::abs(min(0.0, (QTotUnitOut - QSensUnitOut)));
            thisMSHeatPump.TotHeatEnergyRate = 0.0;
            thisMSHeatPump.SensHeatEnergyRate = 0.0;
            thisMSHeatPump.LatHeatEnergyRate = 0.0;
        }

        thisMSHeatPump.AuxElecPower = thisMSHeatPump.AuxOnCyclePower * state.dataHVACMultiSpdHP->SaveCompressorPLR +
                                                 thisMSHeatPump.AuxOffCyclePower * (1.0 - state.dataHVACMultiSpdHP->SaveCompressorPLR);
        Real64 locFanElecPower = 0.0;
        locFanElecPower = Fans::GetFanPower(state, thisMSHeatPump.FanNum);
        thisMSHeatPump.ElecPower = locFanElecPower + state.dataHVACGlobal->DXElecCoolingPower + state.dataHVACGlobal->DXElecHeatingPower +
                                              state.dataHVACGlobal->ElecHeatingCoilPower + state.dataHVACGlobal->SuppHeatingCoilPower +
                                              thisMSHeatPump.AuxElecPower;
    }

    //******************************************************************************

    void GetMSHeatPumpInput(EnergyPlusData &state)
    {
        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Lixing Gu, FSEC
        //       DATE WRITTEN:    July 2007
        //       MODIFIED         na
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        //  This routine will get the input required by the multispeed heat pump model

        // Using/Aliasing
        using BranchNodeConnections::TestCompSet;
        using DataHVACGlobals::FanType_SimpleConstVolume;
        using DataHVACGlobals::FanType_SimpleOnOff;
        using DataSizing::AutoSize;
        using Fans::GetFanIndex;
        using Fans::GetFanInletNode;
        using Fans::GetFanOutletNode;
        using Fans::GetFanType;
        using Fans::GetFanVolFlow;
        using FluidProperties::FindGlycol;

        using BranchNodeConnections::SetUpCompSets;
        using DXCoils::GetDXCoilIndex;
        using NodeInputManager::GetOnlySingleNode;
        using DXCoils::GetDXCoilNumberOfSpeeds;
        using HeatingCoils::GetCoilIndex;
        using HeatingCoils::GetCoilInletNode;
        using HeatingCoils::GetCoilOutletNode;
        using HeatingCoils::GetHeatingCoilIndex;
        using HeatingCoils::GetHeatingCoilNumberOfStages;
        using WaterCoils::GetCoilMaxWaterFlowRate;
        using WaterCoils::GetCoilWaterInletNode;
        using SteamCoils::GetCoilAirOutletNode;
        using SteamCoils::GetCoilSteamInletNode;
        using SteamCoils::GetSteamCoilIndex;
        using DXCoils::SetMSHPDXCoilHeatRecoveryFlag;
        using FluidProperties::GetSatDensityRefrig;

        // Locals
        // PARAMETERS
        static constexpr std::string_view RoutineName("GetMSHeatPumpInput: "); // include trailing blank space
        static constexpr std::string_view RoutineNameNoColon("GetMSHeatPumpInput");

        // LOCAL VARIABLES
        int NumAlphas;                 // Number of elements in the alpha array
        int NumNumbers;                // Number of Numbers for each GetObjectItem call
        int IOStatus;                  // Used in GetObjectItem
        bool ErrorsFound(false);       // True when input errors found
        bool IsNotOK;                  // Flag to verify name
        bool AirNodeFound;             // True when an air node is found
        bool AirLoopFound;             // True when an air loop is found
        int FanType;                   // Fan type
        bool Found;                    // Flag to find autosize
        int HeatingCoilInletNode;      // Heating coil inlet node number
        int HeatingCoilOutletNode;     // Heating coil outlet node number
        int CoolingCoilInletNode;      // Cooling coil inlet node number
        int CoolingCoilOutletNode;     // Cooling coil outlet node number
        int SuppHeatCoilInletNode;     // Supplemental heating coil inlet node number
        int SuppHeatCoilOutletNode;    // Supplemental heating coil outlet node number
        bool LocalError;               // Local error flag
        Array1D_string Alphas;         // Alpha input items for object
        Array1D_string cAlphaFields;   // Alpha field names
        Array1D_string cNumericFields; // Numeric field names
        Array1D<Real64> Numbers;       // Numeric input items for object
        Array1D_bool lAlphaBlanks;     // Logical array, alpha field input BLANK = .TRUE.
        Array1D_bool lNumericBlanks;   // Logical array, numeric field input BLANK = .TRUE.
        int MaxNums(0);                // Maximum number of numeric input fields
        int MaxAlphas(0);              // Maximum number of alpha input fields
        int TotalArgs(0);              // Total number of alpha and numeric arguments (max) for a
        //  certain object in the input file
        bool errFlag;
        int SteamIndex;      // steam coil steam inlet density
        Real64 SteamDensity; // density of steam at 100C

        if (state.dataHVACMultiSpdHP->MSHeatPump.allocated()) return;

        state.dataHVACMultiSpdHP->CurrentModuleObject =
            "AirLoopHVAC:UnitaryHeatPump:AirToAir:MultiSpeed"; // Object type for getting and error messages

        state.dataInputProcessing->inputProcessor->getObjectDefMaxArgs(
            state, state.dataHVACMultiSpdHP->CurrentModuleObject, TotalArgs, NumAlphas, NumNumbers);
        MaxNums = max(MaxNums, NumNumbers);
        MaxAlphas = max(MaxAlphas, NumAlphas);

        Alphas.allocate(MaxAlphas);
        cAlphaFields.allocate(MaxAlphas);
        Numbers.dimension(MaxNums, 0.0);
        cNumericFields.allocate(MaxNums);
        lAlphaBlanks.dimension(MaxAlphas, true);
        lNumericBlanks.dimension(MaxNums, true);

        state.dataHVACMultiSpdHP->NumMSHeatPumps =
            state.dataInputProcessing->inputProcessor->getNumObjectsFound(state, state.dataHVACMultiSpdHP->CurrentModuleObject);

        if (state.dataHVACMultiSpdHP->NumMSHeatPumps <= 0) {
            ShowSevereError(state, format("No {} objects specified in input file.", state.dataHVACMultiSpdHP->CurrentModuleObject));
            ErrorsFound = true;
        }

        // ALLOCATE ARRAYS
        state.dataHVACMultiSpdHP->MSHeatPump.allocate(state.dataHVACMultiSpdHP->NumMSHeatPumps);
        state.dataHVACMultiSpdHP->MSHeatPumpReport.allocate(state.dataHVACMultiSpdHP->NumMSHeatPumps);
        state.dataHVACMultiSpdHP->CheckEquipName.dimension(state.dataHVACMultiSpdHP->NumMSHeatPumps, true);

        // Load arrays with reformulated electric EIR chiller data
        for (int MSHPNum = 1; MSHPNum <= state.dataHVACMultiSpdHP->NumMSHeatPumps; ++MSHPNum) {

            auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHPNum);

            HeatingCoilInletNode = 0;
            HeatingCoilOutletNode = 0;
            CoolingCoilInletNode = 0;
            CoolingCoilOutletNode = 0;
            SuppHeatCoilInletNode = 0;
            SuppHeatCoilOutletNode = 0;

            state.dataInputProcessing->inputProcessor->getObjectItem(state,
                                                                     state.dataHVACMultiSpdHP->CurrentModuleObject,
                                                                     MSHPNum,
                                                                     Alphas,
                                                                     NumAlphas,
                                                                     Numbers,
                                                                     NumNumbers,
                                                                     IOStatus,
                                                                     lNumericBlanks,
                                                                     lAlphaBlanks,
                                                                     cAlphaFields,
                                                                     cNumericFields);

            thisMSHeatPump.Name = Alphas(1);
            if (lAlphaBlanks(2)) {
                thisMSHeatPump.AvaiSchedPtr = DataGlobalConstants::ScheduleAlwaysOn;
            } else {
                thisMSHeatPump.AvaiSchedPtr = GetScheduleIndex(state, Alphas(2));
                if (thisMSHeatPump.AvaiSchedPtr == 0) {
                    ShowSevereError(state,
                                    format("{}, \"{}\" {} not found: {}",
                                           state.dataHVACMultiSpdHP->CurrentModuleObject,
                                           thisMSHeatPump.Name,
                                           cAlphaFields(2),
                                           Alphas(2)));
                    ErrorsFound = true;
                }
            }

            thisMSHeatPump.AirInletNodeName = Alphas(3);
            thisMSHeatPump.AirOutletNodeName = Alphas(4);
            thisMSHeatPump.AirInletNodeNum = GetOnlySingleNode(state,
                                                                    Alphas(3),
                                                                    ErrorsFound,
                                                                    DataLoopNode::ConnectionObjectType::AirLoopHVACUnitaryHeatPumpAirToAirMultiSpeed,
                                                                    Alphas(1),
                                                                    DataLoopNode::NodeFluidType::Air,
                                                                    DataLoopNode::ConnectionType::Inlet,
                                                                    NodeInputManager::CompFluidStream::Primary,
                                                                    ObjectIsParent);

            thisMSHeatPump.AirOutletNodeNum = GetOnlySingleNode(state,
                                                                     Alphas(4),
                                                                     ErrorsFound,
                                                                     DataLoopNode::ConnectionObjectType::AirLoopHVACUnitaryHeatPumpAirToAirMultiSpeed,
                                                                     Alphas(1),
                                                                     DataLoopNode::NodeFluidType::Air,
                                                                     DataLoopNode::ConnectionType::Outlet,
                                                                     NodeInputManager::CompFluidStream::Primary,
                                                                     ObjectIsParent);

            TestCompSet(state, state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1), Alphas(3), Alphas(4), "Air Nodes");

            // Get the Controlling Zone or Location of the engine driven heat pump Thermostat
            thisMSHeatPump.ControlZoneNum = UtilityRoutines::FindItemInList(Alphas(5), state.dataHeatBal->Zone);
            thisMSHeatPump.ControlZoneName = Alphas(5);
            if (thisMSHeatPump.ControlZoneNum == 0) {
                ShowSevereError(state,
                                format("{}, \"{}\" {} not found: {}",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cAlphaFields(5),
                                       thisMSHeatPump.ControlZoneName));
                ErrorsFound = true;
            }

            // Get the node number for the zone with the thermostat
            if (thisMSHeatPump.ControlZoneNum > 0) {
                AirNodeFound = false;
                AirLoopFound = false;
                int ControlledZoneNum = thisMSHeatPump.ControlZoneNum;
                // Find the controlled zone number for the specified thermostat location
                thisMSHeatPump.NodeNumOfControlledZone = state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).ZoneNode;
                // Determine if furnace is on air loop served by the thermostat location specified
                for (int zoneInNode = 1; zoneInNode <= state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).NumInletNodes; ++zoneInNode) {
                    int AirLoopNumber = state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).InletNodeAirLoopNum(zoneInNode);
                    if (AirLoopNumber > 0) {
                        for (int BranchNum = 1; BranchNum <= state.dataAirSystemsData->PrimaryAirSystems(AirLoopNumber).NumBranches; ++BranchNum) {
                            for (int CompNum = 1;
                                 CompNum <= state.dataAirSystemsData->PrimaryAirSystems(AirLoopNumber).Branch(BranchNum).TotalComponents;
                                 ++CompNum) {
                                if (!UtilityRoutines::SameString(
                                        state.dataAirSystemsData->PrimaryAirSystems(AirLoopNumber).Branch(BranchNum).Comp(CompNum).Name,
                                        thisMSHeatPump.Name) ||
                                    !UtilityRoutines::SameString(
                                        state.dataAirSystemsData->PrimaryAirSystems(AirLoopNumber).Branch(BranchNum).Comp(CompNum).TypeOf,
                                        state.dataHVACMultiSpdHP->CurrentModuleObject))
                                    continue;
                                AirLoopFound = true;
                                thisMSHeatPump.AirLoopNumber = AirLoopNumber;
                                break;
                            }
                            thisMSHeatPump.ZoneInletNode = state.dataZoneEquip->ZoneEquipConfig(ControlledZoneNum).InletNode(zoneInNode);
                            if (AirLoopFound) break;
                        }
                        for (int TstatZoneNum = 1; TstatZoneNum <= state.dataZoneCtrls->NumTempControlledZones; ++TstatZoneNum) {
                            if (state.dataZoneCtrls->TempControlledZone(TstatZoneNum).ActualZoneNum != thisMSHeatPump.ControlZoneNum) continue;
                            AirNodeFound = true;
                        }
                        for (int TstatZoneNum = 1; TstatZoneNum <= state.dataZoneCtrls->NumComfortControlledZones; ++TstatZoneNum) {
                            if (state.dataZoneCtrls->ComfortControlledZone(TstatZoneNum).ActualZoneNum != thisMSHeatPump.ControlZoneNum)
                                continue;
                            AirNodeFound = true;
                        }
                        for (int TstatZoneNum = 1; TstatZoneNum <= state.dataZoneTempPredictorCorrector->NumStageCtrZone; ++TstatZoneNum) {
                            if (state.dataZoneCtrls->StageControlledZone(TstatZoneNum).ActualZoneNum != thisMSHeatPump.ControlZoneNum) continue;
                            AirNodeFound = true;
                        }
                    }
                    if (AirLoopFound) break;
                }
                if (!AirNodeFound) {
                    ShowSevereError(state,
                                    format("Did not find Air Node ({}), {} = \"\"{}",
                                           cAlphaFields(5),
                                           state.dataHVACMultiSpdHP->CurrentModuleObject,
                                           thisMSHeatPump.Name));
                    ShowContinueError(state, format("Specified {} = {}", cAlphaFields(5), Alphas(5)));
                    ErrorsFound = true;
                }
                if (!AirLoopFound) {
                    ShowSevereError(state,
                                    format("Did not find correct AirLoopHVAC for {} = {}",
                                           state.dataHVACMultiSpdHP->CurrentModuleObject,
                                           thisMSHeatPump.Name));
                    ShowContinueError(state, format("The {} = {} is not served by this Primary Air Loop equipment.", cAlphaFields(5), Alphas(5)));
                    ErrorsFound = true;
                }
            }

            // Get supply fan data
            thisMSHeatPump.FanName = Alphas(7);
            if (UtilityRoutines::SameString(Alphas(6), "Fan:OnOff") || UtilityRoutines::SameString(Alphas(6), "Fan:ConstantVolume")) {
                if (UtilityRoutines::SameString(Alphas(6), "Fan:OnOff")) {
                    thisMSHeatPump.FanType = FanType_SimpleOnOff;
                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Fan:OnOff",
                                  thisMSHeatPump.FanName,
                                  "UNDEFINED",
                                  "UNDEFINED");
                    thisMSHeatPump.FanInletNode = GetFanInletNode(state, "Fan:OnOff", thisMSHeatPump.FanName, ErrorsFound);
                    thisMSHeatPump.FanOutletNode = GetFanOutletNode(state, "Fan:OnOff", thisMSHeatPump.FanName, ErrorsFound);
                } else {
                    thisMSHeatPump.FanType = FanType_SimpleConstVolume;
                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Fan:ConstantVolume",
                                  thisMSHeatPump.FanName,
                                  "UNDEFINED",
                                  "UNDEFINED");
                    thisMSHeatPump.FanInletNode = GetFanInletNode(state, "Fan:ConstantVolume", thisMSHeatPump.FanName, ErrorsFound);
                    thisMSHeatPump.FanOutletNode = GetFanOutletNode(state, "Fan:ConstantVolume", thisMSHeatPump.FanName, ErrorsFound);
                }
                GetFanIndex(state, Alphas(7), thisMSHeatPump.FanNum, ErrorsFound, state.dataHVACMultiSpdHP->CurrentModuleObject);
                GetFanType(state, Alphas(7), FanType, ErrorsFound);
                if (FanType != thisMSHeatPump.FanType) {
                    ShowSevereError(state,
                                    format("{}, \"{}\", {} and {} do not match in Fan objects.",
                                           state.dataHVACMultiSpdHP->CurrentModuleObject,
                                           thisMSHeatPump.Name,
                                           cAlphaFields(6),
                                           cAlphaFields(7)));
                    ShowContinueError(state, format("The entered {} = {} and {} = {}", cAlphaFields(7), Alphas(7), cAlphaFields(6), Alphas(6)));
                    ErrorsFound = true;
                }
            } else {
                ShowSevereError(state,
                                format("{}, \"{}\", {} is not allowed = {}",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cAlphaFields(6),
                                       Alphas(6)));
                ShowContinueError(state, "Valid choices are Fan:OnOff or Fan:ConstantVolume");
                ErrorsFound = true;
            }

            // Get supply fan placement data
            if (UtilityRoutines::SameString(Alphas(8), "BlowThrough") || UtilityRoutines::SameString(Alphas(8), "DrawThrough")) {
                if (UtilityRoutines::SameString(Alphas(8), "BlowThrough")) {
                    thisMSHeatPump.FanPlaceType = BlowThru;
                } else {
                    thisMSHeatPump.FanPlaceType = DrawThru;
                }
            } else {
                ShowSevereError(state,
                                format("{}, \"{}\", {} is not allowed = {}",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cAlphaFields(8),
                                       Alphas(8)));
                ShowContinueError(state, "Valid choices are BlowThrough or DrawThrough");
                ErrorsFound = true;
            }

            thisMSHeatPump.FanSchedule = Alphas(9);
            thisMSHeatPump.FanSchedPtr = GetScheduleIndex(state, Alphas(9));
            if (thisMSHeatPump.FanSchedPtr == 0) {
                ShowSevereError(state,
                                format("{}, \"{}\" {} not found: {}",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cAlphaFields(9),
                                       Alphas(9)));
                ErrorsFound = true;
            }

            if (thisMSHeatPump.FanSchedPtr > 0 && thisMSHeatPump.FanType == FanType_SimpleConstVolume) {
                if (!CheckScheduleValueMinMax(state, thisMSHeatPump.FanSchedPtr, ">", 0.0, "<=", 1.0)) {
                    ShowSevereError(state, format("{} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state,
                                      format("{} must be continuous (fan operating mode schedule values > 0) for {} = Fan:ConstantVolume.",
                                             cAlphaFields(9),
                                             cAlphaFields(6)));
                    ShowContinueError(state, format("Error found in {} = {}", cAlphaFields(9), Alphas(9)));
                    ShowContinueError(state, "schedule values must be (>0., <=1.)");
                    ErrorsFound = true;
                }
            }

            if (UtilityRoutines::SameString(Alphas(10), "Coil:Heating:DX:MultiSpeed")) {
                thisMSHeatPump.HeatCoilType = MultiSpeedHeatingCoil;
                thisMSHeatPump.HeatCoilNum =
                    state.dataInputProcessing->inputProcessor->getObjectItemNum(state, "Coil:Heating:DX:MultiSpeed", Alphas(11));
                thisMSHeatPump.DXHeatCoilName = Alphas(11);
                if (thisMSHeatPump.HeatCoilNum <= 0) {
                    ShowSevereError(state, format("Configuration error in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ShowContinueError(state, format("{} \"{}\" not found.", cAlphaFields(11), Alphas(11)));
                    ShowContinueError(state, format("{} must be Coil:Heating:DX:MultiSpeed ", cAlphaFields(10)));
                    ShowFatalError(state,
                                   format("{}Errors found in getting {} input. Preceding condition(s) causes termination.",
                                          RoutineName,
                                          state.dataHVACMultiSpdHP->CurrentModuleObject));
                    ErrorsFound = true;
                }
                LocalError = false;
                GetDXCoilIndex(
                    state, thisMSHeatPump.DXHeatCoilName, thisMSHeatPump.DXHeatCoilIndex, LocalError, "Coil:Heating:DX:MultiSpeed");
                if (LocalError) {
                    ShowSevereError(state, format("The index of {} is not found \"{}\"", cAlphaFields(11), Alphas(11)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                HeatingCoilInletNode = DXCoils::GetCoilInletNode(state, Alphas(10), Alphas(11), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The inlet node number of {} is not found \"{}\"", cAlphaFields(11), Alphas(11)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                HeatingCoilOutletNode = DXCoils::GetCoilOutletNode(state, Alphas(10), Alphas(11), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The outlet node number of {} is not found \"{}\"", cAlphaFields(11), Alphas(11)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                thisMSHeatPump.MinOATCompressorHeating = DXCoils::GetMinOATCompressor(state, thisMSHeatPump.DXHeatCoilIndex, LocalError);
                if (LocalError) {
                    ShowContinueError(state,
                                      format("...for heating coil. Occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    LocalError = false;
                }
                SetUpCompSets(state,
                              state.dataHVACMultiSpdHP->CurrentModuleObject,
                              thisMSHeatPump.Name,
                              "Coil:Heating:DX:MultiSpeed",
                              thisMSHeatPump.DXHeatCoilName,
                              "UNDEFINED",
                              "UNDEFINED");
            } else if (UtilityRoutines::SameString(Alphas(10), "Coil:Heating:Electric:MultiStage") ||
                       UtilityRoutines::SameString(Alphas(10), "Coil:Heating:Gas:MultiStage")) {

                if (UtilityRoutines::SameString(Alphas(10), "Coil:Heating:Electric:MultiStage")) {
                    thisMSHeatPump.HeatCoilType = Coil_HeatingElectric_MultiStage;
                    thisMSHeatPump.HeatCoilNum =
                        state.dataInputProcessing->inputProcessor->getObjectItemNum(state, "Coil:Heating:Electric:MultiStage", Alphas(11));
                    if (thisMSHeatPump.HeatCoilNum <= 0) {
                        ShowSevereError(state, format("Configuration error in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                        ShowContinueError(state, format("{} \"{}\" not found.", cAlphaFields(11), Alphas(11)));
                        ShowContinueError(state, format("{} must be Coil:Heating:Electric:MultiStage ", cAlphaFields(10)));
                        ShowFatalError(state,
                                       format("{}Errors found in getting {} input. Preceding condition(s) causes termination.",
                                              RoutineName,
                                              state.dataHVACMultiSpdHP->CurrentModuleObject));
                        ErrorsFound = true;
                    }
                } else {
                    thisMSHeatPump.HeatCoilType = Coil_HeatingGas_MultiStage;
                    thisMSHeatPump.HeatCoilNum =
                        state.dataInputProcessing->inputProcessor->getObjectItemNum(state, "Coil:Heating:Gas:MultiStage", Alphas(11));
                    if (thisMSHeatPump.HeatCoilNum <= 0) {
                        ShowSevereError(state, format("Configuration error in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                        ShowContinueError(state, format("{} \"{}\" not found.", cAlphaFields(11), Alphas(11)));
                        ShowContinueError(state, format("{} must be Coil:Heating:Gas:MultiStage ", cAlphaFields(10)));
                        ShowFatalError(state,
                                       format("{}Errors found in getting {} input. Preceding condition(s) causes termination.",
                                              RoutineName,
                                              state.dataHVACMultiSpdHP->CurrentModuleObject));
                        ErrorsFound = true;
                    }
                }
                thisMSHeatPump.HeatCoilName = Alphas(11);
                LocalError = false;
                if (UtilityRoutines::SameString(Alphas(10), "Coil:Heating:Electric:MultiStage")) {
                    GetCoilIndex(state, thisMSHeatPump.HeatCoilName, thisMSHeatPump.HeatCoilIndex, LocalError);
                } else {
                    GetCoilIndex(state, thisMSHeatPump.HeatCoilName, thisMSHeatPump.HeatCoilIndex, LocalError);
                }
                if (LocalError) {
                    ShowSevereError(state, format("The index of {} is not found \"{}\"", cAlphaFields(11), Alphas(11)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                HeatingCoilInletNode = GetCoilInletNode(state, Alphas(10), Alphas(11), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The inlet node number of {} is not found \"{}\"", cAlphaFields(11), Alphas(11)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                HeatingCoilOutletNode = GetCoilOutletNode(state, Alphas(10), Alphas(11), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The outlet node number of {} is not found \"{}\"", cAlphaFields(11), Alphas(11)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                if (UtilityRoutines::SameString(Alphas(10), "Coil:Heating:Electric:MultiStage")) {
                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Coil:Heating:Electric:MultiStage",
                                  thisMSHeatPump.HeatCoilName,
                                  "UNDEFINED",
                                  "UNDEFINED");
                } else {
                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Coil:Heating:Gas:MultiStage",
                                  thisMSHeatPump.HeatCoilName,
                                  "UNDEFINED",
                                  "UNDEFINED");
                }
            } else if (UtilityRoutines::SameString(Alphas(10), "Coil:Heating:Water")) {
                thisMSHeatPump.HeatCoilType = Coil_HeatingWater;
                ValidateComponent(state, Alphas(10), Alphas(11), IsNotOK, state.dataHVACMultiSpdHP->CurrentModuleObject);
                if (IsNotOK) {
                    ShowContinueError(state, format("...occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                } else { // mine data from heating coil object

                    thisMSHeatPump.HeatCoilName = Alphas(11);
                    // Get the Heating Coil water Inlet or control Node number
                    errFlag = false;
                    thisMSHeatPump.CoilControlNode =
                        GetCoilWaterInletNode(state, "Coil:Heating:Water", thisMSHeatPump.HeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the ReHeat Coil hot water max volume flow rate
                    errFlag = false;
                    thisMSHeatPump.MaxCoilFluidFlow =
                        GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.HeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the lemental Heating Coil Inlet Node
                    errFlag = false;
                    HeatingCoilInletNode = WaterCoils::GetCoilInletNode(state, "Coil:Heating:Water", thisMSHeatPump.HeatCoilName, errFlag);
                    thisMSHeatPump.CoilAirInletNode = HeatingCoilInletNode;
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the lemental Heating Coil Outlet Node
                    errFlag = false;
                    HeatingCoilOutletNode = WaterCoils::GetCoilOutletNode(state, "Coil:Heating:Water", thisMSHeatPump.HeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }
                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Coil:Heating:Water",
                                  thisMSHeatPump.HeatCoilName,
                                  state.dataLoopNodes->NodeID(HeatingCoilInletNode),
                                  state.dataLoopNodes->NodeID(HeatingCoilOutletNode));
                }
            } else if (UtilityRoutines::SameString(Alphas(10), "Coil:Heating:Steam")) {
                thisMSHeatPump.HeatCoilType = Coil_HeatingSteam;
                ValidateComponent(state, Alphas(10), Alphas(11), IsNotOK, state.dataHVACMultiSpdHP->CurrentModuleObject);
                if (IsNotOK) {
                    ShowContinueError(state, format("...occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ErrorsFound = true;
                } else { // mine data from heating coil object

                    thisMSHeatPump.HeatCoilName = Alphas(11);
                    errFlag = false;
                    thisMSHeatPump.HeatCoilNum = GetSteamCoilIndex(state, Alphas(10), thisMSHeatPump.HeatCoilName, errFlag);
                    if (thisMSHeatPump.HeatCoilNum == 0) {
                        ShowSevereError(state,
                                        format("{} illegal {} = {}",
                                               state.dataHVACMultiSpdHP->CurrentModuleObject,
                                               cAlphaFields(10),
                                               thisMSHeatPump.HeatCoilName));
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the lemental Heating Coil steam inlet node number
                    errFlag = false;
                    thisMSHeatPump.CoilControlNode =
                        GetCoilAirOutletNode(state, "Coil:Heating:Steam", thisMSHeatPump.HeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the lemental Heating Coil steam max volume flow rate
                    thisMSHeatPump.MaxCoilFluidFlow = SteamCoils::GetCoilMaxSteamFlowRate(state, thisMSHeatPump.HeatCoilNum, errFlag);
                    if (thisMSHeatPump.MaxCoilFluidFlow > 0.0) {
                        SteamIndex = 0; // Function GetSatDensityRefrig will look up steam index if 0 is passed
                        SteamDensity =
                            GetSatDensityRefrig(state, fluidNameSteam, state.dataHVACMultiSpdHP->TempSteamIn, 1.0, SteamIndex, RoutineNameNoColon);
                        thisMSHeatPump.MaxCoilFluidFlow *= SteamDensity;
                    }

                    // Get the lemental Heating Coil Inlet Node
                    errFlag = false;
                    HeatingCoilInletNode =
                        SteamCoils::GetCoilAirInletNode(state, thisMSHeatPump.HeatCoilNum, thisMSHeatPump.HeatCoilName, errFlag);
                    thisMSHeatPump.CoilAirInletNode = HeatingCoilInletNode;
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the lemental Heating Coil Outlet Node
                    errFlag = false;
                    HeatingCoilOutletNode = GetCoilAirOutletNode(state, thisMSHeatPump.HeatCoilNum, thisMSHeatPump.HeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Coil:Heating:Steam",
                                  thisMSHeatPump.HeatCoilName,
                                  state.dataLoopNodes->NodeID(HeatingCoilInletNode),
                                  state.dataLoopNodes->NodeID(HeatingCoilOutletNode));
                }
            } else {
                ShowSevereError(state,
                                format("The allowed {} are Coil:Heating:DX:MultiSpeed, Coil:Heating:Electric:MultiStage, and "
                                       "Coil:Heating:Gas:MultiStage  in {} \"{}\"",
                                       cAlphaFields(10),
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       Alphas(1)));
                ShowContinueError(state, format("The entered {} = \"{}\".", cAlphaFields(10), Alphas(10)));
                ErrorsFound = true;
            }

            // thisMSHeatPump.MinOATCompressor = Numbers(1); // deprecated, now uses coil MinOAT inputs

            if (UtilityRoutines::SameString(Alphas(12), "Coil:Cooling:DX:MultiSpeed")) {
                thisMSHeatPump.CoolCoilType = MultiSpeedCoolingCoil;
                thisMSHeatPump.DXCoolCoilName = Alphas(13);
                if (state.dataInputProcessing->inputProcessor->getObjectItemNum(state, "Coil:Cooling:DX:MultiSpeed", Alphas(13)) <= 0) {
                    ShowSevereError(state, format("Configuration error in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ShowContinueError(state, format("{} \"{}\" not found.", cAlphaFields(13), Alphas(13)));
                    ShowContinueError(state, format("{} must be Coil:Cooling:DX:MultiSpeed ", cAlphaFields(12)));
                    ShowFatalError(state,
                                   format("{}Errors found in getting {} input. Preceding condition(s) causes termination.",
                                          RoutineName,
                                          state.dataHVACMultiSpdHP->CurrentModuleObject));
                    ErrorsFound = true;
                }
                LocalError = false;
                GetDXCoilIndex(
                    state, thisMSHeatPump.DXCoolCoilName, thisMSHeatPump.DXCoolCoilIndex, LocalError, "Coil:Cooling:DX:MultiSpeed");
                if (LocalError) {
                    ShowSevereError(state, format("The index of {} is not found \"{}\"", cAlphaFields(13), Alphas(13)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                CoolingCoilInletNode = DXCoils::GetCoilInletNode(state, Alphas(12), Alphas(13), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The inlet node number of {} is not found \"{}\"", cAlphaFields(13), Alphas(13)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                CoolingCoilOutletNode = DXCoils::GetCoilOutletNode(state, Alphas(12), Alphas(13), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The outlet node number of {} is not found \"{}\"", cAlphaFields(13), Alphas(13)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                thisMSHeatPump.MinOATCompressorCooling = DXCoils::GetMinOATCompressor(state, thisMSHeatPump.DXCoolCoilIndex, LocalError);
                if (LocalError) {
                    ShowContinueError(state,
                                      format("...for cooling coil. Occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    LocalError = false;
                }
            } else {
                ShowSevereError(state,
                                format("The allowed {} is Coil:Cooling:DX:MultiSpeed in {} \"{}\"",
                                       cAlphaFields(12),
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       Alphas(1)));
                ShowContinueError(state, format("The entered {} = \"{}\".", cAlphaFields(12), Alphas(12)));
                ErrorsFound = true;
            }
            SetUpCompSets(state,
                          state.dataHVACMultiSpdHP->CurrentModuleObject,
                          thisMSHeatPump.Name,
                          "Coil:Cooling:DX:MultiSpeed",
                          thisMSHeatPump.DXCoolCoilName,
                          "UNDEFINED",
                          "UNDEFINED");

            // Get supplemental heating coil data
            thisMSHeatPump.SuppHeatCoilName = Alphas(15);
            if (UtilityRoutines::SameString(Alphas(14), "Coil:Heating:Fuel")) {
                thisMSHeatPump.SuppHeatCoilType = SuppHeatingCoilGas;
                errFlag = false;
                thisMSHeatPump.SuppHeatCoilNum = GetHeatingCoilIndex(state, "Coil:Heating:Fuel", Alphas(15), errFlag);
                if (thisMSHeatPump.SuppHeatCoilNum <= 0 || errFlag) {
                    ShowContinueError(state, format("Configuration error in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ShowContinueError(state, format("{} of type Coil:Heating:Fuel \"{}\" not found.", cAlphaFields(15), Alphas(15)));
                    ErrorsFound = true;
                }

                // Get the Supplemental Heating Coil Node Numbers
                LocalError = false;
                SuppHeatCoilInletNode = HeatingCoils::GetCoilInletNode(state, Alphas(14), Alphas(15), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The inlet node number of {} is not found \"{}\"", cAlphaFields(15), Alphas(15)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                SuppHeatCoilOutletNode = HeatingCoils::GetCoilOutletNode(state, Alphas(14), Alphas(15), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The outlet node number of {} is not found \"{}\"", cAlphaFields(15), Alphas(15)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }

                // Get supplemental heating coil capacity to see if it is autosize
                thisMSHeatPump.DesignSuppHeatingCapacity = HeatingCoils::GetCoilCapacity(state, Alphas(14), Alphas(15), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The capacity {} is not found \"{}\"", cAlphaFields(15), Alphas(15)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                SetUpCompSets(state,
                              state.dataHVACMultiSpdHP->CurrentModuleObject,
                              thisMSHeatPump.Name,
                              "Coil:Heating:Fuel",
                              thisMSHeatPump.SuppHeatCoilName,
                              "UNDEFINED",
                              "UNDEFINED");
            }
            if (UtilityRoutines::SameString(Alphas(14), "Coil:Heating:Electric")) {
                thisMSHeatPump.SuppHeatCoilType = SuppHeatingCoilElec;
                errFlag = false;
                thisMSHeatPump.SuppHeatCoilNum = GetHeatingCoilIndex(state, "Coil:Heating:Electric", Alphas(15), errFlag);
                if (thisMSHeatPump.SuppHeatCoilNum <= 0 || errFlag) {
                    ShowContinueError(state, format("Configuration error in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ShowContinueError(state, format("{} of type Coil:Heating:Electric \"{}\" not found.", cAlphaFields(15), Alphas(15)));
                    ErrorsFound = true;
                }

                // Get the Supplemental Heating Coil Node Numbers
                LocalError = false;
                SuppHeatCoilInletNode = HeatingCoils::GetCoilInletNode(state, Alphas(14), Alphas(15), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The inlet node number of {} is not found \"{}\"", cAlphaFields(15), Alphas(15)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }
                SuppHeatCoilOutletNode = HeatingCoils::GetCoilOutletNode(state, Alphas(14), Alphas(15), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The outlet node number of {} is not found \"{}\"", cAlphaFields(15), Alphas(15)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }

                // Get supplemental heating coil capacity to see if it is autosize
                thisMSHeatPump.DesignSuppHeatingCapacity = HeatingCoils::GetCoilCapacity(state, Alphas(14), Alphas(15), LocalError);
                if (LocalError) {
                    ShowSevereError(state, format("The capacity {} is not found \"{}\"", cAlphaFields(15), Alphas(15)));
                    ShowContinueError(state, format("...occurs in {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                    LocalError = false;
                }

                SetUpCompSets(state,
                              state.dataHVACMultiSpdHP->CurrentModuleObject,
                              thisMSHeatPump.Name,
                              "Coil:Heating:Electric",
                              thisMSHeatPump.SuppHeatCoilName,
                              "UNDEFINED",
                              "UNDEFINED");
            }

            if (UtilityRoutines::SameString(Alphas(14), "Coil:Heating:Water")) {
                thisMSHeatPump.SuppHeatCoilType = Coil_HeatingWater;
                ValidateComponent(state, Alphas(14), thisMSHeatPump.SuppHeatCoilName, IsNotOK, state.dataHVACMultiSpdHP->CurrentModuleObject);
                if (IsNotOK) {
                    ShowContinueError(state, format("...occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1)));
                    ErrorsFound = true;
                } else { // mine data from heating coil object

                    // Get the Heating Coil water Inlet or control Node number
                    errFlag = false;
                    thisMSHeatPump.SuppCoilControlNode =
                        GetCoilWaterInletNode(state, "Coil:Heating:Water", thisMSHeatPump.SuppHeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the ReHeat Coil hot water max volume flow rate
                    errFlag = false;
                    thisMSHeatPump.MaxSuppCoilFluidFlow =
                        GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.SuppHeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the Supplemental Heating Coil Inlet Node
                    errFlag = false;
                    SuppHeatCoilInletNode = WaterCoils::GetCoilInletNode(state, "Coil:Heating:Water", thisMSHeatPump.SuppHeatCoilName, errFlag);
                    thisMSHeatPump.SuppCoilAirInletNode = SuppHeatCoilInletNode;
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the Supplemental Heating Coil Outlet Node
                    errFlag = false;
                    SuppHeatCoilOutletNode = WaterCoils::GetCoilOutletNode(state, "Coil:Heating:Water", thisMSHeatPump.SuppHeatCoilName, errFlag);
                    thisMSHeatPump.SuppCoilAirOutletNode = SuppHeatCoilOutletNode;
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }
                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Coil:Heating:Water",
                                  thisMSHeatPump.SuppHeatCoilName,
                                  state.dataLoopNodes->NodeID(SuppHeatCoilInletNode),
                                  state.dataLoopNodes->NodeID(SuppHeatCoilOutletNode));
                }
            }
            if (UtilityRoutines::SameString(Alphas(14), "Coil:Heating:Steam")) {
                thisMSHeatPump.SuppHeatCoilType = Coil_HeatingSteam;
                ValidateComponent(state, Alphas(14), thisMSHeatPump.SuppHeatCoilName, IsNotOK, state.dataHVACMultiSpdHP->CurrentModuleObject);
                if (IsNotOK) {
                    ShowContinueError(state, format("...occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ErrorsFound = true;
                } else { // mine data from heating coil object

                    errFlag = false;
                    thisMSHeatPump.SuppHeatCoilNum = GetSteamCoilIndex(state, Alphas(14), thisMSHeatPump.SuppHeatCoilName, errFlag);
                    if (thisMSHeatPump.SuppHeatCoilNum == 0) {
                        ShowSevereError(state,
                                        format("{} illegal {} = {}",
                                               state.dataHVACMultiSpdHP->CurrentModuleObject,
                                               cAlphaFields(14),
                                               thisMSHeatPump.SuppHeatCoilName));
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the Supplemental Heating Coil steam inlet node number
                    errFlag = false;
                    thisMSHeatPump.SuppCoilControlNode =
                        GetCoilAirOutletNode(state, "Coil:Heating:Steam", thisMSHeatPump.SuppHeatCoilName, errFlag);
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the Supplemental Heating Coil steam max volume flow rate
                    thisMSHeatPump.MaxSuppCoilFluidFlow = SteamCoils::GetCoilMaxSteamFlowRate(state, thisMSHeatPump.SuppHeatCoilNum, errFlag);
                    if (thisMSHeatPump.MaxSuppCoilFluidFlow > 0.0) {
                        SteamIndex = 0; // Function GetSatDensityRefrig will look up steam index if 0 is passed
                        SteamDensity =
                            GetSatDensityRefrig(state, fluidNameSteam, state.dataHVACMultiSpdHP->TempSteamIn, 1.0, SteamIndex, RoutineNameNoColon);
                        thisMSHeatPump.MaxSuppCoilFluidFlow *= SteamDensity;
                    }

                    // Get the Supplemental Heating Coil Inlet Node
                    errFlag = false;
                    SuppHeatCoilInletNode =
                        SteamCoils::GetCoilAirInletNode(state, thisMSHeatPump.SuppHeatCoilNum, thisMSHeatPump.SuppHeatCoilName, errFlag);
                    thisMSHeatPump.SuppCoilAirInletNode = SuppHeatCoilInletNode;
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    // Get the Supplemental Heating Coil Outlet Node
                    errFlag = false;
                    SuppHeatCoilOutletNode =
                        GetCoilAirOutletNode(state, thisMSHeatPump.SuppHeatCoilNum, thisMSHeatPump.SuppHeatCoilName, errFlag);
                    thisMSHeatPump.SuppCoilAirOutletNode = SuppHeatCoilOutletNode;
                    if (errFlag) {
                        ShowContinueError(state,
                                          format("Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                        ErrorsFound = true;
                    }

                    SetUpCompSets(state,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject,
                                  thisMSHeatPump.Name,
                                  "Coil:Heating:Steam",
                                  thisMSHeatPump.SuppHeatCoilName,
                                  state.dataLoopNodes->NodeID(SuppHeatCoilInletNode),
                                  state.dataLoopNodes->NodeID(SuppHeatCoilOutletNode));
                }
            }

            if (thisMSHeatPump.SuppHeatCoilType == 0) {
                ShowSevereError(state,
                                format("{}, \"{}\", {} is not allowed = {}",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cAlphaFields(14),
                                       Alphas(14)));
                ShowContinueError(state, "Valid choices are Coil:Heating:Fuel,Coil:Heating:Electric,Coil:Heating:Steam,or Coil:Heating:Water");
                ErrorsFound = true;
            }

            thisMSHeatPump.SuppMaxAirTemp = Numbers(2);
            thisMSHeatPump.SuppMaxOATemp = Numbers(3);
            if (thisMSHeatPump.SuppMaxOATemp > 21.0) {
                ShowSevereError(state,
                                format("{}, \"{}\", {} is greater than 21.0",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cNumericFields(3)));
                ShowContinueError(state, format("The input value is {:.2R}", Numbers(3)));
                ErrorsFound = true;
            }

            thisMSHeatPump.AuxOnCyclePower = Numbers(4);
            thisMSHeatPump.AuxOffCyclePower = Numbers(5);
            if (thisMSHeatPump.AuxOnCyclePower < 0.0) {
                ShowSevereError(state,
                                format("{}, \"{}\", A negative value for {} is not allowed ",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cNumericFields(4)));
                ErrorsFound = true;
            }
            if (thisMSHeatPump.AuxOffCyclePower < 0.0) {
                ShowSevereError(state,
                                format("{}, \"{}\", A negative value for {} is not allowed ",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cNumericFields(5)));
                ErrorsFound = true;
            }

            // Heat recovery
            thisMSHeatPump.DesignHeatRecFlowRate = Numbers(6);
            if (thisMSHeatPump.DesignHeatRecFlowRate > 0.0) {
                thisMSHeatPump.HeatRecActive = true;
                thisMSHeatPump.DesignHeatRecMassFlowRate =
                    RhoH2O(DataGlobalConstants::HWInitConvTemp) * thisMSHeatPump.DesignHeatRecFlowRate;
                thisMSHeatPump.HeatRecInletNodeNum =
                    GetOnlySingleNode(state,
                                      Alphas(16),
                                      ErrorsFound,
                                      DataLoopNode::ConnectionObjectType::AirLoopHVACUnitaryHeatPumpAirToAirMultiSpeed,
                                      Alphas(1),
                                      DataLoopNode::NodeFluidType::Water,
                                      DataLoopNode::ConnectionType::Inlet,
                                      NodeInputManager::CompFluidStream::Tertiary,
                                      ObjectIsNotParent);
                if (thisMSHeatPump.HeatRecInletNodeNum == 0) {
                    ShowSevereError(
                        state,
                        format("{}, \"{}\", Missing {}.", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name, cAlphaFields(16)));
                    ErrorsFound = true;
                }
                thisMSHeatPump.HeatRecOutletNodeNum =
                    GetOnlySingleNode(state,
                                      Alphas(17),
                                      ErrorsFound,
                                      DataLoopNode::ConnectionObjectType::AirLoopHVACUnitaryHeatPumpAirToAirMultiSpeed,
                                      Alphas(1),
                                      DataLoopNode::NodeFluidType::Water,
                                      DataLoopNode::ConnectionType::Outlet,
                                      NodeInputManager::CompFluidStream::Tertiary,
                                      ObjectIsNotParent);
                if (thisMSHeatPump.HeatRecOutletNodeNum == 0) {
                    ShowSevereError(
                        state,
                        format("{}, \"{}\", Missing {}.", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name, cAlphaFields(17)));
                    ErrorsFound = true;
                }
                TestCompSet(state, state.dataHVACMultiSpdHP->CurrentModuleObject, Alphas(1), Alphas(16), Alphas(17), "MSHP Heat receovery Nodes");
                SetMSHPDXCoilHeatRecoveryFlag(state, thisMSHeatPump.DXCoolCoilIndex);
                if (thisMSHeatPump.DXHeatCoilIndex > 0) {
                    SetMSHPDXCoilHeatRecoveryFlag(state, thisMSHeatPump.DXHeatCoilIndex);
                }
            } else {
                thisMSHeatPump.HeatRecActive = false;
                thisMSHeatPump.DesignHeatRecMassFlowRate = 0.0;
                thisMSHeatPump.HeatRecInletNodeNum = 0;
                thisMSHeatPump.HeatRecOutletNodeNum = 0;
                if (!lAlphaBlanks(16) || !lAlphaBlanks(17)) {
                    ShowWarningError(state,
                                     format("Since {} = 0.0, heat recovery is inactive for {} = {}",
                                            cNumericFields(6),
                                            state.dataHVACMultiSpdHP->CurrentModuleObject,
                                            Alphas(1)));
                    ShowContinueError(state, format("However, {} or {} was specified.", cAlphaFields(16), cAlphaFields(17)));
                }
            }
            thisMSHeatPump.MaxHeatRecOutletTemp = Numbers(7);
            if (thisMSHeatPump.MaxHeatRecOutletTemp < 0.0) {
                ShowSevereError(state,
                                format("{}, \"{}\", The value for {} is below 0.0",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cNumericFields(7)));
                ErrorsFound = true;
            }
            if (thisMSHeatPump.MaxHeatRecOutletTemp > 100.0) {
                ShowSevereError(state,
                                format("{}, \"{}\", The value for {} is above 100.0",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cNumericFields(7)));
                ErrorsFound = true;
            }

            thisMSHeatPump.IdleVolumeAirRate = Numbers(8);
            if (thisMSHeatPump.IdleVolumeAirRate < 0.0 && thisMSHeatPump.IdleVolumeAirRate != AutoSize) {
                ShowSevereError(state,
                                format("{}, \"{}\", {} cannot be less than zero.",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       thisMSHeatPump.Name,
                                       cNumericFields(8)));
                ErrorsFound = true;
            }

            //     AirFlowControl only valid if fan opmode = ContFanCycCoil
            if (thisMSHeatPump.IdleVolumeAirRate == 0.0) {
                thisMSHeatPump.AirFlowControl = AirflowControl::UseCompressorOnFlow;
            } else {
                thisMSHeatPump.AirFlowControl = AirflowControl::UseCompressorOffFlow;
            }

            //   Initialize last mode of compressor operation
            thisMSHeatPump.LastMode = ModeOfOperation::HeatingMode;

            thisMSHeatPump.NumOfSpeedHeating = Numbers(9);
            if (thisMSHeatPump.NumOfSpeedHeating < 2 || thisMSHeatPump.NumOfSpeedHeating > 4) {
                if (thisMSHeatPump.HeatCoilType == MultiSpeedHeatingCoil) {
                    ShowSevereError(state,
                                    format("{}, The maximum {} is 4, and the minimum number is 2",
                                           state.dataHVACMultiSpdHP->CurrentModuleObject,
                                           cNumericFields(9)));
                    ShowContinueError(state, format("The input value is {:.0R}", Numbers(9)));
                    ErrorsFound = true;
                }
            }
            thisMSHeatPump.NumOfSpeedCooling = Numbers(10);
            if (thisMSHeatPump.NumOfSpeedCooling < 2 || thisMSHeatPump.NumOfSpeedCooling > 4) {
                ShowSevereError(state,
                                format("{}, The maximum {} is 4, and the minimum number is 2",
                                       state.dataHVACMultiSpdHP->CurrentModuleObject,
                                       cNumericFields(10)));
                ShowContinueError(state, format("The input value is {:.0R}", Numbers(10)));
                ErrorsFound = true;
            }

            // Generate a dynamic array for heating
            if (thisMSHeatPump.NumOfSpeedHeating > 0) {
                thisMSHeatPump.HeatMassFlowRate.allocate(thisMSHeatPump.NumOfSpeedHeating);
                thisMSHeatPump.HeatVolumeFlowRate.allocate(thisMSHeatPump.NumOfSpeedHeating);
                thisMSHeatPump.HeatingSpeedRatio.allocate(thisMSHeatPump.NumOfSpeedHeating);
                thisMSHeatPump.HeatingSpeedRatio = 1.0;
                for (int i = 1; i <= thisMSHeatPump.NumOfSpeedHeating; ++i) {
                    thisMSHeatPump.HeatVolumeFlowRate(i) = Numbers(10 + i);
                    if (thisMSHeatPump.HeatCoilType == MultiSpeedHeatingCoil) {
                        if (thisMSHeatPump.HeatVolumeFlowRate(i) <= 0.0 && thisMSHeatPump.HeatVolumeFlowRate(i) != AutoSize) {
                            ShowSevereError(state,
                                            format("{}, \"{}\", {} must be greater than zero.",
                                                   state.dataHVACMultiSpdHP->CurrentModuleObject,
                                                   thisMSHeatPump.Name,
                                                   cNumericFields(10 + i)));
                            ErrorsFound = true;
                        }
                    }
                }
                // Ensure flow rate at high speed should be greater or equal to the flow rate at low speed
                for (int i = 2; i <= thisMSHeatPump.NumOfSpeedHeating; ++i) {
                    if (thisMSHeatPump.HeatVolumeFlowRate(i) == AutoSize) continue;
                    Found = false;
		    int j;
                    for (j = i - 1; j >= 1; --j) {
                        if (thisMSHeatPump.HeatVolumeFlowRate(i) != AutoSize) {
                            Found = true;
                            break;
                        }
                    }
                    if (Found) {
                        if (thisMSHeatPump.HeatVolumeFlowRate(i) < thisMSHeatPump.HeatVolumeFlowRate(j)) {
                            ShowSevereError(state,
                                            format("{}, \"{}\", {}",
                                                   state.dataHVACMultiSpdHP->CurrentModuleObject,
                                                   thisMSHeatPump.Name,
                                                   cNumericFields(10 + i)));
                            ShowContinueError(state, format(" cannot be less than {}", cNumericFields(10 + j)));
                            ErrorsFound = true;
                        }
                    }
                }
            }

            if (state.dataGlobal->DoCoilDirectSolutions) {
                int MaxNumber = std::max(thisMSHeatPump.NumOfSpeedCooling, thisMSHeatPump.NumOfSpeedHeating);
                thisMSHeatPump.FullOutput.allocate(MaxNumber);
                DXCoils::DisableLatentDegradation(state, thisMSHeatPump.DXCoolCoilIndex);
            }
            // Generate a dynamic array for cooling
            if (thisMSHeatPump.NumOfSpeedCooling > 0) {
                thisMSHeatPump.CoolMassFlowRate.allocate(thisMSHeatPump.NumOfSpeedCooling);
                thisMSHeatPump.CoolVolumeFlowRate.allocate(thisMSHeatPump.NumOfSpeedCooling);
                thisMSHeatPump.CoolingSpeedRatio.allocate(thisMSHeatPump.NumOfSpeedCooling);
                thisMSHeatPump.CoolingSpeedRatio = 1.0;
                for (int i = 1; i <= thisMSHeatPump.NumOfSpeedCooling; ++i) {
                    thisMSHeatPump.CoolVolumeFlowRate(i) = Numbers(14 + i);
                    if (thisMSHeatPump.CoolVolumeFlowRate(i) <= 0.0 && thisMSHeatPump.CoolVolumeFlowRate(i) != AutoSize) {
                        ShowSevereError(state,
                                        format("{}, \"{}\", {} must be greater than zero.",
                                               state.dataHVACMultiSpdHP->CurrentModuleObject,
                                               thisMSHeatPump.Name,
                                               cNumericFields(14 + i)));
                        ErrorsFound = true;
                    }
                }
                // Ensure flow rate at high speed should be greater or equal to the flow rate at low speed
                for (int i = 2; i <= thisMSHeatPump.NumOfSpeedCooling; ++i) {
                    if (thisMSHeatPump.CoolVolumeFlowRate(i) == AutoSize) continue;
                    Found = false;
		    int j;
                    for (j = i - 1; j >= 1; --j) {
                        if (thisMSHeatPump.CoolVolumeFlowRate(i) != AutoSize) {
                            Found = true;
                            break;
                        }
                    }
                    if (Found) {
                        if (thisMSHeatPump.CoolVolumeFlowRate(i) < thisMSHeatPump.CoolVolumeFlowRate(j)) {
                            ShowSevereError(state,
                                            format("{}, \"{}\", {}",
                                                   state.dataHVACMultiSpdHP->CurrentModuleObject,
                                                   thisMSHeatPump.Name,
                                                   cNumericFields(14 + i)));
                            ShowContinueError(state, format(" cannot be less than {}", cNumericFields(14 + j)));
                            ErrorsFound = true;
                        }
                    }
                }
            }

            // Check node integrity
            if (thisMSHeatPump.FanPlaceType == BlowThru) {
                if (thisMSHeatPump.FanInletNode != thisMSHeatPump.AirInletNodeNum) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(
                        state, format("When a blow through fan is specified, the fan inlet node name must be the same as the {}", cAlphaFields(3)));
                    ShowContinueError(state,
                                      format("...Fan inlet node name           = {}", state.dataLoopNodes->NodeID(thisMSHeatPump.FanInletNode)));
                    ShowContinueError(state, format("...{} = {}", cAlphaFields(3), state.dataLoopNodes->NodeID(thisMSHeatPump.AirInletNodeNum)));
                    ErrorsFound = true;
                }
                if (thisMSHeatPump.FanOutletNode != CoolingCoilInletNode) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(
                        state,
                        "When a blow through fan is specified, the fan outlet node name must be the same as the cooling coil inlet node name.");
                    ShowContinueError(state,
                                      format("...Fan outlet node name         = {}", state.dataLoopNodes->NodeID(thisMSHeatPump.FanOutletNode)));
                    ShowContinueError(state, format("...Cooling coil inlet node name = {}", state.dataLoopNodes->NodeID(CoolingCoilInletNode)));
                    ErrorsFound = true;
                }
                if (CoolingCoilOutletNode != HeatingCoilInletNode) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state, "The cooling coil outlet node name must be the same as the heating coil inlet node name.");
                    ShowContinueError(state, format("...Cooling coil outlet node name = {}", state.dataLoopNodes->NodeID(CoolingCoilOutletNode)));
                    ShowContinueError(state, format("...Heating coil inlet node name  = {}", state.dataLoopNodes->NodeID(HeatingCoilInletNode)));
                    ErrorsFound = true;
                }
                if (HeatingCoilOutletNode != SuppHeatCoilInletNode) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state,
                                      "When a blow through fan is specified, the heating coil outlet node name must be the same as the reheat coil "
                                      "inlet node name.");
                    ShowContinueError(state, format("...Heating coil outlet node name = {}", state.dataLoopNodes->NodeID(HeatingCoilOutletNode)));
                    ShowContinueError(state, format("...Reheat coil inlet node name   = {}", state.dataLoopNodes->NodeID(SuppHeatCoilInletNode)));
                    ErrorsFound = true;
                }
                if (SuppHeatCoilOutletNode != thisMSHeatPump.AirOutletNodeNum) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state, format("The supplemental heating coil outlet node name must be the same as the {}", cAlphaFields(4)));
                    ShowContinueError(
                        state, format("...Supplemental heating coil outlet node name   = {}", state.dataLoopNodes->NodeID(SuppHeatCoilOutletNode)));
                    ShowContinueError(state,
                                      format("...{} = {}", cAlphaFields(4), state.dataLoopNodes->NodeID(thisMSHeatPump.AirOutletNodeNum)));
                    ErrorsFound = true;
                }
            } else {
                if (CoolingCoilInletNode != thisMSHeatPump.AirInletNodeNum) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(
                        state,
                        format("When a draw through fan is specified, the cooling coil inlet node name must be the same as the {}", cAlphaFields(3)));
                    ShowContinueError(state, format("...Cooling coil inlet node name  = {}", state.dataLoopNodes->NodeID(CoolingCoilInletNode)));
                    ShowContinueError(state, format("...{} = {}", cAlphaFields(3), state.dataLoopNodes->NodeID(thisMSHeatPump.AirInletNodeNum)));
                    ErrorsFound = true;
                }
                if (CoolingCoilOutletNode != HeatingCoilInletNode) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state, "The cooling coil outlet node name must be the same as the heating coil inlet node name.");
                    ShowContinueError(state, format("...Cooling coil outlet node name = {}", state.dataLoopNodes->NodeID(CoolingCoilOutletNode)));
                    ShowContinueError(state, format("...Heating coil inlet node name  = {}", state.dataLoopNodes->NodeID(HeatingCoilInletNode)));
                    ErrorsFound = true;
                }
                if (HeatingCoilOutletNode != thisMSHeatPump.FanInletNode) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(
                        state,
                        "When a draw through fan is specified, the heating coil outlet node name must be the same as the fan inlet node name.");
                    ShowContinueError(state, format("...Heating coil outlet node name = {}", state.dataLoopNodes->NodeID(HeatingCoilOutletNode)));
                    ShowContinueError(state,
                                      format("...Fan inlet node name           = {}", state.dataLoopNodes->NodeID(thisMSHeatPump.FanInletNode)));
                    ErrorsFound = true;
                }
                if (thisMSHeatPump.FanOutletNode != SuppHeatCoilInletNode) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(
                        state, "When a draw through fan is specified, the fan outlet node name must be the same as the reheat coil inlet node name.");
                    ShowContinueError(state,
                                      format("...Fan outlet node name        = {}", state.dataLoopNodes->NodeID(thisMSHeatPump.FanOutletNode)));
                    ShowContinueError(state, format("...Reheat coil inlet node name = {}", state.dataLoopNodes->NodeID(SuppHeatCoilInletNode)));
                    ErrorsFound = true;
                }
                if (SuppHeatCoilOutletNode != thisMSHeatPump.AirOutletNodeNum) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state, format("The reheat coil outlet node name must be the same as the {}", cAlphaFields(4)));
                    ShowContinueError(state, format("...Reheat coil outlet node name   = {}", state.dataLoopNodes->NodeID(SuppHeatCoilOutletNode)));
                    ShowContinueError(state,
                                      format("...{} = {}", cAlphaFields(4), state.dataLoopNodes->NodeID(thisMSHeatPump.AirOutletNodeNum)));
                    ErrorsFound = true;
                }
            }

            // Ensure the numbers of speeds defined in the parent object are equal to the numbers defined in coil objects
            if (thisMSHeatPump.HeatCoilType == MultiSpeedHeatingCoil) {
                
                if (thisMSHeatPump.NumOfSpeedHeating != GetDXCoilNumberOfSpeeds(state, Alphas(10), Alphas(11), ErrorsFound)) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(
                        state, format("The {} is not equal to the number defined in {} = {}", cNumericFields(9), cAlphaFields(11), Alphas(11)));
                    ErrorsFound = true;
                }
            } else if (thisMSHeatPump.HeatCoilType == Coil_HeatingElectric_MultiStage ||
                       thisMSHeatPump.HeatCoilType == Coil_HeatingGas_MultiStage) {

                if (thisMSHeatPump.NumOfSpeedHeating != GetHeatingCoilNumberOfStages(state, Alphas(10), Alphas(11), ErrorsFound)) {
                    ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(
                        state, format("The {} is not equal to the number defined in {} = {}", cNumericFields(9), cAlphaFields(11), Alphas(11)));
                    ErrorsFound = true;
                }
            }

            if (thisMSHeatPump.NumOfSpeedCooling != GetDXCoilNumberOfSpeeds(state, Alphas(12), Alphas(13), ErrorsFound)) {
                ShowSevereError(state, format("For {} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                ShowContinueError(state,
                                  format("The {} is not equal to the number defined in {} = {}", cNumericFields(10), cAlphaFields(13), Alphas(13)));
                ErrorsFound = true;
            }
        }

        if (ErrorsFound) {
            ShowFatalError(state,
                           format("{}Errors found in getting {} input.  Preceding condition(s) causes termination.",
                                  RoutineName,
                                  state.dataHVACMultiSpdHP->CurrentModuleObject));
        }
        // End of multispeed heat pump

        for (int MSHPNum = 1; MSHPNum <= state.dataHVACMultiSpdHP->NumMSHeatPumps; ++MSHPNum) {
		
            auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHPNum);
            // Setup Report Variables for MSHP Equipment
            SetupOutputVariable(state,
                                "Unitary System Ancillary Electricity Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.AuxElecPower,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Cooling Ancillary Electricity Energy",
                                OutputProcessor::Unit::J,
                                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHPNum).AuxElecCoolConsumption,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Summed,
                                thisMSHeatPump.Name,
                                _,
                                "Electricity",
                                "Cooling",
                                _,
                                "System");
            SetupOutputVariable(state,
                                "Unitary System Heating Ancillary Electricity Energy",
                                OutputProcessor::Unit::J,
                                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHPNum).AuxElecHeatConsumption,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Summed,
                                thisMSHeatPump.Name,
                                _,
                                "Electricity",
                                "Heating",
                                _,
                                "System");
            SetupOutputVariable(state,
                                "Unitary System Fan Part Load Ratio",
                                OutputProcessor::Unit::None,
                                thisMSHeatPump.FanPartLoadRatio,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Compressor Part Load Ratio",
                                OutputProcessor::Unit::None,
                                thisMSHeatPump.CompPartLoadRatio,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Electricity Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.ElecPower,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Electricity Energy",
                                OutputProcessor::Unit::J,
                                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHPNum).ElecPowerConsumption,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Summed,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System DX Coil Cycling Ratio",
                                OutputProcessor::Unit::None,
                                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHPNum).CycRatio,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System DX Coil Speed Ratio",
                                OutputProcessor::Unit::None,
                                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHPNum).SpeedRatio,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System DX Coil Speed Level",
                                OutputProcessor::Unit::None,
                                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHPNum).SpeedNum,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Total Cooling Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.TotCoolEnergyRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Total Heating Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.TotHeatEnergyRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Sensible Cooling Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.SensCoolEnergyRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Sensible Heating Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.SensHeatEnergyRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Latent Cooling Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.LatCoolEnergyRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            SetupOutputVariable(state,
                                "Unitary System Latent Heating Rate",
                                OutputProcessor::Unit::W,
                                thisMSHeatPump.LatHeatEnergyRate,
                                OutputProcessor::SOVTimeStepType::System,
                                OutputProcessor::SOVStoreType::Average,
                                thisMSHeatPump.Name);
            if (thisMSHeatPump.HeatRecActive) {
                SetupOutputVariable(state,
                                    "Unitary System Heat Recovery Rate",
                                    OutputProcessor::Unit::W,
                                    thisMSHeatPump.HeatRecoveryRate,
                                    OutputProcessor::SOVTimeStepType::System,
                                    OutputProcessor::SOVStoreType::Average,
                                    thisMSHeatPump.Name);
                SetupOutputVariable(state,
                                    "Unitary System Heat Recovery Inlet Temperature",
                                    OutputProcessor::Unit::C,
                                    thisMSHeatPump.HeatRecoveryInletTemp,
                                    OutputProcessor::SOVTimeStepType::System,
                                    OutputProcessor::SOVStoreType::Average,
                                    thisMSHeatPump.Name);
                SetupOutputVariable(state,
                                    "Unitary System Heat Recovery Outlet Temperature",
                                    OutputProcessor::Unit::C,
                                    thisMSHeatPump.HeatRecoveryOutletTemp,
                                    OutputProcessor::SOVTimeStepType::System,
                                    OutputProcessor::SOVStoreType::Average,
                                    thisMSHeatPump.Name);
                SetupOutputVariable(state,
                                    "Unitary System Heat Recovery Fluid Mass Flow Rate",
                                    OutputProcessor::Unit::kg_s,
                                    thisMSHeatPump.HeatRecoveryMassFlowRate,
                                    OutputProcessor::SOVTimeStepType::System,
                                    OutputProcessor::SOVStoreType::Average,
                                    thisMSHeatPump.Name);
                SetupOutputVariable(state,
                                    "Unitary System Heat Recovery Energy",
                                    OutputProcessor::Unit::J,
                                    state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHPNum).HeatRecoveryEnergy,
                                    OutputProcessor::SOVTimeStepType::System,
                                    OutputProcessor::SOVStoreType::Summed,
                                    thisMSHeatPump.Name);
            }
        }
        if (state.dataGlobal->AnyEnergyManagementSystemInModel) {
            for (int coilNUM = 1; coilNUM <= int(state.dataHVACMultiSpdHP->MSHeatPump.size()); ++coilNUM) {
                SetupEMSActuator(state,
                                 "Coil Speed Control",
                                 state.dataHVACMultiSpdHP->MSHeatPump(coilNUM).Name,
                                 "Unitary System DX Coil Speed Value",
                                 "[]",
                                 state.dataHVACMultiSpdHP->MSHeatPump(coilNUM).EMSOverrideCoilSpeedNumOn,
                                 state.dataHVACMultiSpdHP->MSHeatPump(coilNUM).EMSOverrideCoilSpeedNumValue);
            }
        }
    }

    //******************************************************************************

    void InitMSHeatPump(EnergyPlusData &state,
                        int const MSHeatPumpNum,       // Engine driven heat pump number
                        bool const FirstHVACIteration, // TRUE if first HVAC iteration
                        int const AirLoopNum,          // air loop index
                        Real64 &QZnReq,                // Heating/Cooling load for all served zones
                        Real64 &OnOffAirFlowRatio      // Ratio of compressor ON airflow to average airflow over timestep
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Lixing Gu, FSEC
        //       DATE WRITTEN:    July 2007
        //       MODIFIED         Bereket Nigusse, June 2010 - added a procedure to calculate supply air flow fraction
        //                        through controlled zone
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for initializations of the multispeed heat pump (MSHP) components.

        // METHODOLOGY EMPLOYED:
        // Uses the status flags to trigger initializations. The MSHP system is simulated with no load (coils off) to
        // determine the outlet temperature. A setpoint temperature is calculated on FirstHVACIteration = TRUE.
        // Once the setpoint is calculated, the inlet mass flow rate on FirstHVACIteration = FALSE is used to
        // determine the bypass fraction. The simulation converges quickly on mass flow rate. If the zone
        // temperatures float in the deadband, additional iterations are required to converge on mass flow rate.

        // Using/Aliasing
        using DataSizing::AutoSize;
        using Fans::GetFanIndex;
        using Fans::GetFanVolFlow;
        using FluidProperties::GetDensityGlycol;
        using FluidProperties::GetSatDensityRefrig;

        using PlantUtilities::InitComponentNodes;
        using PlantUtilities::ScanPlantLoopsForObject;
        using PlantUtilities::SetComponentFlowRate;
        using Psychrometrics::PsyRhoAirFnPbTdbW;
        using ScheduleManager::GetCurrentScheduleValue;
        using SteamCoils::SimulateSteamCoilComponents;
        using DXCoils::GetDXCoilAvailSchPtr;
        using WaterCoils::GetCoilMaxWaterFlowRate;
        using WaterCoils::SimulateWaterCoilComponents;

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        static constexpr std::string_view RoutineName("InitMSHeatPump");
        int ZoneInNode; // Zone inlet node number in the controlled zone for MSHP
        Real64 RhoAir;  // Air density at InNode

        Real64 QSensUnitOut; // Output of MSHP system with coils off
        Real64 PartLoadFrac; // Part-load ratio
        int ZoneNum;
        Real64 DeltaMassRate;  // Difference of mass flow rate between inlet node and system outlet node

        int ZoneInletNodeNum(0);                          // zone inlet nodes node number
        Real64 SumOfMassFlowRateMax(0.0);                 // the sum of mass flow rates at inlet to zones in an airloop
        Real64 CntrlZoneTerminalUnitMassFlowRateMax(0.0); // Maximum mass flow rate through controlled zone terminal unit
        bool errFlag;
        Real64 rho;    // local fluid density
        Real64 MdotHR; // local temporary for heat recovery fluid mass flow rate (kg/s)
        Real64 ZoneLoadToCoolSPSequenced;
        Real64 ZoneLoadToHeatSPSequenced;

        bool ErrorsFound(false);        // flag returned from mining call
        int SteamIndex(0);              // index of steam quality for steam heating coil
        Real64 mdot(0.0);               // local temporary for mass flow rate (kg/s)
        Real64 SteamDensity(0.0);       // density of steam at 100C, used for steam heating coils
        Real64 CoilMaxVolFlowRate(0.0); // coil fluid maximum volume flow rate
        Real64 QActual(0.0);            // coil actual capacity
        int CoilAvailSchPtr(0);         // DX coil availability schedule pointer

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        int InNode = thisMSHeatPump.AirInletNodeNum;
        int OutNode = thisMSHeatPump.AirOutletNodeNum;
        int NumOfSpeedCooling = thisMSHeatPump.NumOfSpeedCooling;
        int NumOfSpeedHeating = thisMSHeatPump.NumOfSpeedHeating;

        ++state.dataHVACMultiSpdHP->AirLoopPass;
        if (state.dataHVACMultiSpdHP->AirLoopPass > 2) state.dataHVACMultiSpdHP->AirLoopPass = 1;

        if (thisMSHeatPump.MyPlantScantFlag && allocated(state.dataPlnt->PlantLoop)) {
            if (thisMSHeatPump.HeatRecActive) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisMSHeatPump.Name,
                                        DataPlant::PlantEquipmentType::MultiSpeedHeatPumpRecovery,
                                        thisMSHeatPump.HRPlantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
                if (errFlag) {
                    ShowFatalError(state, "InitMSHeatPump: Program terminated for previous conditions.");
                }

                thisMSHeatPump.MyPlantScantFlag = false;
            } else {
                thisMSHeatPump.MyPlantScantFlag = false;
            }
            if (thisMSHeatPump.HeatCoilType == Coil_HeatingWater) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisMSHeatPump.HeatCoilName,
                                        DataPlant::PlantEquipmentType::CoilWaterSimpleHeating,
                                        thisMSHeatPump.plantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
                if (errFlag) {
                    ShowFatalError(state, "InitMSHeatPump: Program terminated for previous conditions.");
                }
                thisMSHeatPump.MaxCoilFluidFlow =
                    GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.HeatCoilName, ErrorsFound);

                if (thisMSHeatPump.MaxCoilFluidFlow > 0.0) {
                    rho = GetDensityGlycol(state,
                                           state.dataPlnt->PlantLoop(thisMSHeatPump.plantLoc.loopNum).FluidName,
                                           DataGlobalConstants::HWInitConvTemp,
                                           state.dataPlnt->PlantLoop(thisMSHeatPump.plantLoc.loopNum).FluidIndex,
                                           RoutineName);
                    thisMSHeatPump.MaxCoilFluidFlow =
                        GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.HeatCoilName, ErrorsFound) * rho;
                }
                // fill outlet node for coil
                thisMSHeatPump.CoilOutletNode =
                    DataPlant::CompData::getPlantComponent(state, thisMSHeatPump.plantLoc).NodeNumOut;
                thisMSHeatPump.MyPlantScantFlag = false;

            } else if (thisMSHeatPump.HeatCoilType == Coil_HeatingSteam) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisMSHeatPump.HeatCoilName,
                                        DataPlant::PlantEquipmentType::CoilSteamAirHeating,
                                        thisMSHeatPump.plantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
                if (errFlag) {
                    ShowFatalError(state, "InitMSHeatPump: Program terminated for previous conditions.");
                }
                thisMSHeatPump.MaxCoilFluidFlow = SteamCoils::GetCoilMaxSteamFlowRate(state, thisMSHeatPump.HeatCoilNum, ErrorsFound);
                if (thisMSHeatPump.MaxCoilFluidFlow > 0.0) {
                    SteamIndex =
                        0; // Function GetSatDensityRefrig will look up steam index if 0 is passed // TODO: Why do you want to re-look this up?
                    SteamDensity = GetSatDensityRefrig(state, fluidNameSteam, state.dataHVACMultiSpdHP->TempSteamIn, 1.0, SteamIndex, RoutineName);
                    thisMSHeatPump.MaxCoilFluidFlow *= SteamDensity;
                }

                // fill outlet node for coil
                thisMSHeatPump.CoilOutletNode =
                    DataPlant::CompData::getPlantComponent(state, thisMSHeatPump.plantLoc).NodeNumOut;
                thisMSHeatPump.MyPlantScantFlag = false;
            }
            if (thisMSHeatPump.SuppHeatCoilType == Coil_HeatingWater) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisMSHeatPump.SuppHeatCoilName,
                                        DataPlant::PlantEquipmentType::CoilWaterSimpleHeating,
                                        thisMSHeatPump.SuppPlantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
                if (errFlag) {
                    ShowFatalError(state, "InitMSHeatPump: Program terminated for previous conditions.");
                }
                thisMSHeatPump.MaxSuppCoilFluidFlow =
                    GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.SuppHeatCoilName, ErrorsFound);

                if (thisMSHeatPump.MaxSuppCoilFluidFlow > 0.0) {
                    rho = GetDensityGlycol(state,
                                           state.dataPlnt->PlantLoop(thisMSHeatPump.SuppPlantLoc.loopNum).FluidName,
                                           DataGlobalConstants::HWInitConvTemp,
                                           state.dataPlnt->PlantLoop(thisMSHeatPump.SuppPlantLoc.loopNum).FluidIndex,
                                           RoutineName);
                    thisMSHeatPump.MaxSuppCoilFluidFlow =
                        GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.SuppHeatCoilName, ErrorsFound) * rho;
                }
                // fill outlet node for coil
                thisMSHeatPump.SuppCoilOutletNode =
                    DataPlant::CompData::getPlantComponent(state, thisMSHeatPump.SuppPlantLoc).NodeNumOut;
                thisMSHeatPump.MyPlantScantFlag = false;

            } else if (thisMSHeatPump.SuppHeatCoilType == Coil_HeatingSteam) {
                errFlag = false;
                ScanPlantLoopsForObject(state,
                                        thisMSHeatPump.SuppHeatCoilName,
                                        DataPlant::PlantEquipmentType::CoilSteamAirHeating,
                                        thisMSHeatPump.SuppPlantLoc,
                                        errFlag,
                                        _,
                                        _,
                                        _,
                                        _,
                                        _);
                if (errFlag) {
                    ShowFatalError(state, "InitMSHeatPump: Program terminated for previous conditions.");
                }
                thisMSHeatPump.MaxSuppCoilFluidFlow =
                    SteamCoils::GetCoilMaxSteamFlowRate(state, thisMSHeatPump.SuppHeatCoilNum, ErrorsFound);
                if (thisMSHeatPump.MaxSuppCoilFluidFlow > 0.0) {
                    SteamIndex = 0; // Function GetSatDensityRefrig will look up steam index if 0 is passed
                    SteamDensity = GetSatDensityRefrig(state, fluidNameSteam, state.dataHVACMultiSpdHP->TempSteamIn, 1.0, SteamIndex, RoutineName);
                    thisMSHeatPump.MaxSuppCoilFluidFlow *= SteamDensity;
                }

                // fill outlet node for coil
                thisMSHeatPump.SuppCoilOutletNode =
                    DataPlant::CompData::getPlantComponent(state, thisMSHeatPump.SuppPlantLoc).NodeNumOut;
                thisMSHeatPump.MyPlantScantFlag = false;
            }
        } else if (thisMSHeatPump.MyPlantScantFlag && !state.dataGlobal->AnyPlantInModel) {
            thisMSHeatPump.MyPlantScantFlag = false;
        }

        if (!state.dataGlobal->SysSizingCalc && thisMSHeatPump.MySizeFlag) {
            GetFanVolFlow(state, thisMSHeatPump.FanNum, thisMSHeatPump.FanVolFlow);
            SizeMSHeatPump(state, MSHeatPumpNum);
            thisMSHeatPump.FlowFraction = 1.0;
            thisMSHeatPump.MySizeFlag = false;
            // Pass the fan cycling schedule index up to the air loop. Set the air loop unitary system flag.
            state.dataAirLoop->AirLoopControlInfo(AirLoopNum).CycFanSchedPtr = thisMSHeatPump.FanSchedPtr;
            state.dataAirLoop->AirLoopControlInfo(AirLoopNum).UnitarySys = true;
            state.dataAirLoop->AirLoopControlInfo(AirLoopNum).UnitarySysSimulating =
                false; // affects child coil sizing by allowing coil to size itself instead of parent telling coil what size to use
            state.dataAirLoop->AirLoopControlInfo(AirLoopNum).FanOpMode = thisMSHeatPump.OpMode;
        }

        if (allocated(state.dataZoneEquip->ZoneEquipConfig) && thisMSHeatPump.MyCheckFlag) {
            int zoneNum = thisMSHeatPump.ControlZoneNum;
            int zoneInlet = thisMSHeatPump.ZoneInletNode;
            int coolingPriority = 0;
            int heatingPriority = 0;
            // setup furnace zone equipment sequence information based on finding matching air terminal
            if (state.dataZoneEquip->ZoneEquipConfig(zoneNum).EquipListIndex > 0) {
                state.dataZoneEquip->ZoneEquipList(state.dataZoneEquip->ZoneEquipConfig(zoneNum).EquipListIndex)
                    .getPrioritiesForInletNode(state, zoneInlet, coolingPriority, heatingPriority);
                thisMSHeatPump.ZoneSequenceCoolingNum = coolingPriority;
                thisMSHeatPump.ZoneSequenceHeatingNum = heatingPriority;
            }
            thisMSHeatPump.MyCheckFlag = false;
            if (thisMSHeatPump.ZoneSequenceCoolingNum == 0 || thisMSHeatPump.ZoneSequenceHeatingNum == 0) {
                ShowSevereError(state,
                                format("AirLoopHVAC:UnitaryHeatPump:AirToAir:MultiSpeed, \"{}\": Airloop air terminal in the zone equipment list for "
                                       "zone = {} not found or is not allowed Zone Equipment Cooling or Heating Sequence = 0.",
                                       thisMSHeatPump.Name,
                                       thisMSHeatPump.ControlZoneName));
                ShowFatalError(state,
                               "Subroutine InitMSHeatPump: Errors found in getting AirLoopHVAC:UnitaryHeatPump:AirToAir:MultiSpeed input.  Preceding "
                               "condition(s) causes termination.");
            }
        }

        // Find the number of zones (zone Inlet Nodes) attached to an air loop from the air loop number
        int NumAirLoopZones =
            state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).NumZonesCooled + state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).NumZonesHeated;
        if (allocated(state.dataAirLoop->AirToZoneNodeInfo) && thisMSHeatPump.MyFlowFracFlag) {
            state.dataHVACMultiSpdHP->FlowFracFlagReady = true;
            for (int ZoneInSysIndex = 1; ZoneInSysIndex <= NumAirLoopZones; ++ZoneInSysIndex) {
                // zone inlet nodes for cooling
                if (state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).NumZonesCooled > 0) {
                    if (state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).TermUnitCoolInletNodes(ZoneInSysIndex) == -999) {
                        // the data structure for the zones inlet nodes has not been filled
                        state.dataHVACMultiSpdHP->FlowFracFlagReady = false;
                    }
                }
                // zone inlet nodes for heating
                if (state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).NumZonesHeated > 0) {
                    if (state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).TermUnitHeatInletNodes(ZoneInSysIndex) == -999) {
                        // the data structure for the zones inlet nodes has not been filled
                        state.dataHVACMultiSpdHP->FlowFracFlagReady = false;
                    }
                }
            }
        }
        if (allocated(state.dataAirLoop->AirToZoneNodeInfo) && state.dataHVACMultiSpdHP->FlowFracFlagReady) {
            SumOfMassFlowRateMax = 0.0; // initialize the sum of the maximum flows
            for (int ZoneInSysIndex = 1; ZoneInSysIndex <= NumAirLoopZones; ++ZoneInSysIndex) {
                ZoneInletNodeNum = state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).TermUnitCoolInletNodes(ZoneInSysIndex);
                SumOfMassFlowRateMax += state.dataLoopNodes->Node(ZoneInletNodeNum).MassFlowRateMax;
                if (state.dataAirLoop->AirToZoneNodeInfo(AirLoopNum).CoolCtrlZoneNums(ZoneInSysIndex) == thisMSHeatPump.ControlZoneNum) {
                    CntrlZoneTerminalUnitMassFlowRateMax = state.dataLoopNodes->Node(ZoneInletNodeNum).MassFlowRateMax;
                }
            }
            if (SumOfMassFlowRateMax != 0.0 && thisMSHeatPump.MyFlowFracFlag) {
                if (CntrlZoneTerminalUnitMassFlowRateMax >= SmallAirVolFlow) {
                    thisMSHeatPump.FlowFraction = CntrlZoneTerminalUnitMassFlowRateMax / SumOfMassFlowRateMax;
                } else {
                    ShowSevereError(state, format("{} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state, " The Fraction of Supply Air Flow That Goes Through the Controlling Zone is set to 1.");
                }
                BaseSizer::reportSizerOutput(state,
                                             state.dataHVACMultiSpdHP->CurrentModuleObject,
                                             thisMSHeatPump.Name,
                                             "Fraction of Supply Air Flow That Goes Through the Controlling Zone",
                                             thisMSHeatPump.FlowFraction);
                thisMSHeatPump.MyFlowFracFlag = false;
            }
        }

        // Do the Begin Environment initializations
        if (state.dataGlobal->BeginEnvrnFlag && thisMSHeatPump.MyEnvrnFlag) {
            RhoAir = state.dataEnvrn->StdRhoAir;
            // set the mass flow rates from the input volume flow rates
            for (int i = 1; i <= NumOfSpeedCooling; ++i) {
                thisMSHeatPump.CoolMassFlowRate(i) = RhoAir * thisMSHeatPump.CoolVolumeFlowRate(i);
            }
            for (int i = 1; i <= NumOfSpeedHeating; ++i) {
                thisMSHeatPump.HeatMassFlowRate(i) = RhoAir * thisMSHeatPump.HeatVolumeFlowRate(i);
            }
            thisMSHeatPump.IdleMassFlowRate = RhoAir * thisMSHeatPump.IdleVolumeAirRate;
            // set the node max and min mass flow rates
            state.dataLoopNodes->Node(InNode).MassFlowRateMax =
                max(thisMSHeatPump.CoolMassFlowRate(NumOfSpeedCooling), thisMSHeatPump.HeatMassFlowRate(NumOfSpeedHeating));
            state.dataLoopNodes->Node(InNode).MassFlowRateMaxAvail =
                max(thisMSHeatPump.CoolMassFlowRate(NumOfSpeedCooling), thisMSHeatPump.HeatMassFlowRate(NumOfSpeedHeating));
            state.dataLoopNodes->Node(InNode).MassFlowRateMin = 0.0;
            state.dataLoopNodes->Node(InNode).MassFlowRateMinAvail = 0.0;
            state.dataLoopNodes->Node(OutNode) = state.dataLoopNodes->Node(InNode);
            thisMSHeatPump.LoadLoss = 0.0;

            if ((thisMSHeatPump.HeatRecActive) && (!thisMSHeatPump.MyPlantScantFlag)) {

                rho = GetDensityGlycol(state,
                                       state.dataPlnt->PlantLoop(thisMSHeatPump.HRPlantLoc.loopNum).FluidName,
                                       DataGlobalConstants::HWInitConvTemp,
                                       state.dataPlnt->PlantLoop(thisMSHeatPump.HRPlantLoc.loopNum).FluidIndex,
                                       RoutineName);

                thisMSHeatPump.DesignHeatRecMassFlowRate = thisMSHeatPump.DesignHeatRecFlowRate * rho;

                InitComponentNodes(state,
                                   0.0,
                                   thisMSHeatPump.DesignHeatRecMassFlowRate,
                                   thisMSHeatPump.HeatRecInletNodeNum,
                                   thisMSHeatPump.HeatRecOutletNodeNum);
            }
            if (thisMSHeatPump.CoilControlNode > 0) {
                if (thisMSHeatPump.MaxCoilFluidFlow == AutoSize) {
                    if (thisMSHeatPump.HeatCoilType == Coil_HeatingWater) {
                        SimulateWaterCoilComponents(
                            state, thisMSHeatPump.HeatCoilName, FirstHVACIteration, thisMSHeatPump.HeatCoilNum);

                        CoilMaxVolFlowRate =
                            GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.HeatCoilName, ErrorsFound);
                        if (CoilMaxVolFlowRate != AutoSize) {
                            rho = GetDensityGlycol(state,
                                                   state.dataPlnt->PlantLoop(thisMSHeatPump.plantLoc.loopNum).FluidName,
                                                   DataGlobalConstants::HWInitConvTemp,
                                                   state.dataPlnt->PlantLoop(thisMSHeatPump.plantLoc.loopNum).FluidIndex,
                                                   RoutineName);
                            thisMSHeatPump.MaxCoilFluidFlow = CoilMaxVolFlowRate * rho;
                        }
                        InitComponentNodes(state,
                                           0.0,
                                           thisMSHeatPump.MaxCoilFluidFlow,
                                           thisMSHeatPump.CoilControlNode,
                                           thisMSHeatPump.CoilOutletNode);
                    }
                    if (thisMSHeatPump.HeatCoilType == Coil_HeatingSteam) {

                        SimulateSteamCoilComponents(state,
                                                    thisMSHeatPump.HeatCoilName,
                                                    FirstHVACIteration,
                                                    thisMSHeatPump.HeatCoilNum,
                                                    1.0,
                                                    QActual); // QCoilReq, simulate any load > 0 to get max capacity of steam coil
                        CoilMaxVolFlowRate = SteamCoils::GetCoilMaxSteamFlowRate(state, thisMSHeatPump.HeatCoilNum, ErrorsFound);

                        if (CoilMaxVolFlowRate != AutoSize) {
                            SteamIndex = 0; // Function GetSatDensityRefrig will look up steam index if 0 is passed
                            SteamDensity =
                                GetSatDensityRefrig(state, fluidNameSteam, state.dataHVACMultiSpdHP->TempSteamIn, 1.0, SteamIndex, RoutineName);
                            thisMSHeatPump.MaxCoilFluidFlow = CoilMaxVolFlowRate * SteamDensity;
                        }
                        InitComponentNodes(state,
                                           0.0,
                                           thisMSHeatPump.MaxCoilFluidFlow,
                                           thisMSHeatPump.CoilControlNode,
                                           thisMSHeatPump.CoilOutletNode);
                    }
                }
            }
            if (thisMSHeatPump.SuppCoilControlNode > 0) {
                if (thisMSHeatPump.MaxSuppCoilFluidFlow == AutoSize) {
                    if (thisMSHeatPump.SuppHeatCoilType == Coil_HeatingWater) {
                        SimulateWaterCoilComponents(
                            state, thisMSHeatPump.SuppHeatCoilName, FirstHVACIteration, thisMSHeatPump.SuppHeatCoilNum);

                        CoilMaxVolFlowRate =
                            GetCoilMaxWaterFlowRate(state, "Coil:Heating:Water", thisMSHeatPump.SuppHeatCoilName, ErrorsFound);
                        if (CoilMaxVolFlowRate != AutoSize) {
                            rho = GetDensityGlycol(state,
                                                   state.dataPlnt->PlantLoop(thisMSHeatPump.SuppPlantLoc.loopNum).FluidName,
                                                   DataGlobalConstants::HWInitConvTemp,
                                                   state.dataPlnt->PlantLoop(thisMSHeatPump.SuppPlantLoc.loopNum).FluidIndex,
                                                   RoutineName);
                            thisMSHeatPump.MaxSuppCoilFluidFlow = CoilMaxVolFlowRate * rho;
                        }
                        InitComponentNodes(state,
                                           0.0,
                                           thisMSHeatPump.MaxSuppCoilFluidFlow,
                                           thisMSHeatPump.SuppCoilControlNode,
                                           thisMSHeatPump.SuppCoilOutletNode);
                    }
                    if (thisMSHeatPump.SuppHeatCoilType == Coil_HeatingSteam) {

                        SimulateSteamCoilComponents(state,
                                                    thisMSHeatPump.SuppHeatCoilName,
                                                    FirstHVACIteration,
                                                    thisMSHeatPump.SuppHeatCoilNum,
                                                    1.0,
                                                    QActual); // QCoilReq, simulate any load > 0 to get max capacity of steam coil
                        CoilMaxVolFlowRate = SteamCoils::GetCoilMaxSteamFlowRate(state, thisMSHeatPump.SuppHeatCoilNum, ErrorsFound);

                        if (CoilMaxVolFlowRate != AutoSize) {
                            SteamIndex = 0; // Function GetSatDensityRefrig will look up steam index if 0 is passed
                            SteamDensity =
                                GetSatDensityRefrig(state, fluidNameSteam, state.dataHVACMultiSpdHP->TempSteamIn, 1.0, SteamIndex, RoutineName);
                            thisMSHeatPump.MaxSuppCoilFluidFlow = CoilMaxVolFlowRate * SteamDensity;
                        }
                        InitComponentNodes(state,
                                           0.0,
                                           thisMSHeatPump.MaxSuppCoilFluidFlow,
                                           thisMSHeatPump.SuppCoilControlNode,
                                           thisMSHeatPump.SuppCoilOutletNode);
                    }
                }
            }
            thisMSHeatPump.MyEnvrnFlag = false;
        } // end one time inits

        if (!state.dataGlobal->BeginEnvrnFlag) {
            thisMSHeatPump.MyEnvrnFlag = true;
        }

        // IF MSHP system was not autosized and the fan is autosized, check that fan volumetric flow rate is greater than MSHP flow rates
        if (!state.dataGlobal->DoingSizing && thisMSHeatPump.CheckFanFlow) {
            state.dataHVACMultiSpdHP->CurrentModuleObject = "AirLoopHVAC:UnitaryHeatPump:AirToAir:MultiSpeed";
            GetFanVolFlow(state, thisMSHeatPump.FanNum, thisMSHeatPump.FanVolFlow);
            if (thisMSHeatPump.FanVolFlow != AutoSize) {
                //     Check fan versus system supply air flow rates
                if (thisMSHeatPump.FanVolFlow < thisMSHeatPump.CoolVolumeFlowRate(NumOfSpeedCooling)) {
                    ShowWarningError(state,
                                     format("{} - air flow rate = {:.7T} in fan object {} is less than the MSHP system air flow rate when cooling is "
                                            "required ({:.7T}).",
                                            state.dataHVACMultiSpdHP->CurrentModuleObject,
                                            thisMSHeatPump.FanVolFlow,
                                            thisMSHeatPump.FanName,
                                            thisMSHeatPump.CoolVolumeFlowRate(NumOfSpeedCooling)));
                    ShowContinueError(
                        state, " The MSHP system flow rate when cooling is required is reset to the fan flow rate and the simulation continues.");
                    ShowContinueError(state,
                                      format(" Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    thisMSHeatPump.CoolVolumeFlowRate(NumOfSpeedCooling) = thisMSHeatPump.FanVolFlow;
                    // Check flow rates in other speeds and ensure flow rates are not above the max flow rate
                    for (int i = NumOfSpeedCooling - 1; i >= 1; --i) {
                        if (thisMSHeatPump.CoolVolumeFlowRate(i) > thisMSHeatPump.CoolVolumeFlowRate(i + 1)) {
                            ShowContinueError(state,
                                              format(" The MSHP system flow rate when cooling is required is reset to the flow rate at higher speed "
                                                     "and the simulation continues at Speed{}.",
                                                     i));
                            ShowContinueError(
                                state, format(" Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                            thisMSHeatPump.CoolVolumeFlowRate(i) = thisMSHeatPump.CoolVolumeFlowRate(i + 1);
                        }
                    }
                }
                if (thisMSHeatPump.FanVolFlow < thisMSHeatPump.HeatVolumeFlowRate(NumOfSpeedHeating)) {
                    ShowWarningError(state,
                                     format("{} - air flow rate = {:.7T} in fan object {} is less than the MSHP system air flow rate when heating is "
                                            "required ({:.7T}).",
                                            state.dataHVACMultiSpdHP->CurrentModuleObject,
                                            thisMSHeatPump.FanVolFlow,
                                            thisMSHeatPump.FanName,
                                            thisMSHeatPump.HeatVolumeFlowRate(NumOfSpeedHeating)));
                    ShowContinueError(
                        state, " The MSHP system flow rate when heating is required is reset to the fan flow rate and the simulation continues.");
                    ShowContinueError(state,
                                      format(" Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    thisMSHeatPump.HeatVolumeFlowRate(NumOfSpeedHeating) = thisMSHeatPump.FanVolFlow;
                    for (int i = NumOfSpeedHeating - 1; i >= 1; --i) {
                        if (thisMSHeatPump.HeatVolumeFlowRate(i) > thisMSHeatPump.HeatVolumeFlowRate(i + 1)) {
                            ShowContinueError(state,
                                              format(" The MSHP system flow rate when heating is required is reset to the flow rate at higher speed "
                                                     "and the simulation continues at Speed{}.",
                                                     i));
                            ShowContinueError(
                                state,
                                format(" Occurs in {} system = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                            thisMSHeatPump.HeatVolumeFlowRate(i) = thisMSHeatPump.HeatVolumeFlowRate(i + 1);
                        }
                    }
                }
                if (thisMSHeatPump.FanVolFlow < thisMSHeatPump.IdleVolumeAirRate &&
                    thisMSHeatPump.IdleVolumeAirRate != 0.0) {
                    ShowWarningError(state,
                                     format("{} - air flow rate = {:.7T} in fan object {} is less than the MSHP system air flow rate when no heating "
                                            "or cooling is needed ({:.7T}).",
                                            state.dataHVACMultiSpdHP->CurrentModuleObject,
                                            thisMSHeatPump.FanVolFlow,
                                            thisMSHeatPump.FanName,
                                            thisMSHeatPump.IdleVolumeAirRate));
                    ShowContinueError(state,
                                      " The MSHP system flow rate when no heating or cooling is needed is reset to the fan flow rate and the "
                                      "simulation continues.");
                    ShowContinueError(state,
                                      format(" Occurs in {} = {}", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    thisMSHeatPump.IdleVolumeAirRate = thisMSHeatPump.FanVolFlow;
                }
                RhoAir = state.dataEnvrn->StdRhoAir;
                // set the mass flow rates from the reset volume flow rates
                for (int i = 1; i <= NumOfSpeedCooling; ++i) {
                    thisMSHeatPump.CoolMassFlowRate(i) = RhoAir * thisMSHeatPump.CoolVolumeFlowRate(i);
                    if (thisMSHeatPump.FanVolFlow > 0.0) {
                        thisMSHeatPump.CoolingSpeedRatio(i) =
                            thisMSHeatPump.CoolVolumeFlowRate(i) / thisMSHeatPump.FanVolFlow;
                    }
                }
                for (int i = 1; i <= NumOfSpeedHeating; ++i) {
                    thisMSHeatPump.HeatMassFlowRate(i) = RhoAir * thisMSHeatPump.HeatVolumeFlowRate(i);
                    if (thisMSHeatPump.FanVolFlow > 0.0) {
                        thisMSHeatPump.HeatingSpeedRatio(i) =
                            thisMSHeatPump.HeatVolumeFlowRate(i) / thisMSHeatPump.FanVolFlow;
                    }
                }
                thisMSHeatPump.IdleMassFlowRate = RhoAir * thisMSHeatPump.IdleVolumeAirRate;
                if (thisMSHeatPump.FanVolFlow > 0.0) {
                    thisMSHeatPump.IdleSpeedRatio = thisMSHeatPump.IdleVolumeAirRate / thisMSHeatPump.FanVolFlow;
                }
                // set the node max and min mass flow rates based on reset volume flow rates
                state.dataLoopNodes->Node(InNode).MassFlowRateMax =
                    max(thisMSHeatPump.CoolMassFlowRate(NumOfSpeedCooling), thisMSHeatPump.HeatMassFlowRate(NumOfSpeedHeating));
                state.dataLoopNodes->Node(InNode).MassFlowRateMaxAvail =
                    max(thisMSHeatPump.CoolMassFlowRate(NumOfSpeedCooling), thisMSHeatPump.HeatMassFlowRate(NumOfSpeedHeating));
                state.dataLoopNodes->Node(InNode).MassFlowRateMin = 0.0;
                state.dataLoopNodes->Node(InNode).MassFlowRateMinAvail = 0.0;
                state.dataLoopNodes->Node(OutNode) = state.dataLoopNodes->Node(InNode);
                thisMSHeatPump.CheckFanFlow = false;
            }
        }

        if (thisMSHeatPump.FanSchedPtr > 0) {
            if (GetCurrentScheduleValue(state, thisMSHeatPump.FanSchedPtr) == 0.0) {
                thisMSHeatPump.OpMode = CycFanCycCoil;
            } else {
                thisMSHeatPump.OpMode = ContFanCycCoil;
            }
        }

        // Calcuate air distribution losses
        if (!FirstHVACIteration && state.dataHVACMultiSpdHP->AirLoopPass == 1) {
            ZoneInNode = thisMSHeatPump.ZoneInletNode;
            DeltaMassRate = state.dataLoopNodes->Node(OutNode).MassFlowRate -
                            state.dataLoopNodes->Node(ZoneInNode).MassFlowRate / thisMSHeatPump.FlowFraction;
            if (DeltaMassRate < 0.0) DeltaMassRate = 0.0;
            Real64 MassFlowRate(0.0);        // parent mass flow rate
            Real64 LatentOutput(0.0);        // latent output rate
            Real64 TotalOutput(0.0);         // total output rate
            Real64 SensibleOutputDelta(0.0); // delta sensible output rate
            Real64 LatentOutputDelta(0.0);   // delta latent output rate
            Real64 TotalOutputDelta(0.0);    // delta total output rate
            MassFlowRate = state.dataLoopNodes->Node(ZoneInNode).MassFlowRate / thisMSHeatPump.FlowFraction;
            Real64 MinHumRat = state.dataLoopNodes->Node(ZoneInNode).HumRat;
            if (state.dataLoopNodes->Node(OutNode).Temp < state.dataLoopNodes->Node(thisMSHeatPump.NodeNumOfControlledZone).Temp)
                MinHumRat = state.dataLoopNodes->Node(OutNode).HumRat;
            CalcZoneSensibleLatentOutput(MassFlowRate,
                                         state.dataLoopNodes->Node(OutNode).Temp,
                                         MinHumRat,
                                         state.dataLoopNodes->Node(ZoneInNode).Temp,
                                         MinHumRat,
                                         thisMSHeatPump.LoadLoss,
                                         LatentOutput,
                                         TotalOutput);
            CalcZoneSensibleLatentOutput(DeltaMassRate,
                                         state.dataLoopNodes->Node(OutNode).Temp,
                                         MinHumRat,
                                         state.dataLoopNodes->Node(thisMSHeatPump.NodeNumOfControlledZone).Temp,
                                         MinHumRat,
                                         SensibleOutputDelta,
                                         LatentOutputDelta,
                                         TotalOutputDelta);
            thisMSHeatPump.LoadLoss = thisMSHeatPump.LoadLoss + SensibleOutputDelta;
            if (std::abs(thisMSHeatPump.LoadLoss) < 1.0e-6) thisMSHeatPump.LoadLoss = 0.0;
        }

        // Returns load only for zones requesting cooling (heating). If in deadband, Qzoneload = 0.
        ZoneNum = thisMSHeatPump.ControlZoneNum;
        if ((thisMSHeatPump.ZoneSequenceCoolingNum > 0) && (thisMSHeatPump.ZoneSequenceHeatingNum > 0)) {
            ZoneLoadToCoolSPSequenced = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(thisMSHeatPump.ControlZoneNum)
                                            .SequencedOutputRequiredToCoolingSP(thisMSHeatPump.ZoneSequenceCoolingNum);
            ZoneLoadToHeatSPSequenced = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(thisMSHeatPump.ControlZoneNum)
                                            .SequencedOutputRequiredToHeatingSP(thisMSHeatPump.ZoneSequenceHeatingNum);
            if (ZoneLoadToHeatSPSequenced > SmallLoad && ZoneLoadToCoolSPSequenced > SmallLoad) {
                QZnReq = ZoneLoadToHeatSPSequenced;
            } else if (ZoneLoadToHeatSPSequenced < (-1.0 * SmallLoad) && ZoneLoadToCoolSPSequenced < (-1.0 * SmallLoad)) {
                QZnReq = ZoneLoadToCoolSPSequenced;
            } else if (ZoneLoadToHeatSPSequenced <= (-1.0 * SmallLoad) && ZoneLoadToCoolSPSequenced >= SmallLoad) {
                QZnReq = 0.0;
            } else {
                QZnReq = 0.0; // Autodesk:Init Case added to prevent use of uninitialized value (occurred in MultiSpeedACFurnace example)
            }
            QZnReq /= thisMSHeatPump.FlowFraction;
        } else {
            QZnReq = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(ZoneNum).RemainingOutputRequired / thisMSHeatPump.FlowFraction;
        }
        if (state.dataZoneEnergyDemand->CurDeadBandOrSetback(ZoneNum)) QZnReq = 0.0;

        if (QZnReq > SmallLoad) {
            thisMSHeatPump.HeatCoolMode = ModeOfOperation::HeatingMode;
        } else if (QZnReq < (-1.0 * SmallLoad)) {
            thisMSHeatPump.HeatCoolMode = ModeOfOperation::CoolingMode;
        } else {
            thisMSHeatPump.HeatCoolMode = ModeOfOperation::Invalid;
        }

        // Determine the staged status
        if (allocated(state.dataZoneCtrls->StageZoneLogic)) {
            if (state.dataZoneCtrls->StageZoneLogic(ZoneNum)) {
                thisMSHeatPump.Staged = true;
                thisMSHeatPump.StageNum = state.dataZoneEnergyDemand->ZoneSysEnergyDemand(ZoneNum).StageNum;
            } else {
                if (thisMSHeatPump.MyStagedFlag) {
                    ShowWarningError(state,
                                     "ZoneControl:Thermostat:StagedDualSetpoint is found, but is not applied to this "
                                     "AirLoopHVAC:UnitaryHeatPump:AirToAir:MultiSpeed object = ");
                    ShowContinueError(state, format("{}. Please make correction. Simulation continues...", thisMSHeatPump.Name));
                    thisMSHeatPump.MyStagedFlag = false;
                }
            }
        }
        // Set the inlet node mass flow rate
        if (thisMSHeatPump.OpMode == ContFanCycCoil) {
            // constant fan mode
            if (QZnReq > SmallLoad && !state.dataZoneEnergyDemand->CurDeadBandOrSetback(ZoneNum)) {
                state.dataHVACMultiSpdHP->CompOnMassFlow = thisMSHeatPump.HeatMassFlowRate(1);
                state.dataHVACMultiSpdHP->CompOnFlowRatio = thisMSHeatPump.HeatingSpeedRatio(1);
                thisMSHeatPump.LastMode = ModeOfOperation::HeatingMode;
            } else if (QZnReq < (-1.0 * SmallLoad) && !state.dataZoneEnergyDemand->CurDeadBandOrSetback(ZoneNum)) {
                state.dataHVACMultiSpdHP->CompOnMassFlow = thisMSHeatPump.CoolMassFlowRate(1);
                state.dataHVACMultiSpdHP->CompOnFlowRatio = thisMSHeatPump.CoolingSpeedRatio(1);
                thisMSHeatPump.LastMode = ModeOfOperation::CoolingMode;
            } else {
                state.dataHVACMultiSpdHP->CompOnMassFlow = thisMSHeatPump.IdleMassFlowRate;
                state.dataHVACMultiSpdHP->CompOnFlowRatio = thisMSHeatPump.IdleSpeedRatio;
            }
            state.dataHVACMultiSpdHP->CompOffMassFlow = thisMSHeatPump.IdleMassFlowRate;
            state.dataHVACMultiSpdHP->CompOffFlowRatio = thisMSHeatPump.IdleSpeedRatio;
        } else {
            // cycling fan mode
            if (QZnReq > SmallLoad && !state.dataZoneEnergyDemand->CurDeadBandOrSetback(ZoneNum)) {
                state.dataHVACMultiSpdHP->CompOnMassFlow = thisMSHeatPump.HeatMassFlowRate(1);
                state.dataHVACMultiSpdHP->CompOnFlowRatio = thisMSHeatPump.HeatingSpeedRatio(1);
            } else if (QZnReq < (-1.0 * SmallLoad) && !state.dataZoneEnergyDemand->CurDeadBandOrSetback(ZoneNum)) {
                state.dataHVACMultiSpdHP->CompOnMassFlow = thisMSHeatPump.CoolMassFlowRate(1);
                state.dataHVACMultiSpdHP->CompOnFlowRatio = thisMSHeatPump.CoolingSpeedRatio(1);
            } else {
                state.dataHVACMultiSpdHP->CompOnMassFlow = 0.0;
                state.dataHVACMultiSpdHP->CompOnFlowRatio = 0.0;
            }
            state.dataHVACMultiSpdHP->CompOffMassFlow = 0.0;
            state.dataHVACMultiSpdHP->CompOffFlowRatio = 0.0;
        }

        // Set the inlet node mass flow rate
        if (GetCurrentScheduleValue(state, thisMSHeatPump.AvaiSchedPtr) > 0.0 && state.dataHVACMultiSpdHP->CompOnMassFlow != 0.0) {
            OnOffAirFlowRatio = 1.0;
            if (FirstHVACIteration) {
                state.dataLoopNodes->Node(thisMSHeatPump.AirInletNodeNum).MassFlowRate = state.dataHVACMultiSpdHP->CompOnMassFlow;
                PartLoadFrac = 0.0;
            } else {
                if (thisMSHeatPump.HeatCoolMode != ModeOfOperation::Invalid) {
                    PartLoadFrac = 1.0;
                } else {
                    PartLoadFrac = 0.0;
                }
            }
        } else {
            PartLoadFrac = 0.0;
            state.dataLoopNodes->Node(InNode).MassFlowRate = 0.0;
            state.dataLoopNodes->Node(OutNode).MassFlowRate = 0.0;
            state.dataLoopNodes->Node(OutNode).MassFlowRateMaxAvail = 0.0;
            OnOffAirFlowRatio = 1.0;
        }

        // Check availability of DX coils
        if (GetCurrentScheduleValue(state, thisMSHeatPump.AvaiSchedPtr) > 0.0) {
            if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::CoolingMode) {
                CoilAvailSchPtr = GetDXCoilAvailSchPtr( // TODO: Why isn't this stored on the struct?
                    state,
                    "Coil:Cooling:DX:MultiSpeed",
                    thisMSHeatPump.DXCoolCoilName,
                    ErrorsFound,
                    thisMSHeatPump.DXCoolCoilIndex);
                if (ErrorsFound) {
                    ShowFatalError(state, "InitMSHeatPump, The previous error causes termination.");
                }
                if (GetCurrentScheduleValue(state, CoilAvailSchPtr) == 0.0) {
                    if (thisMSHeatPump.CoolCountAvail == 0) {
                        ++thisMSHeatPump.CoolCountAvail;
                        ShowWarningError(
                            state,
                            format("{} is ready to perform cooling, but its DX cooling coil = {} is not available at Available Schedule = {}.",
                                   thisMSHeatPump.Name,
                                   thisMSHeatPump.DXCoolCoilName,
                                   GetScheduleName(state, CoilAvailSchPtr)));
                        ShowContinueErrorTimeStamp(state,
                                                   format("Availability schedule returned={:.1R}", GetCurrentScheduleValue(state, CoilAvailSchPtr)));
                    } else {
                        ++thisMSHeatPump.CoolCountAvail;
                        ShowRecurringWarningErrorAtEnd(state,
                                                       thisMSHeatPump.Name + ": Cooling coil is still not available ...",
                                                       thisMSHeatPump.CoolIndexAvail,
                                                       GetCurrentScheduleValue(state, CoilAvailSchPtr),
                                                       GetCurrentScheduleValue(state, CoilAvailSchPtr));
                    }
                }
            }
            if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::HeatingMode &&
                thisMSHeatPump.HeatCoilType == MultiSpeedHeatingCoil) {
                CoilAvailSchPtr = GetDXCoilAvailSchPtr(state,
                                                       "Coil:Heating:DX:MultiSpeed",
                                                       thisMSHeatPump.DXHeatCoilName,
                                                       ErrorsFound,
                                                       thisMSHeatPump.DXHeatCoilIndex);
                if (ErrorsFound) {
                    ShowFatalError(state, "InitMSHeatPump, The previous error causes termination.");
                }
                if (GetCurrentScheduleValue(state, CoilAvailSchPtr) == 0.0) {
                    if (thisMSHeatPump.HeatCountAvail == 0) {
                        ++thisMSHeatPump.HeatCountAvail;
                        ShowWarningError(
                            state,
                            format("{} is ready to perform heating, but its DX heating coil = {} is not available at Available Schedule = {}.",
                                   thisMSHeatPump.Name,
                                   thisMSHeatPump.DXCoolCoilName,
                                   GetScheduleName(state, CoilAvailSchPtr)));
                        ShowContinueErrorTimeStamp(state,
                                                   format("Availability schedule returned={:.1R}", GetCurrentScheduleValue(state, CoilAvailSchPtr)));
                    } else {
                        ++thisMSHeatPump.HeatCountAvail;
                        ShowRecurringWarningErrorAtEnd(state,
                                                       thisMSHeatPump.Name + ": Heating coil is still not available ...",
                                                       thisMSHeatPump.HeatIndexAvail,
                                                       GetCurrentScheduleValue(state, CoilAvailSchPtr),
                                                       GetCurrentScheduleValue(state, CoilAvailSchPtr));
                    }
                }
            }
        }

        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).CycRatio = 0.0;
        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedRatio = 0.0;
        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedNum = 0;

        CalcMSHeatPump(state,
                       MSHeatPumpNum,
                       FirstHVACIteration,
                       DataHVACGlobals::CompressorOperation::On,
                       1,
                       0.0,
                       PartLoadFrac,
                       QSensUnitOut,
                       QZnReq,
                       OnOffAirFlowRatio,
                       state.dataHVACMultiSpdHP->SupHeaterLoad);

        thisMSHeatPump.TotHeatEnergyRate = 0.0;
        thisMSHeatPump.SensHeatEnergyRate = 0.0;
        thisMSHeatPump.LatHeatEnergyRate = 0.0;
        thisMSHeatPump.TotCoolEnergyRate = 0.0;
        thisMSHeatPump.SensCoolEnergyRate = 0.0;
        thisMSHeatPump.LatCoolEnergyRate = 0.0;
        
        // If unit is scheduled OFF, setpoint is equal to inlet node temperature.
        //!!LKL Discrepancy with < 0
        if (GetCurrentScheduleValue(state, thisMSHeatPump.AvaiSchedPtr) == 0.0) {
            state.dataLoopNodes->Node(OutNode).Temp = state.dataLoopNodes->Node(InNode).Temp;
            return;
        }

        if ((thisMSHeatPump.HeatCoolMode == ModeOfOperation::Invalid && thisMSHeatPump.OpMode == CycFanCycCoil) ||
            state.dataHVACMultiSpdHP->CompOnMassFlow == 0.0) {
            QZnReq = 0.0;
            PartLoadFrac = 0.0;
            state.dataLoopNodes->Node(InNode).MassFlowRate = 0.0;
            state.dataLoopNodes->Node(OutNode).MassFlowRateMaxAvail = 0.0;
        }
        thisMSHeatPump.LoadMet = 0.0;
        SetAverageAirFlow(state, MSHeatPumpNum, PartLoadFrac, OnOffAirFlowRatio);

        // Init maximum available Heat Recovery flow rate
        if ((thisMSHeatPump.HeatRecActive) && (!thisMSHeatPump.MyPlantScantFlag)) {
            if (PartLoadFrac > 0.0) {
                if (FirstHVACIteration) {
                    MdotHR = thisMSHeatPump.DesignHeatRecMassFlowRate;
                } else {
                    if (thisMSHeatPump.HeatRecoveryMassFlowRate > 0.0) {
                        MdotHR = thisMSHeatPump.HeatRecoveryMassFlowRate;
                    } else {
                        MdotHR = thisMSHeatPump.DesignHeatRecMassFlowRate;
                    }
                }
            } else {
                MdotHR = 0.0;
            }

            SetComponentFlowRate(state,
                                 MdotHR,
                                 thisMSHeatPump.HeatRecInletNodeNum,
                                 thisMSHeatPump.HeatRecOutletNodeNum,
                                 thisMSHeatPump.HRPlantLoc);
        }

        // get operating capacity of water and steam coil
        if (FirstHVACIteration) {
            if (thisMSHeatPump.HeatCoilType == Coil_HeatingWater) {
                //     set air-side and steam-side mass flow rates
                state.dataLoopNodes->Node(thisMSHeatPump.CoilAirInletNode).MassFlowRate = state.dataHVACMultiSpdHP->CompOnMassFlow;
                mdot = thisMSHeatPump.MaxCoilFluidFlow;
                SetComponentFlowRate(state,
                                     mdot,
                                     thisMSHeatPump.CoilControlNode,
                                     thisMSHeatPump.CoilOutletNode,
                                     thisMSHeatPump.plantLoc);
                //     simulate water coil to find operating capacity
                SimulateWaterCoilComponents(
                    state, thisMSHeatPump.HeatCoilName, FirstHVACIteration, thisMSHeatPump.HeatCoilNum, QActual);
            } // from IF(thisMSHeatPump%HeatCoilType == Coil_HeatingWater) THEN

            if (thisMSHeatPump.HeatCoilType == Coil_HeatingSteam) {

                //     set air-side and steam-side mass flow rates
                state.dataLoopNodes->Node(thisMSHeatPump.CoilAirInletNode).MassFlowRate = state.dataHVACMultiSpdHP->CompOnMassFlow;
                mdot = thisMSHeatPump.MaxCoilFluidFlow;
                SetComponentFlowRate(state,
                                     mdot,
                                     thisMSHeatPump.CoilControlNode,
                                     thisMSHeatPump.CoilOutletNode,
                                     thisMSHeatPump.plantLoc);

                //     simulate steam coil to find operating capacity
                SimulateSteamCoilComponents(state,
                                            thisMSHeatPump.HeatCoilName,
                                            FirstHVACIteration,
                                            thisMSHeatPump.HeatCoilNum,
                                            1.0,
                                            QActual); // QCoilReq, simulate any load > 0 to get max capacity of steam coil

            } // from IF(thisMSHeatPump%HeatCoilType == Coil_HeatingSteam) THEN
            if (thisMSHeatPump.SuppHeatCoilType == Coil_HeatingWater) {
                //     set air-side and steam-side mass flow rates
                state.dataLoopNodes->Node(thisMSHeatPump.SuppCoilAirInletNode).MassFlowRate = state.dataHVACMultiSpdHP->CompOnMassFlow;
                mdot = thisMSHeatPump.MaxSuppCoilFluidFlow;
                SetComponentFlowRate(state,
                                     mdot,
                                     thisMSHeatPump.SuppCoilControlNode,
                                     thisMSHeatPump.SuppCoilOutletNode,
                                     thisMSHeatPump.SuppPlantLoc);
                //     simulate water coil to find operating capacity
                SimulateWaterCoilComponents(
                    state, thisMSHeatPump.SuppHeatCoilName, FirstHVACIteration, thisMSHeatPump.SuppHeatCoilNum, QActual);
                thisMSHeatPump.DesignSuppHeatingCapacity = QActual;

            } // from IF(thisMSHeatPump%SuppHeatCoilType == Coil_HeatingWater) THEN

            if (thisMSHeatPump.SuppHeatCoilType == Coil_HeatingSteam) {

                //     set air-side and steam-side mass flow rates
                state.dataLoopNodes->Node(thisMSHeatPump.SuppCoilAirInletNode).MassFlowRate = state.dataHVACMultiSpdHP->CompOnMassFlow;
                mdot = thisMSHeatPump.MaxSuppCoilFluidFlow;
                SetComponentFlowRate(state,
                                     mdot,
                                     thisMSHeatPump.SuppCoilControlNode,
                                     thisMSHeatPump.SuppCoilOutletNode,
                                     thisMSHeatPump.SuppPlantLoc);

                //     simulate steam coil to find operating capacity
                SimulateSteamCoilComponents(state,
                                            thisMSHeatPump.SuppHeatCoilName,
                                            FirstHVACIteration,
                                            thisMSHeatPump.SuppHeatCoilNum,
                                            1.0,
                                            QActual); // QCoilReq, simulate any load > 0 to get max capacity of steam coil
                thisMSHeatPump.DesignSuppHeatingCapacity =
                    SteamCoils::GetCoilCapacity(state, "Coil:Heating:Steam", thisMSHeatPump.SuppHeatCoilName, ErrorsFound);

            } // from IF(thisMSHeatPump%SuppHeatCoilType == Coil_HeatingSteam) THEN
        }     // from IF( FirstHVACIteration ) THEN
    }

    //******************************************************************************

    void SizeMSHeatPump(EnergyPlusData &state, int const MSHeatPumpNum) // Engine driven heat pump number
    {
        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Lixing Gu, FSEC
        //       DATE WRITTEN:    June 2007
        //       MODIFIED         na
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine is for sizing multispeed heat pump airflow rates and flow fraction.

        // Using/Aliasing
        using namespace DataSizing;

        using PlantUtilities::RegisterPlantCompDesignFlow;

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int NumOfSpeedCooling; // Number of speeds for cooling
        int NumOfSpeedHeating; // Number of speeds for heating

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);
	auto &thisPrimaryAirSystem = state.dataAirSystemsData->PrimaryAirSystems(state.dataSize->CurSysNum);
	
        if (state.dataSize->CurSysNum > 0 && state.dataSize->CurOASysNum == 0) {
            if (thisMSHeatPump.FanType == DataHVACGlobals::FanType_SystemModelObject) {
                thisPrimaryAirSystem.supFanVecIndex = thisMSHeatPump.FanNum;
                thisPrimaryAirSystem.supFanModelType = DataAirSystems::ObjectVectorOOFanSystemModel;
            } else {
                thisPrimaryAirSystem.SupFanNum = thisMSHeatPump.FanNum;
                thisPrimaryAirSystem.supFanModelType = DataAirSystems::StructArrayLegacyFanModels;
            }
            if (thisMSHeatPump.FanPlaceType == BlowThru) {
                thisPrimaryAirSystem.supFanLocation = DataAirSystems::FanPlacement::BlowThru;
            } else if (thisMSHeatPump.FanPlaceType == DrawThru) {
                thisPrimaryAirSystem.supFanLocation = DataAirSystems::FanPlacement::DrawThru;
            }
        }

        NumOfSpeedCooling = thisMSHeatPump.NumOfSpeedCooling;
        NumOfSpeedHeating = thisMSHeatPump.NumOfSpeedHeating;

        for (int i = NumOfSpeedCooling; i >= 1; --i) {

            if (thisMSHeatPump.CoolVolumeFlowRate(i) == AutoSize) {
                if (state.dataSize->CurSysNum > 0) {
                    if (i == NumOfSpeedCooling) {
                        CheckSysSizing(state, state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name);
                        thisMSHeatPump.CoolVolumeFlowRate(i) = state.dataSize->FinalSysSizing(state.dataSize->CurSysNum).DesMainVolFlow;
                        if (thisMSHeatPump.FanVolFlow < thisMSHeatPump.CoolVolumeFlowRate(i) &&
                            thisMSHeatPump.FanVolFlow != AutoSize) {
                            thisMSHeatPump.CoolVolumeFlowRate(i) = thisMSHeatPump.FanVolFlow;
                            ShowWarningError(state,
                                             format("{} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                            ShowContinueError(state,
                                              "The supply air flow rate at high speed is less than the autosized value for the supply air flow rate "
                                              "in cooling mode. Consider autosizing the fan for this simulation.");
                            ShowContinueError(
                                state,
                                "The air flow rate at high speed in cooling mode is reset to the supply air flow rate and the simulation continues.");
                        }
                    } else {
                        thisMSHeatPump.CoolVolumeFlowRate(i) =
                            thisMSHeatPump.CoolVolumeFlowRate(NumOfSpeedCooling) * i / NumOfSpeedCooling;
                    }
                    if (thisMSHeatPump.CoolVolumeFlowRate(i) < SmallAirVolFlow) {
                        thisMSHeatPump.CoolVolumeFlowRate = 0.0;
                    }
                    // Ensure the flow rate at lower speed has to be less or equal to the flow rate at higher speed
                    if (i != NumOfSpeedCooling) {
                        if (thisMSHeatPump.CoolVolumeFlowRate(i) > thisMSHeatPump.CoolVolumeFlowRate(i + 1)) {
                            thisMSHeatPump.CoolVolumeFlowRate(i) = thisMSHeatPump.CoolVolumeFlowRate(i + 1);
                        }
                    }
                    BaseSizer::reportSizerOutput(state,
                                                 state.dataHVACMultiSpdHP->CurrentModuleObject,
                                                 thisMSHeatPump.Name,
                                                 format("Speed {} Supply Air Flow Rate During Cooling Operation [m3/s]", i),
                                                 thisMSHeatPump.CoolVolumeFlowRate(i));
                }
            }
        }

        for (int i = NumOfSpeedHeating; i >= 1; --i) {
            if (thisMSHeatPump.HeatVolumeFlowRate(i) == AutoSize) {
                if (state.dataSize->CurSysNum > 0) {
                    if (i == NumOfSpeedHeating) {
                        CheckSysSizing(state, state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name);
                        thisMSHeatPump.HeatVolumeFlowRate(i) = state.dataSize->FinalSysSizing(state.dataSize->CurSysNum).DesMainVolFlow;
                        if (thisMSHeatPump.FanVolFlow < thisMSHeatPump.HeatVolumeFlowRate(i) &&
                            thisMSHeatPump.FanVolFlow != AutoSize) {
                            thisMSHeatPump.HeatVolumeFlowRate(i) = thisMSHeatPump.FanVolFlow;
                            ShowWarningError(state,
                                             format("{} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                            ShowContinueError(state,
                                              "The supply air flow rate at high speed is less than the autosized value for the maximum air flow rate "
                                              "in heating mode. Consider autosizing the fan for this simulation.");
                            ShowContinueError(state,
                                              "The maximum air flow rate at high speed in heating mode is reset to the supply air flow rate and the "
                                              "simulation continues.");
                        }
                    } else {
                        thisMSHeatPump.HeatVolumeFlowRate(i) =
                            thisMSHeatPump.HeatVolumeFlowRate(NumOfSpeedHeating) * i / NumOfSpeedHeating;
                    }
                    if (thisMSHeatPump.HeatVolumeFlowRate(i) < SmallAirVolFlow) {
                        thisMSHeatPump.HeatVolumeFlowRate(i) = 0.0;
                    }
                    // Ensure the flow rate at lower speed has to be less or equal to the flow rate at higher speed
                    if (i != NumOfSpeedHeating) {
                        if (thisMSHeatPump.HeatVolumeFlowRate(i) > thisMSHeatPump.HeatVolumeFlowRate(i + 1)) {
                            thisMSHeatPump.HeatVolumeFlowRate(i) = thisMSHeatPump.HeatVolumeFlowRate(i + 1);
                        }
                    }
                    BaseSizer::reportSizerOutput(state,
                                                 state.dataHVACMultiSpdHP->CurrentModuleObject,
                                                 thisMSHeatPump.Name,
                                                 format("Speed{}Supply Air Flow Rate During Heating Operation [m3/s]", i),
                                                 thisMSHeatPump.HeatVolumeFlowRate(i));
                }
            }
        }

        if (thisMSHeatPump.IdleVolumeAirRate == AutoSize) {
            if (state.dataSize->CurSysNum > 0) {
                CheckSysSizing(state, state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name);
                thisMSHeatPump.IdleVolumeAirRate = state.dataSize->FinalSysSizing(state.dataSize->CurSysNum).DesMainVolFlow;
                if (thisMSHeatPump.FanVolFlow < thisMSHeatPump.IdleVolumeAirRate &&
                    thisMSHeatPump.FanVolFlow != AutoSize) {
                    thisMSHeatPump.IdleVolumeAirRate = thisMSHeatPump.FanVolFlow;
                    ShowWarningError(state, format("{} \"{}\"", state.dataHVACMultiSpdHP->CurrentModuleObject, thisMSHeatPump.Name));
                    ShowContinueError(state,
                                      "The supply air flow rate is less than the autosized value for the maximum air flow rate when no heating or "
                                      "cooling is needed. Consider autosizing the fan for this simulation.");
                    ShowContinueError(state,
                                      "The maximum air flow rate when no heating or cooling is needed is reset to the supply air flow rate and the "
                                      "simulation continues.");
                }
                if (thisMSHeatPump.IdleVolumeAirRate < SmallAirVolFlow) {
                    thisMSHeatPump.IdleVolumeAirRate = 0.0;
                }

                BaseSizer::reportSizerOutput(state,
                                             state.dataHVACMultiSpdHP->CurrentModuleObject,
                                             thisMSHeatPump.Name,
                                             "Supply Air Flow Rate When No Cooling or Heating is Needed [m3/s]",
                                             thisMSHeatPump.IdleVolumeAirRate);
            }
        }

        if (thisMSHeatPump.SuppMaxAirTemp == AutoSize) {
            if (state.dataSize->CurSysNum > 0) {
                if (thisMSHeatPump.SuppHeatCoilType == 1) { // Gas
                    CheckZoneSizing(state, "Coil:Heating:Fuel", thisMSHeatPump.Name);
                } else {
                    CheckZoneSizing(state, "Coil:Heating:Electric", thisMSHeatPump.Name);
                }
                thisMSHeatPump.SuppMaxAirTemp = state.dataSize->FinalSysSizing(state.dataSize->CurSysNum).HeatSupTemp;
                BaseSizer::reportSizerOutput(state,
                                             state.dataHVACMultiSpdHP->CurrentModuleObject,
                                             thisMSHeatPump.Name,
                                             "Maximum Supply Air Temperature from Supplemental Heater [C]",
                                             thisMSHeatPump.SuppMaxAirTemp);
            }
        }

        if (thisMSHeatPump.DesignSuppHeatingCapacity == AutoSize) {
            if (state.dataSize->CurSysNum > 0) {
                if (thisMSHeatPump.SuppHeatCoilType == 1) { // Gas
                    CheckSysSizing(state, "Coil:Heating:Fuel", thisMSHeatPump.Name);
                } else {
                    CheckSysSizing(state, "Coil:Heating:Electric", thisMSHeatPump.Name);
                }
                thisMSHeatPump.DesignSuppHeatingCapacity = state.dataSize->FinalSysSizing(state.dataSize->CurSysNum).HeatCap;
            } else {
                thisMSHeatPump.DesignSuppHeatingCapacity = 0.0;
            }
            BaseSizer::reportSizerOutput(state,
                                         state.dataHVACMultiSpdHP->CurrentModuleObject,
                                         thisMSHeatPump.Name,
                                         "Supplemental Heating Coil Nominal Capacity [W]",
                                         thisMSHeatPump.DesignSuppHeatingCapacity);
        }
        state.dataSize->SuppHeatCap = thisMSHeatPump.DesignSuppHeatingCapacity;

        if (thisMSHeatPump.HeatRecActive) {
            RegisterPlantCompDesignFlow(state, thisMSHeatPump.HeatRecInletNodeNum, thisMSHeatPump.DesignHeatRecFlowRate);
        }
    }

    //******************************************************************************

    void ControlMSHPOutputEMS(EnergyPlusData &state,
                              int const MSHeatPumpNum,                                 // Unit index of engine driven heat pump
                              bool const FirstHVACIteration,                           // flag for 1st HVAC iteration in the time step
                              DataHVACGlobals::CompressorOperation const CompressorOp, // compressor operation; 1=on, 0=off
                              int const OpMode,                                        // operating mode: CycFanCycCoil | ContFanCycCoil
                              Real64 const QZnReq,                                     // cooling or heating output needed by zone [W]
                              Real64 const SpeedVal,                                   // continuous speed value
                              int &SpeedNum,                                           // discrete speed level
                              Real64 &SpeedRatio,                                      // unit speed ratio for DX coils
                              Real64 &PartLoadFrac,                                    // unit part load fraction
                              Real64 &OnOffAirFlowRatio, // ratio of compressor ON airflow to AVERAGE airflow over timestep
                              Real64 &SupHeaterLoad      // Supplemental heater load [W]

    )
    {
        OnOffAirFlowRatio = 0.0;
        SupHeaterLoad = 0.0;

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        // Get EMS output
        SpeedNum = ceil(SpeedVal);
        bool useMaxedSpeed = false;
        std::string useMaxedSpeedCoilName;
        if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::HeatingMode) {
            if (SpeedNum > thisMSHeatPump.NumOfSpeedHeating) {
                SpeedNum = thisMSHeatPump.NumOfSpeedHeating;
                useMaxedSpeed = true;
                useMaxedSpeedCoilName = thisMSHeatPump.DXHeatCoilName;
            }
        } else if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::CoolingMode) {
            if (SpeedNum > thisMSHeatPump.NumOfSpeedCooling) {
                SpeedNum = thisMSHeatPump.NumOfSpeedCooling;
                useMaxedSpeed = true;
                useMaxedSpeedCoilName = thisMSHeatPump.DXCoolCoilName;
            }
        }
        if (useMaxedSpeed) {
            thisMSHeatPump.CoilSpeedErrIndex++;
            ShowRecurringWarningErrorAtEnd(state,
                                           "Wrong coil speed EMS override value, for unit=\"" + useMaxedSpeedCoilName +
                                               "\". Exceeding maximum coil speed level. Speed level is set to the maximum coil speed level allowed.",
                                           thisMSHeatPump.CoilSpeedErrIndex,
                                           SpeedVal,
                                           SpeedVal,
                                           _,
                                           "",
                                           "");
        }
        // Calculate TempOutput
        Real64 TempOutput = 0.0; // unit output when iteration limit exceeded [W]

        if (SpeedNum == 1) {
            SpeedRatio = 0.0;
            if (useMaxedSpeed || floor(SpeedVal) == SpeedVal) {
                PartLoadFrac = 1;
            } else {
                PartLoadFrac = SpeedVal - floor(SpeedVal);
            }
            CalcMSHeatPump(state,
                           MSHeatPumpNum,
                           FirstHVACIteration,
                           CompressorOp,
                           SpeedNum,
                           SpeedRatio,
                           PartLoadFrac,
                           TempOutput,
                           QZnReq,
                           OnOffAirFlowRatio,
                           SupHeaterLoad);
        } else {
            PartLoadFrac = 0.0;
            if (useMaxedSpeed || floor(SpeedVal) == SpeedVal) {
                SpeedRatio = 1;
            } else {
                SpeedRatio = SpeedVal - floor(SpeedVal);
            }
            CalcMSHeatPump(state,
                           MSHeatPumpNum,
                           FirstHVACIteration,
                           CompressorOp,
                           SpeedNum,
                           SpeedRatio,
                           PartLoadFrac,
                           TempOutput,
                           QZnReq,
                           OnOffAirFlowRatio,
                           SupHeaterLoad);
        }

        ControlMSHPSupHeater(state,
                             MSHeatPumpNum,
                             FirstHVACIteration,
                             CompressorOp,
                             OpMode,
                             QZnReq,
                             TempOutput,
                             SpeedNum,
                             SpeedRatio,
                             PartLoadFrac,
                             OnOffAirFlowRatio,
                             SupHeaterLoad);
        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).CycRatio = PartLoadFrac;
        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedRatio = SpeedRatio;
        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedNum = SpeedNum;
    }

    void ControlMSHPSupHeater(EnergyPlusData &state,
                              int const MSHeatPumpNum,                                 // Unit index of engine driven heat pump
                              bool const FirstHVACIteration,                           // flag for 1st HVAC iteration in the time step
                              DataHVACGlobals::CompressorOperation const CompressorOp, // compressor operation; 1=on, 0=off
                              int const OpMode,                                        // operating mode: CycFanCycCoil | ContFanCycCoil
                              Real64 const QZnReq,                                     // cooling or heating output needed by zone [W]
                              int const EMSOutput,                                     // unit full output when compressor is operating [W]vvvv
                              int const SpeedNum,                                      // Speed number
                              Real64 SpeedRatio,                                       // unit speed ratio for DX coils
                              Real64 PartLoadFrac,                                     // unit part load fraction
                              Real64 OnOffAirFlowRatio, // ratio of compressor ON airflow to AVERAGE airflow over timestep
                              Real64 &SupHeaterLoad     // Supplemental heater load [W]

    )
    {
        // if the DX heating coil cannot meet the load, trim with supplemental heater
        // occurs with constant fan mode when compressor is on or off
        // occurs with cycling fan mode when compressor PLR is equal to 1
        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        if ((QZnReq > SmallLoad && QZnReq > EMSOutput)) {
            Real64 TempOutput;
            if (state.dataEnvrn->OutDryBulbTemp <= thisMSHeatPump.SuppMaxAirTemp) {
                SupHeaterLoad = QZnReq - EMSOutput;
            } else {
                SupHeaterLoad = 0.0;
            }
            CalcMSHeatPump(state,
                           MSHeatPumpNum,
                           FirstHVACIteration,
                           CompressorOp,
                           SpeedNum,
                           SpeedRatio,
                           PartLoadFrac,
                           TempOutput,
                           QZnReq,
                           OnOffAirFlowRatio,
                           SupHeaterLoad);
        }

        // check the outlet of the supplemental heater to be lower than the maximum supplemental heater supply air temperature
        if (state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).Temp > thisMSHeatPump.SuppMaxAirTemp &&
            SupHeaterLoad > 0.0) {

            //   If the supply air temperature is to high, turn off the supplemental heater to recalculate the outlet temperature
            SupHeaterLoad = 0.0;
            Real64 QCoilActual; // coil load actually delivered returned to calling component
            CalcNonDXHeatingCoils(state, MSHeatPumpNum, FirstHVACIteration, SupHeaterLoad, OpMode, QCoilActual);

            //   If the outlet temperature is below the maximum supplemental heater supply air temperature, reduce the load passed to
            //   the supplemental heater, otherwise leave the supplemental heater off. If the supplemental heater is to be turned on,
            //   use the outlet conditions when the supplemental heater was off (CALL above) as the inlet conditions for the calculation
            //   of supplemental heater load to just meet the maximum supply air temperature from the supplemental heater.
            if (state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).Temp < thisMSHeatPump.SuppMaxAirTemp) {
                Real64 CpAir = Psychrometrics::PsyCpAirFnW(state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).HumRat);
                SupHeaterLoad =
                    state.dataLoopNodes->Node(thisMSHeatPump.AirInletNodeNum).MassFlowRate * CpAir *
                    (thisMSHeatPump.SuppMaxAirTemp - state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).Temp);

            } else {
                SupHeaterLoad = 0.0;
            }
        }
    }

    void ControlMSHPOutput(EnergyPlusData &state,
                           int const MSHeatPumpNum,                                 // Unit index of engine driven heat pump
                           bool const FirstHVACIteration,                           // flag for 1st HVAC iteration in the time step
                           DataHVACGlobals::CompressorOperation const CompressorOp, // compressor operation; 1=on, 0=off
                           int const OpMode,                                        // operating mode: CycFanCycCoil | ContFanCycCoil
                           Real64 const QZnReq,                                     // cooling or heating output needed by zone [W]
                           int const ZoneNum [[maybe_unused]],                      // Index to zone number
                           int &SpeedNum,                                           // Speed number
                           Real64 &SpeedRatio,                                      // unit speed ratio for DX coils
                           Real64 &PartLoadFrac,                                    // unit part load fraction
                           Real64 &OnOffAirFlowRatio,                               // ratio of compressor ON airflow to AVERAGE airflow over timestep
                           Real64 &SupHeaterLoad                                    // Supplemental heater load [W]
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Lixing Gu
        //       DATE WRITTEN   June 2007
        //       MODIFIED       na
        //       RE-ENGINEERED  Revised for multispeed heat pump use based on ControlPTHPOutput

        // PURPOSE OF THIS SUBROUTINE:
        // Determine the part load fraction at low speed, and speed ratio at high speed for this time step.

        // METHODOLOGY EMPLOYED:
        // Use RegulaFalsi technique to iterate on part-load ratio until convergence is achieved.

        using General::SolveRoot;
        using Psychrometrics::PsyCpAirFnW;

        // SUBROUTINE PARAMETER DEFINITIONS:
        int constexpr MaxIte(500); // maximum number of iterations

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 FullOutput;         // unit full output when compressor is operating [W]
        Real64 LowOutput;          // unit full output at low speed [W]
        Real64 TempOutput;         // unit output when iteration limit exceeded [W]
        Real64 NoCompOutput;       // output when no active compressor [W]
        Real64 ErrorToler;         // error tolerance
        int SolFla;                // Flag of RegulaFalsi solver
        Real64 CpAir;              // air specific heat
        Real64 OutsideDryBulbTemp; // Outside air temperature at external node height
        Real64 QCoilActual;        // coil load actually delivered returned to calling component

        SupHeaterLoad = 0.0;
        PartLoadFrac = 0.0;
        SpeedRatio = 0.0;
        SpeedNum = 1;

        OutsideDryBulbTemp = state.dataEnvrn->OutDryBulbTemp;

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        //!!LKL Discrepancy with < 0
        if (GetCurrentScheduleValue(state, thisMSHeatPump.AvaiSchedPtr) == 0.0) return;

        // Get result when DX coil is off
        CalcMSHeatPump(state,
                       MSHeatPumpNum,
                       FirstHVACIteration,
                       CompressorOp,
                       SpeedNum,
                       SpeedRatio,
                       PartLoadFrac,
                       NoCompOutput,
                       QZnReq,
                       OnOffAirFlowRatio,
                       SupHeaterLoad);

        // If cooling and NoCompOutput < QZnReq, the coil needs to be off
        // If heating and NoCompOutput > QZnReq, the coil needs to be off
        if ((QZnReq < (-1.0 * SmallLoad) && NoCompOutput < QZnReq) || (QZnReq > SmallLoad && NoCompOutput > QZnReq) ||
            std::abs(QZnReq) <= SmallLoad) {
            return;
        }

        // Get full load result
        PartLoadFrac = 1.0;
        SpeedRatio = 1.0;
        if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::HeatingMode) {
            SpeedNum = thisMSHeatPump.NumOfSpeedHeating;
            if (thisMSHeatPump.Staged && std::abs(thisMSHeatPump.StageNum) < SpeedNum) {
                SpeedNum = std::abs(thisMSHeatPump.StageNum);
                if (SpeedNum == 1) SpeedRatio = 0.0;
            }
        }
        if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::CoolingMode) {
            SpeedNum = thisMSHeatPump.NumOfSpeedCooling;
            if (thisMSHeatPump.Staged && std::abs(thisMSHeatPump.StageNum) < SpeedNum) {
                SpeedNum = std::abs(thisMSHeatPump.StageNum);
                if (SpeedNum == 1) SpeedRatio = 0.0;
            }
        }

        CalcMSHeatPump(state,
                       MSHeatPumpNum,
                       FirstHVACIteration,
                       CompressorOp,
                       SpeedNum,
                       SpeedRatio,
                       PartLoadFrac,
                       FullOutput,
                       QZnReq,
                       OnOffAirFlowRatio,
                       SupHeaterLoad);

        if (QZnReq < (-1.0 * SmallLoad)) {
            // Since we are cooling, we expect FullOutput to be < 0 and FullOutput < NoCompOutput
            // Check that this is the case; if not set PartLoadFrac = 0.0 (off) and return
            if (FullOutput >= 0.0 || FullOutput >= NoCompOutput) {
                PartLoadFrac = 0.0;
                SpeedRatio = 0.0;
                SpeedNum = 0;
                return;
            }
            //  ! If the QZnReq <= FullOutput the unit needs to run full out
            if (QZnReq <= FullOutput) {
                PartLoadFrac = 1.0;
                SpeedRatio = 1.0;
                if (thisMSHeatPump.Staged && SpeedNum == 1) SpeedRatio = 0.0;
                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).CycRatio = PartLoadFrac;
                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedRatio = SpeedRatio;
                state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedNum = SpeedNum;
                return;
            }
            ErrorToler = 0.001; // Error tolerance for convergence from input deck
        } else {
            // Since we are heating, we expect FullOutput to be > 0 and FullOutput > NoCompOutput
            // Check that this is the case; if not set PartLoadFrac = 0.0 (off)
            if (FullOutput <= 0.0 || FullOutput <= NoCompOutput) {
                PartLoadFrac = 0.0;
                SpeedRatio = 0.0;
                // may need supplemental heating so don't return in heating mode
            }
            if (QZnReq >= FullOutput) {
                PartLoadFrac = 1.0;
                SpeedRatio = 1.0;
                // may need supplemental heating so don't return in heating mode
            }
            ErrorToler = 0.001; // Error tolerance for convergence from input deck
        }

        // Direct solution
        if (state.dataGlobal->DoCoilDirectSolutions && !thisMSHeatPump.Staged) {
            Real64 TempOutput0 = 0.0;
            thisMSHeatPump.FullOutput = 0.0;

            // heating
            if (QZnReq > SmallLoad && QZnReq < FullOutput) {
                CalcMSHeatPump(
                    state, MSHeatPumpNum, FirstHVACIteration, CompressorOp, 1, 0.0, 0.0, TempOutput0, QZnReq, OnOffAirFlowRatio, SupHeaterLoad);

                for (int i = 1; i <= thisMSHeatPump.NumOfSpeedHeating; ++i) {
                    if (i == 1) {
                        CalcMSHeatPump(state,
                                       MSHeatPumpNum,
                                       FirstHVACIteration,
                                       CompressorOp,
                                       i,
                                       0.0,
                                       1.0,
                                       thisMSHeatPump.FullOutput(i),
                                       QZnReq,
                                       OnOffAirFlowRatio,
                                       SupHeaterLoad);
                        if (QZnReq <= thisMSHeatPump.FullOutput(i)) {
                            SpeedNum = i;
                            PartLoadFrac = (QZnReq - TempOutput0) / (thisMSHeatPump.FullOutput(i) - TempOutput0);
                            CalcMSHeatPump(state,
                                           MSHeatPumpNum,
                                           FirstHVACIteration,
                                           CompressorOp,
                                           i,
                                           0.0,
                                           PartLoadFrac,
                                           TempOutput,
                                           QZnReq,
                                           OnOffAirFlowRatio,
                                           SupHeaterLoad);
                            break;
                        }
                    } else {
                        CalcMSHeatPump(state,
                                       MSHeatPumpNum,
                                       FirstHVACIteration,
                                       CompressorOp,
                                       i,
                                       1.0,
                                       1.0,
                                       thisMSHeatPump.FullOutput(i),
                                       QZnReq,
                                       OnOffAirFlowRatio,
                                       SupHeaterLoad);
                        if (QZnReq <= thisMSHeatPump.FullOutput(i)) {
                            SpeedNum = i;
                            PartLoadFrac = 1.0;
                            SpeedRatio = (QZnReq - thisMSHeatPump.FullOutput(i - 1)) /
                                         (thisMSHeatPump.FullOutput(i) - thisMSHeatPump.FullOutput(i - 1));
                            CalcMSHeatPump(state,
                                           MSHeatPumpNum,
                                           FirstHVACIteration,
                                           CompressorOp,
                                           i,
                                           SpeedRatio,
                                           1.0,
                                           TempOutput,
                                           QZnReq,
                                           OnOffAirFlowRatio,
                                           SupHeaterLoad);
                            break;
                        }
                    }
                }
            }

            // Coolling
            if (QZnReq < (-1.0 * SmallLoad) && QZnReq > FullOutput) {
                CalcMSHeatPump(
                    state, MSHeatPumpNum, FirstHVACIteration, CompressorOp, 1, 0.0, 0.0, TempOutput0, QZnReq, OnOffAirFlowRatio, SupHeaterLoad);
                for (int i = 1; i <= thisMSHeatPump.NumOfSpeedCooling; ++i) {
                    if (i == 1) {
                        CalcMSHeatPump(state,
                                       MSHeatPumpNum,
                                       FirstHVACIteration,
                                       CompressorOp,
                                       i,
                                       0.0,
                                       1.0,
                                       thisMSHeatPump.FullOutput(i),
                                       QZnReq,
                                       OnOffAirFlowRatio,
                                       SupHeaterLoad);
                        if (QZnReq >= thisMSHeatPump.FullOutput(i)) {
                            SpeedNum = i;
                            PartLoadFrac = (QZnReq - TempOutput0) / (thisMSHeatPump.FullOutput(i) - TempOutput0);
                            CalcMSHeatPump(state,
                                           MSHeatPumpNum,
                                           FirstHVACIteration,
                                           CompressorOp,
                                           i,
                                           0.0,
                                           PartLoadFrac,
                                           TempOutput,
                                           QZnReq,
                                           OnOffAirFlowRatio,
                                           SupHeaterLoad);
                            break;
                        }
                    } else {
                        CalcMSHeatPump(state,
                                       MSHeatPumpNum,
                                       FirstHVACIteration,
                                       CompressorOp,
                                       i,
                                       1.0,
                                       1.0,
                                       thisMSHeatPump.FullOutput(i),
                                       QZnReq,
                                       OnOffAirFlowRatio,
                                       SupHeaterLoad);
                        if (QZnReq >= thisMSHeatPump.FullOutput(i)) {
                            SpeedNum = i;
                            PartLoadFrac = 1.0;
                            SpeedRatio = (QZnReq - thisMSHeatPump.FullOutput(i - 1)) /
                                         (thisMSHeatPump.FullOutput(i) - thisMSHeatPump.FullOutput(i - 1));
                            CalcMSHeatPump(state,
                                           MSHeatPumpNum,
                                           FirstHVACIteration,
                                           CompressorOp,
                                           i,
                                           SpeedRatio,
                                           1.0,
                                           TempOutput,
                                           QZnReq,
                                           OnOffAirFlowRatio,
                                           SupHeaterLoad);
                            break;
                        }
                    }
                }
            }
        } else {
            // Calculate the part load fraction
            if (((QZnReq > SmallLoad && QZnReq < FullOutput) || (QZnReq < (-1.0 * SmallLoad) && QZnReq > FullOutput)) &&
                (!thisMSHeatPump.Staged)) {
                // Check whether the low speed coil can meet the load or not
                CalcMSHeatPump(
                    state, MSHeatPumpNum, FirstHVACIteration, CompressorOp, 1, 0.0, 1.0, LowOutput, QZnReq, OnOffAirFlowRatio, SupHeaterLoad);
                if ((QZnReq > 0.0 && QZnReq <= LowOutput) || (QZnReq < 0.0 && QZnReq >= LowOutput)) {
                    SpeedRatio = 0.0;
                    SpeedNum = 1;
                    auto f = [&state, MSHeatPumpNum, FirstHVACIteration, QZnReq, OnOffAirFlowRatio, SupHeaterLoad, CompressorOp](
                                 Real64 const PartLoadFrac) {
                        //  Calculates residual function ((ActualOutput - QZnReq)/QZnReq); MSHP output depends on PLR which is being varied to zero
                        //  the residual.
                        Real64 ActualOutput; // delivered capacity of MSHP
                        Real64 tmpAirFlowRatio = OnOffAirFlowRatio;
                        Real64 tmpHeaterLoad = SupHeaterLoad;
                        CalcMSHeatPump(state,
                                       MSHeatPumpNum,
                                       FirstHVACIteration,
                                       CompressorOp,
                                       1,
                                       0.0,
                                       PartLoadFrac,
                                       ActualOutput,
                                       QZnReq,
                                       tmpAirFlowRatio,
                                       tmpHeaterLoad);
                        return (ActualOutput - QZnReq) / QZnReq;
                    };
                    SolveRoot(state, ErrorToler, MaxIte, SolFla, PartLoadFrac, f, 0.0, 1.0);
                    if (SolFla == -1) {
                        if (!state.dataGlobal->WarmupFlag) {
                            if (state.dataHVACMultiSpdHP->ErrCountCyc == 0) {
                                ++state.dataHVACMultiSpdHP->ErrCountCyc; // TODO: Why is the error count shared among all heat pump units?
                                ShowWarningError(state,
                                                 format("Iteration limit exceeded calculating DX unit cycling ratio, for unit={}",
                                                        thisMSHeatPump.Name));
                                ShowContinueErrorTimeStamp(state, format("Cycling ratio returned={:.2R}", PartLoadFrac));
                            } else {
                                ++state.dataHVACMultiSpdHP->ErrCountCyc;
                                ShowRecurringWarningErrorAtEnd(
                                    state,
                                    thisMSHeatPump.Name +
                                        "\": Iteration limit warning exceeding calculating DX unit cycling ratio  continues...",
                                    thisMSHeatPump.ErrIndexCyc,
                                    PartLoadFrac,
                                    PartLoadFrac);
                            }
                        }
                    } else if (SolFla == -2) {
                        ShowFatalError(state,
                                       format("DX unit cycling ratio calculation failed: cycling limits exceeded, for unit={}",
                                              thisMSHeatPump.DXCoolCoilName));
                    }
                } else {
                    // Check to see which speed to meet the load
                    PartLoadFrac = 1.0;
                    SpeedRatio = 1.0;
                    if (QZnReq < (-1.0 * SmallLoad)) { // Cooling
                        for (int i = 2; i <= thisMSHeatPump.NumOfSpeedCooling; ++i) {
                            CalcMSHeatPump(state,
                                           MSHeatPumpNum,
                                           FirstHVACIteration,
                                           CompressorOp,
                                           i,
                                           SpeedRatio,
                                           PartLoadFrac,
                                           TempOutput,
                                           QZnReq,
                                           OnOffAirFlowRatio,
                                           SupHeaterLoad);
                            if (QZnReq >= TempOutput) {
                                SpeedNum = i;
                                break;
                            }
                        }
                    } else {
                        for (int i = 2; i <= thisMSHeatPump.NumOfSpeedHeating; ++i) {
                            CalcMSHeatPump(state,
                                           MSHeatPumpNum,
                                           FirstHVACIteration,
                                           CompressorOp,
                                           i,
                                           SpeedRatio,
                                           PartLoadFrac,
                                           TempOutput,
                                           QZnReq,
                                           OnOffAirFlowRatio,
                                           SupHeaterLoad);
                            if (QZnReq <= TempOutput) {
                                SpeedNum = i;
                                break;
                            }
                        }
                    }
                    auto f = [&state, OnOffAirFlowRatio, SupHeaterLoad, MSHeatPumpNum, FirstHVACIteration, CompressorOp, SpeedNum, QZnReq](
                                 Real64 const SpeedRatio) {
                        //  Calculates residual function ((ActualOutput - QZnReq)/QZnReq) MSHP output depends on PLR which is being varied to zero the
                        //  residual.
                        Real64 localAirFlowRatio = OnOffAirFlowRatio;
                        Real64 localHeaterLoad = SupHeaterLoad;
                        Real64 ActualOutput;
                        CalcMSHeatPump(state,
                                       MSHeatPumpNum,
                                       FirstHVACIteration,
                                       CompressorOp,
                                       SpeedNum,
                                       SpeedRatio,
                                       1.0,
                                       ActualOutput,
                                       QZnReq,
                                       localAirFlowRatio,
                                       localHeaterLoad);
                        return (ActualOutput - QZnReq) / QZnReq;
                    };
                    SolveRoot(state, ErrorToler, MaxIte, SolFla, SpeedRatio, f, 0.0, 1.0);
                    if (SolFla == -1) {
                        if (!state.dataGlobal->WarmupFlag) {
                            if (state.dataHVACMultiSpdHP->ErrCountVar == 0) {
                                ++state.dataHVACMultiSpdHP->ErrCountVar;
                                ShowWarningError(
                                    state,
                                    format("Iteration limit exceeded calculating DX unit speed ratio, for unit={}", thisMSHeatPump.Name));
                                ShowContinueErrorTimeStamp(state, format("Speed ratio returned=[{:.2R}], Speed number ={}", SpeedRatio, SpeedNum));
                            } else {
                                ++state.dataHVACMultiSpdHP->ErrCountVar;
                                ShowRecurringWarningErrorAtEnd(
                                    state,
                                    thisMSHeatPump.Name +
                                        "\": Iteration limit warning exceeding calculating DX unit speed ratio continues...",
                                    thisMSHeatPump.ErrIndexVar,
                                    SpeedRatio,
                                    SpeedRatio);
                            }
                        }
                    } else if (SolFla == -2) {
                        ShowFatalError(state,
                                       format("DX unit compressor speed calculation failed: speed limits exceeded, for unit={}",
                                              thisMSHeatPump.DXCoolCoilName));
                    }
                }
            } else {
                // Staged thermostat performance
                if (thisMSHeatPump.StageNum != 0) {
                    SpeedNum = std::abs(thisMSHeatPump.StageNum);
                    if (SpeedNum == 1) {
                        CalcMSHeatPump(
                            state, MSHeatPumpNum, FirstHVACIteration, CompressorOp, 1, 0.0, 1.0, LowOutput, QZnReq, OnOffAirFlowRatio, SupHeaterLoad);
                        SpeedRatio = 0.0;
                        if ((QZnReq > 0.0 && QZnReq <= LowOutput) || (QZnReq < 0.0 && QZnReq >= LowOutput)) {
                            auto f = [&state, MSHeatPumpNum, FirstHVACIteration, QZnReq, OnOffAirFlowRatio, SupHeaterLoad, CompressorOp](
                                         Real64 const PartLoadFrac) {
                                //  Calculates residual function ((ActualOutput - QZnReq)/QZnReq); MSHP output depends on PLR which is being varied to
                                //  zero the residual.
                                Real64 ActualOutput; // delivered capacity of MSHP
                                Real64 tmpAirFlowRatio = OnOffAirFlowRatio;
                                Real64 tmpHeaterLoad = SupHeaterLoad;
                                CalcMSHeatPump(state,
                                               MSHeatPumpNum,
                                               FirstHVACIteration,
                                               CompressorOp,
                                               1,
                                               0.0,
                                               PartLoadFrac,
                                               ActualOutput,
                                               QZnReq,
                                               tmpAirFlowRatio,
                                               tmpHeaterLoad);
                                return (ActualOutput - QZnReq) / QZnReq;
                            };
                            SolveRoot(state, ErrorToler, MaxIte, SolFla, PartLoadFrac, f, 0.0, 1.0);
                            if (SolFla == -1) {
                                if (!state.dataGlobal->WarmupFlag) {
                                    if (state.dataHVACMultiSpdHP->ErrCountCyc == 0) {
                                        ++state.dataHVACMultiSpdHP->ErrCountCyc;
                                        ShowWarningError(state,
                                                         format("Iteration limit exceeded calculating DX unit cycling ratio, for unit={}",
                                                                thisMSHeatPump.Name));
                                        ShowContinueErrorTimeStamp(state, format("Cycling ratio returned={:.2R}", PartLoadFrac));
                                    } else {
                                        ++state.dataHVACMultiSpdHP->ErrCountCyc;
                                        ShowRecurringWarningErrorAtEnd(
                                            state,
                                            thisMSHeatPump.Name +
                                                "\": Iteration limit warning exceeding calculating DX unit cycling ratio  continues...",
                                            thisMSHeatPump.ErrIndexCyc,
                                            PartLoadFrac,
                                            PartLoadFrac);
                                    }
                                }
                            } else if (SolFla == -2) {
                                ShowFatalError(state,
                                               format("DX unit cycling ratio calculation failed: cycling limits exceeded, for unit={}",
                                                      thisMSHeatPump.DXCoolCoilName));
                            }
                        } else {
                            FullOutput = LowOutput;
                            PartLoadFrac = 1.0;
                        }
                    } else {
                        if (thisMSHeatPump.StageNum < 0) {
                            SpeedNum = min(thisMSHeatPump.NumOfSpeedCooling, std::abs(thisMSHeatPump.StageNum));
                        } else {
                            SpeedNum = min(thisMSHeatPump.NumOfSpeedHeating, std::abs(thisMSHeatPump.StageNum));
                        }
                        CalcMSHeatPump(state,
                                       MSHeatPumpNum,
                                       FirstHVACIteration,
                                       CompressorOp,
                                       SpeedNum,
                                       0.0,
                                       1.0,
                                       LowOutput,
                                       QZnReq,
                                       OnOffAirFlowRatio,
                                       SupHeaterLoad);
                        if ((QZnReq > 0.0 && QZnReq >= LowOutput) || (QZnReq < 0.0 && QZnReq <= LowOutput)) {
                            CalcMSHeatPump(state,
                                           MSHeatPumpNum,
                                           FirstHVACIteration,
                                           CompressorOp,
                                           SpeedNum,
                                           1.0,
                                           1.0,
                                           FullOutput,
                                           QZnReq,
                                           OnOffAirFlowRatio,
                                           SupHeaterLoad);
                            if ((QZnReq > 0.0 && QZnReq <= FullOutput) || (QZnReq < 0.0 && QZnReq >= FullOutput)) {
                                auto f = // (THIS_AUTO_OK)
                                    [&state, OnOffAirFlowRatio, SupHeaterLoad, MSHeatPumpNum, FirstHVACIteration, CompressorOp, SpeedNum, QZnReq](
                                        Real64 const SpeedRatio) {
                                        //  Calculates residual function ((ActualOutput - QZnReq)/QZnReq) MSHP output depends on PLR which is being
                                        //  varied to zero the residual.
                                        Real64 localAirFlowRatio = OnOffAirFlowRatio;
                                        Real64 localHeaterLoad = SupHeaterLoad;
                                        Real64 ActualOutput;
                                        CalcMSHeatPump(state,
                                                       MSHeatPumpNum,
                                                       FirstHVACIteration,
                                                       CompressorOp,
                                                       SpeedNum,
                                                       SpeedRatio,
                                                       1.0,
                                                       ActualOutput,
                                                       QZnReq,
                                                       localAirFlowRatio,
                                                       localHeaterLoad);
                                        return (ActualOutput - QZnReq) / QZnReq;
                                    };
                                SolveRoot(state, ErrorToler, MaxIte, SolFla, SpeedRatio, f, 0.0, 1.0);
                                if (SolFla == -1) {
                                    if (!state.dataGlobal->WarmupFlag) {
                                        if (state.dataHVACMultiSpdHP->ErrCountVar == 0) {
                                            ++state.dataHVACMultiSpdHP->ErrCountVar;
                                            ShowWarningError(state,
                                                             format("Iteration limit exceeded calculating DX unit speed ratio, for unit={}",
                                                                    thisMSHeatPump.Name));
                                            ShowContinueErrorTimeStamp(
                                                state, format("Speed ratio returned=[{:.2R}], Speed number ={}", SpeedRatio, SpeedNum));
                                        } else {
                                            ++state.dataHVACMultiSpdHP->ErrCountVar;
                                            ShowRecurringWarningErrorAtEnd(
                                                state,
                                                thisMSHeatPump.Name +
                                                    "\": Iteration limit warning exceeding calculating DX unit speed ratio continues...",
                                                thisMSHeatPump.ErrIndexVar,
                                                SpeedRatio,
                                                SpeedRatio);
                                        }
                                    }
                                } else if (SolFla == -2) {
                                    ShowFatalError(state,
                                                   format("DX unit compressor speed calculation failed: speed limits exceeded, for unit={}",
                                                          thisMSHeatPump.DXCoolCoilName));
                                }
                            } else {
                                SpeedRatio = 1.0;
                            }
                        } else { // lowOutput provides a larger capacity than needed
                            SpeedRatio = 0.0;
                        }
                    }
                }
            }
        }

        // if the DX heating coil cannot meet the load, trim with supplemental heater
        // occurs with constant fan mode when compressor is on or off
        // occurs with cycling fan mode when compressor PLR is equal to 1
        if ((QZnReq > SmallLoad && QZnReq > FullOutput)) {
            PartLoadFrac = 1.0;
            SpeedRatio = 1.0;
            if (thisMSHeatPump.Staged && SpeedNum == 1) SpeedRatio = 0.0;
            if (OutsideDryBulbTemp <= thisMSHeatPump.SuppMaxAirTemp) {
                SupHeaterLoad = QZnReq - FullOutput;
            } else {
                SupHeaterLoad = 0.0;
            }
            CalcMSHeatPump(state,
                           MSHeatPumpNum,
                           FirstHVACIteration,
                           CompressorOp,
                           SpeedNum,
                           SpeedRatio,
                           PartLoadFrac,
                           TempOutput,
                           QZnReq,
                           OnOffAirFlowRatio,
                           SupHeaterLoad);
        }

        // check the outlet of the supplemental heater to be lower than the maximum supplemental heater supply air temperature
        if (state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).Temp > thisMSHeatPump.SuppMaxAirTemp &&
            SupHeaterLoad > 0.0) {

            //   If the supply air temperature is to high, turn off the supplemental heater to recalculate the outlet temperature
            SupHeaterLoad = 0.0;
            CalcNonDXHeatingCoils(state, MSHeatPumpNum, FirstHVACIteration, SupHeaterLoad, OpMode, QCoilActual);

            //   If the outlet temperature is below the maximum supplemental heater supply air temperature, reduce the load passed to
            //   the supplemental heater, otherwise leave the supplemental heater off. If the supplemental heater is to be turned on,
            //   use the outlet conditions when the supplemental heater was off (CALL above) as the inlet conditions for the calculation
            //   of supplemental heater load to just meet the maximum supply air temperature from the supplemental heater.
            if (state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).Temp < thisMSHeatPump.SuppMaxAirTemp) {
                CpAir = PsyCpAirFnW(state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).HumRat);
                SupHeaterLoad =
                    state.dataLoopNodes->Node(thisMSHeatPump.AirInletNodeNum).MassFlowRate * CpAir *
                    (thisMSHeatPump.SuppMaxAirTemp - state.dataLoopNodes->Node(thisMSHeatPump.AirOutletNodeNum).Temp);

            } else {
                SupHeaterLoad = 0.0;
            }
        }

        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).CycRatio = PartLoadFrac;
        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedRatio = SpeedRatio;
        state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).SpeedNum = SpeedNum;
    }

    //******************************************************************************

    void CalcMSHeatPump(EnergyPlusData &state,
                        int const MSHeatPumpNum,                                 // Engine driven heat pump number
                        bool const FirstHVACIteration,                           // Flag for 1st HVAC iteration
                        DataHVACGlobals::CompressorOperation const CompressorOp, // Compressor on/off; 1=on, 0=off
                        int const SpeedNum,                                      // Speed number
                        Real64 const SpeedRatio,                                 // Compressor speed ratio
                        Real64 const PartLoadFrac,                               // Compressor part load fraction
                        Real64 &LoadMet,                                         // Load met by unit (W)
                        Real64 const QZnReq,                                     // Zone load (W)
                        Real64 &OnOffAirFlowRatio,                               // Ratio of compressor ON airflow to AVERAGE airflow over timestep
                        Real64 &SupHeaterLoad                                    // supplemental heater load (W)
    )
    {
        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Lixing Gu, FSEC
        //       DATE WRITTEN:    June 2007
        //       MODIFIED         na
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        //  This routine will calcultes MSHP performance based on given system load

        // METHODOLOGY EMPLOYED:
        // na

        // REFERENCES: na

        // Using/Aliasing
        using DXCoils::SimDXCoilMultiSpeed;
        using Fans::SimulateFanComponents;
        using HeatingCoils::SimulateHeatingCoilComponents;

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int OutletNode;            // MSHP air outlet node
        int InletNode;             // MSHP air inlet node
        Real64 OutsideDryBulbTemp; // Outdoor dry bulb temperature [C]
        Real64 AirMassFlow;        // Air mass flow rate [kg/s]
        int FanInletNode;          // MSHP air outlet node
        int FanOutletNode;         // MSHP air inlet node
        Real64 SavePartloadRatio;
        Real64 SaveSpeedRatio;
        Real64 QCoilActual;  // coil load actually delivered returned to calling component
        Real64 MinWaterFlow; // minimum water flow rate
        Real64 ErrorToler;   // supplemental heating coil convergence tollerance

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        OutletNode = thisMSHeatPump.AirOutletNodeNum;
        InletNode = thisMSHeatPump.AirInletNodeNum;
        if (thisMSHeatPump.DXHeatCoilIndex > 0) {
            if (state.dataDXCoils->DXCoil(thisMSHeatPump.DXHeatCoilIndex).IsSecondaryDXCoilInZone) {
                OutsideDryBulbTemp = state.dataZoneTempPredictorCorrector
                                         ->zoneHeatBalance(state.dataDXCoils->DXCoil(thisMSHeatPump.DXHeatCoilIndex).SecZonePtr)
                                         .ZT;
            } else {
                OutsideDryBulbTemp = state.dataEnvrn->OutDryBulbTemp;
            }
        } else if (thisMSHeatPump.DXCoolCoilIndex > 0) {
            if (state.dataDXCoils->DXCoil(thisMSHeatPump.DXCoolCoilIndex).IsSecondaryDXCoilInZone) {
                OutsideDryBulbTemp = state.dataZoneTempPredictorCorrector
                                         ->zoneHeatBalance(state.dataDXCoils->DXCoil(thisMSHeatPump.DXCoolCoilIndex).SecZonePtr)
                                         .ZT;
            } else {
                OutsideDryBulbTemp = state.dataEnvrn->OutDryBulbTemp;
            }
        } else {
            OutsideDryBulbTemp = state.dataEnvrn->OutDryBulbTemp;
        }
        FanOutletNode = thisMSHeatPump.FanOutletNode;
        FanInletNode = thisMSHeatPump.FanInletNode;

        state.dataHVACMultiSpdHP->SaveCompressorPLR = 0.0;
        SavePartloadRatio = 0.0;
        MinWaterFlow = 0.0;
        ErrorToler = 0.001;
        // Set inlet air mass flow rate based on PLR and compressor on/off air flow rates
        SetAverageAirFlow(state, MSHeatPumpNum, PartLoadFrac, OnOffAirFlowRatio, SpeedNum, SpeedRatio);

        AirMassFlow = state.dataLoopNodes->Node(InletNode).MassFlowRate;
        // if blow through, simulate fan then coils
        if (thisMSHeatPump.FanPlaceType == BlowThru) {
            SimulateFanComponents(state,
                                  thisMSHeatPump.FanName,
                                  FirstHVACIteration,
                                  thisMSHeatPump.FanNum,
                                  state.dataHVACMultiSpdHP->FanSpeedRatio);
            if (QZnReq < (-1.0 * SmallLoad)) {
                if (OutsideDryBulbTemp > thisMSHeatPump.MinOATCompressorCooling) {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXCoolCoilName,
                                        SpeedRatio,
                                        PartLoadFrac,
                                        thisMSHeatPump.DXCoolCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                    SavePartloadRatio = PartLoadFrac;
                    SaveSpeedRatio = SpeedRatio;
                } else {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXCoolCoilName,
                                        0.0,
                                        0.0,
                                        thisMSHeatPump.DXCoolCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                }
                state.dataHVACMultiSpdHP->SaveCompressorPLR = state.dataDXCoils->DXCoilPartLoadRatio(thisMSHeatPump.DXCoolCoilIndex);
            } else {
                SimDXCoilMultiSpeed(state,
                                    thisMSHeatPump.DXCoolCoilName,
                                    0.0,
                                    0.0,
                                    thisMSHeatPump.DXCoolCoilIndex,
                                    SpeedNum,
                                    thisMSHeatPump.OpMode,
                                    CompressorOp);
            }
            if (thisMSHeatPump.HeatCoilType == MultiSpeedHeatingCoil) {
                if (QZnReq > SmallLoad) {
                    if (OutsideDryBulbTemp > thisMSHeatPump.MinOATCompressorHeating) {
                        SimDXCoilMultiSpeed(state,
                                            thisMSHeatPump.DXHeatCoilName,
                                            SpeedRatio,
                                            PartLoadFrac,
                                            thisMSHeatPump.DXHeatCoilIndex,
                                            SpeedNum,
                                            thisMSHeatPump.OpMode,
                                            CompressorOp);
                        SavePartloadRatio = PartLoadFrac;
                        SaveSpeedRatio = SpeedRatio;
                    } else {
                        SimDXCoilMultiSpeed(state,
                                            thisMSHeatPump.DXHeatCoilName,
                                            0.0,
                                            0.0,
                                            thisMSHeatPump.DXHeatCoilIndex,
                                            SpeedNum,
                                            thisMSHeatPump.OpMode,
                                            CompressorOp);
                    }
                    state.dataHVACMultiSpdHP->SaveCompressorPLR = state.dataDXCoils->DXCoilPartLoadRatio(thisMSHeatPump.DXHeatCoilIndex);
                } else {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXHeatCoilName,
                                        0.0,
                                        0.0,
                                        thisMSHeatPump.DXHeatCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                }
            } else if (thisMSHeatPump.HeatCoilType == Coil_HeatingElectric_MultiStage ||
                       thisMSHeatPump.HeatCoilType == Coil_HeatingGas_MultiStage) {
                if (QZnReq > SmallLoad) {
                    SimulateHeatingCoilComponents(state,
                                                  thisMSHeatPump.HeatCoilName,
                                                  FirstHVACIteration,
                                                  _,
                                                  0,
                                                  _,
                                                  _,
                                                  thisMSHeatPump.OpMode,
                                                  PartLoadFrac,
                                                  SpeedNum,
                                                  SpeedRatio);
                } else {
                    SimulateHeatingCoilComponents(state,
                                                  thisMSHeatPump.HeatCoilName,
                                                  FirstHVACIteration,
                                                  _,
                                                  0,
                                                  _,
                                                  _,
                                                  thisMSHeatPump.OpMode,
                                                  0.0,
                                                  SpeedNum,
                                                  0.0);
                }
            } else {
                CalcNonDXHeatingCoils(state, MSHeatPumpNum, FirstHVACIteration, QZnReq, thisMSHeatPump.OpMode, QCoilActual, PartLoadFrac);
            }
            // Call twice to ensure the fan outlet conditions are updated
            SimulateFanComponents(state,
                                  thisMSHeatPump.FanName,
                                  FirstHVACIteration,
                                  thisMSHeatPump.FanNum,
                                  state.dataHVACMultiSpdHP->FanSpeedRatio);
            if (QZnReq < (-1.0 * SmallLoad)) {
                if (OutsideDryBulbTemp > thisMSHeatPump.MinOATCompressorCooling) {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXCoolCoilName,
                                        SpeedRatio,
                                        PartLoadFrac,
                                        thisMSHeatPump.DXCoolCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                    SavePartloadRatio = PartLoadFrac;
                    SaveSpeedRatio = SpeedRatio;
                } else {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXCoolCoilName,
                                        0.0,
                                        0.0,
                                        thisMSHeatPump.DXCoolCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                }
                state.dataHVACMultiSpdHP->SaveCompressorPLR = state.dataDXCoils->DXCoilPartLoadRatio(thisMSHeatPump.DXCoolCoilIndex);
            } else {
                SimDXCoilMultiSpeed(state,
                                    thisMSHeatPump.DXCoolCoilName,
                                    0.0,
                                    0.0,
                                    thisMSHeatPump.DXCoolCoilIndex,
                                    SpeedNum,
                                    thisMSHeatPump.OpMode,
                                    CompressorOp);
            }
            if (thisMSHeatPump.HeatCoilType == MultiSpeedHeatingCoil) {
                if (QZnReq > SmallLoad) {
                    if (OutsideDryBulbTemp > thisMSHeatPump.MinOATCompressorHeating) {
                        SimDXCoilMultiSpeed(state,
                                            thisMSHeatPump.DXHeatCoilName,
                                            SpeedRatio,
                                            PartLoadFrac,
                                            thisMSHeatPump.DXHeatCoilIndex,
                                            SpeedNum,
                                            thisMSHeatPump.OpMode,
                                            CompressorOp);
                        SavePartloadRatio = PartLoadFrac;
                        SaveSpeedRatio = SpeedRatio;
                    } else {
                        SimDXCoilMultiSpeed(state,
                                            thisMSHeatPump.DXHeatCoilName,
                                            0.0,
                                            0.0,
                                            thisMSHeatPump.DXHeatCoilIndex,
                                            SpeedNum,
                                            thisMSHeatPump.OpMode,
                                            CompressorOp);
                    }
                    state.dataHVACMultiSpdHP->SaveCompressorPLR = state.dataDXCoils->DXCoilPartLoadRatio(thisMSHeatPump.DXHeatCoilIndex);
                } else {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXHeatCoilName,
                                        0.0,
                                        0.0,
                                        thisMSHeatPump.DXHeatCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                }
            } else if (thisMSHeatPump.HeatCoilType == Coil_HeatingElectric_MultiStage ||
                       thisMSHeatPump.HeatCoilType == Coil_HeatingGas_MultiStage) {
                if (QZnReq > SmallLoad) {
                    SimulateHeatingCoilComponents(state,
                                                  thisMSHeatPump.HeatCoilName,
                                                  FirstHVACIteration,
                                                  _,
                                                  0,
                                                  _,
                                                  _,
                                                  thisMSHeatPump.OpMode,
                                                  PartLoadFrac,
                                                  SpeedNum,
                                                  SpeedRatio);
                } else {
                    SimulateHeatingCoilComponents(state,
                                                  thisMSHeatPump.HeatCoilName,
                                                  FirstHVACIteration,
                                                  _,
                                                  0,
                                                  _,
                                                  _,
                                                  thisMSHeatPump.OpMode,
                                                  0.0,
                                                  SpeedNum,
                                                  0.0);
                }
            } else {
                CalcNonDXHeatingCoils(state, MSHeatPumpNum, FirstHVACIteration, QZnReq, thisMSHeatPump.OpMode, QCoilActual, PartLoadFrac);
            }
            //  Simulate supplemental heating coil for blow through fan
            if (thisMSHeatPump.SuppHeatCoilNum > 0) {
                CalcNonDXHeatingCoils(state, MSHeatPumpNum, FirstHVACIteration, SupHeaterLoad, thisMSHeatPump.OpMode, QCoilActual);
            }
        } else { // otherwise simulate DX coils then fan then supplemental heater
            if (QZnReq < (-1.0 * SmallLoad)) {
                if (OutsideDryBulbTemp > thisMSHeatPump.MinOATCompressorCooling) {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXCoolCoilName,
                                        SpeedRatio,
                                        PartLoadFrac,
                                        thisMSHeatPump.DXCoolCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                    SavePartloadRatio = PartLoadFrac;
                    SaveSpeedRatio = SpeedRatio;
                } else {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXCoolCoilName,
                                        0.0,
                                        0.0,
                                        thisMSHeatPump.DXCoolCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                }
                state.dataHVACMultiSpdHP->SaveCompressorPLR = state.dataDXCoils->DXCoilPartLoadRatio(thisMSHeatPump.DXCoolCoilIndex);
            } else {
                SimDXCoilMultiSpeed(state,
                                    thisMSHeatPump.DXCoolCoilName,
                                    0.0,
                                    0.0,
                                    thisMSHeatPump.DXCoolCoilIndex,
                                    SpeedNum,
                                    thisMSHeatPump.OpMode,
                                    CompressorOp);
            }
            if (thisMSHeatPump.HeatCoilType == MultiSpeedHeatingCoil) {
                if (QZnReq > SmallLoad) {
                    if (OutsideDryBulbTemp > thisMSHeatPump.MinOATCompressorHeating) {
                        SimDXCoilMultiSpeed(state,
                                            thisMSHeatPump.DXHeatCoilName,
                                            SpeedRatio,
                                            PartLoadFrac,
                                            thisMSHeatPump.DXHeatCoilIndex,
                                            SpeedNum,
                                            thisMSHeatPump.OpMode,
                                            CompressorOp);
                        SavePartloadRatio = PartLoadFrac;
                        SaveSpeedRatio = SpeedRatio;
                    } else {
                        SimDXCoilMultiSpeed(state,
                                            thisMSHeatPump.DXHeatCoilName,
                                            0.0,
                                            0.0,
                                            thisMSHeatPump.DXHeatCoilIndex,
                                            SpeedNum,
                                            thisMSHeatPump.OpMode,
                                            CompressorOp);
                    }
                    state.dataHVACMultiSpdHP->SaveCompressorPLR = state.dataDXCoils->DXCoilPartLoadRatio(thisMSHeatPump.DXHeatCoilIndex);
                } else {
                    SimDXCoilMultiSpeed(state,
                                        thisMSHeatPump.DXHeatCoilName,
                                        0.0,
                                        0.0,
                                        thisMSHeatPump.DXHeatCoilIndex,
                                        SpeedNum,
                                        thisMSHeatPump.OpMode,
                                        CompressorOp);
                }
            } else if (thisMSHeatPump.HeatCoilType == Coil_HeatingElectric_MultiStage ||
                       thisMSHeatPump.HeatCoilType == Coil_HeatingGas_MultiStage) {
                if (QZnReq > SmallLoad) {
                    SimulateHeatingCoilComponents(state,
                                                  thisMSHeatPump.HeatCoilName,
                                                  FirstHVACIteration,
                                                  _,
                                                  0,
                                                  _,
                                                  _,
                                                  thisMSHeatPump.OpMode,
                                                  PartLoadFrac,
                                                  SpeedNum,
                                                  SpeedRatio);
                } else {
                    SimulateHeatingCoilComponents(state,
                                                  thisMSHeatPump.HeatCoilName,
                                                  FirstHVACIteration,
                                                  _,
                                                  0,
                                                  _,
                                                  _,
                                                  thisMSHeatPump.OpMode,
                                                  0.0,
                                                  SpeedNum,
                                                  0.0);
                }
            } else {
                CalcNonDXHeatingCoils(state, MSHeatPumpNum, FirstHVACIteration, QZnReq, thisMSHeatPump.OpMode, QCoilActual, PartLoadFrac);
            }
            SimulateFanComponents(state,
                                  thisMSHeatPump.FanName,
                                  FirstHVACIteration,
                                  thisMSHeatPump.FanNum,
                                  state.dataHVACMultiSpdHP->FanSpeedRatio);
            //  Simulate supplemental heating coil for draw through fan
            if (thisMSHeatPump.SuppHeatCoilNum > 0) {
                CalcNonDXHeatingCoils(state, MSHeatPumpNum, FirstHVACIteration, SupHeaterLoad, thisMSHeatPump.OpMode, QCoilActual);
            }
        }

        // calculate sensible load met
        Real64 SensibleOutput(0.0); // sensible output rate
        // calculate sensible load met using delta enthalpy at a constant (minimum) humidity ratio)
        Real64 MinHumRat = state.dataLoopNodes->Node(thisMSHeatPump.NodeNumOfControlledZone).HumRat;
        if (state.dataLoopNodes->Node(OutletNode).Temp < state.dataLoopNodes->Node(thisMSHeatPump.NodeNumOfControlledZone).Temp)
            MinHumRat = state.dataLoopNodes->Node(OutletNode).HumRat;
        SensibleOutput = AirMassFlow *
                         Psychrometrics::PsyDeltaHSenFnTdb2W2Tdb1W1(state.dataLoopNodes->Node(OutletNode).Temp,
                                                                    MinHumRat,
                                                                    state.dataLoopNodes->Node(thisMSHeatPump.NodeNumOfControlledZone).Temp,
                                                                    MinHumRat);
        LoadMet = SensibleOutput - thisMSHeatPump.LoadLoss;

        thisMSHeatPump.LoadMet = LoadMet;
    }

    void UpdateMSHeatPump(EnergyPlusData &state, int const MSHeatPumpNum) // Engine driven heat pump number
    {
        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Lixing Gu, FSEC
        //       DATE WRITTEN:    June 2007
        //       MODIFIED         na
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        //  This routine will update MSHP performance and calculate heat recovery rate and crankcase heater power

        // METHODOLOGY EMPLOYED:
        // na

        // Calculate heat recovery
        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        if (thisMSHeatPump.HeatRecActive) {
            MSHPHeatRecovery(state, MSHeatPumpNum);
        }

        if (state.afn->distribution_simulated) {
            auto &thisAirLoopAFNInfo = state.dataAirLoop->AirLoopAFNInfo(thisMSHeatPump.AirLoopNumber);
            thisAirLoopAFNInfo.LoopSystemOnMassFlowrate = state.dataHVACMultiSpdHP->CompOnMassFlow;
            thisAirLoopAFNInfo.LoopSystemOffMassFlowrate = state.dataHVACMultiSpdHP->CompOffMassFlow;
            thisAirLoopAFNInfo.LoopFanOperationMode = thisMSHeatPump.OpMode;
            thisAirLoopAFNInfo.LoopOnOffFanPartLoadRatio = thisMSHeatPump.FanPartLoadRatio;
            thisAirLoopAFNInfo.LoopCompCycRatio = state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum).CycRatio;
        }
    }

    //******************************************************************************

    void ReportMSHeatPump(EnergyPlusData &state, int const MSHeatPumpNum) // Engine driven heat pump number
    {
        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Lixing Gu, FSEC
        //       DATE WRITTEN:    June 2007
        //       MODIFIED         na
        //       RE-ENGINEERED    na

        // PURPOSE OF THIS SUBROUTINE:
        //  This routine will write values to output variables in MSHP

        // METHODOLOGY EMPLOYED:

        // REFERENCES: na

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);
        auto &thisMSHeatPumpReport = state.dataHVACMultiSpdHP->MSHeatPumpReport(MSHeatPumpNum);

        Real64 TimeStepSysSec = state.dataHVACGlobal->TimeStepSys * DataGlobalConstants::SecInHour;
        thisMSHeatPumpReport.ElecPowerConsumption = thisMSHeatPump.ElecPower * TimeStepSysSec; // + &
        thisMSHeatPumpReport.HeatRecoveryEnergy = thisMSHeatPump.HeatRecoveryRate * TimeStepSysSec;

        thisMSHeatPumpReport.AuxElecHeatConsumption = 0.0;
        thisMSHeatPumpReport.AuxElecCoolConsumption = 0.0;

        thisMSHeatPump.AuxElecPower = thisMSHeatPump.AuxOnCyclePower * state.dataHVACMultiSpdHP->SaveCompressorPLR +
                                                 thisMSHeatPump.AuxOffCyclePower * (1.0 - state.dataHVACMultiSpdHP->SaveCompressorPLR);
        if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::CoolingMode) {
            thisMSHeatPumpReport.AuxElecCoolConsumption =
                thisMSHeatPump.AuxOnCyclePower * state.dataHVACMultiSpdHP->SaveCompressorPLR * TimeStepSysSec;
        }
        else if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::HeatingMode) {
            thisMSHeatPumpReport.AuxElecHeatConsumption =
                thisMSHeatPump.AuxOnCyclePower * state.dataHVACMultiSpdHP->SaveCompressorPLR * TimeStepSysSec;
        }

        if (thisMSHeatPump.LastMode == ModeOfOperation::HeatingMode) {
            thisMSHeatPumpReport.AuxElecHeatConsumption +=
                thisMSHeatPump.AuxOffCyclePower * (1.0 - state.dataHVACMultiSpdHP->SaveCompressorPLR) * TimeStepSysSec;
        } else {
            thisMSHeatPumpReport.AuxElecCoolConsumption +=
                thisMSHeatPump.AuxOffCyclePower * (1.0 - state.dataHVACMultiSpdHP->SaveCompressorPLR) * TimeStepSysSec;
        }

        if (thisMSHeatPump.FirstPass) {
            if (!state.dataGlobal->SysSizingCalc) {
                DataSizing::resetHVACSizingGlobals(
                    state, state.dataSize->CurZoneEqNum, state.dataSize->CurSysNum, thisMSHeatPump.FirstPass);
            }
        }

        // reset to 1 in case blow through fan configuration (fan resets to 1, but for blow thru fans coil sets back down < 1)
        state.dataHVACGlobal->OnOffFanPartLoadFraction = 1.0;
    }

    void MSHPHeatRecovery(EnergyPlusData &state, int const MSHeatPumpNum) // Number of the current electric MSHP being simulated
    {
        // SUBROUTINE INFORMATION:
        //       AUTHOR:          Lixing Gu
        //       DATE WRITTEN:    June 2007
        //       MODIFIED:        na
        //       RE-ENGINEERED    Revised to calculate MSHP heat recovery rate based on EIR Chiller heat recovery subroutine
        // PURPOSE OF THIS SUBROUTINE:
        //  Calculate the heat recovered from MSHP

        // Using/Aliasing
        using FluidProperties::GetSpecificHeatGlycol;
        using PlantUtilities::SafeCopyPlantNode;

        // SUBROUTINE PARAMETER DEFINITIONS:
        static constexpr std::string_view RoutineName("MSHPHeatRecovery");

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 HeatRecOutletTemp;   // Heat reclaim outlet temp [C]
	Real64 CpHeatRec;

        // Begin routine
	auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);
        int HeatRecInNode = thisMSHeatPump.HeatRecInletNodeNum;
        int HeatRecOutNode = thisMSHeatPump.HeatRecOutletNodeNum;

        // Inlet node to the heat recovery heat exchanger
        Real64 HeatRecInletTemp = state.dataLoopNodes->Node(HeatRecInNode).Temp;

        // Set heat recovery mass flow rates
        Real64 HeatRecMassFlowRate = state.dataLoopNodes->Node(HeatRecInNode).MassFlowRate;

        Real64 QHeatRec = state.dataHVACGlobal->MSHPWasteHeat;

        if (HeatRecMassFlowRate > 0.0) {

            CpHeatRec =
                GetSpecificHeatGlycol(state,
                                      state.dataPlnt->PlantLoop(thisMSHeatPump.HRPlantLoc.loopNum).FluidName,
                                      HeatRecInletTemp,
                                      state.dataPlnt->PlantLoop(thisMSHeatPump.HRPlantLoc.loopNum).FluidIndex,
                                      RoutineName);

            HeatRecOutletTemp = QHeatRec / (HeatRecMassFlowRate * CpHeatRec) + HeatRecInletTemp;
            if (HeatRecOutletTemp > thisMSHeatPump.MaxHeatRecOutletTemp) {
                HeatRecOutletTemp = max(HeatRecInletTemp, thisMSHeatPump.MaxHeatRecOutletTemp);
                QHeatRec = HeatRecMassFlowRate * CpHeatRec * (HeatRecOutletTemp - HeatRecInletTemp);
            }
        } else {
            HeatRecOutletTemp = HeatRecInletTemp;
            QHeatRec = 0.0;
        }

        SafeCopyPlantNode(state, HeatRecInNode, HeatRecOutNode);
        // changed outputs
        state.dataLoopNodes->Node(HeatRecOutNode).Temp = HeatRecOutletTemp;

        thisMSHeatPump.HeatRecoveryRate = QHeatRec;
        thisMSHeatPump.HeatRecoveryInletTemp = HeatRecInletTemp;
        thisMSHeatPump.HeatRecoveryOutletTemp = HeatRecOutletTemp;
        thisMSHeatPump.HeatRecoveryMassFlowRate = HeatRecMassFlowRate;
    }

    void SetAverageAirFlow(EnergyPlusData &state,
                           int const MSHeatPumpNum,                     // Unit index
                           Real64 const PartLoadRatio,                  // unit part load ratio
                           Real64 &OnOffAirFlowRatio,                   // ratio of compressor ON airflow to average airflow over timestep
                           ObjexxFCL::Optional_int_const SpeedNum,      // Speed number
                           ObjexxFCL::Optional<Real64 const> SpeedRatio // Speed ratio
    )
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Lixing
        //       DATE WRITTEN   June 2007
        //       MODIFIED       na
        //       RE-ENGINEERED  Resived to meet requirements of multispeed heat pump based on the same subroutine
        //                      in PTHP module

        // PURPOSE OF THIS SUBROUTINE:
        // Set the average air mass flow rates using the part load fraction of the heat pump for this time step
        // Set OnOffAirFlowRatio to be used by DX coils

        // Using/Aliasing
        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        int InletNode;              // inlet node number for PTHPNum
        Real64 AverageUnitMassFlow; // average supply air mass flow rate over time step

        state.dataHVACGlobal->MSHPMassFlowRateLow = 0.0;  // Mass flow rate at low speed
        state.dataHVACGlobal->MSHPMassFlowRateHigh = 0.0; // Mass flow rate at high speed

        if (!state.dataZoneEnergyDemand->CurDeadBandOrSetback(thisMSHeatPump.ControlZoneNum) &&
            present(SpeedNum)) {
            if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::HeatingMode) {
                if (SpeedNum == 1) {
                    state.dataHVACMultiSpdHP->CompOnMassFlow = thisMSHeatPump.HeatMassFlowRate(SpeedNum);
                    state.dataHVACMultiSpdHP->CompOnFlowRatio = thisMSHeatPump.HeatingSpeedRatio(SpeedNum);
                    state.dataHVACGlobal->MSHPMassFlowRateLow = thisMSHeatPump.HeatMassFlowRate(1);
                    state.dataHVACGlobal->MSHPMassFlowRateHigh = thisMSHeatPump.HeatMassFlowRate(1);
                } else if (SpeedNum > 1) {
                    state.dataHVACMultiSpdHP->CompOnMassFlow =
                        SpeedRatio * thisMSHeatPump.HeatMassFlowRate(SpeedNum) +
                        (1.0 - SpeedRatio) * thisMSHeatPump.HeatMassFlowRate(SpeedNum - 1);
                    state.dataHVACMultiSpdHP->CompOnFlowRatio =
                        SpeedRatio * thisMSHeatPump.HeatingSpeedRatio(SpeedNum) +
                        (1.0 - SpeedRatio) * thisMSHeatPump.HeatingSpeedRatio(SpeedNum - 1);
                    state.dataHVACGlobal->MSHPMassFlowRateLow = thisMSHeatPump.HeatMassFlowRate(SpeedNum - 1);
                    state.dataHVACGlobal->MSHPMassFlowRateHigh = thisMSHeatPump.HeatMassFlowRate(SpeedNum);
                }
            } else if (thisMSHeatPump.HeatCoolMode == ModeOfOperation::CoolingMode) {
                if (SpeedNum == 1) {
                    state.dataHVACMultiSpdHP->CompOnMassFlow = thisMSHeatPump.CoolMassFlowRate(SpeedNum);
                    state.dataHVACMultiSpdHP->CompOnFlowRatio = thisMSHeatPump.CoolingSpeedRatio(SpeedNum);
                    state.dataHVACGlobal->MSHPMassFlowRateLow = thisMSHeatPump.CoolMassFlowRate(1);
                    state.dataHVACGlobal->MSHPMassFlowRateHigh = thisMSHeatPump.CoolMassFlowRate(1);
                } else if (SpeedNum > 1) {
                    state.dataHVACMultiSpdHP->CompOnMassFlow =
                        SpeedRatio * thisMSHeatPump.CoolMassFlowRate(SpeedNum) +
                        (1.0 - SpeedRatio) * thisMSHeatPump.CoolMassFlowRate(SpeedNum - 1);
                    state.dataHVACMultiSpdHP->CompOnFlowRatio =
                        SpeedRatio * thisMSHeatPump.CoolingSpeedRatio(SpeedNum) +
                        (1.0 - SpeedRatio) * thisMSHeatPump.CoolingSpeedRatio(SpeedNum - 1);
                    state.dataHVACGlobal->MSHPMassFlowRateLow = thisMSHeatPump.CoolMassFlowRate(SpeedNum - 1);
                    state.dataHVACGlobal->MSHPMassFlowRateHigh = thisMSHeatPump.CoolMassFlowRate(SpeedNum);
                }
            }
        }
        InletNode = thisMSHeatPump.AirInletNodeNum;

        // Set up fan flow rate during compressor off time
        if (thisMSHeatPump.OpMode == ContFanCycCoil && present(SpeedNum)) {
            if (thisMSHeatPump.AirFlowControl == AirflowControl::UseCompressorOnFlow &&
                state.dataHVACMultiSpdHP->CompOnMassFlow > 0.0) {
                if (thisMSHeatPump.LastMode == ModeOfOperation::HeatingMode) {
                    state.dataHVACMultiSpdHP->CompOffMassFlow = thisMSHeatPump.HeatMassFlowRate(SpeedNum);
                    state.dataHVACMultiSpdHP->CompOffFlowRatio = thisMSHeatPump.HeatingSpeedRatio(SpeedNum);
                } else {
                    state.dataHVACMultiSpdHP->CompOffMassFlow = thisMSHeatPump.CoolMassFlowRate(SpeedNum);
                    state.dataHVACMultiSpdHP->CompOffFlowRatio = thisMSHeatPump.CoolingSpeedRatio(SpeedNum);
                }
            }
        }

        if (present(SpeedNum)) {
            if (SpeedNum > 1) {
                AverageUnitMassFlow = state.dataHVACMultiSpdHP->CompOnMassFlow;
                state.dataHVACMultiSpdHP->FanSpeedRatio = state.dataHVACMultiSpdHP->CompOnFlowRatio;
            } else {
                AverageUnitMassFlow =
                    (PartLoadRatio * state.dataHVACMultiSpdHP->CompOnMassFlow) + ((1 - PartLoadRatio) * state.dataHVACMultiSpdHP->CompOffMassFlow);
                if (state.dataHVACMultiSpdHP->CompOffFlowRatio > 0.0) {
                    state.dataHVACMultiSpdHP->FanSpeedRatio = (PartLoadRatio * state.dataHVACMultiSpdHP->CompOnFlowRatio) +
                                                              ((1 - PartLoadRatio) * state.dataHVACMultiSpdHP->CompOffFlowRatio);
                } else {
                    state.dataHVACMultiSpdHP->FanSpeedRatio = state.dataHVACMultiSpdHP->CompOnFlowRatio;
                }
            }
        } else {
            AverageUnitMassFlow =
                (PartLoadRatio * state.dataHVACMultiSpdHP->CompOnMassFlow) + ((1 - PartLoadRatio) * state.dataHVACMultiSpdHP->CompOffMassFlow);
            if (state.dataHVACMultiSpdHP->CompOffFlowRatio > 0.0) {
                state.dataHVACMultiSpdHP->FanSpeedRatio =
                    (PartLoadRatio * state.dataHVACMultiSpdHP->CompOnFlowRatio) + ((1 - PartLoadRatio) * state.dataHVACMultiSpdHP->CompOffFlowRatio);
            } else {
                state.dataHVACMultiSpdHP->FanSpeedRatio = state.dataHVACMultiSpdHP->CompOnFlowRatio;
            }
        }

        //!!LKL Discrepancy with > 0
        if (GetCurrentScheduleValue(state, thisMSHeatPump.AvaiSchedPtr) == 0.0) {
            state.dataLoopNodes->Node(InletNode).MassFlowRate = 0.0;
            OnOffAirFlowRatio = 0.0;
        } else {
            state.dataLoopNodes->Node(InletNode).MassFlowRate = AverageUnitMassFlow;
            state.dataLoopNodes->Node(InletNode).MassFlowRateMaxAvail = AverageUnitMassFlow;
            if (AverageUnitMassFlow > 0.0) {
                OnOffAirFlowRatio = state.dataHVACMultiSpdHP->CompOnMassFlow / AverageUnitMassFlow;
            } else {
                OnOffAirFlowRatio = 0.0;
            }
        }
    }

    void CalcNonDXHeatingCoils(EnergyPlusData &state,
                               int const MSHeatPumpNum,       // multispeed heatpump index
                               bool const FirstHVACIteration, // flag for first HVAC iteration in the time step
                               Real64 const HeatingLoad,      // supplemental coil load to be met by unit (watts)
                               int const FanMode,             // fan operation mode
                               Real64 &HeatCoilLoadmet,       // Heating Load Met
                               ObjexxFCL::Optional<Real64 const> PartLoadFrac)
    {

        // SUBROUTINE INFORMATION:
        //       AUTHOR         Bereket Nigusse, FSEC/UCF
        //       DATE WRITTEN   January 2012
        //       MODIFIED       na
        //       RE-ENGINEERED  na

        // PURPOSE OF THIS SUBROUTINE:
        // This subroutine simulates the four non dx heating coil types: Gas, Electric, hot water and steam.

        // METHODOLOGY EMPLOYED:
        // Simply calls the different heating coil component.  The hot water flow rate matching the coil load
        // is calculated iteratively.

        // REFERENCES:
        // na

        // USE STATEMENTS:

        // Using/Aliasing
        using DataHVACGlobals::SmallLoad;

        using General::SolveRoot;
        using HeatingCoils::SimulateHeatingCoilComponents;
        using PlantUtilities::SetComponentFlowRate;
        using SteamCoils::SimulateSteamCoilComponents;
        using WaterCoils::SimulateWaterCoilComponents;

        // Locals
        static constexpr std::string_view CurrentModuleObject("AirLoopHVAC:UnitaryHeatPump:AirToAir:MultiSpeed");

        // SUBROUTINE ARGUMENT DEFINITIONS:

        // SUBROUTINE PARAMETER DEFINITIONS:
        Real64 constexpr ErrTolerance(0.001); // convergence limit for hotwater coil
        int constexpr SolveMaxIter(50);

        // INTERFACE BLOCK SPECIFICATIONS
        // na

        // DERIVED TYPE DEFINITIONS
        // na

        // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        Real64 QCoilActual;     // actual heating load met
        Real64 mdot;            // heating coil steam or hot water mass flow rate
        Real64 MinWaterFlow;    // coil minimum hot water mass flow rate, kg/s
        Real64 MaxHotWaterFlow; // coil maximum hot water mass flow rate, kg/s
        Real64 HotWaterMdot;    // actual hot water mass flow rate
        int SolFlag;

        int HeatCoilType;
        int HeatCoilNum;
        Real64 MaxCoilFluidFlow;
        Real64 SteamCoilHeatingLoad;
        int CoilControlNode;
        int CoilOutletNode;
        PlantLocation plantLoc{};

        QCoilActual = 0.0;

        auto &thisMSHeatPump = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);

        if (present(PartLoadFrac)) {
            HeatCoilType = thisMSHeatPump.HeatCoilType;
            state.dataHVACMultiSpdHP->HeatCoilName = thisMSHeatPump.HeatCoilName;
            HeatCoilNum = thisMSHeatPump.HeatCoilNum;
            MaxCoilFluidFlow = thisMSHeatPump.MaxCoilFluidFlow;
            CoilControlNode = thisMSHeatPump.CoilControlNode;
            CoilOutletNode = thisMSHeatPump.CoilOutletNode;
            plantLoc = thisMSHeatPump.plantLoc;
        } else {
            HeatCoilType = thisMSHeatPump.SuppHeatCoilType;
            state.dataHVACMultiSpdHP->HeatCoilName = thisMSHeatPump.SuppHeatCoilName;
            HeatCoilNum = thisMSHeatPump.SuppHeatCoilNum;
            MaxCoilFluidFlow = thisMSHeatPump.MaxSuppCoilFluidFlow;
            CoilControlNode = thisMSHeatPump.SuppCoilControlNode;
            CoilOutletNode = thisMSHeatPump.SuppCoilOutletNode;
            plantLoc = thisMSHeatPump.SuppPlantLoc;
        }

        thisMSHeatPump.HotWaterPlantLoc = plantLoc;
        thisMSHeatPump.HotWaterCoilControlNode = CoilControlNode;
        thisMSHeatPump.HotWaterCoilOutletNode = CoilOutletNode;
        thisMSHeatPump.HotWaterCoilName = state.dataHVACMultiSpdHP->HeatCoilName;
        thisMSHeatPump.HotWaterCoilNum = HeatCoilNum;

        if (HeatingLoad > SmallLoad) {

            switch (HeatCoilType) {
            case SuppHeatingCoilGas:
            case SuppHeatingCoilElec: {
                SimulateHeatingCoilComponents(
                    state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatingLoad, HeatCoilNum, QCoilActual, true, FanMode);
            } break;
            case Coil_HeatingWater: {
                if (present(PartLoadFrac)) {
                    MaxHotWaterFlow = MaxCoilFluidFlow * PartLoadFrac;
                    SetComponentFlowRate(state, MaxHotWaterFlow, CoilControlNode, CoilOutletNode, plantLoc);
                    SimulateWaterCoilComponents(state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatCoilNum, QCoilActual, FanMode);
                } else {
                    MaxHotWaterFlow = MaxCoilFluidFlow;
                    SetComponentFlowRate(state, MaxHotWaterFlow, CoilControlNode, CoilOutletNode, plantLoc);
                    SimulateWaterCoilComponents(state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatCoilNum, QCoilActual, FanMode);
                    if (QCoilActual > (HeatingLoad + SmallLoad)) {
                        // control water flow to obtain output matching HeatingLoad
                        SolFlag = 0;
                        MinWaterFlow = 0.0;
                        auto f = [&state, MSHeatPumpNum, FirstHVACIteration, HeatingLoad](Real64 const HWFlow) {
                            // Calculates residual function (QCoilActual - SupHeatCoilLoad) / SupHeatCoilLoad
                            // coil actual output depends on the hot water flow rate which is varied to minimize the residual.
                            Real64 targetHeatingCoilLoad = HeatingLoad;
                            Real64 calcHeatingCoilLoad = targetHeatingCoilLoad;
                            Real64 mdot = HWFlow;
                            auto &hp = state.dataHVACMultiSpdHP->MSHeatPump(MSHeatPumpNum);
                            PlantUtilities::SetComponentFlowRate(
                                state, mdot, hp.HotWaterCoilControlNode, hp.HotWaterCoilOutletNode, hp.HotWaterPlantLoc);
                            // simulate the hot water supplemental heating coil
                            WaterCoils::SimulateWaterCoilComponents(
                                state, hp.HotWaterCoilName, FirstHVACIteration, hp.HotWaterCoilNum, calcHeatingCoilLoad, hp.OpMode);
                            if (targetHeatingCoilLoad != 0.0) {
                                return (calcHeatingCoilLoad - targetHeatingCoilLoad) / targetHeatingCoilLoad;
                            } else { // Autodesk:Return Condition added to assure return value is set
                                return 0.0;
                            }
                        };
                        SolveRoot(state, ErrTolerance, SolveMaxIter, SolFlag, HotWaterMdot, f, MinWaterFlow, MaxHotWaterFlow);
                        if (SolFlag == -1) {
                            if (thisMSHeatPump.HotWaterCoilMaxIterIndex == 0) {
                                ShowWarningMessage(state,
                                                   format("CalcNonDXHeatingCoils: Hot water coil control failed for {}=\"{}\"",
                                                          CurrentModuleObject,
                                                          thisMSHeatPump.Name));
                                ShowContinueErrorTimeStamp(state, "");
                                ShowContinueError(state,
                                                  format("  Iteration limit [{}] exceeded in calculating hot water mass flow rate", SolveMaxIter));
                            }
                            ShowRecurringWarningErrorAtEnd(
                                state,
                                format("CalcNonDXHeatingCoils: Hot water coil control failed (iteration limit [{}]) for {}=\"{}",
                                       SolveMaxIter,
                                       CurrentModuleObject,
                                       thisMSHeatPump.Name),
                                thisMSHeatPump.HotWaterCoilMaxIterIndex);
                        } else if (SolFlag == -2) {
                            if (thisMSHeatPump.HotWaterCoilMaxIterIndex2 == 0) {
                                ShowWarningMessage(state,
                                                   format("CalcNonDXHeatingCoils: Hot water coil control failed (maximum flow limits) for {}=\"{}\"",
                                                          CurrentModuleObject,
                                                          thisMSHeatPump.Name));
                                ShowContinueErrorTimeStamp(state, "");
                                ShowContinueError(state, "...Bad hot water maximum flow rate limits");
                                ShowContinueError(state, format("...Given minimum water flow rate={:.3R} kg/s", MinWaterFlow));
                                ShowContinueError(state, format("...Given maximum water flow rate={:.3R} kg/s", MaxHotWaterFlow));
                            }
                            ShowRecurringWarningErrorAtEnd(state,
                                                           "CalcNonDXHeatingCoils: Hot water coil control failed (flow limits) for " +
                                                               std::string{CurrentModuleObject} + "=\"" + thisMSHeatPump.Name + "\"",
                                                           thisMSHeatPump.HotWaterCoilMaxIterIndex2,
                                                           MaxHotWaterFlow,
                                                           MinWaterFlow,
                                                           _,
                                                           "[kg/s]",
                                                           "[kg/s]");
                        }
                        // simulate hot water supplemental heating coil
                        SimulateWaterCoilComponents(
                            state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatCoilNum, QCoilActual, FanMode);
                    }
                }
            } break;
            case Coil_HeatingSteam: {
                if (present(PartLoadFrac)) {
                    mdot = thisMSHeatPump.MaxCoilFluidFlow * PartLoadFrac;
                    SteamCoilHeatingLoad = HeatingLoad * PartLoadFrac;
                } else {
                    mdot = thisMSHeatPump.MaxCoilFluidFlow;
                    SteamCoilHeatingLoad = HeatingLoad;
                }
                SetComponentFlowRate(state, mdot, CoilControlNode, CoilOutletNode, plantLoc);
                // simulate steam supplemental heating coil
                SimulateSteamCoilComponents(
                    state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatCoilNum, SteamCoilHeatingLoad, QCoilActual, FanMode);
            } break;
            default:
                break;
            }

        } else { // end of IF (HeatingLoad > SmallLoad) THEN

            switch (HeatCoilType) {
            case SuppHeatingCoilGas:
            case SuppHeatingCoilElec: {
                SimulateHeatingCoilComponents(
                    state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatingLoad, HeatCoilNum, QCoilActual, true, FanMode);
            } break;
            case Coil_HeatingWater: {
                mdot = 0.0;
                SetComponentFlowRate(state, mdot, CoilControlNode, CoilOutletNode, plantLoc);
                SimulateWaterCoilComponents(state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatCoilNum, QCoilActual, FanMode);
            } break;
            case Coil_HeatingSteam: {
                mdot = 0.0;
                SetComponentFlowRate(state, mdot, CoilControlNode, CoilOutletNode, plantLoc);
                // simulate the steam supplemental heating coil
                SimulateSteamCoilComponents(
                    state, state.dataHVACMultiSpdHP->HeatCoilName, FirstHVACIteration, HeatCoilNum, HeatingLoad, QCoilActual, FanMode);
            } break;
            default:
                break;
            }
        }
        HeatCoilLoadmet = QCoilActual;
    }

} // namespace HVACMultiSpeedHeatPump

} // namespace EnergyPlus

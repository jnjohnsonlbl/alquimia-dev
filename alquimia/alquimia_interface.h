/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
** through Lawrence Berkeley National Laboratory (subject to receipt of any 
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
** 
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
**
** If you have questions about your rights to use or distribute this software, 
** please contact Berkeley Lab's Technology Transfer and Intellectual Property 
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
** 
** NOTICE.  This software was developed under funding from the U.S. Department 
** of Energy.  As such, the U.S. Government has been granted for itself and 
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
** license in the Software to reproduce, prepare derivative works, and perform 
** publicly and display publicly.  Beginning five (5) years after the date 
** permission to assert copyright is obtained from the U.S. Department of Energy, 
** and subject to any subsequent five (5) year renewals, the U.S. Government is 
** granted for itself and others acting on its behalf a paid-up, nonexclusive, 
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly, 
** and to permit others to do so.
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

#ifndef ALQUIMIA_INTERFACE_H_
#define ALQUIMIA_INTERFACE_H_

/*******************************************************************************
 **
 ** C implementation of the alquimia interface.
 **
 ******************************************************************************/

#include "alquimia/alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus 

  // NOTE(bja): AlquimiaData is a convenience container, and not a data
  // structure that should be passed between the driver and
  // engine. Conditions are not included here because they are short-lived
  // data structures. You need one of these per thread. 
  typedef struct {
    void* engine_state;
    AlquimiaSizes sizes;
    AlquimiaEngineFunctionality functionality;
    AlquimiaState state;
    AlquimiaProperties properties;
    AlquimiaAuxiliaryData aux_data;
    AlquimiaProblemMetaData meta_data;
    AlquimiaAuxiliaryOutputData aux_output;
  } AlquimiaData;

  // NOTE(bja): The alquimia interface should contain nothing but
  // function pointers. In order to thread chemistry, we need just one
  // interface object.
  typedef struct {
    // Reads data from the input file into containers, initializing memory
    // where needed (includes reading the database, swapping the basis, etc.)
    void (*Setup)(
        const char* input_filename,
        bool hands_off,
        void* pft_engine_state,
        AlquimiaSizes* sizes,
        AlquimiaEngineFunctionality* functionality,
        AlquimiaEngineStatus* status);

    // Gracefully shuts down the engine, cleaning memory 
    void (*Shutdown)(
      void* pft_engine_state,
      AlquimiaEngineStatus* status);

    // Applies a geogeochemical condition to enforce boundary/initial conditions. 
    // Called once for each IC/BC. 
    void (*ProcessCondition)(
        void* pft_engine_state,
        AlquimiaGeochemicalCondition* condition,
        AlquimiaProperties* props,
        AlquimiaState* state,
        AlquimiaAuxiliaryData* aux_data,
        AlquimiaEngineStatus* status);

    // Takes one reaction step in operator split mode 
    void (*ReactionStepOperatorSplit)(
        void* pft_engine_state,
        double delta_t,
        AlquimiaProperties* props,
        AlquimiaState* state,
        AlquimiaAuxiliaryData* aux_data,
        AlquimiaEngineStatus* status);
    
    // Computes the Jacobian matrix and the residual vector for a reaction
    // at a single site, given the properties and state. The equilibrium (J_eq) 
    // and kinetic (J_kin) Jacobians are Nc x Nc dense matrices, where Nc is 
    // the number of concentrations; J_eq and J_kin are stored in column-major 
    // order. The kinetic contribution to the residual will be stored in the 
    // (Nc x 1) vector R_kin.
    //
    // (NOTE: This does not currently support surface complexation. For that, 
    //  NOTE: we'll need stuff mentioned in the PFlotran interface.)
    //
    // UNITS: 
    //   J_eq: [kg H2O / L H2O], or [molarity / molality]
    //   J_kin [kg H2O / sec]
    //   R_kin [mol / sec]
    void (*ComputeJacobianAndResidual)(
        void* pft_engine_state,
        AlquimiaProperties* props,
        AlquimiaState* state,
        AlquimiaAuxiliaryData* aux_data,
        double* J_eq,
        double* J_kin,
        double* R_kin,
        AlquimiaEngineStatus* status);
    
    // Retrieves selected geochemical data for output, storing it in 
    // aux_out.
    void (*GetAuxiliaryOutput)(
        void* pft_engine_state,
        AlquimiaProperties* props,
        AlquimiaState* state,
        AlquimiaAuxiliaryData* aux_data,
        AlquimiaAuxiliaryOutputData* aux_out,
        AlquimiaEngineStatus* status);
    
    // Retrieves metadata for the problem, storing it in meta_data.
    void (*GetProblemMetaData)(
        void* pft_engine_state,
        AlquimiaProblemMetaData* meta_data,
        AlquimiaEngineStatus* status);
  } AlquimiaInterface;


  // Creates an Alquimia interface for the chemistry engine identified by 
  // engine_name.
  void CreateAlquimiaInterface(const char* const engine_name,
                               AlquimiaInterface* interface,
                               AlquimiaEngineStatus* status);

#ifdef __cplusplus
}
#endif // __cplusplus 

#endif  // ALQUIMIA_INTERFACE_H_ 

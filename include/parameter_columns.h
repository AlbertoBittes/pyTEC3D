/**
 * @file parameter_columns.h
 * @brief Column layout of the material and real-constant rows cached from the Python model.
 */
#pragma once

/**
 * @brief Column indices into the cached material property row.
 * @details MUR_{X,Y,Z} = relative permeability; RHO_{X,Y,Z} = resistivity; M_{X,Y,Z} = magnetization.
 */
enum MaterialColumnIndex : unsigned int
{
    MAT_MUR_X,
    MAT_MUR_Y,
    MAT_MUR_Z,
    MAT_RHO_X,
    MAT_RHO_Y,
    MAT_RHO_Z,
    MAT_M_X,
    MAT_M_Y,
    MAT_M_Z,
    MAT_NUM_COLS
};

/**
 * @brief Column indices into the cached real-constants row.
 * @details THK = thickness; CS = cross section; TURNS = number of turns; VOL = coil volume;
 * RAD = coil radius; T_{X,Y,Z} = current-direction unit vector; R = resistance; S = symmetry factor.
 */
enum ConstantColumnIndex : unsigned int
{
    CONST_THK = 0,
    CONST_CS = 1,
    CONST_TURNS = 2,
    CONST_VOL = 3,
    CONST_RAD = 4,
    CONST_T_X = 5,
    CONST_T_Y = 6,
    CONST_T_Z = 7,
    CONST_R = 8,
    CONST_S = 9,
    CONST_NUM_COLS = 10
};

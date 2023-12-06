
// The data provided for this two-shaft engine is at a DESIGN POINT at cruise conditions.
// h = 35000 (ft)
// M_f = 0.78;
// pi_f = 1.5;              FAN
// pi_LP = 2.5;             BOOSTER
// pi_HP = 18;              HIGH PRESSURE COMPRESSOR
// T_04 = 1500 (K)          TET
// pi_LPT = 11.13;          LOW PRESSURE TURBINE

#include <iostream>
#include <iomanip>
#include "Design.h"
#include "Off_Design.h"

int main() {
    std::cout << "DESIGN POINT CALCULATION " << std::endl;
    std::cout << "============================== " << std::endl;
    std::cout << "      " << std::endl;

    Design core_nozzle;
    core_nozzle.callDesignCoreChoke();

    Design bypass_nozzle;
    bypass_nozzle.callDesignBypassChoke();

    std::cout << "DESIGN - CORE " << std::endl;
    std::cout << "============================== " << std::endl;
    std::cout << "h (m) = " << Design().h << std::setprecision(8) << std::endl;
    //std::cout << "h (m) = " << Design().t << std::setprecision(8) << std::endl;
    std::cout << "T_a = " << Design().T_a << std::endl;
    std::cout << "P_a = " << Design().P_a << std::endl;
    std::cout << "V_f = " << Design().V_f << std::endl;
    std::cout << "T_02 = " << Design().T_02 << std::endl;
    std::cout << "P_02 = " << Design().P_02 << std::endl;
    std::cout << "T_013 = " << Design().T_013 << std::endl;
    std::cout << "P_013 = " << Design().P_013 << std::endl;
    std::cout << "T_023 = " << Design().T_023 << std::endl;
    std::cout << "P_023 = " << Design().P_023 << std::endl;
    std::cout << "T_03 = " << Design().T_03 << std::endl;
    std::cout << "P_03 = " << Design().P_03 << std::endl;
    std::cout << "Overall Pressure Ratio - Compressor = " << Design().pi_c << std::endl;
    std::cout << "T_04 = " << Design().T_04 << std::endl;
    std::cout << "P_04 = " << Design().P_04 << std::endl;
    std::cout << "Fuel-air ratio = " << Design().f << std::endl;
    std::cout << "mbar_4 = " << Design().m4_bar << std::endl;       // = mbar_9 when choked
    std::cout << "T_045 = " << Design().T_045 << std::endl;
    std::cout << "P_045 = " << Design().P_045 << std::endl;
    std::cout << "mbar_45 = " << Design().m45_bar << std::endl;     // = mbar_4
    std::cout << "T_05 = " << Design().T_05 << std::endl;
    std::cout << "P_05 = " << Design().P_05 << std::endl;
    std::cout << "T_9 = " << Design().T9 << std::endl;
    std::cout << "P_9 = " << Design().P9 << std::endl;
    std::cout << "V_9 = " << Design().V9 << std::endl;
    std::cout << "mbar_9 = " << Design().m9_bar << std::endl;
    std::cout << "Mass low per unit area - HP (kg/sm^2) = " << Design().D_HP_D << std::endl;
    std::cout << "Mass low per unit area - LP (kg/sm^2) = " << Design().D_LP_D << std::endl;
    std::cout << "Mass low per unit area - CORE (kg/sm^2) = " << Design().D_CoreN_D << std::endl;
    std::cout << "Area ratio - CORE [A45/A9] = " << Design().phi2 << std::endl;
    std::cout << "T_04/T_02 = " << Design().T_ratio << std::endl;

    std::cout << "\nDESIGN - BYPASS " << std::endl;
    std::cout << "============================== " << std::endl;
    std::cout << "T_019 = " << Design().T_019 << std::endl;
    std::cout << "P_019 = " << Design().P_019 << std::endl;
    std::cout << "T_19 = " << Design().T19 << std::endl;
    std::cout << "P_19 = " << Design().P19 << std::endl;
    std::cout << "V_19 = " << Design().V19 << std::endl;
    std::cout << "mbar_19 = " << Design().m19_bar << std::endl;
    std::cout << "BPR = " << Design().B_D << std::endl;
    std::cout << "pi_LPT = " << Design().piLPT << std::endl;
    std::cout << "pi_HPC = " << Design().piHPC << std::endl;
    std::cout << "pi_LPC = " << Design().piLPC << std::endl;
    std::cout << "pi_OPR = " << Design().piOPR << std::endl;
    std::cout << "Mass low per unit area - BYPASS (kg/sm^2) = " << Design().D_BypN_D << std::endl;
    std::cout << "Area ratio - BYPASS [A19/A4] = " << Design().phi1 << std::endl;

    std::cout << "\nSpecific thrust (F_s) = " << Design().Fs << std::endl;
    std::cout << "SFC = " << Design().sfc << std::endl;
    std::cout << "Net thrust (F_n) = " << Design().Fn << std::endl;
    std::cout << "Fuel-air ratio = " << Design().f << std::endl;
    std::cout << "Propulsive efficiency = " << Design().eta_p << std::endl;
    std::cout << "Thermal efficiency = " << Design().eta_t << std::endl;
    std::cout << "Overall efficiency = " << Design().eta_o << std::endl;

    std::cout << "mbar_2 = " << Design().m2_bar << std::endl;
    std::cout << "mbar_23 = " << Design().m23_bar << std::endl;
    std::cout << "mdot_9 = " << Design().mdot9 << std::endl;
    std::cout << "mdot_19 = " << Design().mdot19 << std::endl;
    std::cout << "mdot_e = " << Design().mdote << std::endl;

    std::cout << "\nArea 4 = " << Design().A4 << std::endl;
    std::cout << "Area 45 = " << Design().A45 << std::endl;
    std::cout << "Area 9 = " << Design().A9 << std::endl;
    std::cout << "Area 19 = " << Design().A19 << std::endl;
    //Design core_nozzle;
    //core_nozzle.callDesignCoreChoke();

    //Design bypass_nozzle;
    //bypass_nozzle.callDesignBypassChoke();

    std::cout << "\nOFF-DESIGN CALCULATION " << std::endl;
    std::cout << "============================== " << std::endl;
    std::cout << "      " << std::endl;
    Off_Design OffDes;
    OffDes.callOffDesignFunction();


    return 0;
}

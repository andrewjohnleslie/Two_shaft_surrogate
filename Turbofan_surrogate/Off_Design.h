//
// Created by Andrew on 21/11/2023.
//

#ifndef TURBOFAN_MANOEUVRES_SIMPLE_OFF_DESIGN_H
#define TURBOFAN_MANOEUVRES_SIMPLE_OFF_DESIGN_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <iomanip>
#include "Design.h"

class Off_Design : public Design{

public:
    int steps = 5/0.1;   // 50 rows


    // Set off-design conditions (... km, ... m/s)
    //double h_ft_OD = 0;                   // h in ft
    //double h_OD = h_ft_OD * 0.3048;              // CONVERTS ALTITUDE TO M
    double h_OD = 20000;                    // h in m

    // USE THIS BLOCK OF CODE WHEN   [h <= 11km]   - TROPOSPHERIC CONDITIONS
    //double t_OD = (6356577 * h_OD)/(6356577 + h_OD);          //Geo-potential height
    //double Ta_OD = T_std - (6.5 * (t_OD/1000));
    //double Ta_OD = T_std - (0.0065 * h_OD);
    //double Pa_OD = P_std * (pow( (1 - (0.0065 * (h_OD/T_std))), 5.2561));

    // USE THIS BLOCK OF CODE WHEN   [h > 11km <= 20km]   - STRATOSPHERIC CONDITIONS
    double P_trop = 22625.04;
    double g = 9.81;
    double Ta_OD = 216.65;
    double ind = ((-1 * g)/(R * Ta_OD)) * (h_OD - 11000);
    double Pa_OD = P_trop * exp(ind);

//public:
    //GUESS THE INITIAL VALUES FOR PI_FAN
    //
    //double pif_min = 1.01;
    //double pif_max = 2.5;
    //double pif_OD = 1.4;

    //double ERR = 100;
    //double count = 0;
    //double m9_bar_OD = 1;
    //double m19_bar_OD = 1;

    void myOffDesignFunction() {

        std::ofstream outputFile("output.csv");
        if (outputFile.is_open()) {

    for(int i = 0; i <= steps; ++i) {
        //double Mf_OD = 0.9;
        double Mf_OD = i * 0.1;

        double a_OD = pow((gamma * R * Ta_OD), 0.5);
        double Vf_OD = a_OD * Mf_OD;

        //Member variables
        double T02_OD = Ta_OD * (1 + ((gamma - 1) / 2) * pow(Mf_OD, 2));
        double P02_OD = Pa_OD * pow((1 + ((gamma - 1) / 2) * pow(Mf_OD, 2)), gamma / (gamma - 1));
        double T04_OD = 1305;
        //double B_OD = 13.377544; //10;
        double Tratio_OD = T04_OD / T02_OD;

        //double ERR = 100;
        //double count = 0;

//public:
        //GUESS THE INITIAL VALUES FOR PI_FAN
        //
        double pif_min = 1.01;
        double pif_max = 2.5;
        double pif_OD = 1.4;

        double ERR = 100;
        double count = 0;
        double m9_bar_OD = 1;
        double m19_bar_OD = 1;

        //void myOffDesignFunction() {
        do {

            //double count = 0;


            double P013_OD = P02_OD * pif_OD;
            double T013_OD = T02_OD * pow(pif_OD, ((gamma - 1) / (e_c * gamma)));
            double D_BypN_OD = D_BypN_D * sqrt(T_013) * P013_OD / (P_013 * sqrt(T013_OD));

            double T023_OD = T02_OD + kt * (T013_OD - T02_OD);  //first t02 wa at od
            double P023_OD = P02_OD * pow((T023_OD / T02_OD), ((e_c * gamma) / (gamma - 1)));
            double P03_OD = P023_OD * pow(((cpe / cp) * (T04_OD / T023_OD) * kHP + 1), ((e_c * gamma) / (gamma - 1)));
            double T03_OD = T023_OD * pow((P03_OD / P023_OD), ((gamma - 1) / (e_c * gamma)));
            double piLPC_OD = P023_OD / P02_OD;
            double piHPC_OD = P03_OD / P023_OD;
            double piOPR_OD = piLPC_OD * piHPC_OD;

            double P04_OD = P03_OD;
            double T045_OD = T04_OD * (1 - kHP);
            double P045_OD = P04_OD * pow((T045_OD / T04_OD), (gamma_e / (e_t * (gamma_e - 1))));
            double D_HP_ODa = D_HP_D * sqrt(T_04) * P04_OD / (P_04 * sqrt(T04_OD));
            double B_OD = phi1 * D_BypN_OD / D_HP_ODa;                                                          // Bypass ratio - OFF-DESIGN
            double T05_OD = T045_OD * (1 - ((cp / cpe) * (T02_OD / T045_OD) * (B_OD + kt) *
                                            (pow(pif_OD, (gamma - 1) / (e_c * gamma)) - 1)));
            double P05_OD = P045_OD * pow((T05_OD / T045_OD), (gamma_e / (e_t * (gamma_e - 1))));
            double piLPT_OD = (P045_OD / P05_OD);
            double piHPT_OD = (P04_OD / P045_OD);

            double P9_OD = Pa_OD;
            double T9_OD = T05_OD * pow((P9_OD / P05_OD), ((gamma_e - 1) / gamma_e));
            double V9_OD = sqrt(2 * cpe * (T05_OD - T9_OD));

            double y_core_OD = P05_OD / Pa_OD;                        // OFF DESIGN: Core nozzle pressure ratio
            double T09_OD = T05_OD;

            //double m9_bar_OD = 1;

            if (y_core_OD <= Pcrit_core) {                   // If P05/Pa < P05/Pcrit (Core Nozzle PR < Critical PR)
                std::cout << std::endl;
                std::cout << "OFF DESIGN: Core nozzle is not choked! \n" << std::endl;

                double *m9OD = &m9_bar_OD;
                *m9OD = (gamma_e / (gamma_e - 1)) * sqrt(2 * (pow((P05_OD / P9_OD), (-2 / gamma_e)) -
                                                              pow((P05_OD / P9_OD),
                                                                  (-(gamma_e + 1) / gamma_e))));     // Unchoked
            } else {
                std::cout << std::endl;
                std::cout << "OFF DESIGN: Core nozzle is choked! \n " << std::endl;

                double *m9OD = &m9_bar_OD;
                *m9OD = m4_bar;        // Choked
            }

            double P019_OD = P013_OD;
            double T019_OD = T013_OD;
            double P19_OD = Pa_OD;
            double T19_OD = T019_OD * pow((P19_OD / P019_OD), ((gamma - 1) / gamma));
            double V19_OD = sqrt(2 * cp * (T019_OD - T19_OD));

            double y_bypass_OD = P019_OD / Pa_OD;             // OFF DESIGN: Bypass nozzle pressure ratio

            if (y_bypass_OD <= Pcrit_bypass) {                   // If P05/Pa < P05/Pcrit (Core Nozzle PR < Critical PR)
                std::cout << std::endl;
                std::cout << "OFF DESIGN: Bypass nozzle is not choked! \n" << std::endl;
                double *m19OD = &m19_bar_OD;
                *m19OD = (gamma / (gamma - 1)) * sqrt(2 * (pow((P019_OD / P19_OD), (-2 / gamma)) -
                                                           pow((P019_OD / P19_OD),
                                                               (-(gamma + 1) / gamma))));      // Unchoked
            } else {
                std::cout << std::endl;
                std::cout << "OFF DESIGN: Bypass nozzle is choked! \n " << std::endl;
                double *m19OD = &m19_bar_OD;
                *m19OD = (gamma / (sqrt(gamma - 1))) *
                         (pow(0.5 * (gamma + 1), ((-1 * (gamma + 1)) / (2 * (gamma - 1)))));         // Choked
            }

            double D_CoreN_OD = (m9_bar_OD * P05_OD) / sqrt(cpe * T05_OD);
            double D_LP_OD = D_CoreN_OD / phi2;
            double D_HP_ODb = D_HP_D * D_LP_OD / D_LP_D;


            // The variables D_HP_ODb & D_HP_ODa should be of equal value
            double *error = &ERR;
            *error = 100 * (D_HP_ODa - D_HP_ODb) / D_HP_ODa;

            std::cout << "\nBEFORE " << std::endl;
            std::cout << "ERROR = " << ERR << std::endl;

            if (!std::isnan(ERR)) {
                if (ERR < 0) {
                    double *pi_fan_min = &pif_min;
                    *pi_fan_min = pif_OD;
                    double *pi_fan_OD = &pif_OD;
                    *pi_fan_OD = pif_OD + 0.5 * (pif_max - pif_OD);

                    std::cout << "pi_fan_min = " << pif_min << std::endl;
                    std::cout << "pi_fan_max = " << pif_max << std::endl;
                    std::cout << "pi_fan_OD = " << pif_OD << std::endl;
                } else {
                    double *pi_fan_max = &pif_max;
                    *pi_fan_max = pif_OD;
                    double *pi_fan_OD = &pif_OD;
                    *pi_fan_OD = pif_min + 0.5 * (pif_OD - pif_min);

                    std::cout << "pi_fan_min = " << pif_min << std::endl;
                    std::cout << "pi_fan_max = " << pif_max << std::endl;
                    std::cout << "pi_fan_OD = " << pif_OD << std::endl;
                }
            } else {
                double *pi_fan_max = &pif_max;
                *pi_fan_max = pif_OD;
                double *pi_fan_OD = &pif_OD;
                *pi_fan_OD = pif_min + 0.5 * (pif_max - pif_min);

                std::cout << "pi_fan_min = " << pif_min << std::endl;
                std::cout << "pi_fan_max = " << pif_max << std::endl;
                std::cout << "pi_fan_OD = " << pif_OD << std::endl;
            }


            count++;
            std::cout << "count = " << count << "\n" << std::endl;

            if (std::isnan(ERR)) {
                double *error2 = &ERR;
                *error2 = 4;
            }

            double m_core_rat = D_HP_ODa / D_HP_D;
            double m_byp_rat = D_BypN_OD / D_BypN_D;
            double m2_bar_OD = (sqrt(cp * T02_OD) / P02_OD) *
                               ((m4_bar * (1 / D_HP_ODa) * m_core_rat * P03_OD / (sqrt(cpe * T04_OD))) +
                                (m19_bar_OD * (1 / D_BypN_OD) * m_byp_rat * P013_OD / sqrt(cp * T013_OD)));
            double m23_bar_OD = (sqrt(cp * T023_OD) * m4_bar * (1 / D_HP_ODa) * m_core_rat * P03_OD) /
                                (sqrt(cpe * T04_OD) * P023_OD);
            double m2cor_OD_ND = m2_bar_OD / m2_bar;
            double m23cor_OD_ND = m23_bar_OD / m23_bar;

            double f_OD = (cpe * (T04_OD - 298) - cp * (T03_OD - 298)) / (Qr - cpe * (T04_OD - 298));
            double Fs_OD = (1 / (1 + B_OD)) * (V9_OD - Vf_OD) + (B_OD / (1 + B_OD)) * (V19_OD - Vf_OD);
            double Fs_rat = Fs_OD / Fs;
            double FsG_OD = (1 / (1 + B_OD)) * V9_OD + (B_OD / (1 + B_OD)) * V19_OD;
            double SFC_OD = f_OD / (Fs_OD * (1 + B_OD));
            double SFC_rat = SFC_OD / sfc;
            double eta_p_OD = (Fs_OD * Vf_OD) / (Fs_OD * Vf_OD + 0.5 * (1 / (1 + B_OD)) * pow((V9_OD - Vf_OD), 2) +
                                                 0.5 * (B_OD / (1 + B_OD)) * pow((V19_OD - Vf_OD), 2));
            double eta_t_OD = (Fs_OD * Vf_OD + 0.5 * (1 / (1 + B_OD)) * pow((V9_OD - Vf_OD), 2) +
                               0.5 * (B_OD / (1 + B_OD)) * pow((V19_OD - Vf_OD), 2)) / ((f_OD * Qr) / (1 + B_OD));
            double eta_o_OD = (Fs_OD * Vf_OD) / ((f_OD * Qr) / (1 + B_OD));

            double Fn_rat = (D_CoreN_OD / D_CoreN_D) * ((B_OD * (V19_OD - Vf_OD) + (V9_OD - Vf_OD))) /
                            (B_D * (V19 - V_f) + (V9 - V_f));
            double Fg_rat = (D_CoreN_OD / D_CoreN_D) * (B_OD * V19_OD + V9_OD) / ((B_D * V19) + V9);

            double mdot9_OD = D_CoreN_OD * A9;
            double mdot19_OD = D_BypN_OD * A19;
            double mdot_e = mdot9_OD + mdot19_OD;
            double mdot_f_OD = mdot9_OD * f_OD;
            double mdot_a_OD = mdot9_OD - mdot_f_OD;
            double f_OD_new = mdot_f_OD/(0.5 * mdot_a_OD);
            double phi = f_OD_new/0.0676;
            double F_n = Fs_OD * mdot_e;

            /*
            std::cout << Mf_OD << "," << Vf_OD << "," << std::setprecision(10) << Ta_OD << "," << Pa_OD << "," << T02_OD << "," << P02_OD << ","
                         << T023_OD << "," << P023_OD << "," << T03_OD << "," << P03_OD
                         << "," << piLPC_OD << "," << piHPC_OD << "," << piOPR_OD << "," << T04_OD << "," << P04_OD << ","
                         << Tratio_OD << "," << T045_OD << "," << P045_OD << ","
                         << T05_OD << "," << P05_OD << "," << piLPT_OD << "," << piHPT_OD << "," << T09_OD << ","
                         << P9_OD << "," << T9_OD << "," << m9_bar_OD << "," << V9_OD << ","
                         << Pcrit_core << "," << y_core_OD << "," << T013_OD << "," << P013_OD << "," << pif_OD << ","
                         << T019_OD << "," << P019_OD << ","
                         << T19_OD << "," << P19_OD << "," << V19_OD << "," << m19_bar_OD << "," << Pcrit_bypass << "," << y_bypass_OD << ","
                         << F_n << ","  << std::setprecision(10) << SFC_OD << "," << Fs_OD << std::setprecision(8) << "," << f_OD << "," << B_OD << ","
                         << eta_p_OD << "," << eta_t_OD << "," << eta_o_OD << ","
                         << m2_bar_OD << "," << m23_bar_OD << "," << D_HP_ODa << "," << D_HP_ODb << "," << D_LP_OD << "," << D_CoreN_OD << "," << D_BypN_OD << ","
                         << m_core_rat << "," << m_byp_rat << ","
                         << mdot9_OD << "," << mdot19_OD << std::endl;
            */

            if (std::abs(ERR) < 1e-6) {
                std::cout << "Off-design: Flight Mach number = " << Mf_OD << std::endl;
                std::cout << "Off-design: Flight speed = " << Vf_OD << std::endl;
                std::cout << "Off-design: Ta = " << Ta_OD << std::endl;
                std::cout << "Off-design: Pa = " << Pa_OD << std::endl;

                std::cout << "\n OFF-DESIGN - CORE " << std::endl;
                std::cout << "============================ " << std::endl;
                std::cout << "Off-design: h (m) = " << h_OD << std::endl;
                std::cout << "Off-design: T_02 = " << T02_OD << std::endl;
                std::cout << "Off-design: P_02 = " << P02_OD << std::endl;
                std::cout << "Off-design: T_013 = " << T013_OD << std::endl;
                std::cout << "Off-design: P_013 = " << P013_OD << std::endl;
                std::cout << "Off-design: T_023 = " << T023_OD << std::endl;
                std::cout << "Off-design: P_023 = " << P023_OD << std::endl;
                std::cout << "Off-design: T_03 = " << T03_OD << std::endl;
                std::cout << "Off-design: P_03 = " << P03_OD << std::endl;
                std::cout << "Off-design: T_04 = " << T04_OD << std::endl;
                std::cout << "Off-design: P_04 = " << P04_OD << std::endl;
                std::cout << "Off-design: T_045 = " << T045_OD << std::endl;
                std::cout << "Off-design: P_045 = " << P045_OD << std::endl;
                std::cout << "Off-design: T_05 = " << T05_OD << std::endl;
                std::cout << "Off-design: P_05 = " << P05_OD << std::endl;
                std::cout << "Off-design: T_09 = " << T09_OD << std::endl;
                std::cout << "Off-design: P_9 = " << P9_OD << std::endl;
                std::cout << "Off-design: T_9 = " << T9_OD << std::endl;
                std::cout << "Off-design: V_9 = " << V9_OD << std::endl;

                std::cout << "\nOFF-DESIGN - BYPASS " << std::endl;
                std::cout << "============================ " << std::endl;
                std::cout << "Off-design: T_019 = " << T019_OD << std::endl;
                std::cout << "Off-design: P_019 = " << P019_OD << std::endl;
                std::cout << "Off-design: P_19 = " << P19_OD << std::endl;
                std::cout << "Off-design: T_19 = " << T19_OD << std::endl;
                std::cout << "Off-design: V_19 = " << V19_OD << std::endl;
                std::cout << "Off-design: BPR = " << B_OD << std::endl;

                std::cout << "\nOff-design: Fan PR = " << pif_OD << std::endl;
                std::cout << "Off-design: LPC PR = " << piLPC_OD << std::endl;
                std::cout << "Off-design: HPC PR = " << piHPC_OD << std::endl;
                std::cout << "Off-design: OPR = " << piOPR_OD << std::endl;
                std::cout << "Off-design: LPT PR = " << piLPT_OD << std::endl;
                std::cout << "Off-design: HPT PR = " << piHPT_OD << std::endl;
                std::cout << "Off-design: T_04/T_02 = " << Tratio_OD << std::endl;
                std::cout << "Off-design: Nozzle PR [CORE] = " << y_core_OD << std::endl;
                std::cout << "Off-design: Nozzle PR [BYPASS] = " << y_bypass_OD << std::endl;

                std::cout << "\nOff-design: mass flow per unit area - HP (kg/sm^2) [1] = " << D_HP_ODa << std::endl;
                std::cout << "Off-design: mass flow per unit area - HP (kg/sm^2) [2] = " << D_HP_ODb << std::endl;
                std::cout << "Off-design: mass flow per unit area - LP (kg/sm^2) = " << D_LP_OD << std::endl;
                std::cout << "Off-design: mass flow per unit area - CORE (kg/sm^2) = " << D_CoreN_OD << std::endl;
                std::cout << "Off-design: mass flow per unit area - BYPASS (kg/sm^2) = " << D_BypN_OD << std::endl;
                std::cout << "Off-design: m9_bar = " << m9_bar_OD << std::endl;
                std::cout << "Off-design: m19_bar = " << m19_bar_OD << std::endl;
                std::cout << "Off-design: Core mass flow ratio = " << m_core_rat << std::endl;
                std::cout << "Off-design: Bypass mass flow ratio = " << m_byp_rat << std::endl;
                std::cout << "Off-design: m2_bar = " << m2_bar_OD << std::endl;
                std::cout << "Off-design: m23_bar = " << m23_bar_OD << std::endl;
                std::cout << "Off-design: m2_bar ratio = " << m2cor_OD_ND << std::endl;
                std::cout << "Off-design: m23_bar ratio = " << m23cor_OD_ND << std::endl;
                std::cout << "Off-design: mdot_9 = " << mdot9_OD << std::endl;
                std::cout << "Off-design: mdot_19 = " << mdot19_OD << std::endl;

                std::cout << "\nOff-design: Fuel-air ratio = " << f_OD << std::endl;
                std::cout << "Off-design: Specific thrust (F_s) = " << Fs_OD << std::endl;
                std::cout << "Off-design: SFC = " << SFC_OD << std::endl;
                std::cout << "Off-design: Net thrust (F_n) = " << F_n << std::endl;
                std::cout << "Off-design: Specific thrust ratio = " << Fs_rat << std::endl;
                std::cout << "Off-design: Net thrust ratio = " << Fn_rat << std::endl;
                std::cout << "Off-design: Gross thrust ratio = " << Fg_rat << std::endl;
                std::cout << "Off-design: SFC ratio = " << SFC_rat << std::endl;
                std::cout << "Off-design: Propulsive efficiency = " << eta_p_OD << std::endl;
                std::cout << "Off-design: Thermal efficiency = " << eta_t_OD << std::endl;
                std::cout << "Off-design: Overall efficiency = " << eta_o_OD << std::endl;

                outputFile << Mf_OD << "," << Vf_OD << "," << std::setprecision(10) << Ta_OD << "," << Pa_OD << ","
                           << T02_OD << "," << P02_OD << ","
                           << T023_OD << "," << P023_OD << "," << T03_OD << "," << P03_OD
                           << "," << piLPC_OD << "," << piHPC_OD << "," << piOPR_OD << "," << T04_OD << ","
                           << P04_OD << ","
                           << Tratio_OD << "," << T045_OD << "," << P045_OD << ","
                           << T05_OD << "," << P05_OD << "," << piLPT_OD << "," << piHPT_OD << "," << T09_OD << ","
                           << P9_OD << "," << T9_OD << "," << V9_OD << ","
                           << Pcrit_core << "," << y_core_OD << "," << T013_OD << "," << P013_OD << "," << pif_OD
                           << ","
                           << T019_OD << "," << P019_OD << ","
                           << T19_OD << "," << P19_OD << "," << V19_OD << "," << Pcrit_bypass << "," << y_bypass_OD
                           << ","
                           << F_n << "," << std::setprecision(10) << SFC_OD << "," << Fs_OD << std::setprecision(8)
                           << "," << f_OD << "," << B_OD << ","
                           << eta_p_OD << "," << eta_t_OD << "," << eta_o_OD << "," << SFC_rat << "," << Fs_rat
                           << ","
                           << m2_bar_OD << "," << m23_bar_OD << "," << D_HP_ODa << "," << D_HP_ODb << "," << D_LP_OD
                           << "," << D_CoreN_OD << "," << D_BypN_OD << ","
                           << m_core_rat << "," << m_byp_rat << "," << m2cor_OD_ND << "," << m23cor_OD_ND << ","
                           << mdot9_OD << "," << mdot19_OD << "," << mdot_f_OD << "," << mdot_a_OD << "," << f_OD_new
                           << "," << phi << std::endl;

            }


        } while (std::abs(ERR) > 1e-6);
    }

                outputFile.close();
                std::cout << "\nData exported to output.csv" << std::endl;
        } else {
            //std::cerr << "\nError opening file." << std::endl;
        }
    };
    //};

    void callOffDesignFunction(){
        myOffDesignFunction();
    }

    Off_Design() = default;
    ~Off_Design() = default;
};
#endif //TURBOFAN_MANOEUVRES_SIMPLE_OFF_DESIGN_H

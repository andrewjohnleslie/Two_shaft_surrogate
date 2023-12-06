//
// Created by Andrew on 21/11/2023.
//

#ifndef TURBOFAN_MANOEUVRES_SIMPLE_DESIGN_H
#define TURBOFAN_MANOEUVRES_SIMPLE_DESIGN_H

#include <iostream>
#include <cmath>

class Design{

public:

    //Constants
    double gamma = 1.4;
    double gamma_e = 1.3;
    double cp = 1005;
    double R = 287;
    //double cpe = 1243.7;                           // USE THIS FOR MATLAB
    double cpe = (gamma_e * R)/(gamma_e - 1);        // cpe = 1244 J/kgK     USE THIS FOR EXCEL
    double Qr = 43000000;


    //AMBIENT CONDITIONS - DESIGN (CRUISE: h = 35,000 ft [10.668km], Mf = 0.7)
    double P_std = 101325;
    double T_std = 288.15;

    //Design Conditions - Ambient
    double h_ft = 0; //35000;                   // h in ft
    double h = h_ft * 0.3048;              // CONVERTS ALTITUDE TO KM
    //double h = 12000;                    // h in m

    double T_a = T_std - (0.0065 * h);
    //double t = (6356577 * h)/(6356577 + h);
    //double T_a = T_std - (6.5 * (t/1000));
    double P_a = P_std * (pow( (1 - (0.0065 * (h/T_std))), 5.2561));

    double P_trop = 22625.04;
    double g = 9.81;

    //double T_a = 216.65;
    double ind = ((-1 * g)/(R * T_a)) * (h - 11000);
    //double P_a = P_trop * exp(ind);

    double M_f = 0;
    double a = pow((gamma * R * T_a), 0.5);
    double V_f = a * M_f;

    // DESIGN PARAMETERS
    //========================
    double pi_f = 1.3;  //1.5;
    double pi_LP = 2.4;  //2.5;
    double pi_HP = 6.666666667;  //18;
    double T_04 = 1305;   //1500;
    double pi_LPT = 1.94972107;  //11.13;
    double m_in = 185.065687;   //392;

    /*
    double pi_f = 1.6731324; //1.6527886;  //sheet 6 & sheet 5
    double pi_LP = 3.1280177; //3.0517064;
    double pi_HP = 19.942067; //19.862484;
    double T_04 = 1500;
    double pi_LPT = 10.720268; //10.555515;
    double m_in = 251.247356; //304.971953;
*/

// INTAKE
// ====================
// ====================
    double T_02 = T_a * (1 + ((gamma - 1) / 2) * pow(M_f, 2));
    double P_02 = P_a * pow((1 + ((gamma - 1) / 2) * pow(M_f, 2)),gamma/(gamma - 1));
    double rho = P_02 / (R * T_02);


// COMPRESSOR
// ===================
// ===================

    double e_c = 0.9;

    //FAN
    double P_013 = pi_f * P_02;
    double T_013 = T_02 * pow(pi_f, (gamma - 1)/(e_c * gamma));

    //LPC
    double P_023 = pi_LP * P_02;
    double T_023 = T_02 * pow(pi_LP, (gamma - 1)/(e_c * gamma));

    //HPC
    double P_03 = pi_HP * P_023;
    double T_03 = T_023 * pow(pi_HP, (gamma - 1)/(e_c * gamma));

    double kt = (T_023 - T_02)/(T_013 - T_02);
    double pi_c = P_03/P_02;


// COMBUSTOR
// ====================
// ====================

    double P_04 = P_03;                 // Pressure loss is negligible

    double f = (cpe * (T_04 - 298) - cp * (T_03 - 298))/(Qr - cpe * (T_04 - 298));

    double m4_bar = gamma_e/pow(gamma_e - 1, 0.5) * pow(0.5 * (gamma_e + 1), (-1 * (gamma_e + 1)/(2 * (gamma_e - 1))));

    double piLPC = P_023/P_02;
    double piHPC = P_03/P_023;
    double piOPR = piLPC * piHPC;

// TURBINE
// ====================
// ====================

    double e_t = 0.9;
    double T_ratio = T_04/T_02;

    //HPT
    double T_045 = T_04 - ((cp/cpe) * (T_03 - T_023));
    double P_045 = P_04 * (pow((T_045/T_04), (gamma_e/(e_t * (gamma_e - 1)))));         //Issues arise from here - values slightly diverge ~20Pa out
    double m45_bar = m4_bar;

    //LPT
    double P_05 = P_045 / pi_LPT;
    double T_05 = T_045 * pow((1/pi_LPT), ((e_t * (gamma_e - 1))/gamma_e));

    double x = P_05/P_04;                       // Turbine pressure ratio

    double piLPT = (P_045/P_05);
    double piHPT = (P_04/P_045);

    double kHP = 1 - (T_045/T_04);

// NOZZLE
// ====================
// ====================

    //CORE
    double Pcrit_core = pow(0.5 * (gamma_e + 1), (gamma_e/(gamma_e - 1)));
    double P9 = P_a;
    double T9 = T_05 * pow((P9/P_05), ((gamma_e - 1)/gamma_e));
    double V9 = pow(2 * cpe * (T_05 - T9), 0.5);

    double y_core = P_05/P_a;                        // Nozzle pressure ratio
    double T_09 = T_05;


    //double m9_bar = 1.389;

    void DesignCoreChoke() {
        if (y_core < Pcrit_core){                   // If P05/Pa < P05/Pcrit (Core Nozzle PR < Critical PR)
            std::cout << std::endl;
            std::cout << "Core nozzle is not choked! \n" << std::endl;

            //double *m9D = &m9_bar;
            //      *m9D = (gamma_e/(gamma_e - 1)) * pow(2 * (pow( P_05/P9, -2/gamma_e) - pow(P_05/P9, (-1 * (gamma_e + 1)/gamma_e))), 0.5);      // Unchoked
        }else{
            std::cout << std::endl;
            std::cout << "Core nozzle is choked! \n " << std::endl;

            // double *m9D =  &m9_bar;
            //       *m9D = m4_bar;
        }
    }

    void callDesignCoreChoke(){
        DesignCoreChoke();
    }

    //double m9 = m9_bar;
    //double m9_bar = (gamma_e/(gamma_e - 1)) * pow(2 * (pow( P_05/P9, -2/gamma_e) - pow(P_05/P9, (-1 * (gamma_e + 1)/gamma_e))), 0.5);      // Unchoked
    double m9_bar = m4_bar;                                                                                                                                // Choked

    //BYPASS
    double Pcrit_bypass = pow(0.5 * (gamma + 1), (gamma/(gamma - 1)));
    double P_019 = P_013;       //Pressure loss is negligible
    double T_019 = T_013;
    double P19 = P_a;
    double T19 = T_019 * pow((P19/P_019), ((gamma - 1)/gamma));
    double V19 = pow(2 * cp * (T_019 - T19), 0.5);

    double y_bypass = P_019/P_a;

    void DesignBypassChoke() {
        if (y_bypass < Pcrit_bypass){                   // If P05/Pa < P05/Pcrit (Core Nozzle PR < Critical PR)
            std::cout << std::endl;
            std::cout << "Bypass nozzle is not choked! \n" << std::endl;
        }else{
            std::cout << std::endl;
            std::cout << "Bypass nozzle is choked! \n " << std::endl;
        }
    }

    void callDesignBypassChoke(){
        DesignBypassChoke();
    }

    double m19_bar = (gamma/(gamma - 1)) * pow(2 * (pow((P_019/P19), (-2/gamma)) - pow((P_019/P19), (-1 * (gamma + 1)/gamma))), 0.5);       //Unchoked
    //double m19_bar = (gamma/pow((gamma - 1), 0.5)) * pow(0.5 * (gamma + 1), (-1 * (gamma + 1)/(2 * (gamma - 1))));                //Choked

    double D_HP_D = (m4_bar * P_04)/pow(cpe * T_04, 0.5);                  //as choked - mass flow rate per unit area (kg/sm^2)
    double D_LP_D = (m45_bar * P_045)/pow(cpe * T_045, 0.5);               //as choked
    double D_BypN_D = (m19_bar * P_019)/pow(cp * T_019, 0.5);
    double D_CoreN_D = (m9_bar * P_05)/pow(cpe * T_05, 0.5);
    double D_HP = 1/D_HP_D;                                                      //as choked - (m^2 s/kg)
    double D_LP = 1/D_LP_D;
    double D_BypN = 1/D_BypN_D;
    double D_CoreN = 1/D_CoreN_D;

    double m2_bar = (sqrt(cp * T_02) / P_02) * ((m4_bar * (1/D_HP_D) * P_03/(sqrt(cpe * T_04))) + (m19_bar * (1/D_BypN_D) * P_013/sqrt(cp * T_013)));
    double m23_bar = (sqrt(cp * T_023) * m4_bar * (1/D_HP_D) * P_03)/(sqrt(cpe * T_04) * P_023);


    double B_D = ((cpe * (T_045 - T_05)) - (cp * (T_023 - T_02)))/(cp * (T_013 - T_02));            // Bypass Ratio - DESIGN
    double phi1 = B_D * D_HP_D/D_BypN_D;            // Area ratio [A19/A4]
    double phi2 = D_CoreN_D/D_LP_D;                 // Area ratio [A45/A9]
    double A4 = D_HP * m_in * (1/(B_D + 1));
    double A45 = D_LP * m_in * (1/(B_D + 1));
    double A9 = D_CoreN * m_in * (1/(B_D + 1));
    double A19 = D_BypN * m_in * (B_D/(B_D + 1));

    double mdot9 = D_CoreN_D * A9;
    double mdot19 = D_BypN_D * A19;
    double mdote = mdot9 + mdot19;

    // PERFORMANCE
    double Fs = (1/(1 + B_D)) * (V9 - V_f) + (B_D/(1 + B_D)) * (V19 - V_f);
    double Fn = mdote * Fs;
    double sfc = f/(Fs * (1 + B_D));

    double eta_p = (Fs * V_f)/(Fs * V_f + (0.5 * (1/(1 + B_D))) * pow((V9 - V_f), 2) + 0.5 * (B_D/(1 + B_D)) * pow((V19 - V_f),2));
    double eta_t = (Fs * V_f + (0.5 * (1/(1 + B_D))) * pow((V9 - V_f), 2) + (0.5 * (B_D/(1 + B_D))) * pow((V19 - V_f), 2))/((f * Qr)/(1 + B_D));
    double eta_o = (Fs * V_f)/((f * Qr)/(1 + B_D));
    //double eta_o = eta_p * eta_t;

    Design() = default;
    ~Design() = default;
};

#endif //TURBOFAN_MANOEUVRES_SIMPLE_DESIGN_H

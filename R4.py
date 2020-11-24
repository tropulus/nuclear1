import math
from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS) # m/kg/sec/°C/bar/W


class Reaktor:

    def __init__(self, termiskEffekt, drifttryck, n_bransleelement, branslevikt, P, anrikning, N_Th232_N_U238):

        # Reaktorparametrar
        self.termiskEffekt = termiskEffekt  # [W]
        self.drifttryck = drifttryck  # [MPa]
        self.n_bransleelement = n_bransleelement  # Antal bransleelement
        self.branslevikt = branslevikt  # Vikt per bransleelement [kg]
        self.bransleelement = 157
        self.fuel_w = self.branslevikt * self.bransleelement
        self.radie = 0.41  # cm
        self.B1_th = 0.0082768
        self.fuel_T = 800  # K
        self.vu_vm = 3.02  # volymförhållande fr. test.py
        self.xi = 0.91  # KSU s.82
        self.stavar = 264
        self.langd = 3.42 * 100  # uträknat med data från specifikationen\bransleelement densitet & n_stavar
        self.thermCon = 0.024  # W/cm*K
        self.tot_kyl_flow = 142223  # kg/s
        self.sek_kontakt = 3.02165418695  # sekunder som kylvattnet kommer ha kontakt med samma stav

        self.P = P  # Ickeläckagefaktor
        self.anrikning = anrikning  # Anrikningsgrad
        self.p_U = 1
        self.p_Th = 1
        self.N_Th232_N_U238 = N_Th232_N_U238  # Förhållandet mellan thorium-232 och uran-238
        self.rho_UO2 = 10.97  # Densitet uranoxid g*cm^-3
        self.rho_ThO2 = 10  # Densitet Thouriumoxid g*cm^-3
        self.nu = 2.43  # Antalet neutroner som frigörs vid fission
        self.epsilon = 1.06  # Snabba fissionsfaktorn
        self.k = 1
        self.reak = 0
        self.l = 0.0001
        self.timeStep=1E-3
        self.vatten_temp = 282.7  # kylvatten temp (t_in)
        self.U_vatten = 14.3  # W/cm/K
        self.c_p_UO2 = 0.4 * 1000  # J/(kgּK)
        self.c_p_H2O = 0.419 * 1000  # J/(kgּK)
        self.rho_set_w = steamTable.rho_pt(self.drifttryck*10, self.vatten_temp)


        # Antal atomkärnor
        self.N_Pa233 = 0
        self.N_Pu239 = 0
        self.N_U233 = 0
        self.N_U235 = self.calc_atom_karnor(self.rho_UO2, 235 + 2*16)*self.anrikning
        self.N_U238 = self.calc_atom_karnor(self.rho_UO2, 235 + 2 * 16) * (1-self.anrikning) * (1-self.N_Th232_N_U238)
        self.N_Th232 = self.calc_atom_karnor(self.rho_UO2, 232 + 2 * 16) * (1-self.anrikning) * (self.N_Th232_N_U238)

        # Halveringstider
        self.halveringstid_Pa233 = 26.975 * 24 * 3600  # [s]

        # Tvärsnitt
        self.barn = 1E-24
        self.sig_235_f = 576*self.barn
        self.sig_235_g = 97.6 * self.barn
        self.sig_239_f = 801 * self.barn
        self.sig_239_g = 281 * self.barn
        self.sig_238_g = 2.6 * self.barn
        self.sig_238_f = 1.76E-5*self.barn
        self.sig_232_a = 7.3*self.barn
        self.sig_233_f = 514*self.barn
        self.sig_233_g = 42*self.barn
        self.sig_233_a = self.sig_233_f + self.sig_233_g
        self.sig_235_a = self.sig_235_f + self.sig_235_g
        self.sig_238_a = self.sig_238_f + self.sig_238_g
        self.sig_239_a = self.sig_239_f + self.sig_239_g
        self.sigma_tot_w = 44 * self.barn  # dubbelkolla den här (för vatten)

    def calc_konversion(self):  # vi kommer ha en för båda
        self.c_U = (self.sig_238_a*self.N_U238+self.sig_232_a*self.N_Th232)/(self.sig_235_a*self.N_U235+self.sig_233_a*self.N_U233) \
                 + ((self.sig_235_f*self.N_U235 + self.sig_233_f*self.N_U233)/(self.sig_235_a*self.N_U235 + self.sig_233_a*self.N_U233))\
                 * self.epsilon*self.nu*self.P*(1-self.p)
        self.c_Th = self.c_U

    def calc_p(self):  # Beräkning av resonaspassagefaktor
        B1_U = 6.1 * 10 ** -3 + 0.94 * 10 ** -2 / (self.radie * self.rho_UO2)
        B1_Th = self.B1_th  # uträknat med integralförhållande
        sigma_300K_U = (3.0 + 39.6 / math.sqrt(self.radie * self.rho_UO2)) * 10 ** -24  # cm^2
        S = self.radie * 2 * math.pi * 1  # yta bränslekuts
        m = self.radie ** 2 * math.pi * 1 * self.rho_ThO2  # volym * densitet = massa
        sigma_300K_Th = (6.5 + 15.6 / math.sqrt(S/m)) * 10 ** -24  # cm^2
        sigma_fuel_fuel_T = sigma_300K_U * (1 + B1_U * (math.sqrt(self.fuel_T) - math.sqrt(300)))
        sigma_fuel_T_Th = sigma_300K_Th * (1 + B1_Th * (math.sqrt(self.fuel_T) - math.sqrt(300)))
        self.p_U = math.exp(-(1 - self.anrikning) * self.N_U238 * sigma_fuel_fuel_T * self.vu_vm /(self.xi * self.sigma_tot_w * self.calc_atom_karnor(1, 18)))
        self.p_Th = math.exp(-(1 - self.anrikning) * self.N_Th232 * sigma_fuel_T_Th * self.vu_vm /
                     (self.xi * self.sigma_tot_w * self.calc_atom_karnor(1, 18)))

    def calc_fission(self):
        denominator_f = self.N_U235 * self.sig_235_f + self.N_Pu239*self.sig_239_f + self.N_U233*self.sig_233_f
        chans_235 = self.N_U235 * self.sig_235_f / denominator_f
        chans_233 = self.N_U233 * self.sig_233_f / denominator_f
        fission_235 = self.N_U235 * self.sig_235_f / denominator_f * self.FR * 3600  # fissionerade 235
        fission_233 = self.N_U233 * self.sig_233_f / denominator_f * self.FR * 3600  # fissionerade 235
        fission_239 = (1 - chans_235 - chans_233) * self.FR * 3600  # fissionerade 239

        absorption_235 = self.N_U235 * self.sig_235_a * self.neutronflux
        absorption_233 = self.N_U233 * self.sig_233_a * self.neutronflux

        self.N_Pu239 += absorption_235*self.c_U - fission_239
        self.N_Pa233 += absorption_233*self.c_Th - self.N_Pa233 * math.exp(-3600 / self.halveringstid_Pa233)
        self.N_U233 += - fission_233 + self.N_Pa233 * math.exp(-3600 / self.halveringstid_Pa233)
        self.N_U235 -= fission_235
        self.N_U238 -= fission_235*self.c_U
        self.N_Th232 -= fission_233*self.c_Th

    def calc_eta(self):  # Beräkna snabba fissionsfaktorn
        self.eta = self.N_U235*self.nu/(self.N_U235*self.sig_235_a + self.N_U238*self.sig_238_a)

    def calc_FR(self):  # Beräknar fissionsraten
        self.FR = self.termiskEffekt / (3.2E-11) * self.rho_UO2/\
                  (self.branslevikt*self.n_bransleelement*1000)  # Konveterar vikten till gram, beräknar fissionsraten
        self.neutronflux = self.FR / (self.calc_atom_karnor(
            self.rho_UO2, 233) * self.sig_233_f + self.calc_atom_karnor(self.rho_UO2, 235) * self.sig_235_f)

    def calc_n_phi(self):  # Beräknar neutronflödet
        self.n_phi = self.FR / (self.N_U235*self.sig_235_f)

    def calc_atom_karnor(self, densitet, atommassa):  # Beräknar neutronflödet
        u = 1.66043E-24
        N = densitet/(atommassa*u)
        return N

    def calc_anrikning(self):  # Beräknar anrikning
        self.anrikning = (self.N_U235 + self.N_Pu239 + self.N_U233)/(self.N_U235 + self.N_Pu239 + self.N_U233 +
                                                                     self.N_U238 + self.N_Th232)

    def calc_reaktivitet(self):
        self.reak = (self.k - 1) / self.k

    def calc_effekt(self):
        self.termiskEffekt = self.termiskEffekt*math.exp(self.reak*3600/self.l)  # W

    def calc_lin_heat_rate(self):
        self.lin_Q = self.termiskEffekt/self.bransleelement/self.stavar/self.langd  # W/cm
        self.tempDiff = self.lin_Q/(4 * math.pi * self.thermCon)  # slide 10 F10
        q_water = self.U_vatten * (self.fuel_T - self.tempDiff)
        heat_to_w = q_water * self.sek_kontakt
        m = 0.42114504425  # uträknat värde för vatten kring bränslestaven
        self.t_out = heat_to_w/(m*self.c_p_H2O) + self.vatten_temp
        self.rho_w = steamTable.rho_pt(self.drifttryck*10, (self.t_out + self.vatten_temp)/2)

    def calc_volymf(self):
        r_c = 3.355 / 2  # m, uträknad härdradie
        r_b = self.radie / 100  # m, bränslekutsradie
        area_b = r_b ** 2 * math.pi * 264 * 157
        area_s = r_b ** 2 * math.pi * 48 * 24
        area_m = r_c ** 2 * math.pi - area_b - area_s
        area_m *= self.rho_w/self.rho_set_w
        self.vu_vm = area_m/area_b

    def calcdT_dt(self):
        p_bort = self.lin_Q*self.bransleelement/self.stavar/self.langd
        self.dT_dt = ((self.termiskEffekt - p_bort) * self.timeStep) / (self.c_p_UO2 * self.fuel_w)
        self.fuel_T += self.dT_dt



def main():
    R4 = Reaktor(3292E6, 15.5, 157, 523, 0.97, 0.03, 0.10)
    R4.calc_lin_heat_rate()
    R4.calcdT_dt()
    print(R4.fuel_T)
    # for _ in range(1_000):
    #         R4.calc_konversion()
    #         R4.calc_fission()


if __name__ == "__main__":
   main()

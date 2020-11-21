import math
import numpy as np


class Reaktor:
    def __init__(self, termiskEffekt, drifttryck, n_bransleelement, branslevikt, P, anrikning, N_Th232_N_U238):

        #Reaktorparametrar
        self.termiskEffekt = termiskEffekt #[W]
        self.drifttryck = drifttryck #[MPa]
        self.n_bransleelement = n_bransleelement #Antal bränsleelement
        self.branslevikt = branslevikt #Vikt per bränsleelement [kg]
        self.radie = 0.41  # cm
        self.B1_th = 0.0082768
        self.fuel_T = 800  # K
        self.vu_vm = 3  # volymförhållande (scholars remain divided of the value of this parameter)
        self.xi = 0.92  # dubbelkolla den här

        self.P = P #Ickeläckagefaktor
        self.anrikning = anrikning #Anrikningsgrad
        self.N_Th232_N_U238 = N_Th232_N_U238 #Förhållandet mellan thorium-232 och uran-238
        self.rho_UO2 = 10.97 #Densitet uranoxid g*cm^-3
        self.rho_ThO2 = 10 #Densitet Thouriumoxid g*cm^-3
        self.nu = 2.43 #Antalet neutroner som frigörs vid fission
        self.epsilon = 1.06 #Snabba fissionsfaktorn

        #Antal atomkärnor
        self.N_Pa233 = 0
        self.N_Pu239 = 0
        self.N_U233 = 0
        self.N_U235 = self.calc_atom_karnor(self.rho_UO2, 235 + 2*16)*self.anrikning
        self.N_U238 = self.calc_atom_karnor(self.rho_UO2, 235 + 2 * 16) * (1-self.anrikning) * (1-self.N_Th232_N_U238)
        self.N_Th232 = self.calc_atom_karnor(self.rho_UO2, 232 + 2 * 16) * (1-self.anrikning) * (self.N_Th232_N_U238)

        #Halveringstider
        self.halveringstid_Pa233 = 26.975  # [Dagar]
        
        #Tvärsnitt
        self.barn = 1E-24
        self.sig_235_f = 576*self.barn
        self.sig_238_f = 1.76E-5*self.barn
        self.sig_239_f = 801 * self.barn
        self.sig_235_g = 97.6 * self.barn
        self.sig_238_g = 2.6 * self.barn
        self.sig_235_a = self.sig_235_f + self.sig_235_g
        self.sig_238_a = self.sig_238_f + self.sig_238_g
        self.sig_232_a = 7.3*self.barn
        self.sig_233_f = 514*self.barn
        self.sig_233_g = 42*self.barn
        self.sig_233_a = self.sig_233_f + self.sig_233_g
        self.sigma_tot_w = 44 * self.barn  # dubbelkolla den här (för vatten)


    def calc_konversion(self):
        self.c = (self.sig_238_a*self.N_U238+self.sig_232_a*self.N_Th232)/(self.sig_235_a*self.N_U235+self.sig_233_a*self.N_U233) \
                 + ((self.sig_235_f*self.N_U235 + self.sig_233_f*self.N_U233)/(self.sig_235_a*self.N_U235 + self.sig_233_a*self.N_U233))\
                 *self.epsilon*self.nu*self.P*(1-self.p)

    def calc_p(self): #Beräkning av resonaspassagefaktor
        B1_U = 6.1 * 10 ** -3 + 0.94 * 10 ** -2 / (self.radie * self.rho_UO2)
        B1_Th = self.B1_th  # uträknat med integralförhållande
        sigma_300K_U = (3.0 + 39.6 / math.sqrt(self.radie * self.rho_UO2)) * 10 ** -24  # cm^2
        S = self.radie * 2 * math.pi * 1  # yta bränslekuts
        m = self.radie ** 2 * math.pi * 1 * self.rho_ThO2  # volym * densitet = massa
        sigma_300K_Th = (6.5 + 15.6 / math.sqrt(S/m)) * 10 ** -24  # cm^2
        sigma_fuel_T_U = sigma_300K_U * (1 + B1_U * (math.sqrt(self.fuel_T) - math.sqrt(300)))
        sigma_fuel_T_U = sigma_300K_Th * (1 + B1_Th * (math.sqrt(self.fuel_T) - math.sqrt(300)))
        p_U = math.exp(-(1 - self.anrikning) * self.N_U238 * sigma_fuel_T_U * self.vu_vm /
                     (self.xi * self.sigma_tot_w * self.calc_atom_karnor(1, 18)))
        p_Th = math.exp(-(1 - self.anrikning) * self.N_Th232 * sigma_fuel_T_U * self.vu_vm /
                     (self.xi * self.sigma_tot_w * self.calc_atom_karnor(1, 18)))
        return p_U, p_Th



    def calc_fission(self):
        chans_235 = self.N_U235*self.sig_235_f/(self.N_U235*self.sig_235_f + self.N_Pu239*self.sig_239_f)
        fission_235 = self.N_U235*self.sig_235_f / (self.N_U235*self.sig_235_f + self.N_Pu239*self.sig_239_f) * self.FR * 3600  # fissionerade 235
        fission_239 = (1 - chans_235) * self.FR * 3600  # fissionerade 239

        self.N_Pu239 += fission_235*self.c - fission_239
        self.N_U235 -= fission_235
        self.N_U238 -= fission_235*self.c

    def calc_eta(self): #Beräkna snabba fissionsfaktorn
        self.eta = self.N_U235*self.nu/(self.N_U235*self.sig_235_a + self.N_U238*self.sig_238_a)

    def calc_FR(self): #Beräknar fissionsraten
        self.FR = self.termiskEffekt / (3.2E-11) * self.rho_UO2/(self.branslevikt*self.n_bransleelement*1000)  #Konveterar vikten till gram, beräknar fissionsraten

    def calc_n_phi(self): #Beräknar neutronflödet
        self.n_phi = self.FR / (self.N_U235*self.sig_235_f)

    def calc_atom_karnor(self, densitet, atommassa): #Beräknar neutronflödet
        u = 1.66043E-24
        N = densitet/(atommassa*u)
        return N

    def calc_anrikning(self): #Beräknar anrikning
        self.anrikning = (self.N_U235 + self.N_Pu239 + self.N_U233)/(self.N_U235 + self.N_Pu239 + self.N_U233 + self.N_U238 + self.N_Th232)

def main():
    R4 = Reaktor(3292E6, 15.5, 157, 523, 0.97, 0.03, 0.10)
    R4.calc_FR()
    R4.calc_n_phi()
    R4.calc_eta()
    R4.calc_p()  # uppdatera anrikning först

    for _ in range(1_000):
        R4.calc_konversion()
        R4.calc_fission()

if __name__ == "__main__":
    main()